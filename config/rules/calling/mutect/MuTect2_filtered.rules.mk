############################
# MUTECT Calling Rules
# Release: 0.9.2
# Date: 03/02/2023
# Author: Antony Le Bechec
############################


#########################
# MUTECT2 flagged GATK4 #
#########################

# Parameters
THREADS_GATK4_MUTECT2_FILTERED?=$(THREADS_BY_CALLER)
MINPRUNING_GATK4_MUTECT2_FILTERED?=20
MAXREADS_GATK4_MUTECT2_FILTERED?=1000
DPMIN_MUTECT2_FILTERED?=30
VAF_MUTECT2_FILTERED?=0.01
VAF_MUTECT2_FILTERED_HOM?=0.80

GATK4_MUTECT2_FILTERED_FLAGS_SHARED?=--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter --max-reads-per-alignment-start $(MAXREADS_GATK4_MUTECT2_FILTERED) --dont-use-soft-clipped-bases true --min-pruning $(MINPRUNING_GATK4_MUTECT2_FILTERED) --callable-depth $(DPMIN_MUTECT2_FILTERED) --verbosity ERROR --native-pair-hmm-threads $(THREADS_GATK4_MUTECT2_FILTERED)

%.MuTect2_filtered$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.design.bed.interval_list
	# Calling by MuTect2
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK4) Mutect2 $(GATK4_MUTECT2_FILTERED_FLAGS_SHARED) \
		-R $(GENOME) \
		-I $< \
		-tumor $$(basename $< | cut -d"." -f1) \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-O $@.tmp.unfiltered.vcf;
	# Normalization
	grep "^##" $@.tmp.unfiltered.vcf | sed s/ID=TLOD,Number=A/ID=TLOD,Number=./gi > $@.tmp.unfiltered.TLOD.vcf
	grep "^##" -v $@.tmp.unfiltered.vcf | cut -f1-10 >> $@.tmp.unfiltered.TLOD.vcf
	# HOWARD VAF calculation & Filtration by BCFTOOLS
	 +if (($$($(BCFTOOLS) view -H $@.tmp.unfiltered.TLOD.vcf | wc -l ))); then \
	 	#$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.tmp.unfiltered.TLOD.vcf --output=$@.tmp.unfiltered.TLOD.HOWARD.vcf --calculation="VAF"; \
		$(HOWARD2) $(HOWARD_CONFIG_OPTIONS) --input=$@.tmp.unfiltered.TLOD.vcf --output=$@.tmp.unfiltered.TLOD.HOWARD.vcf --calculation="VAF"; \
	 	$(BCFTOOLS) view -i 'FORMAT/DP>=$(DPMIN_MUTECT2_FILTERED) && FORMAT/VAF>=$(VAF_MUTECT2_FILTERED) && FORMAT/VAF<$(VAF_MUTECT2_FILTERED_HOM)' $@.tmp.unfiltered.TLOD.HOWARD.vcf | $(BCFTOOLS) view -e 'GT="0/0"' -o $@.tmp.cleaned.vcf; \
	 else \
	 	cp $@.tmp.unfiltered.TLOD.vcf $@.tmp.cleaned.vcf; \
	 fi;
	cp $@.tmp.unfiltered.TLOD.vcf $@.tmp.cleaned.vcf;
	# Sorting and contig
	$(JAVA) -jar $(PICARD) SortVcf -I $@.tmp.cleaned.vcf -O $@ -SD $(DICT);
	# Clean
	rm -f $@.tmp* $@.idx


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING: MUTECT2 with GATK4 to identify somatic variants and generate *.MuTect2_filtered.vcf files, with DP and VAF filtration, with FilterMutectCalls filteration by removing variant WT. Parameters GATK4_MUTECT2_FILTERED_FLAGS_SHARED='$(GATK4_MUTECT2_FILTERED_FLAGS_SHARED)', DPMIN_MUTECT2_FILTERED='$(DPMIN_MUTECT2_FILTERED)', VAF_MUTECT2_FILTERED='$(VAF_MUTECT2_FILTERED)', VAF_MUTECT2_FILTERED_HOM='$(VAF_MUTECT2_FILTERED_HOM)' "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:MuTect2_filtered:MuTect2 with GATK4 and FilterMutectCalls - by default with DP, VAF and WildType filtration. "
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
