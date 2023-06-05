############################
# MUTECT Calling Rules
# Release: 0.9.2
# Date: 03/02/2023
# Author: Antony Le Bechec
############################


###########################
# MUTECT2 stringent GATK4 #
###########################

# Parameters
THREADS_GATK4_MUTECT2_STRINGENT?=$(THREADS_BY_CALLER)
MINPRUNING_GATK4_MUTECT2_STRINGENT?=20
MAXREADS_GATK4_MUTECT2_STRINGENT?=1000
DPMIN_MUTECT2_STRINGENT?=30
VAF_MUTECT2_STRINGENT?=0.01
VAF_MUTECT2_STRINGENT_HOM?=0.80

GATK4_MUTECT2_STRINGENT_FLAGS_SHARED?=--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter --max-reads-per-alignment-start $(MAXREADS_GATK4_MUTECT2_STRINGENT) --dont-use-soft-clipped-bases true --min-pruning $(MINPRUNING_GATK4_MUTECT2_STRINGENT) --callable-depth $(DPMIN_MUTECT2_STRINGENT) --verbosity ERROR --native-pair-hmm-threads $(THREADS_GATK4_MUTECT2_STRINGENT)


%.MuTect2_stringent$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.dict %.design.bed.interval_list
	# Calling by MuTect2
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK4) Mutect2 $(GATK4_MUTECT2_STRINGENT_FLAGS_SHARED) \
		-R $$(cat $*.genome) \
		-I $< \
		-tumor $$(basename $< | cut -d"." -f1) \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-O $@.tmp.unfiltered.vcf;
	# Filtration by MuTect2
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK4) FilterMutectCalls \
		-R $$(cat $*.genome) \
		-V $@.tmp.unfiltered.vcf \
		-O $@.tmp.FilterMutectCalls.vcf;
	# Normalization
	grep "^##" $@.tmp.FilterMutectCalls.vcf | sed s/ID=TLOD,Number=A/ID=TLOD,Number=./gi > $@.tmp.FilterMutectCalls.TLOD.vcf
	grep "^##" -v $@.tmp.FilterMutectCalls.vcf | cut -f1-10 >> $@.tmp.FilterMutectCalls.TLOD.vcf
	# Sorting
	$(JAVA) -jar $(PICARD) SortVcf -I $@.tmp.FilterMutectCalls.TLOD.vcf -O $@.tmp.FilterMutectCalls.TLOD.sorted.vcf -SD $$(cat $*.dict);
	# Filtration by BCFTOOLS
	$(BCFTOOLS) view -i 'FORMAT/DP>=$(DPMIN_MUTECT2_STRINGENT) && FORMAT/AF>=$(VAF_MUTECT2_STRINGENT) && FORMAT/AF<$(VAF_MUTECT2_STRINGENT_HOM)' $@.tmp.FilterMutectCalls.TLOD.sorted.vcf | $(BCFTOOLS) view -f PASS  -e 'GT="0/0"' -o $@;
	# Clean
	rm -f $@.tmp* $@.idx



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING: MUTECT2 with GATK4 to identify somatic variants and generate *.MuTect2_stringent.vcf files, with DP and VAF filtration, with FilterMutectCalls filteration by removing variant WT and with quality not PASS. Parameters GATK4_MUTECT2_STRINGENT_FLAGS_SHARED='$(GATK4_MUTECT2_STRINGENT_FLAGS_SHARED)', DPMIN_MUTECT2_STRINGENT='$(DPMIN_MUTECT2_STRINGENT)', VAF_MUTECT2_STRINGENT='$(VAF_MUTECT2_STRINGENT)', VAF_MUTECT2_STRINGENT_HOM='$(VAF_MUTECT2_STRINGENT_HOM)' "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:MuTect2_stringent:MuTect2 with GATK4 and FilterMutectCalls - by default with DP, VAF, WildType and FilterMutectCalls filtration. "
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
