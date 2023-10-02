############################
# MUTECT Calling Rules
# Release: 0.9.2
# Date: 03/02/2023
# Author: Antony Le Bechec
############################


#################
# MUTECT2 GATK4 #
#################

# Parameters
THREADS_GATK4_MUTECT2?=$(THREADS_BY_CALLER)
MINPRUNING_GATK4_MUTECT2?=20
MAXREADS_GATK4_MUTECT2?=1000
DPMIN_MUTECT2?=30

GATK4_MUTECT2_FLAGS_SHARED?=--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter --max-reads-per-alignment-start $(MAXREADS_GATK4_MUTECT2) --dont-use-soft-clipped-bases true --min-pruning $(MINPRUNING_GATK4_MUTECT2) --callable-depth $(DPMIN_MUTECT2) --verbosity ERROR --native-pair-hmm-threads $(THREADS_GATK4_MUTECT2)


%.MuTect2$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.dict %.design.bed.interval_list
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK4) Mutect2 $(GATK4_MUTECT2_FLAGS_SHARED) \
		-R $(GENOME) \
		-I $< \
		-tumor $$(basename $< | cut -d"." -f1) \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-O $@.tmp1;
	# Normalize
	grep "^##" $@.tmp1 | sed s/ID=TLOD,Number=A/ID=TLOD,Number=./gi > $@.tmp.vcf
	grep "^##" -v $@.tmp1 | cut -f1-10 >> $@.tmp.vcf
	# sort
	$(JAVA) -jar $(PICARD) SortVcf -I $@.tmp.vcf -O $@.tmp2.vcf -SD $$(cat $*.dict);
	# DPMIN_MUTECT2
	#$(VCFUTILS) varFilter -d $(DPMIN_MUTECT2) $@.tmp2.vcf > $@.tmp3.vcf;
	$(BCFTOOLS) view  -i 'FORMAT/DP>=$(DPMIN_MUTECT2)' $@.tmp2.vcf > $@.tmp3.vcf;
	# Empty
	if [ ! -e $@.tmp3.vcf ]; then cp $*.empty.vcf $@.tmp3.vcf; fi;
	# Copy
	if [ ! -e $@ ]; then cp $@.tmp3.vcf $@; fi;
	# Clean
	rm -f $@.tmp* $@.idx



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING: MUTECT2 with GATK4 to identify somatic variants and generate *.MuTect2.vcf files, with DP filtration. Parameters GATK4_MUTECT2_FLAGS_SHARED='$(GATK4_MUTECT2_FLAGS_SHARED)', DPMIN_MUTECT2='$(DPMIN_MUTECT2)' "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:MuTect2:MuTect2 with GATK4 - by default with DP filtration. "
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
