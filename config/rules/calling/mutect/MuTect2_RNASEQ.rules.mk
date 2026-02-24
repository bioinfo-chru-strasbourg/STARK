############################
# MUTECT Calling Rules
# Release: 1.0.0
# Date: 19/02/2026
# Author: Antony Le Bechec
############################


# Release note
# 1.0.0-19/02/2026: Create MuTect2_RNASEQ with specific parameters for RNASeq data. No dbSNP because not same contig. Post alignment splitncigar mandatory.

#################
# MUTECT2 GATK4 #
#################

# Parameters
THREADS_GATK4_MUTECT2_RNASEQ?=$(THREADS_BY_CALLER)
MINPRUNING_GATK4_MUTECT2_RNASEQ?=20
MAXREADS_GATK4_MUTECT2_RNASEQ?=1000
DPMIN_MUTECT2_RNASEQ?=30

GATK4_MUTECT2_RNASEQ_FLAGS_SHARED?=--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
	--max-reads-per-alignment-start $(MAXREADS_GATK4_MUTECT2_RNASEQ) \
	--dont-use-soft-clipped-bases true \
	--min-pruning $(MINPRUNING_GATK4_MUTECT2_RNASEQ) \
	--callable-depth $(DPMIN_MUTECT2_RNASEQ) \
	--verbosity ERROR \
	--native-pair-hmm-threads $(THREADS_GATK4_MUTECT2_RNASEQ)


%.MuTect2_RNASEQ$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.design.bed.interval_list
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK4) Mutect2 $(GATK4_MUTECT2_RNASEQ_FLAGS_SHARED) \
		-R $(GENOME) \
		-I $< \
		-tumor $$(basename $< | cut -d"." -f1) \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-O $@.tmp1;
	# Normalize
	grep "^##" $@.tmp1 | sed s/ID=TLOD,Number=A/ID=TLOD,Number=./gi > $@.tmp.vcf
	grep "^##" -v $@.tmp1 | cut -f1-10 >> $@.tmp.vcf
	# sort
	$(JAVA) -jar $(PICARD) SortVcf -I $@.tmp.vcf -O $@.tmp2.vcf -SD $(DICT);
	# DPMIN_MUTECT2_RNASEQ
	$(BCFTOOLS) view  -i 'FORMAT/DP>=$(DPMIN_MUTECT2_RNASEQ)' $@.tmp2.vcf > $@.tmp3.vcf;
	# Empty
	if [ ! -e $@.tmp3.vcf ]; then cp $*.empty.vcf $@.tmp3.vcf; fi;
	# Copy
	if [ ! -e $@ ]; then cp $@.tmp3.vcf $@; fi;
	# Clean
	rm -f $@.tmp* $@.idx



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING: MUTECT2 RNASEQ with GATK4 to identify somatic variants on RNA Seq data and generate *.MuTect2_RNASEQ.vcf files, with DP filtration. Parameters GATK4_MUTECT2_RNASEQ_FLAGS_SHARED='$(GATK4_MUTECT2_RNASEQ_FLAGS_SHARED)', DPMIN_MUTECT2_RNASEQ='$(DPMIN_MUTECT2_RNASEQ)' "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:MuTect2_RNASEQ:MuTect2 RNASEQ with GATK4 - for RNA Seq data with DP filtration. "
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
