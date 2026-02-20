############################
# GATK Calling Rules
# Release: 1.0.0
# Date: 19/02/2026
# Author: Antony Le Bechec
############################


# Release note
# 1.0.0-19/02/2026: Create gatkHC_RNASEQ with specific parameters for RNASeq data. No dbSNP because not same contig. Post alignment splitncigar mandatory.


###########
# gatkHC #
###########

# GATKHC Flags
GATKHC_RNASEQ_FLAGS=--interval-padding $(INTERVAL_PADDING) \
	--native-pair-hmm-threads $(THREADS_BY_CALLER) \
	--min-base-quality-score 17 \
	--min-pruning 4 \
	--max-reads-per-alignment-start 1000 \
	--standard-min-confidence-threshold-for-calling 30 \
	--dont-use-soft-clipped-bases true \
	--recover-dangling-heads true


%.gatkHC_RNASEQ$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.design.bed.interval_list
	$(JAVA) $(JAVA_FLAGS) -XX:ParallelGCThreads=$(THREADS_BY_CALLER) -jar $(GATK4) \
		HaplotypeCaller \
		$(GATKHC_RNASEQ_FLAGS) \
		-R $(GENOME) \
		-I $< \
		-O $@ \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;);
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -e $@ ]; then touch $@; fi;
	-rm -f $@.idx


# RELEASE_COMMENT := "\#\# CALLING GATK RNASEQ '$(MK_RELEASE)': GATK tool identify variants from aligned BAM with parameters: GATKHC_RNASEQ_FLAGS='$(GATKHC_RNASEQ_FLAGS)'"
# RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


RELEASE_COMMENT := "\#\# CALLING GATKHC_RNASEQ identify variants and generate *.gatkHC_RNASEQ.vcf files with parameters: GATKHC_RNASEQ_FLAGS='$(GATKHC_RNASEQ_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "CALLER:gatkHC_RNASEQ:GATK4 Haplotype Caller - for RNASeq:GATKHC_RNASEQ_FLAGS='$(GATKHC_RNASEQ_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
