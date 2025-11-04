############################
# GATK Calling Rules
# Release: 0.9.5
# Date: 31/10/2025
# Author: Antony Le Bechec
############################


# Release note
# 0.9.1beta-10/03/2015: change genome reference location, in the file %.genome
# 0.9.1.1beta-30/04/2015: VariantFiltration correction
# 0.9.2b-23/11/2015: adding gatkUG_CPSGEN_MASTR pipeline
# 0.9.3b-25/01/2016: adding gatkUG_DIAG and gatkHC_DIAG pipelines
# 0.9.3.1b-08/02/2016: adding gatkUG_DIAGGEN pipeline
# 0.9.3.2b-03/03/2016: adding gatkUG_DIAG_PARIS pipeline
# 0.9.3.3b-11/04/2016: change gatkHC pipeline to fit with standard analysis
# 0.9.3.4b-03/05/2016: Modification of gatkUG_HUSHEMATO rule
# 0.9.3.5b-04/05/2016: Rewrite rules and update release information
# 0.9.3.6b-10/11/2017: adding gatkUG_ONCOGENET pipeline
# 0.9.3.7b-22/03/2019: Add --dontUseSoftClippedBases for GATKHC
# 0.9.3.8-29/07/2022: Remove --dontUseSoftClippedBases for GATKHC
# 0.9.3.9-29/07/2022: Add --dontUseSoftClippedBases for GATKHC, add GATKUG_LONG_INDELS and GATKHC_LONG_INDELS
# 0.9.4-03/02/2023: Extract gatk4HC
# 0.9.5-31/10/2025: Change GATK4HC to GATKHC

GATKHC_FLAGS_SHARED?=--dont-use-soft-clipped-bases

###########
# gatkHC #
###########

# GATKHC Flags
GATKHC_FLAGS=$(GATK4HC_FLAGS_SHARED) \
	--interval-padding $(INTERVAL_PADDING) \
	--dbsnp $(VCFDBSNP) \
	--native-pair-hmm-threads $(THREADS_BY_CALLER) \
	--min-base-quality-score 17 \
	--min-pruning 4 \
	--max-reads-per-alignment-start 1000 \
	--standard-min-confidence-threshold-for-calling 30


%.gatkHC$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.design.bed.interval_list
	$(JAVA) $(JAVA_FLAGS) -XX:ParallelGCThreads=$(THREADS_BY_CALLER) -jar $(GATK4) \
		HaplotypeCaller \
		$(GATKHC_FLAGS) \
		-R $(GENOME) \
		-I $< \
		-O $@ \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;);
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -e $@ ]; then touch $@; fi;
	-rm -f $@.idx


RELEASE_COMMENT := "\#\# CALLING GATK '$(MK_RELEASE)': GATK tool identify variants from aligned BAM with shared parameters: GATKHC_FLAGS='$(GATKHC_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


RELEASE_COMMENT := "\#\# CALLING GATK4HC identify variants and generate *.gatkHC.vcf files with parameters: GATKHC_FLAGS='$(GATKHC_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "CALLER:gatkHC:GATK4 Haplotype Caller - by default:GATKHC_FLAGS='$(GATKHC_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
