############################
# GATK Calling Rules
# Release: 0.9.4
# Date: 03/02/2023
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



###########
# gatk4HC #
###########

MINPRUNING_GATK4HC?=4
THREADS_GATK4HC?=$(THREADS_BY_CALLER)
STAND_CALL_CONF_GATK4HC=30
maxreadsperalignmentstart_GATK4HC=1000
GATK4HC_FLAGS_SHARED?=
GATK4HC_FLAGS= --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC) --min-pruning $(MINPRUNING_GATK4HC) --native-pair-hmm-threads $(THREADS_GATK4HC) $(GATK4HC_FLAGS_SHARED) --max-reads-per-alignment-start $(maxreadsperalignmentstart_GATK4HC) --standard-min-confidence-threshold-for-calling $(STAND_CALL_CONF_GATK4HC)

%.gatk4HC$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	#
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK4) HaplotypeCaller \
		$(GATK4HC_FLAGS) \
		-R $$(cat $*.genome) \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-O $@;
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -e $@ ]; then touch $@; fi;
	-rm -f $@.idx




RELEASE_COMMENT := "\#\# CALLING GATK4 '$(MK_RELEASE)': GATK tool identify variants from aligned BAM with shared parameters: GATK4='$(GATK4)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


RELEASE_COMMENT := "\#\# CALLING GATK4HC identify variants and generate *.gatk4HC.vcf files with parameters: GATK4HC_FLAGS='$(GATK4HC_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "CALLER:gatk4HC:GATK4 Haplotype Caller - by default:GATK4HC_FLAGS='$(GATK4HC_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
