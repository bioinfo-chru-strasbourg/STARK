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
# 0.9.4-03/02/2023: Extract gatkHC SOLIDTUMOR



#####################
# gatkHC_SOLIDTUMOR #
#####################

DFRAC_HC_SOLIDTUMOR=1
MINPRUNING_SOLIDTUMOR?=20
THREADS_GATKHC_SOLIDTUMOR?=$(THREADS_BY_CALLER)
maxReadsInRegionPerSample_SOLIDTUMOR=1000
MBQ_HC_SOLIDTUMOR=17
GATKHC_SOLIDTUMOR_FLAGS= -nct $(THREADS_GATKHC_SOLIDTUMOR) -stand_call_conf 10 -dfrac $(DFRAC_HC_SOLIDTUMOR) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample_SOLIDTUMOR) --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC_SOLIDTUMOR) -minPruning $(MINPRUNING_SOLIDTUMOR)  $(GATKHC_FLAGS_SHARED)

%.gatkHC_SOLIDTUMOR$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.design.bed.interval_list 
	$(JAVA8) $(JAVA_FLAGS) -jar $(GATK3) $(GATKHC_SOLIDTUMOR_FLAGS) \
		-T HaplotypeCaller \
		-R $(GENOME) \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@;
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -e $@ ]; then touch $@; fi;
	-rm -f $@.idx



RELEASE_COMMENT := "\#\# CALLING GATKHC SOLIDTUMOR identify variants and generate *.gatkHC_SOLIDTUMOR.vcf files with parameters: GATKHC_SOLIDTUMOR_FLAGS='$(GATKHC_SOLIDTUMOR_FLAGS)', INTERVAL_PADDING='$(INTERVAL_PADDING)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkHC_SOLIDTUMOR:GATK Haplotype Caller - designed for SOLIDTUMOR variant discovery:GATKHC_SOLIDTUMOR_FLAGS='$(GATKHC_SOLIDTUMOR_FLAGS)', INTERVAL_PADDING='$(INTERVAL_PADDING)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
