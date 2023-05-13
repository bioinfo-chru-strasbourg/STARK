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
# 0.9.4-03/02/2023: Extract gatkHC GERMLINE



###################
# gatkHC_GERMLINE #
###################

MINPRUNING_HC_GERMLINE?=4
THREADS_GATKHC_GERMLINE?=$(THREADS_BY_CALLER)
maxReadsInRegionPerSample_HC_GERMLINE=1000
MBQ_HC_GERMLINE=17
STAND_CALL_CONF_HC_GERMLINE=30
DFRAC_HC_GERMLINE=1
DPMIN_HC_GERMLINE=4
GATKHC_GERMLINE_FLAGS= -nct $(THREADS_GATKHC_GERMLINE) -stand_call_conf $(STAND_CALL_CONF_HC_GERMLINE) -dfrac $(DFRAC_HC_GERMLINE) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample_HC_GERMLINE) --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC_GERMLINE) -minPruning $(MINPRUNING_HC_GERMLINE)  $(GATKHC_FLAGS_SHARED)

%.gatkHC_GERMLINE$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list
	$(JAVA8) $(JAVA_FLAGS) -jar $(GATK3) $(GATKHC_GERMLINE_FLAGS) \
		-T HaplotypeCaller \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 				# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_HC_GERMLINE) $@.tmp > $@ 	# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 					# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx



RELEASE_COMMENT := "\#\# CALLING GATKHC GERMLINE identify variants and generate *.gatkHC_GERMLINE.vcf files with parameters: GATKHC_GERMLINE_FLAGS='$(GATKHC_GERMLINE_FLAGS)', INTERVAL_PADDING='$(INTERVAL_PADDING)', DPMIN_HC_GERMLINE='$(DPMIN_HC_GERMLINE)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkHC_GERMLINE:GATK Haplotype Caller - designed for GERMLINE discovery:GATKHC_GERMLINE_FLAGS='$(GATKHC_GERMLINE_FLAGS)', INTERVAL_PADDING='$(INTERVAL_PADDING)', DPMIN_HC_GERMLINE='$(DPMIN_HC_GERMLINE)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
