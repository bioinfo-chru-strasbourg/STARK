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
# 0.9.4-03/02/2023: Extract gatkUG LONG INDELS



######################
# gatkUG_LONG_INDELS #
######################

DFRAC_UG_LONG_INDELS=1
MBQ_UG_LONG_INDELS=17
THREADS_GATKUG_LONG_INDELS?=$(THREADS_BY_CALLER)
INTERVAL_PADDING?=0
DPMIN_UG_LONG_INDELS?=4
GATKUG_INDEL_SIZE_LONG_INDELS?=20
GATKUG_LONG_INDELS_FLAGS= -nct $(THREADS_GATKUG_LONG_INDELS) -glm INDEL \
		-minIndelFrac 0.01 \
		-minIndelCnt 2 \
		-deletions 0.01 \
		-baq OFF \
		-stand_call_conf 10 -dfrac $(DFRAC_UG_LONG_INDELS) --dbsnp $(VCFDBSNP) -mbq $(MBQ_UG_LONG_INDELS) -rf BadCigar -dt NONE \
		-allowPotentiallyMisencodedQuals

%.gatkUG_LONG_INDELS$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list
	$(JAVA8) $(JAVA_FLAGS) -jar $(GATK3) $(GATKUG_LONG_INDELS_FLAGS) \
		-T UnifiedGenotyper \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 								# in case of no vcf creation, to not kill the pipeline
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK4) SelectVariants -V $@.tmp --select-type-to-include INDEL --min-indel-size $(GATKUG_INDEL_SIZE_LONG_INDELS) -O $@.tmp2;
	$(BCFTOOLS) view -i "FORMAT/DP>$(DPMIN_UG_LONG_INDELS)" $@.tmp2 > $@ 		# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 									# in case of error in previous line
	-rm -f $@.tmp* $@.idx



RELEASE_COMMENT := "\#\# CALLING GATKUG LONG INDELS identify long indels and generate *.gatkUG_LONG INDELS.vcf files with parameters: GATKUG_LONG_INDELS_FLAGS='$(GATKUG_LONG_INDELS_FLAGS)', DPMIN_UG_LONG_INDELS='$(DPMIN_UG_LONG_INDELS)', GATKUG_INDEL_SIZE_LONG_INDELS='$(GATKUG_INDEL_SIZE_LONG_INDELS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkUG_LONG_INDELS:GATK Unified Genotyper - designed for LONG INDELS discovery:GATKUG_LONG_INDELS_FLAGS='$(GATKUG_LONG_INDELS_FLAGS)', DPMIN_UG_LONG_INDELS='$(DPMIN_UG_LONG_INDELS)', GATKUG_INDEL_SIZE_LONG_INDELS='$(GATKUG_INDEL_SIZE_LONG_INDELS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
