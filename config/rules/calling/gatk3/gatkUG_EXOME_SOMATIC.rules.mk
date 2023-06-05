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
# 0.9.4-03/02/2023: Extract gatkUG EXOME SOMATIC



########################
# gatkUG_EXOME_SOMATIC #
########################

# GATKUG Flags
GATKUG_THREADS_EXOME_SOMATIC?=1
DFRAC_UG_EXOME_SOMATIC=1
MBQ_UG_EXOME_SOMATIC=17
STAND_CALL_CONF_UG_EXOME_SOMATIC=30
minIndelFrac_UG_EXOME_SOMATIC=0.01
minIndelCnt_UG_EXOME_SOMATIC=2
deletions_UG_EXOME_SOMATIC=0.01
DPMIN_UG_EXOME_SOMATIC=50
GATKUG_EXOME_SOMATIC_FLAGS= -nct $(GATKUG_THREADS_EXOME_SOMATIC) -glm BOTH \
		-minIndelFrac $(minIndelFrac_UG_EXOME_SOMATIC)  \
		-minIndelCnt $(minIndelCnt_UG_EXOME_SOMATIC) \
		-deletions $(deletions_UG_EXOME_SOMATIC) \
		-baq OFF \
		-stand_call_conf $(STAND_CALL_CONF_UG_EXOME_SOMATIC) -dfrac $(DFRAC_UG_EXOME_SOMATIC) --dbsnp $(VCFDBSNP) -mbq $(MBQ_UG_EXOME_SOMATIC) -rf BadCigar -dt NONE

%.gatkUG_EXOME_SOMATIC$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list
	$(JAVA8) $(JAVA_FLAGS) -jar $(GATK3) $(GATKUG_EXOME_SOMATIC_FLAGS) \
		-T UnifiedGenotyper \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 							# in case of no vcf creation, to not kill the pipeline
	$(BCFTOOLS) view -i "FORMAT/DP>$(DPMIN_UG_EXOME_SOMATIC)" $@.tmp > $@ 	# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 								# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx



RELEASE_COMMENT := "\#\# CALLING GATKUG EXOME_SOMATIC identify variants and generate *.gatkUG_EXOME_SOMATIC.vcf files with parameters: GATKUG_EXOME_SOMATIC_FLAGS='$(GATKUG_EXOME_SOMATIC_FLAGS)', INTERVAL_PADDING='$(INTERVAL_PADDING)', DPMIN_UG_EXOME_SOMATIC='$(DPMIN_UG_EXOME_SOMATIC)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkUG_EXOME_SOMATIC:GATK Unified Genotyper - designed for EXOME_SOMATIC discovery:GATKUG_EXOME_SOMATIC_FLAGS='$(GATKUG_EXOME_SOMATIC_FLAGS)', INTERVAL_PADDING='$(INTERVAL_PADDING)', DPMIN_UG_EXOME_SOMATIC='$(DPMIN_UG_EXOME_SOMATIC)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
