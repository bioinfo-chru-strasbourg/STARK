############################
# Variants Filtration Rules
# Release: 0.9.0.1
# Date: 25/05/2021
# Author: Antony Le Bechec
############################
MK_RELEASE="0.9.1.0"
MK_DATE="09/11/2021"

# Release notes:
# 0.9.0.0-14/10/2019: Creation
# 0.9.0.1-25/05/2021: Recheck Expressions
# 0.9.1.0-09/11/2021: Expression variables


#####################
# VariantFiltration #
#####################


# VARIABLES

FILTRATION_VCF_FILTER_EXPRESSION?=--filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "HARD_TO_VALIDATE" --filterExpression "DP == 0" --filterName "VeryVeryLowDepth" --filterExpression "DP > 0 && DP < 10" --filterName "VeryLowDepth" --filterExpression "DP >= 10 && DP < 30" --filterName "LowDepth" --filterExpression "QUAL == 0" --filterName "VeryVeryLowQual" --filterExpression "QUAL > 0 && QUAL < 30.0" --filterName "VeryLowQual" --filterExpression "QUAL >= 30.0 && QUAL < 50.0" --filterName "LowQual" --filterExpression "QD >= 0.0 && QD < 1.5" --filterName "LowQD"
FILTRATION_VCF_GENOTYPE_FILTER_EXPRESSION?=--genotypeFilterExpression "GQ == 0" --genotypeFilterName "VeryVeryLowGQ"  --genotypeFilterExpression "GQ > 0 && GQ < 50.0" --genotypeFilterName "VeryLowGQ"  --genotypeFilterExpression "GQ >= 50.0 && GQ < 90.0" --genotypeFilterName "LowGQ" --genotypeFilterExpression "DP == 0" --genotypeFilterName "VeryVeryLowDP" --genotypeFilterExpression "DP >= 0 && DP < 10" --genotypeFilterName "VeryLowDP"  --genotypeFilterExpression "DP >= 10 && DP < 30" --genotypeFilterName "LowDP"

FILTRATION_VCF_INVALIDATE_PREVIOUS_FILTERS_OPTION?=$(shell if (( $(FILTRATION_VCF_INVALIDATE_PREVIOUS_FILTERS) )); then echo " --invalidatePreviousFilters "; fi )



#%.vcf: %.filtration.vcf %.filtration.vcf.idx %.empty.vcf %.genome
%.vcf: %.filtration.vcf %.empty.vcf %.genome
	if [ ! -s $*.filtration.vcf.idx ]; then rm -f $*.filtration.vcf.idx; fi;
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R $$(cat $*.genome) --variant $< -o $@ $(FILTRATION_VCF_INVALIDATE_PREVIOUS_FILTERS_OPTION) $(FILTRATION_VCF_FILTER_EXPRESSION) $(FILTRATION_VCF_GENOTYPE_FILTER_EXPRESSION)
	if [ ! -e $@ ]; then cp $*.filtration.vcf $@; fi;
	if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	if [ ! -e $@ ]; then touch $@; fi;
	rm -f $<.idx $@.idx $*.filtration.vcf*

#if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	


RELEASE_COMMENT := "\#\# VCF FILTRATION: GATK VariantFiltration to tag FILTER info in VCF file."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_CALLING:filtration:VariantFiltration of VCF."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

#PIPELINES_COMMENT := "POST_ANNOTATION:filtration:VariantFiltration of VCF."
#PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
