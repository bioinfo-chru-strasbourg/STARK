############################
# Variants Filtration Rules
# Release: 0.9
# Date: 14/10/2019
# Author: Antony Le Bechec
############################
MK_RELEASE="0.9"
MK_DATE="14/10/2019"

# Release notes:
# 0.9-14/10/2019: Creation


#####################
# VariantFiltration #
#####################

filterExpression=--filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "HARD_TO_VALIDATE" --filterExpression "DP < 30" --filterName "LowCoverage" --filterExpression "DP < 10" --filterName "VeryLowCoverage" --filterExpression "QUAL == 0 " --filterName "VeryVeryLowQual" --filterExpression "QUAL < 30.0 && QUAL > 0 " --filterName "VeryLowQual" --filterExpression "QUAL > 30.0 && QUAL < 50.0" --filterName "LowQual" --filterExpression "QD < 1.5" --filterName "LowQD"
genotypeFilterExpression=--genotypeFilterExpression "GQ < 50.0" --genotypeFilterName "VeryLowGQ"  --genotypeFilterExpression "GQ >= 50.0 && GQ < 90.0" --genotypeFilterName "LowGQ" --genotypeFilterExpression "DP < 10" --genotypeFilterName "VeryLowDP"  --genotypeFilterExpression "DP < 30" --genotypeFilterName "LowDP"


#--filterExpression "QUAL < 30.0 && QUAL > 0 " --filterName "VeryLowQual" --filterExpression "QUAL == 0 " --filterName "VeryVeryLowQual"

%.vcf: %.filtration.vcf %.filtration.vcf.idx %.empty.vcf %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $< -o $@ $(filterExpression) $(genotypeFilterExpression)
	if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	if [ ! -e $@ ]; then touch $@; fi;
	rm -f $<.idx $@.idx
	#-rm $*.unfiltered.vcf*
	#-rm -f $*.filtration.vcf* $@.idx


RELEASE_COMMENT := "\#\# VCF FILTRATION: GATK VariantFiltration to tag FILTER info in VCF file."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_CALLING:filtration:VariantFiltration of VCF."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "POST_ANNOTATION:filtration:VariantFiltration of VCF."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
