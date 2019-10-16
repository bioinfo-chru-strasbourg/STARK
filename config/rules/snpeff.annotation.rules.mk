############################
# HOWARD Annotation Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9b"
MK_DATE="21/04/2016"

## Release note
# 21/04/2016-V0.9b: Create SNPEFF


#HOWARD
SNPEFF?=$(NGSscripts)
SNPEFF_CONFIG?=$(SNPEFF)/snpEff.config


%.snpeff$(POST_ANNOTATION).vcf: %.vcf
	mkdir -p $@.stats
	$(JAVA) -jar $(SNPEFF) $(ASSEMBLY) $< -v -stats $@.stats/$(@F).stats.html > $@ ;


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# ANNOTATION '$(MK_RELEASE)': snpEff annotates, and generates *.snpeff.vcf file and *.snpeff.vcf.stats folder. Releases: SNPEFF_RELEASE='$(SNPEFF_RELEASE)'."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "ANNOTATOR:snpeff:snpEff annotation"
#MPILEUP_SAMTOOLS_OPTIONS='$(MPILEUP_SAMTOOLS_OPTIONS)',SAMTOOLS_FILTERS='$(SAMTOOLS_DP100_FILTERS)', SAMTOOLS_DP100_BCFTOOLS_FILTERS='$(SAMTOOLS_DP100_BCFTOOLS_FILTERS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
