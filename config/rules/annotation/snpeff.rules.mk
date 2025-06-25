############################
# SnpEff Annotation Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.1.0"
MK_DATE="25/06/2025"

## Release note
# 21/04/2016-V0.9b: Create SNPEFF
# 25/06/2025-V0.9.1.0: STARK release 19 compatibility



%.snpeff$(POST_ANNOTATION).vcf: %.vcf
	mkdir -p $@.stats
	$(JAVA) -jar $(SNPEFF) $(ASSEMBLY) $< -v -stats $@.stats/$(@F).stats.html > $@ ;


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# ANNOTATION '$(MK_RELEASE)': snpEff annotates, and generates *.snpeff.vcf file and *.snpeff.vcf.stats folder. Releases: SNPEFF_RELEASE='$(SNPEFF_RELEASE)'."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "ANNOTATOR:snpeff:snpEff annotation"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
