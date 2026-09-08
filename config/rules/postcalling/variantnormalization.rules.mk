############################
# Variant Normalization Rules
# Release: 1.0.0
# Date: 07/09/2026
# Author: Antony Le Bechec
############################

# Release notes:
# 1.0.0-07/09/2026: Creation, variant normalization

# OPTIONS



# Variant Normalization
######################


# Normalization of variants after calling
%.vcf: %.variantnormalization.vcf
	$(BCFTOOLS) norm -m- --multi-overlaps 0 -f $(GENOME) $< | $(BCFTOOLS) norm --rm-dup exact | $(BCFTOOLS) +setGT -- -t . -n 0 | $(BCFTOOLS) +fixploidy -- | $(BCFTOOLS) +fill-tags -- -t all > $@



RELEASE_COMMENT := "\#\# Variant Normalization: BCFTools variant normalization."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_CALLING:variantnormalization:Variant Normalization of VCF using BCFTools."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
