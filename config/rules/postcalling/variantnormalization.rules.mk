############################
# Variant Normalization Rules
# Release: 1.1.0
# Date: 12/10/2026
# Author: Antony Le Bechec
############################

# Release notes:
# 1.0.0-07/09/2026: Creation, variant normalization
# 1.1.0-12/10/2026: Add update of VCF sequence dictionary
# OPTIONS



# Variant Normalization
######################


# Normalization of variants after calling
%.vcf: %.variantnormalization.vcf
	# Add or replace dictionary of contig in the VCF file
	if $(BCFTOOLS) view -h "$<" | grep -q '^##contig='; then \
		ln -s $< $@.with_contigs; \
	else \
		$(JAVA) $(JAVA_FLAGS_GATK4) -jar $(GATK4) \
			UpdateVCFSequenceDictionary \
			--sequence-dictionary $(DICT) \
			--replace True \
			--variant $< \
			--output $@.with_contigs; \
	fi;
	# Normalize the VCF file with BCFTOOLS
	$(BCFTOOLS) norm -m- --multi-overlaps 0 -f $(GENOME) $@.with_contigs | $(BCFTOOLS) norm --rm-dup exact | $(BCFTOOLS) +setGT -- -t . -n 0 | $(BCFTOOLS) +fixploidy -- | $(BCFTOOLS) +fill-tags -- -t all > $@
	rm -f $@.with_contigs



RELEASE_COMMENT := "\#\# Variant Normalization: BCFTools variant normalization."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_CALLING:variantnormalization:Variant Normalization of VCF using BCFTools."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
