############################
# Variant Normalization Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE=1.1.0
MK_DATE=12/10/2026

# Release notes:
# 1.0.0-07/09/2026: Creation, variant normalization
# 1.1.0-12/10/2026: Add update of VCF sequence dictionary


# Variant Normalization
######################


# Normalization of variants after calling
%.vcf: %.variantnormalization.vcf
	# Add or replace dictionary of contig in the VCF file
	if $(BCFTOOLS) view -h "$<" | grep -q '^##contig='; then \
		ln -s $< $@.with_contigs.vcf; \
	else \
		$(JAVA) $(JAVA_FLAGS_GATK4) -jar $(GATK4) \
			UpdateVCFSequenceDictionary \
			--sequence-dictionary $(DICT) \
			--replace True \
			--variant $< \
			--output $@.with_contigs.vcf; \
	fi;
	# Normalize the VCF file with BCFTOOLS
	$(BCFTOOLS) norm -m- --multi-overlaps 0 -f $(GENOME) $@.with_contigs.vcf | $(BCFTOOLS) norm --rm-dup exact | $(BCFTOOLS) +setGT -- -t . -n 0 | $(BCFTOOLS) +fixploidy -- | $(BCFTOOLS) +fill-tags -- -t all > $@
	rm -f $@.with_contigs.*


# Documentation
$(file >>$(RELEASE_INFOS),## POST_CALLING variantnormalization '$(MK_RELEASE)': Variant Normalization of VCF using BCFTools.)
$(file >>$(PIPELINES_INFOS),POST_CALLING:variantnormalization:Variant Normalization of VCF using BCFTools.)
