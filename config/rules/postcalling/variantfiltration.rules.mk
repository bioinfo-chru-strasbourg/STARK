############################
# GATK4 Rules
# Release: 0.9.2
# Date: 14/13/2026
# Author: Antony Le Bechec
############################

# Release notes:
# 0.9.0-29/07/2022: Creation, variant filtration and variant recalibration
# 0.9.1-02/02/2023: Extract Variant Filtration
# 0.9.2-14/03/2026: Use BCFTools to select variants SNP, INDELS and OTHERS

# OPTIONS

# JAVA flags
JAVA_FLAGS_GATK4_CALLING_STEP?=$(JAVA_FLAGS) -XX:+UseParallelGC -XX:ParallelGCThreads=$(THREADS_BY_CALLER) -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_tribble=false

# Variant Filtration
VARIANTFILTRATION_OPTIONS?=
VARIANTFILTRATION_SNP_FILTER_OPTION?=
VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION?=
VARIANTFILTRATION_INDEL_FILTER_OPTION?=
VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION?=
VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS?=0
VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS_OPTION?=$(shell if (( $(VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS) )); then echo " --invalidate-previous-filters "; fi )



# Variant Filtration
######################

# SNP
%.POST_CALLING_VARIANTFILTRATION_SNP.vcf: %.variantfiltration.vcf
	$(BCFTOOLS) view -v snps,mnps --threads=$(THREADS_BY_CALLER) $< | $(BCFTOOLS) sort -o $@.tmp.SNP.vcf;
	$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
		VariantFiltration \
		-R $(GENOME) \
		-V $@.tmp.SNP.vcf \
		-O $@.tmp.SNP.invalidate.vcf \
		$(VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS_OPTION) \
		--verbosity ERROR; \
	if [ ! -z '$(VARIANTFILTRATION_SNP_FILTER_OPTION)' ] && [ ! -z '$(VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION)' ]; then \
		$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			VariantFiltration \
			-R $(GENOME) \
			-V $@.tmp.SNP.invalidate.vcf \
			-O $@ \
			--create-output-variant-index false \
			$(VARIANTFILTRATION_SNP_FILTER_OPTION) \
			$(VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION) \
			--verbosity ERROR; \
	else \
		cp $@.tmp.SNP.invalidate.vcf $@; \
	fi;
	rm -rf $@.tmp*


# INDEL
%.POST_CALLING_VARIANTFILTRATION_InDel.vcf: %.variantfiltration.vcf
	$(BCFTOOLS) view -v indels --threads=$(THREADS_BY_CALLER) $< | $(BCFTOOLS) sort -o $@.tmp.INDEL.vcf;
	$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
		VariantFiltration \
		-R $(GENOME) \
		-V $@.tmp.INDEL.vcf \
		-O $@.tmp.INDEL.invalidate.vcf \
		$(VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS_OPTION) \
		--verbosity ERROR; \
	if [ ! -z '$(VARIANTFILTRATION_INDEL_FILTER_OPTION)' ] && [ ! -z '$(VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION)' ]; then \
		$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			VariantFiltration \
			-R $(GENOME) \
			-V $@.tmp.INDEL.invalidate.vcf \
			-O $@ \
			--create-output-variant-index false \
			$(VARIANTFILTRATION_INDEL_FILTER_OPTION) \
			$(VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION) \
			--verbosity ERROR; \
	else \
		cp $@.tmp.INDEL.invalidate.vcf $@; \
	fi;
	rm -rf $@.tmp*

# Other varaints
%.POST_CALLING_VARIANTFILTRATION_OTHER_variants.vcf: %.variantfiltration.vcf
	$(BCFTOOLS) view -V snps,mnps,indels --threads=$(THREADS_BY_CALLER) $< | $(BCFTOOLS) sort -o $@;


# MERGE SNP and InDel VCF for Post calling steps. Because of loop in rules
%.vcf: %.POST_CALLING_VARIANTFILTRATION_SNP.vcf %.POST_CALLING_VARIANTFILTRATION_InDel.vcf %.POST_CALLING_VARIANTFILTRATION_OTHER_variants.vcf
	-$(JAVA) $(JAVA_FLAGS_GATK4) -jar $(GATK4) \
		MergeVcfs \
		-I $*.POST_CALLING_VARIANTFILTRATION_SNP.vcf \
		-I $*.POST_CALLING_VARIANTFILTRATION_InDel.vcf \
		-I $*.POST_CALLING_VARIANTFILTRATION_OTHER_variants.vcf \
		--CREATE_INDEX false \
		--SEQUENCE_DICTIONARY $(DICT) \
		-O $@;
	# Clear
	-rm -f $*.POST_CALLING_VARIANTFILTRATION_SNP.vcf* $*.POST_CALLING_VARIANTFILTRATION_InDel.vcf* $*.POST_CALLING_VARIANTFILTRATION_OTHER_variants.vcf*



RELEASE_COMMENT := "\#\# Variant Filtration: GATK4 VariantFiltration."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_CALLING:variantfiltration:VariantFiltration of VCF using GATK databases."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
