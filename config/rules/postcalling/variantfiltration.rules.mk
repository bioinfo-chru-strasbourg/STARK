############################
# GATK4 Rules
# Release: 0.9.1
# Date: 02/02/2023
# Author: Antony Le Bechec
############################

# Release notes:
# 0.9.0-29/07/2022: Creation, variant filtration and variant recalibration
# 0.9.1-02/02/2023: Extract Variant Filtration

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
%.POST_CALLING_SNP.vcf: %.variantfiltration.vcf %.genome
	$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
		SelectVariants \
		-R $$(cat $*.genome) \
		-V $< \
		--select-type-to-include SNP \
		--select-type-to-include MIXED \
		--select-type-to-include MNP \
		--select-type-to-include SYMBOLIC \
		--select-type-to-include NO_VARIATION \
		-O $@.tmp.SNP.vcf;
	$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
		VariantFiltration \
		-R $$(cat $*.genome) \
		-V $@.tmp.SNP.vcf \
		-O $@.tmp.SNP.invalidate.vcf \
		$(VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS_OPTION) \
		--verbosity ERROR; \
	if [ ! -z '$(VARIANTFILTRATION_SNP_FILTER_OPTION)' ] && [ ! -z '$(VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION)' ]; then \
		$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			VariantFiltration \
			-R $$(cat $*.genome) \
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
%.POST_CALLING_InDel.vcf: %.variantfiltration.vcf %.genome
	$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
		SelectVariants \
		-R $$(cat $*.genome) \
		-V $< \
		--select-type-to-include INDEL \
		-O $@.tmp.INDEL.vcf;
	$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
		VariantFiltration \
		-R $$(cat $*.genome) \
		-V $@.tmp.INDEL.vcf \
		-O $@.tmp.INDEL.invalidate.vcf \
		$(VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS_OPTION) \
		--verbosity ERROR; \
	if [ ! -z '$(VARIANTFILTRATION_INDEL_FILTER_OPTION)' ] && [ ! -z '$(VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION)' ]; then \
		$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			VariantFiltration \
			-R $$(cat $*.genome) \
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



RELEASE_COMMENT := "\#\# Variant Filtration: GATK4 VariantFiltration."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_CALLING:variantfiltration:VariantFiltration of VCF using GATK databases."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
