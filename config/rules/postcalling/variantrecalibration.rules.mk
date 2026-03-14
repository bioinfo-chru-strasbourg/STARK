############################
# GATK4 Rules
# Release: 0.9.2
# Date: 14/03/2026
# Author: Antony Le Bechec
############################

# Release notes:
# 0.9.0-29/07/2022: Creation, variant filtration and variant recalibration
# 0.9.1-02/02/2023: Extract Variant Recalibration
# 0.9.2-14/03/2026: Use BCFTools to select variants SNP, INDELS and OTHERS

# OPTIONS

# JAVA flags
JAVA_FLAGS_GATK4_CALLING_STEP?=$(JAVA_FLAGS) -XX:+UseParallelGC -XX:ParallelGCThreads=$(THREADS_BY_CALLER) -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_tribble=false

# Variant Recalibrator
VARIANTRECALIBRATOR_OPTIONS?=
VARIANTRECALIBRATION_SNP_ANNOTATIONS?=
VARIANTRECALIBRATION_SNP_TRANCHES?=
VARIANTRECALIBRATION_INDEL_ANNOTATIONS?=
VARIANTRECALIBRATION_INDEL_TRANCHES?=
VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_OPTION?=
VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION?=
VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_OPTION?=
VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION?=
VARIANTRECALIBRATOR_SNP_OPTIONS?=$(VARIANTRECALIBRATION_SNP_RESOURCES_OPTION) $(VARIANTRECALIBRATION_SNP_ANNOTATIONS) $(VARIANTRECALIBRATION_SNP_TRANCHES)
VARIANTRECALIBRATOR_INDEL_OPTIONS?=$(VARIANTRECALIBRATION_INDEL_RESOURCES_OPTION) $(VARIANTRECALIBRATION_INDEL_ANNOTATIONS) $(VARIANTRECALIBRATION_INDEL_TRANCHES)



# Variant Recalibrator
########################

# SNP
%.POST_CALLING_VARIANTRECALIBRATION_SNP.vcf: %.variantrecalibration.vcf
# 	$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
# 		SelectVariants \
# 		-R $(GENOME) \
# 		-V $< \
# 		--select-type-to-include SNP \
# 		--select-type-to-include MIXED \
# 		--select-type-to-include MNP \
# 		--select-type-to-include SYMBOLIC \
# 		--select-type-to-include NO_VARIATION \
# 		-O $@.tmp.SNP.vcf;
	$(BCFTOOLS) view -v snps,mnps $< -o $@.tmp.SNP.vcf;
	-if (($(VARIANTRECALIBRATION_CHECK))); then \
		if ! $(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			VariantRecalibrator \
			$(VARIANTRECALIBRATOR_OPTIONS) \
			$(VARIANTRECALIBRATOR_SNP_OPTIONS) \
			-R $(GENOME) \
			-V $@.tmp.SNP.vcf \
			--mode SNP \
			--output $@.tmp.SNP.recal.vcf \
			--tranches-file $@.tmp.SNP.tranches 2>/dev/null; \
		then \
			echo "[WARNING] No Recal on SNP due to lack on training variant in the input callset"; \
		fi; \
	fi;
	if (($$(grep '^#' -vc $@.tmp.SNP.recal.vcf))) && (($(VARIANTRECALIBRATION_CHECK))); then \
		$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			ApplyVQSR \
			-R $(GENOME) \
			-V $@.tmp.SNP.vcf \
			-O $@ \
			--recal-file $@.tmp.SNP.recal.vcf \
			--tranches-file $@.tmp.SNP.tranches \
			--truth-sensitivity-filter-level 99.0 \
			--mode SNP; \
	else \
		echo "[WARNING] No ApplyVQRS on SNP due to resources error or lack of variant in the input callset"; \
		$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			VariantFiltration \
			-R $(GENOME) \
			-V $@.tmp.SNP.vcf \
			-O $@.tmp.SNP.invalidate.vcf \
			$(VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS_OPTION) \
			--verbosity ERROR; \
		if [ ! -z '$(VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_OPTION)' ] && [ ! -z '$(VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION)' ]; then \
			$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
				VariantFiltration \
				-R $(GENOME) \
				-V $@.tmp.SNP.invalidate.vcf \
				-O $@ \
				--create-output-variant-index false \
				$(VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_OPTION) \
				$(VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION) \
				--verbosity ERROR; \
		else \
			cp $@.tmp.SNP.invalidate.vcf $@; \
		fi; \
	fi;


# INDEL
%.POST_CALLING_VARIANTRECALIBRATION_InDel.vcf: %.variantrecalibration.vcf
# 	$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
# 		SelectVariants \
# 		-R $(GENOME) \
# 		-V $< \
# 		--select-type-to-include INDEL \
# 		-O $@.tmp.InDel.vcf;
	$(BCFTOOLS) view -v indels $< -o $@.tmp.INDEL.vcf;
	-if (($(VARIANTRECALIBRATION_CHECK))); then \
		if ! $(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			VariantRecalibrator \
			$(VARIANTRECALIBRATOR_OPTIONS) \
			$(VARIANTRECALIBRATOR_INDEL_OPTIONS) \
			-R $(GENOME) \
			-V $@.tmp.INDEL.vcf \
			--mode INDEL \
			--output $@.tmp.INDEL.recal.vcf \
			--tranches-file $@.tmp.INDEL.tranches 2>/dev/null; \
		then \
			echo "[WARNING] No Recal on INDEL due to lack on training variant in the input callset"; \
		fi; \
	fi;
	if (($$(grep '^#' -vc $@.tmp.INDEL.recal.vcf))) && (($(VARIANTRECALIBRATION_CHECK))); then \
		$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			ApplyVQSR \
			-R $(GENOME) \
			-V $@.tmp.INDEL.vcf \
			-O $@ \
			--recal-file $@.tmp.INDEL.recal.vcf \
			--tranches-file $@.tmp.INDEL.tranches \
			--truth-sensitivity-filter-level 99.0 \
			--mode INDEL; \
	else \
		echo "[WARNING] No ApplyVQRS on INDEL due to resources error or lack of variant in the input callset"; \
		$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			VariantFiltration \
			-R $(GENOME) \
			-V $@.tmp.INDEL.vcf \
			-O $@.tmp.INDEL.invalidate.vcf \
			$(VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS_OPTION) \
			--verbosity ERROR; \
		if [ ! -z '$(VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_OPTION)' ] && [ ! -z '$(VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION)' ]; then \
			$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
				VariantFiltration \
				-R $(GENOME) \
				-V $@.tmp.INDEL.invalidate.vcf \
				-O $@ \
				--create-output-variant-index false \
				$(VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_OPTION) \
				$(VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION) \
				--verbosity ERROR; \
		else \
			cp $@.tmp.INDEL.invalidate.vcf $@; \
		fi; \
	fi;


# Other variants
%.POST_CALLING_VARIANTRECALIBRATION_OTHER_variants.vcf: %.variantrecalibration.vcf
# 	$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
# 		SelectVariants \
# 		-R $(GENOME) \
# 		-V $< \
# 		--select-type-to-exclude SNP \
# 		--select-type-to-exclude INDEL \
# 		--select-type-to-exclude MIXED \
# 		--select-type-to-exclude MNP \
# 		--select-type-to-exclude SYMBOLIC \
# 		--select-type-to-exclude NO_VARIATION \
# 		-O $@;
	$(BCFTOOLS) view -V snps,mnps,indels $< -o $@;


# MERGE SNP and InDel VCF for Post calling steps. Because of loop in rules
%.vcf: %.POST_CALLING_VARIANTRECALIBRATION_SNP.vcf %.POST_CALLING_VARIANTRECALIBRATION_InDel.vcf %.POST_CALLING_VARIANTRECALIBRATION_OTHER_variants.vcf
	$(JAVA) $(JAVA_FLAGS_GATK4) -jar $(GATK4) \
		MergeVcfs \
		-I $*.POST_CALLING_VARIANTRECALIBRATION_SNP.vcf \
		-I $*.POST_CALLING_VARIANTRECALIBRATION_InDel.vcf \
		-I $*.POST_CALLING_VARIANTRECALIBRATION_OTHER_variants.vcf \
		--CREATE_INDEX false \
		--SEQUENCE_DICTIONARY $(DICT) \
		-O $@;
	# Clear
	-rm -f $*.POST_CALLING_VARIANTRECALIBRATION_SNP.vcf* $*.POST_CALLING_VARIANTRECALIBRATION_InDel.vcf* $*.POST_CALLING_VARIANTRECALIBRATION_OTHER_variants.vcf*




RELEASE_COMMENT := "\#\# Variant Recalibrator: GATK4 VariantRecalibrator."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_CALLING:variantrecalibration:VariantRecalibrator of VCF using GATK databases."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
