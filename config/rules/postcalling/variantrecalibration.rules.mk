############################
# GATK4 Rules
# Release: 0.9.1
# Date: 02/02/2023
# Author: Antony Le Bechec
############################

# Release notes:
# 0.9.0-29/07/2022: Creation, variant filtration and variant recalibration
# 0.9.1-02/02/2023: Extract Variant Recalibration

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
#%.POST_CALLING_SNP.vcf: %.variantrecalibration.vcf
%.POST_CALLING_SNP.vcf: %.variantrecalibration.vcf
	
	# Original VCF
	$(BCFTOOLS) view -V indels $< -Oz1 -o $@.tmp.SNP.vcf.gz --threads $(THREADS_BY_SAMPLE);
	$(TABIX) $@.tmp.SNP.vcf.gz;
	
	# Minimal VCF
	$(BCFTOOLS) view -G $< | $(BCFTOOLS) annotate -x INFO > $@.tmp.minimal.vcf;
	
	# Select SNP for recalibration
	# Select SNP, MNP, MIXED, SYMBOLIC and NO_VARIATION variants for recalibration
	# Note: --select-type-to-include SYMBOLIC is used to include symbolic variants
	# Note: --lenient is used to allow the use of symbolic variants
	$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
		SelectVariants \
		-R $(GENOME) \
		-V $@.tmp.minimal.vcf \
		--select-type-to-include SNP \
		--select-type-to-include MIXED \
		--select-type-to-include MNP \
		--select-type-to-include SYMBOLIC \
		--select-type-to-include NO_VARIATION \
		-O $@.tmp.SNP.for_calibration.vcf \
		--lenient;
	
	# Recalibrate SNP variants
	# Note: If the input callset does not contain any training variants, the VariantRecalibrator will not run and will output a warning message.
	# Note: The --tranches-file option is used to output the tranches file, which contains the thresholds for filtering variants based on their recalibration scores.
	# Note: The --mode SNP option is used to specify that the recalibration is for SNP variants.
	# Note: The --output option is used to specify the output VCF file for the recalibrated variants.
	# Note: The --recal-file option is used to specify the output VCF file for the recalibrated variants.
	# Note: The --truth-sensitivity-filter-level option is used to specify the sensitivity level for filtering variants based on their recalibration scores.
	# Note: The --create-output-variant-index false option is used to disable the creation of an index for the output VCF file.
	# Note: The --verbosity ERROR option is used to set the verbosity level to ERROR, which means that only error messages
	-if (($(VARIANTRECALIBRATION_CHECK))); then \
		if ! $(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			VariantRecalibrator \
			$(VARIANTRECALIBRATOR_OPTIONS) \
			$(VARIANTRECALIBRATOR_SNP_OPTIONS) \
			-R $(GENOME) \
			-V $@.tmp.SNP.for_calibration.vcf \
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
			-V $@.tmp.minimal.vcf \
			-O $@.tmp.output.vcf.gz \
			--recal-file $@.tmp.SNP.recal.vcf \
			--tranches-file $@.tmp.SNP.tranches \
			--truth-sensitivity-filter-level 99.0 \
			--mode SNP; \
	else \
		echo "[WARNING] No ApplyVQRS on SNP due to resources error or lack of variant in the input callset"; \
		$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			VariantFiltration \
			-R $(GENOME) \
			-V $@.tmp.minimal.vcf \
			-O $@.tmp.SNP.invalidate.vcf \
			$(VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS_OPTION) \
			--verbosity ERROR; \
		if [ ! -z '$(VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_OPTION)' ] && [ ! -z '$(VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION)' ]; then \
			$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
				VariantFiltration \
				-R $(GENOME) \
				-V $@.tmp.SNP.invalidate.vcf \
				-O $@.tmp.output.vcf.gz \
				--create-output-variant-index false \
				$(VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_OPTION) \
				$(VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION) \
				--verbosity ERROR; \
		else \
			cp $@.tmp.SNP.invalidate.vcf $@.tmp.output.vcf.gz; \
		fi; \
	fi;
	
	# Index the output VCF file using tabix
	$(TABIX) $@.tmp.output.vcf.gz;

	# Reannotate the VCF with the original INFO fields
	# Note: The -a option is used to specify the annotation file, which contains the original INFO fields.
	# Note: The -c INFO option is used to specify that the INFO fields should be annotated.
	# Note: The -o option is used to specify the output VCF file.
	# Note: The $@.tmp.output.vcf.gz file is the output VCF file after the ApplyVQSR step, which contains the recalibrated variants.
	# Note: The $@ file is the final output VCF file after the annotation step.
	# $(BCFTOOLS) annotate \
  	# 	-a $< \
  	# 	-c INFO \
  	# 	-o $@ \
  	# 	$@.tmp.output.vcf.gz;
	# rm -f $@;
	# mv $@.tmp.reannotated.vcf $@;
	# Check if VCF annotation file contains variants using bcftools
	if (($(BCFTOOLS) view -H $@.tmp.output.vcf.gz | head -n1 | wc -l)); then \
		echo "[INFO] VCF file contains variants, proceeding with annotation."; \
		$(BCFTOOLS) annotate \
			-a $@.tmp.output.vcf.gz \
			-c CHROM,POS,REF,ALT,QUAL,FILTER \
			-o $@ \
			--threads $(THREADS_BY_SAMPLE) \
			$@.tmp.SNP.vcf.gz; \
	else \
		echo "[WARNING] VCF file does not contain variants, skipping annotation."; \
		cp $< $@; \
	fi;

	# Clear
	rm -rf $@.tmp*


# INDEL
%.POST_CALLING_InDel.vcf: %.variantrecalibration.vcf
	$(BCFTOOLS) view -v indels $< > $@.tmp.InDel.vcf
	$(BCFTOOLS) view -G $< | $(BCFTOOLS) annotate -x INFO > $@.tmp.minimal.vcf;
	$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
		SelectVariants \
		-R $(GENOME) \
		-V $@.tmp.minimal.vcf \
		--select-type-to-include INDEL \
		-O $@.tmp.InDel.for_calibration.vcf \
		--lenient;
	-if (($(VARIANTRECALIBRATION_CHECK))); then \
		if ! $(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			VariantRecalibrator \
			$(VARIANTRECALIBRATOR_OPTIONS) \
			$(VARIANTRECALIBRATOR_INDEL_OPTIONS) \
			-R $(GENOME) \
			-V $@.tmp.InDel.for_calibration.vcf \
			--mode INDEL \
			--output $@.tmp.InDel.recal.vcf \
			--tranches-file $@.tmp.InDel.tranches 2>/dev/null; \
		then \
			echo "[WARNING] No Recal on INDEL due to lack on training variant in the input callset"; \
		fi; \
	fi;
	if (($$(grep '^#' -vc $@.tmp.InDel.recal.vcf))) && (($(VARIANTRECALIBRATION_CHECK))); then \
		$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			ApplyVQSR \
			-R $(GENOME) \
			-V $@.tmp.InDel.vcf \
			-O $@ \
			--recal-file $@.tmp.InDel.recal.vcf \
			--tranches-file $@.tmp.InDel.tranches \
			--truth-sensitivity-filter-level 99.0 \
			--mode INDEL; \
	else \
		echo "[WARNING] No ApplyVQRS on INDEL due to resources error or lack of variant in the input callset"; \
		$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			VariantFiltration \
			-R $(GENOME) \
			-V $@.tmp.InDel.vcf \
			-O $@.tmp.InDel.invalidate.vcf \
			$(VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS_OPTION) \
			--verbosity ERROR; \
		if [ ! -z '$(VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_OPTION)' ] && [ ! -z '$(VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION)' ]; then \
			$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
				VariantFiltration \
				-R $(GENOME) \
				-V $@.tmp.InDel.invalidate.vcf \
				-O $@ \
				--create-output-variant-index false \
				$(VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_OPTION) \
				$(VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION) \
				--verbosity ERROR; \
		else \
			cp $@.tmp.InDel.invalidate.vcf $@; \
		fi; \
	fi;
	# Clear
	rm -rf $@.tmp*



RELEASE_COMMENT := "\#\# Variant Recalibrator: GATK4 VariantRecalibrator."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_CALLING:variantrecalibration:VariantRecalibrator of VCF using GATK databases."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
