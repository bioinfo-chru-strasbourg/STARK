############################
# GATK4 Rules
# Release: 0.9.0
# Date: 28/07/2022
# Author: Antony Le Bechec
############################

# Release notes:
# 0.9.0-29/07/2022: Creation, variant filtration and variant recalibration

# OPTIONS


# JAVA flags

JAVA_FLAGS_GATK4_CALLING_STEP?=$(JAVA_FLAGS) -XX:+UseParallelGC -XX:ParallelGCThreads=$(THREADS_BY_CALLER) -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_tribble=false


JAVA_FLAGS_GATK4_ALIGNMENT_STEP?=$(JAVA_FLAGS) -XX:+UseParallelGC -XX:ParallelGCThreads=$(THREADS_BY_ALIGNER) -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_tribble=false

# Variant Filtration

VARIANTFILTRATION_OPTIONS?=

VARIANTFILTRATION_SNP_FILTER_OPTION?=
VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION?=
VARIANTFILTRATION_INDEL_FILTER_OPTION?=
VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION?=

VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS?=0
VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS_OPTION?=$(shell if (( $(VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS) )); then echo " --invalidate-previous-filters "; fi )


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


# Variant Filtration
######################

# SNP

%.POST_CALLING_SNP.vcf: %.variantfiltration.vcf %.genome
	$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
		SelectVariants \
		-R $$(cat $*.genome) \
		-V $< \
		--select-type-to-include SNP \
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



# Variant Recalibrator
########################

# SNP

%.POST_CALLING_SNP.vcf: %.variantrecalibration.vcf %.genome
	$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
		SelectVariants \
		-R $$(cat $*.genome) \
		-V $< \
		--select-type-to-include SNP \
		-O $@.tmp.SNP.vcf;
	-if (($(VARIANTRECALIBRATION_CHECK))); then \
		if ! $(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			VariantRecalibrator \
			$(VARIANTRECALIBRATOR_OPTIONS) \
			$(VARIANTRECALIBRATOR_SNP_OPTIONS) \
			-R $$(cat $*.genome) \
			-V $@.tmp.SNP.vcf \
			--mode SNP \
			--output $@.tmp.SNP.recal.vcf \
			--tranches-file $@.tmp.SNP.tranches 2>/dev/null; \
		then \
			echo "[WARNING] No Recal on SNP due to lack on training variant in the input callset"; \
		fi; \
	fi;
	if (($$(grep '^#' -vc $@.tmp.SNP.recal.vcf))) && (($(VARIANTRECALIBRATION_CHECK))); then \
		$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			ApplyVQSR \
			-R $$(cat $*.genome) \
			-V $@.tmp.SNP.vcf \
			-O $@ \
			--recal-file $@.tmp.SNP.recal.vcf \
			--tranches-file $@.tmp.SNP.tranches \
			--truth-sensitivity-filter-level 99.0 \
			--mode SNP; \
	else \
		echo "[WARNING] No ApplyVQRS on SNP due to resources error or lack of variant in the input callset"; \
		$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			VariantFiltration \
			-R $$(cat $*.genome) \
			-V $@.tmp.SNP.vcf \
			-O $@.tmp.SNP.invalidate.vcf \
			$(VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS_OPTION) \
			--verbosity ERROR; \
		if [ ! -z '$(VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_OPTION)' ] && [ ! -z '$(VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION)' ]; then \
			$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
				VariantFiltration \
				-R $$(cat $*.genome) \
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

%.POST_CALLING_InDel.vcf: %.variantrecalibration.vcf %.genome
	$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
		SelectVariants \
		-R $$(cat $*.genome) \
		-V $< \
		--select-type-to-include INDEL \
		-O $@.tmp.InDel.vcf;
	-if (($(VARIANTRECALIBRATION_CHECK))); then \
		if ! $(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			VariantRecalibrator \
			$(VARIANTRECALIBRATOR_OPTIONS) \
			$(VARIANTRECALIBRATOR_INDEL_OPTIONS) \
			-R $$(cat $*.genome) \
			-V $@.tmp.InDel.vcf \
			--mode INDEL \
			--output $@.tmp.InDel.recal.vcf \
			--tranches-file $@.tmp.InDel.tranches 2>/dev/null; \
		then \
			echo "[WARNING] No Recal on INDEL due to lack on training variant in the input callset"; \
		fi; \
	fi;
	if (($$(grep '^#' -vc $@.tmp.InDel.recal.vcf))) && (($(VARIANTRECALIBRATION_CHECK))); then \
		$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			ApplyVQSR \
			-R $$(cat $*.genome) \
			-V $@.tmp.InDel.vcf \
			-O $@ \
			--recal-file $@.tmp.InDel.recal.vcf \
			--tranches-file $@.tmp.InDel.tranches \
			--truth-sensitivity-filter-level 99.0 \
			--mode INDEL; \
	else \
		echo "[WARNING] No ApplyVQRS on INDEL due to resources error or lack of variant in the input callset"; \
		$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
			VariantFiltration \
			-R $$(cat $*.genome) \
			-V $@.tmp.InDel.vcf \
			-O $@.tmp.InDel.invalidate.vcf \
			$(VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS_OPTION) \
			--verbosity ERROR; \
		if [ ! -z '$(VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_OPTION)' ] && [ ! -z '$(VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION)' ]; then \
			$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
				VariantFiltration \
				-R $$(cat $*.genome) \
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



# BAM Recalibration
####################
# TODO

# %.bam: %.BQSRSpark.bam %.BQSRSpark.bam.bai %.BQSRSpark.from_manifest.interval_list %.genome %.design.bed #%.recalibration.bam.grp #%.recalibration.bam.bai %.from_manifest.intervals %.recalibration.from_manifest.intervals
# 	# Recalibrate BAM with BQSRPipelineSpark 
# 	$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) \
# 		BQSRPipelineSpark \
# 		-I $< \
# 		-O $@ \
# 		-R $$(cat $*.genome) \
# 		--known-sites $(VCFDBSNP) \
# 		--emit-original-quals true \
# 		--intervals $*.design.bed \
# 		--spark-master local[$(THREADS_BY_ALIGNER)] ;
# 	# clean
# 	#-rm -f $<;
# 	-rm -f $*.BQSRSpark.*;



RELEASE_COMMENT := "\#\# Variant Filtration: GATK4 VariantFiltration."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_CALLING:variantfiltration:VariantFiltration of VCF using GATK databases."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


RELEASE_COMMENT := "\#\# Variant Recalibrator: GATK4 VariantRecalibrator."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_CALLING:variantrecalibration:VariantRecalibrator of VCF using GATK databases."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

