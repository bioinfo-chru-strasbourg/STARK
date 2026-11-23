############################
# STAR-Fusion Calling Rules
# Author: Samuel Nicaise, Thomas Lavaux
############################
# Release
MK_RELEASE=1.0.0
MK_DATE=26/08/2026

# Release note
# 1.0.0-26/08/2026: Initial release of STARFusion calling rules.


##############
# STARFusion #
##############

# Detect fusions based on RNA-Seq data
# STARFusion need raw alignement, without splitNcigar, to work properly. So we need to use bam without splitNcigar directly from STAR alignments (%$(POST_ALIGNMENT).bam).

# Flags
# This variable can be defined directly on Apps
STARFUSION_FLAGS?=

# Filters Flags
# This variable can be defined directly on Apps
STARFUSION_FILTERS_FLAGS?= \
	--exclude-genotypes '' \
	--exclude '' \
	--include '' \
	--filter '' \
	--sort


# Rule
%.STARFusion$(POST_CALLING).vcf: %.star_raw.bam %.star_raw.bam.bai %.empty.vcf
	mkdir -p $*.fusion.reports;
	$(MAMBA) run -p $(STARFUSION_ENV) $(STARFUSION) $(STARFUSION_FLAGS) \
		--chimeric_junction $*.junction \
		--genome_lib_dir $$(dirname $(GENOME_RNA)) \
		--output_dir $*.fusion.reports;
	mv $*.fusion.reports/star-fusion.fusion_predictions.tsv $*.fusion.reports/$$(echo $(@F) | rev | cut -d"." -f4-  | rev).star-fusion.tsv
	mv $*.fusion.reports/star-fusion.fusion_predictions.abridged.tsv $*.fusion.reports/$$(echo $(@F) | rev | cut -d"." -f4-  | rev).star-fusion.abridged.tsv
	# VariantConvert
	# Need to create a specific config for variantconvert with the path to the genome fasta and the assembly name to be able to convert STARFusion output to vcf with correct header (configuration file is in config/variantconvert/GENOME/starfusion.json)
	cp $(VARIANTCONVERT_FOLDER_CONFIG)/$(ASSEMBLY)/starfusion.json $@.variantconvert.config.json
	$(VARIANTCONVERT) config -c $@.variantconvert.config.json --set GENOME.path=$(GENOME_RNA) --fill_genome_header
	# Convert to vcf
	$(VARIANTCONVERT) convert \
		-i $*.fusion.reports/$$(echo $(@F) | rev | cut -d"." -f4-  | rev).star-fusion.abridged.tsv \
		-o $@.tmp.vcf \
		-c $@.variantconvert.config.json;
	# Filters
	bash $(STARK_FOLDER_BIN)/filter_vcf.sh --input $@.tmp.vcf --output $@ --threads $(THREADS_BY_CALLER) --bcftools $(BCFTOOLS) $(STARFUSION_FILTERS_FLAGS);
	# Clean
	-rm $@.tmp* $@.variantconvert.config.json


# Documentation
$(file >>$(RELEASE_INFOS),## CALLER STARFusion '$(MK_RELEASE)': STARFusion caller to detect fusions based on RNA-Seq data, using CTAT, with filtrations. Parameters STARFUSION_FLAGS='$(STARFUSION_FLAGS)', STARFUSION_FILTERS_FLAGS='$(STARFUSION_FILTERS_FLAGS)')
$(file >>$(PIPELINES_INFOS),CALLER:STARFusion:STARFusion - by default with filtration:STARFUSION_FLAGS='$(STARFUSION_FLAGS)', STARFUSION_FILTERS_FLAGS='$(STARFUSION_FILTERS_FLAGS)')
