############################
# STAR-Fusion Calling Rules
# Release: 1.0.0
# Date: 26/08/2026
# Author: Samuel Nicaise, Thomas Lavaux
############################

# Release note
# 1.0.0-26/08/2026: Initial release of STARFusion calling rules.


##############
# STARFusion #
##############

# Detect fusions based on RNA-Seq data
# STARFusion need raw alignement, without splitNcigar, to work properly. So we need to use bam without splitNcigar directly from STAR alignments (%$(POST_ALIGNMENT).bam).

# STARFUSION filters Flags
# This variable can be defined directly on Apps
STARFUSION_FILTERS_FLAGS?= \
	--exclude-genotypes '' \
	--exclude '' \
	--include '' \
	--filter '' \
	--sort


%.STARFusion$(POST_CALLING).vcf: %.star_raw.bam %.star_raw.bam.bai %.empty.vcf
	mkdir -p $*.fusion.reports;
	$(MAMBA) run -p $(STARFUSION_ENV) $(STARFUSION) \
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
		-o $@ \
		-c $@.variantconvert.config.json;
	# Filters
	bash $(STARK_FOLDER_BIN)/filter_vcf.sh --input $@ --output $@ --threads $(THREADS_BY_CALLER) --bcftools $(BCFTOOLS) $(STARFUSION_FILTERS_FLAGS);
	# Clean
	-rm $@.variantconvert.config.json


# ITDSEEK Release documentation

RELEASE_COMMENT := "\#\# CALLING STARFusion: Tool to detect fusions based on RNA-Seq data"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING STARFusion: STARFusion Tool to detect fusions based on RNA-Seq data (using CTAT). Parameters STARFUSION_FILTERS_FLAGS='$(STARFUSION_FILTERS_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:STARFusion:STARFusion - by default: no options, STARFUSION_FILTERS_FLAGS='$(STARFUSION_FILTERS_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
