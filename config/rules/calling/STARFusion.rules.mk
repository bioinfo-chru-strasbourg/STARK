############################
# STAR-Fusion Calling Rules
# Release: 0.9.4.8
# Date: 26/08/2022
# Author: Samuel Nicaise, Thomas Lavaux
############################


# STARFusion need raw alignement, without splitNcigar, to work properly. So we need to use bam without splitNcigar directly from STAR alignments (%$(POST_ALIGNMENT).bam).

#%.STARFusion$(POST_CALLING).vcf: %$(POST_ALIGNMENT).bam %$(POST_ALIGNMENT).bam.bai %.empty.vcf %.junction
%.STARFusion$(POST_CALLING).vcf: %.star_raw.bam %.star_raw.bam.bai %.empty.vcf #%.junction
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
	# Clean
	-rm $@.variantconvert.config.json

# -fi breakpoints \
# -fo vcf \

# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING STARFusion '$(MK_RELEASE)': CTAT Tool to detect fusions based on RNA-Seq data"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:STARFusion:STARFusion - by default: no options"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
