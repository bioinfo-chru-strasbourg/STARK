############################
# STAR-Fusion Calling Rules
# Release: 0.9.4.8
# Date: 26/08/2022
# Author: Samuel Nicaise, Thomas Lavaux
############################

%.STARFusion$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.junction
	mkdir -p $*.fusion.reports;
	STAR-Fusion --chimeric_junction $*.junction --genome_lib_dir $$(cat $*.genome | xargs -0 dirname) --output_dir $*.fusion.reports;
	variantconvert convert -i $*.fusion.reports/star-fusion.fusion_predictions.tsv -o $@ -fi breakpoints -fo vcf -c /STARK/data/variantconvert/configs/config_starfusion_bas.json;

# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING STARFusion '$(MK_RELEASE)': CTAT Tool to detect fusions based on RNA-Seq data"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:STARFusion:STARFusion - by default: no options"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )