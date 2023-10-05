############################
# STAR-Fusion Calling Rules
# Release: 0.9.4.8
# Date: 26/08/2022
# Author: Samuel Nicaise, Thomas Lavaux
############################

# CTAT_DATABASES=$DBFOLDER/CTAT_LIB/current

%.STARFusion$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.junction
	mkdir -p $*.fusion.reports;
	$STARFUSION \
		--chimeric_junction $*.junction \
		--genome_lib_dir $$CTAT_DATABASES/$ASSEMBLY/ \
		--output_dir $*.fusion.reports;
	# rename ...reports/star-fusion.fusion_predictions.tsv to ...reports/<sample>.star-fusion.tsv
	mv $*.fusion.reports/star-fusion.fusion_predictions.tsv $*.fusion.reports/$$(echo $(@F) | rev | cut -d"." -f4-  | rev).star-fusion.tsv
	mv $*.fusion.reports/star-fusion.fusion_predictions.abridged.tsv $*.fusion.reports/$$(echo $(@F) | rev | cut -d"." -f4-  | rev).star-fusion.abridged.tsv
	# convert to vcf
	variantconvert convert \
		-i $*.fusion.reports/$$(echo $(@F) | rev | cut -d"." -f4-  | rev).star-fusion.abridged.tsv \
		-o $@ \
		-fi breakpoints \
		-fo vcf \
		-c $(ASSEMBLY)/starfusion.json;

# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING STARFusion '$(MK_RELEASE)': CTAT Tool to detect fusions based on RNA-Seq data"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:STARFusion:STARFusion - by default: no options"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
