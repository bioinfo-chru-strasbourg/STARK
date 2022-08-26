############################
# STAR-Fusion Calling Rules
# Release: 0.9.4.8
# Date: 26/08/2022
# Author: Samuel Nicaise, Thomas Lavaux
############################

%.STARFusion$(POST_CALLING).SNP.vcf: %.bam %.genome
	STAR-Fusion --chimeric_junction $$(dirname $<)/Chimeric.out.junction --genome_lib_dir $$(cat $*.genome | dirname)

# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING STARFusion '$(MK_RELEASE)': CTAT Tool to detect fusions based on RNA-Seq data"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:STARFusion:STARFusion - by default: no options"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
