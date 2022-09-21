############################
# Arriba Calling Rules
# Release: 0.9.4.8
# Date: 26/08/2022
# Author: Samuel Nicaise, Thomas Lavaux
############################

%.Arriba$(POST_CALLING).vcf: %.star$(POST_ALIGNMENT).bam %.empty.vcf %.genome
	mkdir -p $*.arriba.reports;
	arriba -x $< -g $$(cat $*.genome) -k /STARK/tools/arriba/current/database/known_fusions_hg19_hs37d5_GRCh37_v2.2.1.tsv.gz -b /STARK/tools/arriba/current/database/blacklist_hg19_hs37d5_GRCh37_v2.2.1.tsv.gz -p /STARK/tools/arriba/current/database/protein_domains_hg19_hs37d5_GRCh37_v2.2.1.gff3 -o $*.arriba.reports/arriba.fusions.tsv -O $*.arriba.reports/arriba.fusions.discarded.tsv
	variantconvert convert -i $*.arriba.reports/arriba.fusions.tsv -o $@ -fi breakpoints -fo vcf -c /STARK/data/variantconvert/configs/config_starfusion_bas.json

# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING Arriba '$(MK_RELEASE)': Tool to detect fusions based on RNA-Seq data"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:Arriba:Arriba - by default: using blacklists packaged with Arriba"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
