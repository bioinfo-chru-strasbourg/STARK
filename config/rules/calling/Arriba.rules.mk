############################
# Arriba Calling Rules
# Release: 0.9.4.8
# Date: 26/08/2022
# Author: Samuel Nicaise, Thomas Lavaux
############################

%.Arriba$(POST_CALLING).vcf: %.bam %.empty.vcf %.genome
	mkdir -p $*.arriba.reports;
	arriba \
		-x $< \
		-a $$(cat $*.genome) \
		-g $$(cat $*.genome | xargs -0 dirname)/ref_annot.gtf \
		-k $$(ls $(ARRIBA_DATABASES)/known_fusions_$(ASSEMBLY)_*.tsv.gz) \
		-b $$(ls $(ARRIBA_DATABASES)/blacklist_$(ASSEMBLY)_*.tsv.gz) \
		-p $$(ls $(ARRIBA_DATABASES)/protein_domains_$(ASSEMBLY)_*.gff3) \
		-o $*.arriba.reports/arriba.fusions.tsv \
		-O $*.arriba.reports/arriba.fusions.discarded.tsv;
	# rename ...reports/arriba.fusions.tsv to ...reports/<sample>.arriba.fusions.tsv
	mv $*.arriba.reports/arriba.fusions.tsv $*.arriba.reports/$$(echo $(@F) | rev | cut -d"." -f4-  | rev).arriba.fusions.tsv
	mv $*.arriba.reports/arriba.fusions.discarded.tsv $*.arriba.reports/$$(echo $(@F) | rev | cut -d"." -f4-  | rev).arriba.fusions.discarded.tsv
	# convert to vcf
	variantconvert convert \
		-i $*.arriba.reports/$$(echo $(@F) | rev | cut -d"." -f4-  | rev).arriba.fusions.tsv \
		-o $@ \
		-fi breakpoints \
		-fo vcf \
		-c $(ASSEMBLY)/arriba.json

# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING Arriba '$(MK_RELEASE)': Tool to detect fusions based on RNA-Seq data"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:Arriba:Arriba - by default: using blacklists packaged with Arriba"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )