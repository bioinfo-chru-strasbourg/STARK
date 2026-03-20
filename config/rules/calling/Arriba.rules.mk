############################
# Arriba Calling Rules
# Release: 0.9.4.8
# Date: 26/08/2022
# Author: Samuel Nicaise, Thomas Lavaux
############################

# Arriba need raw alignement, without splitNcigar, to work properly. So we need to use bam without splitNcigar directly from STAR alignments (%$(POST_ALIGNMENT).bam).

%.Arriba$(POST_CALLING).vcf: %.star_raw.bam %.star_raw.bam.bai %.empty.vcf
	mkdir -p $*.arriba.reports;
	$(ARRIBA) \
		-x $< \
		-a $(GENOME) \
		-g $$(dirname $(GENOME_RNA))/ref_annot.gtf \
		-k $$(ls $(ARRIBA_DATABASES)/$(ASSEMBLY)/known_fusions_$(ASSEMBLY)_*.tsv.gz) \
		-b $$(ls $(ARRIBA_DATABASES)/$(ASSEMBLY)/blacklist_$(ASSEMBLY)_*.tsv.gz) \
		-p $$(ls $(ARRIBA_DATABASES)/$(ASSEMBLY)/protein_domains_$(ASSEMBLY)_*.gff3) \
		-o $*.Arriba.reports/arriba.fusions.tsv \
		-O $*.Arriba.reports/arriba.fusions.discarded.tsv;
	mv $*.Arriba.reports/arriba.fusions.tsv $*.Arriba.reports/$$(echo $(@F) | rev | cut -d"." -f4-  | rev).arriba.fusions.tsv
	mv $*.Arriba.reports/arriba.fusions.discarded.tsv $*.Arriba.reports/$$(echo $(@F) | rev | cut -d"." -f4-  | rev).arriba.fusions.discarded.tsv
	# VariantConvert
	# Need to create a specific config for variantconvert with the path to the genome fasta and the assembly name to be able to convert Arriba output to vcf with correct header (configuration file is in config/variantconvert/GENOME/arriba.json)
	cp $(VARIANTCONVERT_FOLDER_CONFIG)/$(ASSEMBLY)/arriba.json $@.variantconvert.config.json
	$(VARIANTCONVERT) config -c $@.variantconvert.config.json --set GENOME.path=$(GENOME_RNA) --fill_genome_header
	# Convert to vcf
	$(VARIANTCONVERT) convert \
		-i $*.Arriba.reports/$$(echo $(@F) | rev | cut -d"." -f4-  | rev).arriba.fusions.tsv \
		-o $@ \
		-c $@.variantconvert.config.json
	# Clean
	-rm $@.variantconvert.config.json

# -g $$(dirname $(GENOME_RNA))/ref_annot.gtf \
# -g $(REFSEQ_GENES_GTF) \

# VARIANTCONVERT_FOLDER_CONFIG

# -fi breakpoints \
# -fo vcf \

# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING Arriba '$(MK_RELEASE)': Tool to detect fusions based on RNA-Seq data"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:Arriba:Arriba - by default: using blacklists packaged with Arriba"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
