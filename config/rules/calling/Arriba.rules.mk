############################
# Arriba Calling Rules
# Release: 1.0.0
# Date: 08/11/2026
# Author: Samuel Nicaise, Thomas Lavaux
############################

# Release note
# 1.0.0-09/11/2026: Initial release of Arriba calling rules.


##########
# Arriba #
##########

# Detect fusions based on RNA-Seq data
# Arriba need raw alignement, without splitNcigar, to work properly. So we need to use bam without splitNcigar directly from STAR alignments (%$(POST_ALIGNMENT).bam).

# ARRIBA filters Flags
# This variable can be defined directly on Apps
ARRIBA_FILTERS_FLAGS?= \
	--exclude-genotypes '' \
	--exclude '' \
	--include '' \
	--filter '' \
	--sort


%.Arriba$(POST_CALLING).vcf: %.star_raw.bam %.star_raw.bam.bai %.empty.vcf
	mkdir -p $*.Arriba.reports;
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
	mv $*.Arriba.reports/arriba.fusions.discarded.tsv $*.Arriba.reports/$$(echo $(@F) | rev | cut -d"." -f4- | rev).arriba.fusions.discarded.tsv
	# VariantConvert
	# Need to create a specific config for variantconvert with the path to the genome fasta and the assembly name to be able to convert Arriba output to vcf with correct header (configuration file is in config/variantconvert/GENOME/arriba.json)
	cp $(VARIANTCONVERT_FOLDER_CONFIG)/$(ASSEMBLY)/arriba.json $@.variantconvert.config.json
	$(VARIANTCONVERT) config -c $@.variantconvert.config.json --set GENOME.path=$(GENOME_RNA) --fill_genome_header
	# Convert to vcf
	$(VARIANTCONVERT) convert \
		-i $*.Arriba.reports/$$(echo $(@F) | rev | cut -d"." -f4- | rev).arriba.fusions.tsv \
		-o $@.tmp.vcf \
		-c $@.variantconvert.config.json
	# Filters
	bash $(STARK_FOLDER_BIN)/filter_vcf.sh --input $@.tmp.vcf --output $@ --threads $(THREADS_BY_CALLER) --bcftools $(BCFTOOLS) $(ARRIBA_FILTERS_FLAGS);
	# Clean
	-rm $@.variantconvert.config.json $@.tmp.vcf


# ITDSEEK Release documentation

RELEASE_COMMENT := "\#\# CALLING Arriba: Tool to detect fusions based on RNA-Seq data (using blacklists packaged with Arriba and CTAT)"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING: Arriba fusion detection, generate *.Arriba.vcf files with with parameters. Parameters ARRIBA_FILTERS_FLAGS='$(ARRIBA_FILTERS_FLAGS)' "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:Arriba:Arriba - by default: using blacklists packaged with Arriba, ARRIBA_FILTERS_FLAGS='$(ARRIBA_FILTERS_FLAGS)' "
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
