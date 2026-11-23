############################
# Arriba Calling Rules
# Author: Samuel Nicaise, Thomas Lavaux
############################
# Release
MK_RELEASE=1.0.0
MK_DATE=08/11/2026

# Release note
# 1.0.0-08/11/2026: Initial release of Arriba calling rules.


##########
# Arriba #
##########

# Detect fusions based on RNA-Seq data
# Arriba need raw alignement, without splitNcigar, to work properly. So we need to use bam without splitNcigar directly from STAR alignments (%$(POST_ALIGNMENT).bam).

# Flags
# This variable can be defined directly on Apps
ARRIBA_FLAGS?=

# Filters Flags
# This variable can be defined directly on Apps
ARRIBA_FILTERS_FLAGS?= \
	--exclude-genotypes '' \
	--exclude '' \
	--include '' \
	--filter '' \
	--sort


# Rule
%.Arriba$(POST_CALLING).vcf: %.star_raw.bam %.star_raw.bam.bai %.empty.vcf
	mkdir -p $*.Arriba.reports;
	$(ARRIBA) \
		-x $< \
		-a $(GENOME) \
		-g $$(dirname $(GENOME_RNA))/ref_annot.gtf \
		-k $$(ls $(ARRIBA_DATABASES)/$(ASSEMBLY)/known_fusions_$(ASSEMBLY)_*.tsv.gz) \
		-b $$(ls $(ARRIBA_DATABASES)/$(ASSEMBLY)/blacklist_$(ASSEMBLY)_*.tsv.gz) \
		-p $$(ls $(ARRIBA_DATABASES)/$(ASSEMBLY)/protein_domains_$(ASSEMBLY)_*.gff3) \
		$(ARRIBA_FLAGS) \
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


# Documentation
$(file >>$(RELEASE_INFOS),## CALLER Arriba '$(MK_RELEASE)': Arriba caller tool identify fusions from RNA-Seq data, generates *.Arriba.vcf files. Parameters ARRIBA_FLAGS='$(ARRIBA_FLAGS)', ARRIBA_FILTERS_FLAGS='$(ARRIBA_FILTERS_FLAGS)')
$(file >>$(PIPELINES_INFOS),CALLER:Arriba:Arriba - by default:ARRIBA_FLAGS='$(ARRIBA_FLAGS)', ARRIBA_FILTERS_FLAGS='$(ARRIBA_FILTERS_FLAGS)')
