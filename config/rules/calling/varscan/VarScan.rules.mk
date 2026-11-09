############################
# VarScan Calling Rules
# Release: 1.0.0
# Date: 04/11/2026
# Author: Antony Le Bechec
############################

# Release notes:
# 0.9.1-24/04/2015: Add 'varscanlowfreqlowcovsnp' calling
# 0.9.2-05/10/2015: Add 'VarScan' calling, default VarScan parameters
# 0.9.3-07/10/2015: Add filters on VarScan vcf
# 0.9.3.1-10/11/2015: Add pipeline VarScanHEMATO
# 0.9.4-10/11/2015: Add Merge SNP and InDel
# 0.9.4.1-03/05/2016: Bug correction, if empty VCF after calling. Add Release infos
# 0.9.4.2-04/05/2016: Bug correction
# 0.9.4.3-17/05/2016: Bug correction
# 0.9.4.4-21/03/2019: Change parameters for mpileup by adding FATBAM Soft Clip to Q0
# 0.9.4.5-27/09/2019: Change FATBAM to CAP tool
# 0.9.4.6-25/05/2021: Change VARSCAN_FILTERS
# 0.9.4.7-28/06/2021: Homogenize VARSCAN thresholds 
# 0.9.4.8-29/07/2022: Homogenize VARSCAN rules path 
# 0.9.5-03/02/2023: Extract VarScan rules
# 0.9.6-12/05/2023: Remove VariantFiltration of GATK
# 1.0.0-04/11/2026: Improve parallelize and filter, with common algorithm


###########
# VARSCAN #
###########

# VARSCAN SNP Flags
# This variable can be defined directly on Apps
VARSCAN_SNP_FLAGS?= \
	--min-var-freq 0.01 \
	--min-reads2 3 \
	--min-coverage 30 \
	--p-value 1e-3 \
	--min-avg-qual 30

# VARSCAN INDEL Flags
# This variable can be defined directly on Apps
VARSCAN_INDEL_FLAGS?= \
	--min-var-freq 0.01 \
	--min-reads2 3 \
	--min-coverage 30 \
	--p-value 1e-3 \
	--min-avg-qual 10

# VARSCAN filters Flags
# This variable can be defined directly on Apps
VARSCAN_FILTERS_FLAGS?= \
	--exclude-genotypes '' \
	--exclude '' \
	--include '' \
	--filter '' \
	--sort


# VARSCAN target rule
%.VarScan$(POST_CALLING).vcf: %.bam %.bam.bai %.bam.mpileup %.empty.vcf %.design.bed %.bam.chroms
	# List of vcf and sample files to be generated per chromosome
	cut -f1 $*.bam.chroms | sed 's|.*|$@.split.&.vcf.gz|' | tr '\n' ' ' > $@.split.vcfs.list;
	$(BCFTOOLS) view $*.empty.vcf | grep "^#CHROM" | cut -f10 > $@.split.sample.txt;
	# Generate bed per chromosome
	awk '{print > ("$@.split."$$1".bed")}' $*.design.bed
	# Split the mpileup by chromosome in a single pass (chroms.txt keeps the mpileup order)
	awk '$$1 != prev { if (prev != "") close(out); out = "$@.split."$$1".mpileup"; prev = $$1 } { print > out }' $*.bam.mpileup;
	# Generate the sub Makefile: one explicit rule per chromosome (empty "all" if no chromosome)
	{ \
	printf 'all: %s\n' "$$(cat $@.split.vcfs.list)"; \
	for c in $$(cut -f1 $*.bam.chroms); do \
		printf '%s:\n' "$@.split.$$c.SNP.vcf"; \
		printf '\tcat %s | $(JAVA) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_SNP_FLAGS) --vcf-sample-list %s > %s\n' "$@.split.$$c.mpileup" "$@.split.sample.txt" "$@.split.$$c.SNP.vcf"; \
		printf '%s:\n' "$@.split.$$c.INDEL.vcf"; \
		printf '\tcat %s | $(JAVA) -jar $(VARSCAN) mpileup2indel --output-vcf 1 $(VARSCAN_INDEL_FLAGS) --vcf-sample-list %s > %s\n' "$@.split.$$c.mpileup" "$@.split.sample.txt" "$@.split.$$c.INDEL.vcf"; \
		printf '%s: %s %s\n' "$@.split.$$c.vcf.gz" "$@.split.$$c.SNP.vcf" "$@.split.$$c.INDEL.vcf"; \
		printf '\t$(JAVA) $(JAVA_FLAGS) -jar $(GATK) MergeVcfs -I %s -I %s --CREATE_INDEX false --SEQUENCE_DICTIONARY $(DICT) -O %s\n' "$@.split.$$c.SNP.vcf" "$@.split.$$c.INDEL.vcf" "$@.split.$$c.vcf.gz" ; \
	done; \
	} > $@.split.sub.mk
	# Add rule to concat per-chromosome VCFs into the final VCF (empty VCF if no chromosome)
	{ \
		printf '%s: %s\n' "$@.tmp.unfiltered.vcf" "$$(cat $@.split.vcfs.list) $*.empty.vcf"; \
		printf '\t$(BCFTOOLS) concat %s -Ov --threads $(THREADS_BY_CALLER) -o %s\n' "$$(cat $@.split.vcfs.list) $*.empty.vcf" "$@.tmp.unfiltered.vcf"; \
	} >> $@.split.sub.mk;
	# Calling all chromosomes in parallel (nothing to do if no chromosome, empty VCF will be used)
	+make -f $@.split.sub.mk $@.tmp.unfiltered.vcf
	# Filters
	bash $(STARK_FOLDER_BIN)/filter_vcf.sh --input $@.tmp.unfiltered.vcf --output $@ --threads $(THREADS_BY_CALLER) --bcftools $(BCFTOOLS) $(VARSCAN_FILTERS_FLAGS)
	# Cleaning
	rm -rf $@.tmp* $@.split* $@.idx $@.tbi


# VARSCAN Release documentation

RELEASE_COMMENT := "\#\# CALLING VARSCAN: VARSCAN tool identify variants from mpileup file, generated from aligned BAM"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING: VarScan identify variants in parallel to identify somatic variants with low VAF and generate *.VarScan.vcf files. Parameters VARSCAN_SNP_FLAGS='$(VARSCAN_SNP_FLAGS)', VARSCAN_INDEL_FLAGS='$(VARSCAN_INDEL_FLAGS)' and VARSCAN_FILTERS_FLAGS='$(VARSCAN_FILTERS_FLAGS)' "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:VarScan:SAMTOOLS/VARSCAN - by default:VARSCAN_SNP_FLAGS='$(VARSCAN_SNP_FLAGS)', VARSCAN_INDEL_FLAGS='$(VARSCAN_INDEL_FLAGS)', VARSCAN_FILTERS_FLAGS='$(VARSCAN_FILTERS_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
