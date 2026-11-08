############################
# MUTECT Calling Rules
# Release: 1.0.0
# Date: 03/11/2026
# Author: Antony Le Bechec
############################

# Release note
# 0.9.2-03/02/2023: Initial release of MuTect calling rules.
# 1.0.0-03/11/2026: Parallelization and filtration


##########
# MUTECT #
##########

# MUTECT Flags
# VCFDBSNP_WITH_GATK variable depends on USE_VCFDBSNP_WITH_GATK setting on Apps
# INTERVAL_PADDING variable defines the padding around intervals for all callers, settings on Apps
# This variable can be defined directly on Apps
MUTECT_FLAGS?= \
	$(VCFDBSNP_WITH_GATK) \
	--interval_padding $(INTERVAL_PADDING) \
	--read_filter BadCigar \
	--baq OFF \
	--downsampling_type NONE \
	--downsample_to_fraction 1 \
	--min_qscore 17 \
	--allow_potentially_misencoded_quality_scores

# MUTECT filters Flags
# This variable can be defined directly on Apps
MUTECT_FILTERS_FLAGS?= \
	--exclude-genotypes '0/0' \
	--exclude 'FILTER="REJECT"' \
	--include 'FORMAT/DP>=30 && FORMAT/FA>=0.01' \
	--filter 'PASS' \
	--sort


# MUTECT target rule
%.MuTect$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.design.bed %.bam.chroms
	# List of vcf files to be generated per chromosome
	cut -f1 $*.bam.chroms | sed 's|.*|$@.split.&.vcf.gz|' | tr '\n' ' ' > $@.split.vcfs.list;
	# Generate bed per chromosome
	awk '{print > ("$@.split."$$1".bed")}' $*.design.bed
	# Generate the sub Makefile: one explicit rule per chromosome (empty "all" if no chromosome)
	{ \
	printf 'all: %s\n' "$$(cat $@.split.vcfs.list)"; \
	for c in $$(cut -f1 $*.bam.chroms); do \
		sample=$$(bcftools view $*.empty.vcf | grep "^#CHROM" | cut -f10); \
		if [ -s $@.split.$$c.bed ]; then \
			L_ARG="$@.split.$$c.bed"; \
		else \
			L_ARG="$$c"; \
		fi; \
		printf '%s:\n' "$@.split.$$c.vcf.gz"; \
		printf '\t$(JAVA7) $(JAVA_FLAGS) -jar $(MUTECT) --analysis_type MuTect $(MUTECT_FLAGS) -L %s --reference_sequence $(GENOME) --input_file:tumor $< --vcf %s -o /dev/null --num_cpu_threads_per_data_thread 1\n' "$$L_ARG" "$@.split.$$c.tmp.vcf.gz"; \
		printf '\t$(BCFTOOLS) view %s -s %s -Oz -o %s\n' "$@.split.$$c.tmp.vcf.gz" "$$sample" "$@.split.$$c.vcf.gz"; \
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
	bash $(STARK_FOLDER_BIN)/filter_vcf.sh --input $@.tmp.unfiltered.vcf --output $@ --threads $(THREADS_BY_CALLER) --bcftools $(BCFTOOLS) $(MUTECT_FILTERS_FLAGS)
	# Cleaning
	rm -rf $@.tmp* $@.split* $@.idx $@.tbi


# MUTECT Release documentation

RELEASE_COMMENT := "\#\# CALLING MuTect Caller '$(MK_RELEASE)': MuTect tool identify variants from aligned BAM"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING: MUTECT to identify somatic variants and generate *.MuTect.vcf files. Parameteres MUTECT_FLAGS='$(MUTECT_FLAGS)' and MUTECT_FILTERS_FLAGS='$(MUTECT_FILTERS_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:MuTect:MuTect - by default, no filtering:MUTECT_FLAGS='$(MUTECT_FLAGS)' and MUTECT_FILTERS_FLAGS='$(MUTECT_FILTERS_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
