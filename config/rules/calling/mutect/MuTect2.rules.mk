############################
# MUTECT Calling Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE=1.2.0
MK_DATE=03/11/2026

# Release note
# 0.9.2-03/02/2023: Initial release of MuTect calling rules.
# 1.0.0-21/10/2026: Parallelize MuTect2 calling.
# 1.1.0-23/10/2026: Improve filter parameters.
# 1.2.0-03/11/2026: Normalize parallelization and filter parameters


###########
# MUTECT2 #
###########

# Flags
# This variable can be defined directly on Apps
MUTECT2_FLAGS?= \
	$(VCFDBSNP_WITH_GATK) \
	--interval-padding $(INTERVAL_PADDING) \
	--min-base-quality-score 17 \
	--max-reads-per-alignment-start 1000 \
	--min-pruning 20 \
	--callable-depth 30 \
	--dont-use-soft-clipped-bases true \
	--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
	--verbosity ERROR

# Filters Flags
# This variable can be defined directly on Apps
MUTECT2_FILTERS_FLAGS?= \
	--exclude-genotypes '0/0' \
	--exclude 'FILTER="FAIL"' \
	--include 'FORMAT/DP>=30 && FORMAT/AF>=0.01' \
	--filter 'PASS' \
	--sort


# Rule
%.MuTect2$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.design.bed %.bam.chroms
	# List of vcf and stats files to be generated per chromosome
	cut -f1 $*.bam.chroms | sed 's|.*|$@.split.&.vcf.gz|' | tr '\n' ' ' > $@.split.vcfs.list;
	cut -f1 $*.bam.chroms | sed 's|.*|-stats $@.split.&.vcf.gz.stats|' | tr '\n' ' ' > $@.split.vcfs.stats.param;
	# Generate bed per chromosome
	awk '{print > ("$@.split."$$1".bed")}' $*.design.bed
	# Generate the sub Makefile: one explicit rule per chromosome (empty "all" if no chromosome)
	{ \
	printf 'all: %s\n' "$$(cat $@.split.vcfs.list)"; \
	for c in $$(cut -f1 $*.bam.chroms); do \
		if [ -s $@.split.$$c.bed ]; then \
			L_ARG="$@.split.$$c.bed"; \
		else \
			L_ARG="$$c"; \
		fi; \
		printf '%s:\n' "$@.split.$$c.vcf.gz"; \
		printf '\t$(JAVA) $(JAVA_FLAGS) -jar $(GATK) Mutect2 $(MUTECT2_FLAGS) -L %s -R $(GENOME) -I $< -O %s --native-pair-hmm-threads 1\n' "$$L_ARG" "$@.split.$$c.vcf.gz"; \
	done; \
	} > $@.split.sub.mk
	# Add rule to concat per-chromosome VCFs into the final VCF (empty VCF if no chromosome)
	# Perform MuTectCalls using stats on the merged per-chromosome VCFs (if any stats are available)
	{ \
		printf '%s: %s\n' "$@.tmp.unfiltered.vcf" "$$(cat $@.split.vcfs.list) $*.empty.vcf"; \
		printf '\t$(BCFTOOLS) concat %s -Ov --threads $(THREADS_BY_CALLER) -o %s\n' "$$(cat $@.split.vcfs.list) $*.empty.vcf" "$@.tmp.merged.vcf"; \
		printf '\tif [ -s %s ]; then \\\n' "$@.split.vcfs.stats.param"; \
		printf '\t\t$(JAVA) $(JAVA_FLAGS) -jar $(GATK) MergeMutectStats %s -O %s && \\\n' "$$(cat $@.split.vcfs.stats.param)" "$@.tmp.merged.vcf.stats"; \
		printf '\t\t$(JAVA) $(JAVA_FLAGS) -jar $(GATK) FilterMutectCalls -V %s -O %s -R $(GENOME); \\\n' "$@.tmp.merged.vcf" "$@.tmp.unfiltered.vcf"; \
		printf '\telse \\\n'; \
		printf '\t\tcp %s %s; \\\n' "$@.tmp.merged.vcf" "$@.tmp.unfiltered.vcf"; \
		printf '\tfi\n'; \
	} >> $@.split.sub.mk;
	# Calling all chromosomes in parallel (nothing to do if no chromosome, empty VCF will be used)
	+make -f $@.split.sub.mk $@.tmp.unfiltered.vcf
	# Filters
	bash $(STARK_FOLDER_BIN)/filter_vcf.sh --input $@.tmp.unfiltered.vcf --output $@ --threads $(THREADS_BY_CALLER) --bcftools $(BCFTOOLS) $(MUTECT2_FILTERS_FLAGS)
	# Cleaning
	rm -rf $@.tmp* $@.split* $@.idx $@.tbi


# Documentation
$(file >>$(RELEASE_INFOS),## CALLER MuTect2 '$(MK_RELEASE)': MuTect2 caller with GATK4 in parallel to identify somatic variants and generate *.MuTect2.vcf files, with MuTect2 and filtrations. Parameters MUTECT2_FLAGS='$(MUTECT2_FLAGS)', MUTECT2_FILTERS_FLAGS='$(MUTECT2_FILTERS_FLAGS)')
$(file >>$(PIPELINES_INFOS),CALLER:MuTect2:MuTect2 with GATK4 - by default with MuTect2 and filtration:MUTECT2_FLAGS='$(MUTECT2_FLAGS)', MUTECT2_FILTERS_FLAGS='$(MUTECT2_FILTERS_FLAGS)')
