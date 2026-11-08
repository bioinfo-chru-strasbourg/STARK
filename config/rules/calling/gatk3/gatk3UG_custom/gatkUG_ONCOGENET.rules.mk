############################
# GATK Calling Rules
# Release: 1.0.0
# Date: 23/10/2026
# Author: Antony Le Bechec
############################


# Release note
# 0.9.1beta-10/03/2015: change genome reference location, in the file %.genome
# 0.9.1.1beta-30/04/2015: VariantFiltration correction
# 0.9.2b-23/11/2015: adding gatkUG_CPSGEN_MASTR pipeline
# 0.9.3b-25/01/2016: adding gatkUG_DIAG and gatkHC_DIAG pipelines
# 0.9.3.1b-08/02/2016: adding gatkUG_DIAGGEN pipeline
# 0.9.3.2b-03/03/2016: adding gatkUG_DIAG_PARIS pipeline
# 0.9.3.3b-11/04/2016: change gatkHC pipeline to fit with standard analysis
# 0.9.3.4b-03/05/2016: Modification of gatkUG_HUSHEMATO rule
# 0.9.3.5b-04/05/2016: Rewrite rules and update release information
# 0.9.3.6b-10/11/2017: adding gatkUG_ONCOGENET pipeline
# 0.9.3.7b-22/03/2019: Add --dontUseSoftClippedBases for GATKHC
# 0.9.3.8-29/07/2022: Remove --dontUseSoftClippedBases for GATKHC
# 0.9.3.9-29/07/2022: Add --dontUseSoftClippedBases for GATKHC, add GATKUG_LONG_INDELS and GATKHC_LONG_INDELS
# 0.9.5.0-31/10/2025: Parallelised GATKUG
# 1.0.0-23/10/2026: Parallelised implementation for gatkUG oncogenetic variant calling


GATKUG_FLAGS_SHARED?=

##########
# gatkUG #
##########


# GATKUG DPMIN filter
DPMIN_GATKUG_ONCOGENET=4

# INTERVAL PADDING for GATKUG
INTERVAL_PADDING?=0
INTERVAL_PADDING_GATKUG_ONCOGENET?=$(INTERVAL_PADDING)

# GATKUG Flags
GATKUG_ONCOGENET_FLAGS=$(GATKUG_FLAGS_SHARED) \
	$(VCFDBSNP_WITH_GATK) \
	--interval_padding $(INTERVAL_PADDING_GATKUG_ONCOGENET) \
	--num_cpu_threads_per_data_thread 1 \
	--genotype_likelihoods_model BOTH \
	--read_filter BadCigar \
	--baq OFF \
	--downsampling_type NONE \
	--min_indel_fraction_per_sample 0.05 \
	--min_indel_count_for_genotyping 5 \
	--max_deletion_fraction 0.05 \
	--standard_min_confidence_threshold_for_calling 10 \
	--downsample_to_fraction 1 \
	--min_base_quality_score 17 \
	--allow_potentially_misencoded_quality_scores


# GATKUG ONCOGENET target rule
%.gatkUG_ONCOGENET$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.design.bed
	# Clean and create split directory
	rm -rf $@.split && mkdir -p $@.split && : > $@.split/chroms.txt
	# Contigs with reads (idxstats failure stops the rule)
	$(SAMTOOLS) idxstats $< > $@.split/idxstats.raw
	awk -F'\t' '$$1!="*" && $$3+$$4>0' $@.split/idxstats.raw > $@.split/idxstats.txt
	# Chromosomes of the design (empty file = empty design)
	awk -F'\t' '!/^#/ && !/^(track|browser)/ && NF{print $$1}' $*.design.bed | sort -u > $@.split/bedchroms.txt
	# Chromosomes to call (BAM header order): design chromosomes with reads, or all chromosomes with reads if the design is empty
	if [ -s $@.split/bedchroms.txt ]; then \
		cut -f1 $@.split/idxstats.txt | grep -Fx -f $@.split/bedchroms.txt > $@.split/chroms.txt || { rc=$$?; [ $$rc -eq 1 ] || exit $$rc; }; \
	else \
		echo "[INFO] GATKUG: empty design for $<. Calling all chromosomes with reads"; \
		cut -f1 $@.split/idxstats.txt > $@.split/chroms.txt; \
	fi
	# Execution order: most reads first (merge order stays chroms.txt)
	awk -F'\t' 'NR==FNR{k[$$1];next} $$1 in k{print $$3+$$4"\t"$$1}' $@.split/chroms.txt $@.split/idxstats.txt \
		| sort -k1,1nr | cut -f2 > $@.split/order.txt
	# Generate the sub Makefile: one explicit rule per chromosome (empty "all" if no chromosome)
	{ \
	printf '.DELETE_ON_ERROR:\n'; \
	printf 'all: %s\n' "$$(sed 's|.*|$@.split/&.vcf.gz|' $@.split/order.txt | tr '\n' ' ')"; \
	for c in $$(cat $@.split/order.txt); do \
		if [ -s $@.split/bedchroms.txt ]; then \
			awk -F'\t' -v c="$$c" '$$1==c' $*.design.bed > $@.split/$$c.bed; \
			L_ARG="$@.split/$$c.bed"; \
		else \
			L_ARG="$$c"; \
		fi; \
		printf '%s:\n' "$@.split/$$c.vcf.gz"; \
		printf '\t$(JAVA8) $(JAVA_FLAGS) -jar $(GATK3) $(GATKUG_ONCOGENET_FLAGS) -L %s -T UnifiedGenotyper -R $(GENOME) -I $< -o %s\n' "$$L_ARG" "$@.split/$$c.vcf.gz"; \
	done; \
	} > $@.split/sub.mk
	# Calling by GATK UnifiedGenotyper, all chromosomes in parallel (nothing to do if no chromosome)
	+make -f $@.split/sub.mk all
	# No chromosome: empty VCF. Otherwise: merge + normalize + DP filter in a single pipe
	if [ ! -s $@.split/chroms.txt ]; then \
		cp $*.empty.vcf $@; \
		echo "[WARNING] GATKUG: no chromosome with reads (and BED lines if a design is given) for $<. Empty VCF"; \
	else \
		set -eo pipefail; \
		sed 's|.*|$@.split/&.vcf.gz|' $@.split/chroms.txt > $@.split/vcfs.list; \
		if [ "$$(grep -c '' $@.split/vcfs.list)" -eq 1 ]; then \
			MERGE="$(BCFTOOLS) view -O u $$(cat $@.split/vcfs.list)"; \
		else \
			MERGE="$(BCFTOOLS) concat -f $@.split/vcfs.list -O u --threads $(THREADS_BY_CALLER)"; \
		fi; \
		$$MERGE \
			| $(BCFTOOLS) norm -m -any -O u --threads $(THREADS_BY_CALLER) \
			| $(BCFTOOLS) view -i "FORMAT/DP>$(DPMIN_GATKUG_ONCOGENET)" --threads $(THREADS_BY_CALLER) -o $@; \
	fi
	# Cleaning
	rm -rf $@.split $@.idx


RELEASE_COMMENT := "\#\# CALLING GATKUG ONCOGENET identify variants and generate *.gatkUG_ONCOGENET.vcf files with parameters: GATKUG_ONCOGENET_FLAGS='$(GATKUG_ONCOGENET_FLAGS)' and DPMIN_GATKUG_ONCOGENET='$(DPMIN_GATKUG_ONCOGENET)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkUG_ONCOGENET:GATK Unified Genotyper - designed for oncogenetic variant discovery:GATKUG_ONCOGENET_FLAGS='$(GATKUG_ONCOGENET_FLAGS)' and DPMIN_GATKUG_ONCOGENET='$(DPMIN_GATKUG_ONCOGENET)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
