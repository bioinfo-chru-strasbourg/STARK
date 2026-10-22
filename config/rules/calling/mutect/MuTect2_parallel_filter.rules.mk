############################
# MUTECT Calling Rules
# Release: 1.0.0
# Date: 21/10/2026
# Author: Antony Le Bechec
############################

# Release note
# 0.9.2-03/02/2023: Initial release of MuTect calling rules.
# 1.0.0-21/10/2026: Update release date and version for MuTect2 parallel filter rules.

#################################
# MUTECT2 GATK4 PARALLEL FILTER #
#################################

# Parameters
THREADS_GATK4_MUTECT2_PARALLEL_FILTER?=1 #$(THREADS_BY_CALLER)
MINPRUNING_GATK4_MUTECT2_PARALLEL_FILTER?=20
MAXREADS_GATK4_MUTECT2_PARALLEL_FILTER?=1000
DPMIN_MUTECT2_PARALLEL_FILTER?=30

GATK4_MUTECT2_PARALLEL_FILTER_FLAGS_SHARED?=--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
	--max-reads-per-alignment-start $(MAXREADS_GATK4_MUTECT2_PARALLEL_FILTER) \
	--dont-use-soft-clipped-bases true \
	--min-pruning $(MINPRUNING_GATK4_MUTECT2_PARALLEL_FILTER) \
	--callable-depth $(DPMIN_MUTECT2_PARALLEL_FILTER) \
	--verbosity ERROR \
	--native-pair-hmm-threads $(THREADS_GATK4_MUTECT2_PARALLEL_FILTER)



%.MuTect2_parallel_filter$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.design.bed
	# Clean and create split directory
	rm -rf $@.split && mkdir -p $@.split
	# Generate idxstats for the BAM file (contigs with reads only)
	$(SAMTOOLS) idxstats $< | awk '$$3>0' > $@.split/idxstats.txt
	# List of chromosomes (BAM header order), restricted to the BED if it has data lines
	if (( $$(grep -v -e "^#" -e "^$$" -c $*.design.bed) )); then \
		grep -v -e "^#" -e "^$$" $*.design.bed | cut -f1 | sort -u > $@.split/bedchroms.txt; \
		cut -f1 $@.split/idxstats.txt | grep -Fx -f $@.split/bedchroms.txt > $@.split/chroms.txt || [ "$$?" -eq 1 ]; \
	else \
		cut -f1 $@.split/idxstats.txt > $@.split/chroms.txt; \
	fi;
	# Execution order: most reads first (merge order stays chroms.txt)
	awk 'NR==FNR{k[$$1];next} $$1 in k{print $$3"\t"$$1}' $@.split/chroms.txt $@.split/idxstats.txt \
		| sort -k1,1nr | cut -f2 > $@.split/order.txt;
	# Generate the sub Makefile: one explicit rule per chromosome (empty "all" if no chromosome)
	{ \
	SAMPLE=$$(basename $< | cut -d"." -f1); \
	printf '.DELETE_ON_ERROR:\n'; \
	printf 'all: %s\n' "$$(sed 's|.*|$@.split/&.vcf|' $@.split/order.txt | tr '\n' ' ')"; \
	for c in $$(cat $@.split/order.txt); do \
		L_OPT="-L $$c"; \
		if [ -s $@.split/bedchroms.txt ]; then \
			awk -F'\t' -v c="$$c" '$$1==c' $*.design.bed > $@.split/$$c.bed; \
			L_OPT="-L $@.split/$$c.bed"; \
		fi; \
		printf '%s:\n' "$@.split/$$c.vcf"; \
		printf '\t$(JAVA) $(JAVA_FLAGS) -jar $(GATK4) Mutect2 $(GATK4_MUTECT2_PARALLEL_FILTER_FLAGS_SHARED) -R $(GENOME) -I $< -tumor %s %s -O %s\n' "$$SAMPLE" "$$L_OPT" "$@.split/$$c.vcf"; \
	done; \
	} > $@.split/sub.mk
	# Calling by MuTect2, all chromosomes in parallel (nothing to do if no chromosome)
	+make -f $@.split/sub.mk all
	# No chromosome: the final output is the empty VCF
	if [ ! -s $@.split/chroms.txt ]; then cp $*.empty.vcf $@; fi
	# Otherwise: merge parts (VCF + stats required by FilterMutectCalls), filter, normalize, sort, filter
	if [ -s $@.split/chroms.txt ]; then \
		set -e; \
		if [ "$$(grep -c '' $@.split/chroms.txt)" -eq 1 ]; then \
			c=$$(cat $@.split/chroms.txt); \
			cp $@.split/$$c.vcf $@.tmp.unfiltered.vcf; \
			cp $@.split/$$c.vcf.stats $@.tmp.unfiltered.vcf.stats; \
		else \
			sed 's|.*|$@.split/&.vcf|' $@.split/chroms.txt > $@.split/vcfs.list; \
			$(BCFTOOLS) concat -f $@.split/vcfs.list -o $@.tmp.unfiltered.vcf; \
			$(JAVA) $(JAVA_FLAGS) -jar $(GATK4) MergeMutectStats \
				$$(sed 's|.*|-stats $@.split/&.vcf.stats|' $@.split/chroms.txt | tr '\n' ' ') \
				-O $@.tmp.unfiltered.vcf.stats; \
		fi; \
		$(JAVA) $(JAVA_FLAGS) -jar $(GATK4) FilterMutectCalls \
			-R $(GENOME) \
			-V $@.tmp.unfiltered.vcf \
			-O $@.tmp.FilterMutectCalls.vcf; \
		grep "^##" $@.tmp.FilterMutectCalls.vcf | sed s/ID=TLOD,Number=A/ID=TLOD,Number=./gi > $@.tmp.FilterMutectCalls.TLOD.vcf; \
		grep "^##" -v $@.tmp.FilterMutectCalls.vcf | cut -f1-10 >> $@.tmp.FilterMutectCalls.TLOD.vcf; \
		$(JAVA) -jar $(PICARD) SortVcf -I $@.tmp.FilterMutectCalls.TLOD.vcf -O $@.tmp.FilterMutectCalls.TLOD.sorted.vcf -SD $(DICT); \
		$(BCFTOOLS) view -i 'FORMAT/DP>=$(DPMIN_MUTECT2_PARALLEL_FILTER)' $@.tmp.FilterMutectCalls.TLOD.sorted.vcf | $(BCFTOOLS) view -f PASS -e 'GT="0/0"' -o $@; \
	fi
	# Clean
	rm -rf $@.split
	rm -f $@.tmp* $@.idx



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING: MUTECT2 with GATK4 in parallel to identify somatic variants and generate *.MuTect2_parallel_filter.vcf files, with MuTect2 and DP filtration. Parameters GATK4_MUTECT2_PARALLEL_FILTER_FLAGS_SHARED='$(GATK4_MUTECT2_PARALLEL_FILTER_FLAGS_SHARED)', DPMIN_MUTECT2_PARALLEL_FILTER='$(DPMIN_MUTECT2_PARALLEL_FILTER)' "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:MuTect2_parallel_filter:MuTect2 with GATK4 - by default with MuTect2 and DP filtration. Parameters GATK4_MUTECT2_PARALLEL_FILTER_FLAGS_SHARED='$(GATK4_MUTECT2_PARALLEL_FILTER_FLAGS_SHARED)', DPMIN_MUTECT2_PARALLEL_FILTER='$(DPMIN_MUTECT2_PARALLEL_FILTER)' "
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
