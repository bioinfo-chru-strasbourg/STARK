############################
# OutLyzer Calling Rules
# Release: 0.9.2
# Date: 03/02/2023
# Author: Antony Le Bechec
############################


#####################
# OUTLYZER Filtered #
#####################

# Parameters
OUTLYZER_FILTERED_EXCLUDE?=
OUTLYZER_FILTERED_INCLUDE?=FORMAT/DP>=30 && FORMAT/VAF>=0.01


%.outLyzer_filtered$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.design.bed
	# Clean and create split directory
	rm -rf $@.split && mkdir -p $@.split && : > $@.split/chroms.txt
	# List of chromosomes (BAM header order): BED chromosomes having reads in the BAM
	if [ -s $*.design.bed ]; then \
		$(SAMTOOLS) idxstats $< | awk '$$3>0' > $@.split/idxstats.txt; \
		grep -v -e "^#" -e "^$$" $*.design.bed | cut -f1 | sort -u > $@.split/bedchroms.txt; \
		cut -f1 $@.split/idxstats.txt | grep -Fx -f $@.split/bedchroms.txt > $@.split/chroms.txt || true; \
	fi;
	# Execution order: most BED bases first (merge order stays chroms.txt)
	awk -F'\t' 'NR==FNR{k[$$1];next} !/^#/ && NF && ($$1 in k){s[$$1]+=$$3-$$2} END{for(c in s) print s[c]"\t"c}' $@.split/chroms.txt $*.design.bed \
		| sort -k1,1nr | cut -f2 > $@.split/order.txt
	# Generate the sub Makefile: one explicit rule per chromosome (empty "all" if no chromosome)
	{ \
	printf '.DELETE_ON_ERROR:\n'; \
	printf 'all: %s\n' "$$(sed 's|.*|$@.split/&.done|' $@.split/order.txt | tr '\n' ' ')"; \
	for c in $$(cat $@.split/order.txt); do \
		awk -F'\t' -v c="$$c" '$$1==c' $*.design.bed > $@.split/$$c.bed; \
		printf '%s:\n' "$@.split/$$c.done"; \
		printf '\tmkdir -p %s && $(PYTHON) $(OUTLYZER) calling -pythonPath=$(PYTHON) -samtools=$(SAMTOOLS) -bed %s -bam $< -ref $(GENOME) -output %s/ -core 1 -verbose 1 && touch %s\n' "$@.split/$$c.out" "$@.split/$$c.bed" "$@.split/$$c.out" "$@.split/$$c.done"; \
	done; \
	} > $@.split/sub.mk
	# Calling by OutLyzer, all chromosomes in parallel (nothing to do if no chromosome)
	+make -f $@.split/sub.mk all
	# List the non-empty VCFs produced by OutLyzer (chromosome order)
	for c in $$(cat $@.split/chroms.txt); do find $@.split/$$c.out -name "*.vcf" -size +0 | sort; done > $@.split/vcfs.list
	# No VCF: empty VCF. Otherwise: merge (single header), normalize, sort + contigs, filter
	if [ ! -s $@.split/vcfs.list ]; then \
		cp $*.empty.vcf $@; \
		if [ -s $*.design.bed ]; then \
			echo "[WARNING] No OutLyzer VCF produced. Empty OutLyzer VCF"; \
		else \
			echo "[ERROR] Empty BED. Empty OutLyzer VCF"; \
		fi; \
	else \
		set -e; \
		grep "^#" "$$(head -n 1 $@.split/vcfs.list)" > $@.tmp.merged.vcf; \
		xargs -r grep -hv "^#" < $@.split/vcfs.list >> $@.tmp.merged.vcf || true; \
		n_parts=$$(xargs -r grep -hv "^#" < $@.split/vcfs.list | wc -l); \
		n_merged=$$(grep -vc "^#" $@.tmp.merged.vcf || true); \
		if [ "$$n_parts" -ne "$$n_merged" ]; then echo "[ERROR] OutLyzer merge: $$n_parts variant lines in parts, $$n_merged in merged VCF" >&2; exit 1; fi; \
		if [ "$$(xargs -r grep -h "^#CHROM" < $@.split/vcfs.list | sort -u | wc -l)" -ne 1 ]; then echo "[WARNING] OutLyzer VCF #CHROM lines differ between parts" >&2; fi; \
		awk -f $(STARK_FOLDER_BIN)/outlyzer_norm.awk $@.tmp.merged.vcf > $@.tmp.normalized.vcf; \
		$(JAVA) -jar $(PICARD) SortVcf -I $@.tmp.normalized.vcf -O $@.tmp.normalized.sorted.vcf -SD $(DICT); \
		CMD="bcftools view -e 'GT=\"0/0\"' --threads=$(THREADS_BY_CALLER) $@.tmp.normalized.sorted.vcf"; \
		if [ -n "$(OUTLYZER_FILTERED_EXCLUDE)" ]; then \
			CMD="$$CMD | bcftools view -e '$(OUTLYZER_FILTERED_EXCLUDE)' --threads=$(THREADS_BY_CALLER)"; \
		fi; \
		if [ -n "$(OUTLYZER_FILTERED_INCLUDE)" ]; then \
			CMD="$$CMD | bcftools view -i '$(OUTLYZER_FILTERED_INCLUDE)' --threads=$(THREADS_BY_CALLER)"; \
		fi; \
		CMD="$$CMD > $@" ; \
		echo $$CMD; \
		eval $$CMD; \
	fi
	# Cleaning
	rm -rf $@.split $@.tmp* $@.idx


RELEASE_CMD := $(shell echo "\#\# OutLyzer: OutLyzer caller to detect variant with low VAF, generate *.outLyzer_filtered.vcf files. VCF file is modified to include GT/Genotype information. Filtered on DP and VAF. DPMIN_OUTLYZER_FILTERED='$(DPMIN_OUTLYZER_FILTERED)', VAF_OUTLYZER_FILTERED='$(VAF_OUTLYZER_FILTERED)' " >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:outLyzer_filtered:OutLyzer - OutLyzer caller, GT info added, filter on DP and VAF:DP and VAF"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
