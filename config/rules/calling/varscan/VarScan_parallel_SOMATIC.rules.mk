############################
# VarScan Calling Rules
# Release: 1.0.0
# Date: 21/10/2026
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
# 1.0.0-21/10/2026: Update release date and version for VarScan parallel rules.



####################
# VARSCAN Standard #
####################

VARSCAN_PARALLEL_SOMATIC_VAF=0.01
VARSCAN_PARALLEL_SOMATIC_ALT=5
VARSCAN_PARALLEL_SOMATIC_DP=50
VARSCAN_PARALLEL_SOMATIC_PVAL=1e-3
VARSCAN_PARALLEL_SOMATIC_BOTH_OPTIONS= --min-var-freq $(VARSCAN_PARALLEL_SOMATIC_VAF) --min-reads2 $(VARSCAN_PARALLEL_SOMATIC_ALT) --min-coverage $(VARSCAN_PARALLEL_SOMATIC_DP) --p-value $(VARSCAN_PARALLEL_SOMATIC_PVAL)
VARSCAN_PARALLEL_SOMATIC_SNP_OPTIONS= $(VARSCAN_PARALLEL_SOMATIC_BOTH_OPTIONS) --min-avg-qual 30
VARSCAN_PARALLEL_SOMATIC_INDEL_OPTIONS= $(VARSCAN_PARALLEL_SOMATIC_BOTH_OPTIONS) --min-avg-qual 10


%.VarScan_parallel_SOMATIC$(POST_CALLING).SNP.vcf: %.bam.mpileup %.empty.vcf
	# Clean and create split directory
	rm -rf $@.split && mkdir -p $@.split && : > $@.split/chroms.txt
	# Split the mpileup by chromosome in a single pass (chroms.txt keeps the mpileup order)
	if [ -s $< ]; then \
		awk -F'\t' -v d="$@.split" '$$1!=p{ if(f!="") close(f); p=$$1; f=d"/"p".mpileup"; print p >> (d"/chroms.txt") } { print > f }' $<; \
		basename $* | cut -d"." -f1 > $@.split/sample.txt; \
	fi;
	# Execution order: biggest chromosomes first (merge order stays chroms.txt)
	for c in $$(cat $@.split/chroms.txt); do echo "$$(stat -c %s $@.split/$$c.mpileup) $$c"; done | sort -k1,1nr | cut -d" " -f2 > $@.split/order.txt
	# Generate the sub Makefile: one explicit rule per chromosome (empty "all" if no chromosome)
	{ \
	printf '.DELETE_ON_ERROR:\n'; \
	printf 'all: %s\n' "$$(sed 's|.*|$@.split/&.vcf|' $@.split/order.txt | tr '\n' ' ')"; \
	for c in $$(cat $@.split/order.txt); do \
		printf '%s:\n' "$@.split/$$c.vcf"; \
		printf '\tcat %s | $(JAVA) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_PARALLEL_SOMATIC_SNP_OPTIONS) --vcf-sample-list %s > %s && rm -f %s\n' "$@.split/$$c.mpileup" "$@.split/sample.txt" "$@.split/$$c.vcf" "$@.split/$$c.mpileup"; \
	done; \
	} > $@.split/sub.mk
	# Calling by VarScan, all chromosomes in parallel (nothing to do if no chromosome)
	+make -f $@.split/sub.mk all
	# Merge: header from the first part having one, then the variant lines in mpileup order
	if [ -s $@.split/chroms.txt ]; then \
		: > $@; \
		for c in $$(cat $@.split/chroms.txt); do \
			if grep -q "^#" $@.split/$$c.vcf; then grep "^#" $@.split/$$c.vcf > $@; break; fi; \
		done; \
		for c in $$(cat $@.split/chroms.txt); do grep -hv "^#" $@.split/$$c.vcf >> $@ || true; done; \
	fi;
	if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	# Remove intermediate files
	rm -rf $@.split
	rm -f $@.idx $@.unfiltered.vcf*


%.VarScan_parallel_SOMATIC$(POST_CALLING).InDel.vcf: %.bam.mpileup %.empty.vcf
	# Clean and create split directory
	rm -rf $@.split && mkdir -p $@.split && : > $@.split/chroms.txt
	# Split the mpileup by chromosome in a single pass (chroms.txt keeps the mpileup order)
	if [ -s $< ]; then \
		awk -F'\t' -v d="$@.split" '$$1!=p{ if(f!="") close(f); p=$$1; f=d"/"p".mpileup"; print p >> (d"/chroms.txt") } { print > f }' $<; \
		basename $* | cut -d"." -f1 > $@.split/sample.txt; \
	fi;
	# Execution order: biggest chromosomes first (merge order stays chroms.txt)
	for c in $$(cat $@.split/chroms.txt); do echo "$$(stat -c %s $@.split/$$c.mpileup) $$c"; done | sort -k1,1nr | cut -d" " -f2 > $@.split/order.txt
	# Generate the sub Makefile: one explicit rule per chromosome (empty "all" if no chromosome)
	{ \
	printf '.DELETE_ON_ERROR:\n'; \
	printf 'all: %s\n' "$$(sed 's|.*|$@.split/&.vcf|' $@.split/order.txt | tr '\n' ' ')"; \
	for c in $$(cat $@.split/order.txt); do \
		printf '%s:\n' "$@.split/$$c.vcf"; \
		printf '\tcat %s | $(JAVA) -jar $(VARSCAN) mpileup2indel --output-vcf 1 $(VARSCAN_PARALLEL_SOMATIC_INDEL_OPTIONS) --vcf-sample-list %s > %s && rm -f %s\n' "$@.split/$$c.mpileup" "$@.split/sample.txt" "$@.split/$$c.vcf" "$@.split/$$c.mpileup"; \
	done; \
	} > $@.split/sub.mk
	# Calling by VarScan, all chromosomes in parallel (nothing to do if no chromosome)
	+make -f $@.split/sub.mk all
	# Merge: header from the first part having one, then the variant lines in mpileup order
	if [ -s $@.split/chroms.txt ]; then \
		: > $@; \
		for c in $$(cat $@.split/chroms.txt); do \
			if grep -q "^#" $@.split/$$c.vcf; then grep "^#" $@.split/$$c.vcf > $@; break; fi; \
		done; \
		for c in $$(cat $@.split/chroms.txt); do grep -hv "^#" $@.split/$$c.vcf >> $@ || true; done; \
	fi;
	if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	# Remove intermediate files
	rm -rf $@.split
	rm -f $@.idx $@.unfiltered.vcf*


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING VARSCAN '$(MK_RELEASE)': VARSCAN tool in parallel identify variants from mpileup file, generated from aligned BAM"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_CMD := $(shell echo "\#\# VARSCAN_parallel_SOMATIC: identify variants and generate *.VarScan_parallel_SOMATIC.vcf files with parameters SAMTOOLS mpileup and VARSCAN standard application 'VARSCAN_PARALLEL_SOMATIC_SNP_OPTIONS: $(VARSCAN_PARALLEL_SOMATIC_SNP_OPTIONS)' 'VARSCAN_PARALLEL_SOMATIC_INDEL_OPTIONS: $(VARSCAN_PARALLEL_SOMATIC_INDEL_OPTIONS)' " >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:VarScan_parallel_SOMATIC:SAMTOOLS/VARSCAN - by default in parallel:VARSCAN_PARALLEL_SOMATIC_SNP_OPTIONS='$(VARSCAN_PARALLEL_SOMATIC_SNP_OPTIONS)', VARSCAN_PARALLEL_SOMATIC_INDEL_OPTIONS='$(VARSCAN_PARALLEL_SOMATIC_INDEL_OPTIONS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
