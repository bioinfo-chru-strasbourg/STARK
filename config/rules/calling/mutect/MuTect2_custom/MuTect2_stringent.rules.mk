############################
# MUTECT Calling Rules
# Release: 1.1.0
# Date: 23/10/2026
# Author: Antony Le Bechec
############################

# Release note
# 0.9.2-03/02/2023: Initial release of MuTect calling rules.
# 1.0.0-21/10/2026: Update release date and version for MuTect2 parallel filter rules.
# 1.1.0-23/10/2026: Minor updates for MuTect2 filtered rules.

#################
# MUTECT2 GATK4 #
#################

# Parameters
GATK4_MUTECT2_STRINGENT_EXCLUDE?=
GATK4_MUTECT2_STRINGENT_INCLUDE?=FORMAT/AF>=0.01 && FORMAT/AF<0.80
GATK4_MUTECT2_STRINGENT_FILTER?=PASS

# Shared flags for MuTect2 parallel stringent calling
GATK4_MUTECT2_STRINGENT_FLAGS_SHARED?=--native-pair-hmm-threads 1 \
	--max-reads-per-alignment-start 1000 \
	--min-pruning 20 \
	--callable-depth 30 \
	--dont-use-soft-clipped-bases true \
	--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
	--verbosity ERROR


%.MuTect2_stringent$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.design.bed
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
		printf '\t$(JAVA) $(JAVA_FLAGS) -jar $(GATK4) Mutect2 $(GATK4_MUTECT2_STRINGENT_FLAGS_SHARED) -R $(GENOME) -I $< -tumor %s %s -O %s\n' "$$SAMPLE" "$$L_OPT" "$@.split/$$c.vcf"; \
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
		if [ -n "$(GATK4_MUTECT2_STRINGENT_FILTER)" ]; then \
			$(JAVA) $(JAVA_FLAGS) -jar $(GATK4) FilterMutectCalls \
				-R $(GENOME) \
				-V $@.tmp.unfiltered.vcf \
				-O $@.tmp.FilterMutectCalls.vcf; \
		else \
			mv $@.tmp.unfiltered.vcf $@.tmp.FilterMutectCalls.vcf; \
		fi; \
		grep "^##" $@.tmp.FilterMutectCalls.vcf | sed s/ID=TLOD,Number=A/ID=TLOD,Number=./gi > $@.tmp.FilterMutectCalls.TLOD.vcf; \
		grep "^##" -v $@.tmp.FilterMutectCalls.vcf | cut -f1-10 >> $@.tmp.FilterMutectCalls.TLOD.vcf; \
		$(JAVA) -jar $(PICARD) SortVcf -I $@.tmp.FilterMutectCalls.TLOD.vcf -O $@.tmp.FilterMutectCalls.TLOD.sorted.vcf -SD $(DICT); \
		CMD="bcftools view -e 'GT=\"0/0\"' --threads=$(THREADS_BY_CALLER) $@.tmp.FilterMutectCalls.TLOD.sorted.vcf"; \
		if [ -n "$(GATK4_MUTECT2_STRINGENT_EXCLUDE)" ]; then \
			CMD="$$CMD | bcftools view -e '$(GATK4_MUTECT2_STRINGENT_EXCLUDE)' --threads=$(THREADS_BY_CALLER)"; \
		fi; \
		if [ -n "$(GATK4_MUTECT2_STRINGENT_INCLUDE)" ]; then \
			CMD="$$CMD | bcftools view -i '$(GATK4_MUTECT2_STRINGENT_INCLUDE)' --threads=$(THREADS_BY_CALLER)"; \
		fi; \
		if [ -n "$(GATK4_MUTECT2_STRINGENT_FILTER)" ]; then \
			CMD="$$CMD | bcftools view -f '$(GATK4_MUTECT2_STRINGENT_FILTER)' --threads=$(THREADS_BY_CALLER)"; \
		fi; \
		CMD="$$CMD > $@" ; \
		echo $$CMD; \
		eval $$CMD; \
	fi
	# Clean
	rm -rf $@.split $@.tmp* $@.idx


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING: MUTECT2 with GATK4 in parallel to identify somatic variants and generate *.MuTect2_filtered.vcf files, with MuTect2 and filtration on DP and VAF. Parameters GATK4_MUTECT2_FILTERED_FLAGS_SHARED='$(GATK4_MUTECT2_FILTERED_FLAGS_SHARED)', GATK4_MUTECT2_FILTERED_EXCLUDE='$(GATK4_MUTECT2_FILTERED_EXCLUDE)', GATK4_MUTECT2_FILTERED_INCLUDE='$(GATK4_MUTECT2_FILTERED_INCLUDE)', GATK4_MUTECT2_FILTERED_FILTER='$(GATK4_MUTECT2_FILTERED_FILTER)' "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:MuTect2_filtered:MuTect2 with GATK4 - by default with MuTect2 and DP filtration. Parameters GATK4_MUTECT2_FILTERED_FLAGS_SHARED='$(GATK4_MUTECT2_FILTERED_FLAGS_SHARED)', GATK4_MUTECT2_FILTERED_EXCLUDE='$(GATK4_MUTECT2_FILTERED_EXCLUDE)', GATK4_MUTECT2_FILTERED_INCLUDE='$(GATK4_MUTECT2_FILTERED_INCLUDE)', GATK4_MUTECT2_FILTERED_FILTER='$(GATK4_MUTECT2_FILTERED_FILTER)' "
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
