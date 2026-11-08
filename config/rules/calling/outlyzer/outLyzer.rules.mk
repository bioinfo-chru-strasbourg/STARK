############################
# OutLyzer Calling Rules
# Release: 1.0.0
# Date: 03/02/2023
# Author: Antony Le Bechec
############################

# Release note
# 1.0.0-03/02/2023: Initial release of OutLyzer calling rules.


############
# OUTLYZER #
############

# OUTLYZER options:
#   -cut CUT              defines into how many parts bed file is divide [3]
#   -t T                  Student t value used in modified Thompson tau technique [0.001]
#   -bal BAL              minimum Forward / Reverse read proportion [0.3]
#   -Q Q                  minimum average Phred Score to be considered as a real mutation (only relevant for SNP) [20]
#   -SDQ SDQ              maximum Standard deviation authorized for average Phred Score [7]
#   -WS WS                Window Size: region (number of bp) around the mutation on which background noise have to be determined [200]
#   -WSmin WSMIN          Window Size Minimum Size: minimum region size (number of bp) required for analysis [10]
#   -x X                  Multiplicative factor that specifies how often the mutation must be above background noise [2]
#   -AS                   Analysis sensitivity: Returns an additional file containing analysis average sensitivity for each line of bed file
#   -FRcor                Forward Reverse Correction: take into account any imbalance in the Forward-Reverse reads distribution in the Forward / Reverse alternative Read Proportion (-bal option)
#   -force FORCE          Force noise background determination to go below a fixed proportion of the Depth [default = disabled; 0 -> 1 ex: "-force 0.2"]
#   -HSM HSM              HotSpot Metrics: Produce sensitivity Threshold for HotSpot positions, in an additional file. Requires formated HotSpot File in argument (see documentation for more details).

# OUTLYZER Flags
# This variable can be defined directly on Apps
OUTLYZER_FLAGS?=-FRcor

# OUTLYZER filters Flags
# This variable can be defined directly on Apps
OUTLYZER_FILTERS_FLAGS?= \
	--exclude-genotypes '' \
	--exclude '' \
	--include 'INFO/DP>=10' \
	--filter '' \
	--sort


# MUTECT2 target rule
%.outLyzer$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.design.bed %.bam.chroms
	# List of vcf and auxiliary files to be generated per chromosome
	cut -f1 $*.bam.chroms | sed 's|.*|$@.split.&.vcf.gz|' | tr '\n' ' ' > $@.split.vcfs.list;
	$(BCFTOOLS) view $*.empty.vcf | grep "^#CHROM" | cut -f10 > $@.split.sample.txt;
	# Generate bed per chromosome
	awk '{print > ("$@.split."$$1".bed")}' $*.design.bed
	# Generate the sub Makefile: one explicit rule per chromosome (empty "all" if no chromosome)
	{ \
	printf 'all: %s\n' "$$(cat $@.split.vcfs.list)"; \
	for c in $$(cut -f1 $*.bam.chroms); do \
		printf '%s:\n' "$@.split.$$c.vcf.gz"; \
		printf '\tmkdir -p %s\n' "$@.split.$$c.vcf.gz.out"; \
		if [ -s $@.split.$$c.bed ]; then \
			L_ARG="-bed $@.split.$$c.bed"; \
			B_ARG="-bam $<"; \
		else \
			printf '\t$(SAMTOOLS) view %s -b -1 -o %s\n' "$<" "$@.split.$$c.vcf.gz.out/$$c.bam"; \
			L_ARG=""; \
			B_ARG="-bam $@.split.$$c.vcf.gz.out/$$c.bam"; \
		fi; \
		printf '\t$(PYTHON) $(OUTLYZER) calling -pythonPath=$(PYTHON) -samtools=$(SAMTOOLS) $(OUTLYZER_FLAGS) %s %s -ref $(GENOME) -output %s/ -core 1 -verbose 0\n' "$$L_ARG" "$$B_ARG" "$@.split.$$c.vcf.gz.out"; \
		printf '\tfind %s -name "*.vcf" | xargs -I {} $(BCFTOOLS) reheader -s %s {} -o %s\n' "$@.split.$$c.vcf.gz.out" "$@.split.sample.txt" "$@.split.$$c.vcf.gz.out.vcf"; \
		printf '\t$(JAVA) $(JAVA_FLAGS) -jar $(GATK) UpdateVCFSequenceDictionary --sequence-dictionary $(DICT) --replace True --variant %s --output %s\n' "$@.split.$$c.vcf.gz.out.vcf" "$@.split.$$c.vcf.gz"; \
	done; \
	} > $@.split.sub.mk
	# Add rule to concat per-chromosome VCFs into the final VCF (empty VCF if no chromosome)
	{ \
		printf '%s: %s\n' "$@.tmp.unfiltered.vcf" "$$(cat $@.split.vcfs.list) $*.empty.vcf"; \
		printf '\t$(BCFTOOLS) concat %s -Ov --threads $(THREADS_BY_CALLER) | awk -f $(STARK_FOLDER_BIN)/outlyzer_norm.awk > %s\n' "$$(cat $@.split.vcfs.list) $*.empty.vcf" "$@.tmp.unfiltered.vcf"; \
	} >> $@.split.sub.mk;
	# Calling all chromosomes in parallel (nothing to do if no chromosome, empty VCF will be used)
	+make -f $@.split.sub.mk $@.tmp.unfiltered.vcf
	# Filters
	bash $(STARK_FOLDER_BIN)/filter_vcf.sh --input $@.tmp.unfiltered.vcf --output $@ --threads $(THREADS_BY_CALLER) --bcftools $(BCFTOOLS) $(OUTLYZER_FILTERS_FLAGS)
	# Cleaning
	rm -rf $@.tmp* $@.split* $@.idx $@.tbi


# OUTLYZER Release documentation

RELEASE_COMMENT := "\#\# CALLING OutLyzer '$(MK_RELEASE)': OutLyzer tool identify variants from aligned BAM"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING: OutLyzer in parallel to identify somatic variants with low VAF and generate *.outLyzer.vcf files. Parameters OUTLYZER_FLAGS='$(OUTLYZER_FLAGS)' and OUTLYZER_FILTERS_FLAGS='$(OUTLYZER_FILTERS_FLAGS)' "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:outLyzer:OutLyzer - OutLyzer caller, GT info added, filter on DP:OUTLYZER_FLAGS='$(OUTLYZER_FLAGS)', OUTLYZER_FILTERS_FLAGS='$(OUTLYZER_FILTERS_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
