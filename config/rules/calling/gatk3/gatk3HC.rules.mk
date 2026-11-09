############################
# GATK Calling Rules
# Release: 1.1.0
# Date: 03/11/2026
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
# 0.9.4-03/02/2023: Extract gatkHC
# 0.9.5-31/10/2025: Change GATKHC to GATK3HC
# 1.1.0-03/11/2026: Normalize variant calling parallelization and parametrization (use %.bam.chroms rule, use filter_vcf.sh for filtering)


########################
# GATK HaplotypeCaller #
########################

# Call germline SNPs and indels via local re-assembly of haplotypes
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php

# The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other. It also makes the HaplotypeCaller much better at calling indels than position-based callers like UnifiedGenotyper.
# In the so-called GVCF mode used for scalable variant calling in DNA sequence data, HaplotypeCaller runs per-sample to generate an intermediate genomic gVCF (gVCF), which can then be used for joint genotyping of multiple samples in a very efficient way, which enables rapid incremental processing of samples as they roll off the sequencer, as well as scaling to very large cohort sizes (e.g. the 92K exomes of ExAC).
# In addition, HaplotypeCaller is able to handle non-diploid organisms as well as pooled experiment data. Note however that the #algorithms used to calculate variant likelihoods is not well suited to extreme allele frequencies (relative to ploidy) so its #use is not recommended for somatic (cancer) variant discovery. For that purpose, use MuTect2 instead.
# Finally, HaplotypeCaller is also able to correctly handle the splice junctions that make RNAseq a challenge for most variant callers.

# GATK3HC Flags
# VCFDBSNP_WITH_GATK variable depends on USE_VCFDBSNP_WITH_GATK setting on Apps
# INTERVAL_PADDING variable defines the padding around intervals for all callers, settings on Apps
# This variable can be defined directly on Apps
GATK3HC_FLAGS?= \
	$(VCFDBSNP_WITH_GATK) \
	--interval_padding $(INTERVAL_PADDING) \
	--baq OFF \
	--read_filter BadCigar \
	--allow_potentially_misencoded_quality_scores \
	--dontUseSoftClippedBases \
	--standard_min_confidence_threshold_for_calling 10 \
	--downsample_to_fraction 1 \
	--maxReadsInRegionPerSample 250 \
	--min_base_quality_score 17 \
	--minPruning 4

# GATK3HC filters Flags
# This variable can be defined directly on Apps
GATK3HC_FILTERS_FLAGS?= \
	--exclude-genotypes '' \
	--exclude '' \
	--include '' \
	--filter '' \
	--sort


# GATK3HC target rule
%.gatk3HC$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.design.bed %.bam.chroms
	# List of vcf files to be generated per chromosome
	cut -f1 $*.bam.chroms | sed 's|.*|$@.split.&.vcf.gz|' | tr '\n' ' ' > $@.split.vcfs.list;
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
		printf '\t$(JAVA8) $(JAVA_FLAGS) -jar $(GATK3) -T HaplotypeCaller $(GATK3HC_FLAGS) -L %s -R $(GENOME) -I $< -o %s --num_cpu_threads_per_data_thread 1\n' "$$L_ARG" "$@.split.$$c.vcf.gz"; \
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
	bash $(STARK_FOLDER_BIN)/filter_vcf.sh --input $@.tmp.unfiltered.vcf --output $@ --threads $(THREADS_BY_CALLER) --bcftools $(BCFTOOLS) $(GATK3HC_FILTERS_FLAGS)
	# Cleaning
	rm -rf $@.tmp* $@.split* $@.idx $@.tbi


# GATK3HC Release documentation

RELEASE_COMMENT := "\#\# CALLING GATK3 Haplotype Caller: GATK3 Haplotype Caller tool identify variants from aligned BAM"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATK3HC identify variants and generate *.gatk3HC.vcf files with parameters: GATK3HC_FLAGS='$(GATK3HC_FLAGS)' and GATK3HC_FILTERS_FLAGS='$(GATK3HC_FILTERS_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:gatk3HC:GATK Haplotype Caller - by default:GATK3HC_FLAGS='$(GATK3HC_FLAGS)' and GATK3HC_FILTERS_FLAGS='$(GATK3HC_FILTERS_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
