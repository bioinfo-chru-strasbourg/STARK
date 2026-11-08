#!/bin/bash
## STARK application LONG_INDELS
# Detect long indels in gene panel sequencing data
# Use bwamem for alignment and GATK HaplotypeCaller, GATK UnifiedGenotyper, or MuTect2 for variant calling
# Disable MuTect2 to avoid calling somatic variants in this analysis
# To combine callers in pipelines, create a specific rule for long indels
# Example: gatkUG_LONG_INDELS would be a rule callinglong indels with GATK UnifiedGenotyper
# Example of pipeline variable: PIPELINES="bwamem.gatkUG.howard bwamem.gatkUG_LONG_INDELS.howard"

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP

# APPLICATION INFOS
#####################
APP_NAME="LONG_INDELS"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect long indels in gene panel sequencing data"
APP_GROUP="GENETIC"
APP_PROJECT="LONG_INDELS"

# ANALYSIS PARAMETERS
#######################

# PIPELINES
PIPELINES="bwamem.gatkHC.howard bwamem.gatkUG.howard bwamem.MuTect2.howard"

# INTERVAL_PADDING / add some padding to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=100


# Callers parameters
######################

# GATKHC Flags
# VCFDBSNP_WITH_GATK variable depends on USE_VCFDBSNP_WITH_GATK setting on Apps
# INTERVAL_PADDING variable defines the padding around intervals for all callers, settings on Apps
# This variable can be defined directly on Apps
export GATKHC_FLAGS=" \
    --interval-padding 100 \
    --min-base-quality-score 17 \
    --min-pruning 4 \
    --max-reads-per-alignment-start 1000 \
    --max-assembly-region-size 500 \
    --max-num-haplotypes-in-population 256 \
    --min-dangling-branch-length 2 \
    --standard-min-confidence-threshold-for-calling 30 \
    --dont-use-soft-clipped-bases true \
    --recover-dangling-heads true"

# GATKHC filters Flags
# This variable can be defined directly on Apps
export GATKHC_FILTERS_FLAGS=" \
    --exclude-genotypes '' \
    --exclude '' \
    --include 'FORMAT/DP>=4 && TYPE=\"indel\" && abs(strlen(REF)-strlen(ALT)) > 10' \
    --filter '' \
    --sort"

# GATKUG Flags
# VCFDBSNP_WITH_GATK variable depends on USE_VCFDBSNP_WITH_GATK setting on Apps
# INTERVAL_PADDING variable defines the padding around intervals for all callers, settings on Apps
# This variable can be defined directly on Apps
export GATKUG_FLAGS=" \
    --interval_padding 100 \
    --genotype_likelihoods_model INDEL \
    --read_filter BadCigar \
    --baq OFF \
    --downsampling_type NONE \
    --min_indel_fraction_per_sample 0.01 \
    --min_indel_count_for_genotyping 2 \
    --max_deletion_fraction 0.01 \
    --standard_min_confidence_threshold_for_calling 10 \
    --downsample_to_fraction 1 \
    --min_base_quality_score 17 \
    --max_alternate_alleles 6 \
    --allow_potentially_misencoded_quality_scores"

# GATKUG filters Flags
# This variable can be defined directly on Apps
export GATKUG_FILTERS_FLAGS=" \
    --exclude-genotypes '0/0' \
    --exclude '' \
    --include 'FORMAT/DP>=4 && TYPE=\"indel\" && abs(strlen(REF)-strlen(ALT)) > 10' \
    --filter '' \
    --sort"

# MUTECT2 Flags
# This variable can be defined directly on Apps
export MUTECT2_FLAGS=" \
	--interval-padding 100 \
	--min-base-quality-score 17 \
	--max-reads-per-alignment-start 1000 \
    --max-assembly-region-size 500 \
	--min-pruning 20 \
	--callable-depth 30 \
    --min-dangling-branch-length 2 \
	--dont-use-soft-clipped-bases true \
	--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
	--verbosity ERROR"

# MUTECT2 filters Flags
# This variable can be defined directly on Apps
export MUTECT2_FILTERS_FLAGS=" \
	--exclude-genotypes '0/0' \
	--exclude 'FILTER=\"FAIL\"' \
	--include 'FORMAT/AF>=0.01 && FORMAT/AF<0.80 && TYPE=\"indel\" && abs(strlen(REF)-strlen(ALT)) > 10' \
	--filter 'PASS' \
	--sort"

