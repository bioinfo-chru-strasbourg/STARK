#!/bin/bash
## STARK application ONCOGENET

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP

# APPLICATION INFOS
#####################
APP_NAME="ONCOGENET"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect germline mutations in oncology gene panel sequencing data"
APP_GROUP="GENETIC"
APP_PROJECT="ONCOGENETIC"

# ANALYSIS PARAMETERS
#######################

# PIPELINES
PIPELINES="bwamem.gatkHC.howard bwamem.gatkUG.howard"


# Callers parameters
######################

# GATKHC Flags
# VCFDBSNP_WITH_GATK variable depends on USE_VCFDBSNP_WITH_GATK setting on Apps
# INTERVAL_PADDING variable defines the padding around intervals for all callers, settings on Apps
# This variable can be defined directly on Apps
export GATKHC_FLAGS=" \
    --interval-padding 100 \
    --min-base-quality-score 17 \
    --min-pruning 20 \
    --max-reads-per-alignment-start 8000 \
    --standard-min-confidence-threshold-for-calling 10 \
    --dont-use-soft-clipped-bases true \
    --recover-dangling-heads true"

# GATKHC filters Flags
# This variable can be defined directly on Apps
export GATKHC_FILTERS_FLAGS=" \
    --exclude-genotypes '' \
    --exclude '' \
    --include 'FORMAT/DP>=4' \
    --filter '' \
    --sort"

# GATKUG Flags
# VCFDBSNP_WITH_GATK variable depends on USE_VCFDBSNP_WITH_GATK setting on Apps
# INTERVAL_PADDING variable defines the padding around intervals for all callers, settings on Apps
# This variable can be defined directly on Apps
export GATKUG_FLAGS=" \
    --interval_padding 100 \
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
    --allow_potentially_misencoded_quality_scores"

# GATKUG filters Flags
# This variable can be defined directly on Apps
export GATKUG_FILTERS_FLAGS=" \
    --exclude-genotypes '0/0' \
    --exclude '' \
    --include 'FORMAT/DP>=4' \
    --filter '' \
    --sort"
