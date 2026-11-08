#!/bin/bash
## STARK application EXOME_SOMATIC

# DEFAULT ENV
######################
source_app SOMATIC,CAPTURE

# APPLICATION INFOS
#####################
APP_NAME="EXOME_SOMATIC"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect somatic mutations in exome sequencing data"
APP_GROUP="SOMATIC"
APP_PROJECT="EXOME"

# ANALYSIS PARAMETERS
#######################

# PIPELINES
#PIPELINES="bwamem.gatkHC_EXOME_SOMATIC.howard bwamem.gatkUG_EXOME_SOMATIC.howard bwamem.outLyzer_filtered.howard bwamem.VarScan_EXOME_SOMATIC.howard bwamem.MuTect2_stringent.howard"
PIPELINES="bwamem.gatkHC.howard bwamem.gatkUG.howard bwamem.outLyzer.howard bwamem.VarScan.howard bwamem.MuTect2.howard"

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
    --standard-min-confidence-threshold-for-calling 30 \
    --dont-use-soft-clipped-bases true \
    --recover-dangling-heads true"

# GATKHC filters Flags
# This variable can be defined directly on Apps
export GATKHC_FILTERS_FLAGS=" \
    --exclude-genotypes '' \
    --exclude '' \
    --include 'FORMAT/DP>=50' \
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
    --min_indel_fraction_per_sample 0.01 \
    --min_indel_count_for_genotyping 2 \
    --max_deletion_fraction 0.01 \
    --standard_min_confidence_threshold_for_calling 30 \
    --downsample_to_fraction 1 \
    --min_base_quality_score 17 \
    --allow_potentially_misencoded_quality_scores"

# GATKUG filters Flags
# This variable can be defined directly on Apps
export GATKUG_FILTERS_FLAGS=" \
    --exclude-genotypes '0/0' \
    --exclude '' \
    --include 'FORMAT/DP>=50' \
    --filter '' \
    --sort"

# OUTLYZER Flags
# This variable can be defined directly on Apps
export OUTLYZER_FLAGS=-FRcor

# OUTLYZER filters Flags
# This variable can be defined directly on Apps
export OUTLYZER_FILTERS_FLAGS=" \
    --exclude-genotypes '' \
	--exclude '' \
	--include 'FORMAT/DP>=30 && FORMAT/VAF>=0.01' \
	--filter '' \
	--sort"


# VARSCAN SNP Flags
# This variable can be defined directly on Apps
export VARSCAN_SNP_FLAGS=" \
	--min-var-freq 0.01 \
	--min-reads2 3 \
	--min-coverage 30 \
	--p-value 1e-3 \
	--min-avg-qual 30"

# VARSCAN INDEL Flags
# This variable can be defined directly on Apps
export VARSCAN_INDEL_FLAGS=" \
	--min-var-freq 0.01 \
	--min-reads2 3 \
	--min-coverage 30 \
	--p-value 1e-3 \
	--min-avg-qual 10"

# VARSCAN filters Flags
# This variable can be defined directly on Apps
export VARSCAN_FILTERS_FLAGS=" \
	--exclude-genotypes '' \
	--exclude '' \
	--include '' \
	--filter '' \
	--sort"

# MUTECT2 Flags
# This variable can be defined directly on Apps
export MUTECT2_FLAGS=" \
	--interval-padding 100 \
	--min-base-quality-score 17 \
	--max-reads-per-alignment-start 1000 \
	--min-pruning 20 \
	--callable-depth 30 \
	--dont-use-soft-clipped-bases true \
	--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
	--verbosity ERROR"

# MUTECT2 filters Flags
# This variable can be defined directly on Apps
export MUTECT2_FILTERS_FLAGS=" \
	--exclude-genotypes '0/0' \
	--exclude 'FILTER=\"FAIL\"' \
	--include 'FORMAT/AF>=0.01 && FORMAT/AF<0.80' \
	--filter 'PASS' \
	--sort"
