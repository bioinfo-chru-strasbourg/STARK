#!/bin/bash
## STARK application RNASEQ

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP

# APPLICATION INFOS
#####################
APP_NAME="RNASEQ"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect somatic mutations in RNA-Seq data"
APP_GROUP="UNKNOWN"
APP_PROJECT="UNKNOWN"

# DEMULTIPLEXING
#######################
#FASTP/UMI options are for Takara. Comment if used with Agilent kit. 
# FASTP_ADDITIONAL_OPTIONS=" --trim_front2=6"
# UMI_LOC="index2"
# UMI_BARCODE_PATTERN="NNNNNNNN"
THREADS_LOADING=2

# ANALYSIS PARAMETERS
#######################

REF="/home1/BAS/DOCKER_STARK_MAIN_FOLDER/data/users/nicaises/rnaseq/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa"

POST_ALIGNMENT_STEPS="sorting splitncigar recalibration compress"

PIPELINES="star.STARFusion.howard star.gatkHC_SOMATIC.howard star_raw.Arriba.howard"

POST_CALLING_MERGING_STEPS="sorting"

# COVERAGE CRITERIA (default "1,30")
# For gene coverage metrics
# the criteria to calculate the percent of bases over $COVERAGE_CRITERIA X (eg 30 for 30X)
COVERAGE_CRITERIA="1,5,10,20,30,50,100,200,300,500,1000"

# COVERAGE DP THRESHOLD (default "30" "100" "1")
# For gene coverage metrics
# the criteria test if genes failed (or just warning) the coverage threshold
SEQUENCING_DEPTH="1" # Sequencing depth threshold
SEQUENCING_COVERAGE_THRESHOLD="1" # Sequencing coverage threshold
MINIMUM_DEPTH="100" # fail DP threshold (default 30X)
EXPECTED_DEPTH="300" # warn DP threshold (default 100X)
DEPTH_COVERAGE_THRESHOLD="1" # threshold percentage of bases over the DP threshold

# HOWARD ANNOTATION/PRIOTITIZATION/TRANSLATION CONFIGURATION

# ANNOTATION
# Default annotation with HOWARD for intermediate VCF (for each caller) used by default with annotation rule "howard"
HOWARD_ANNOTATION="null"
# Default annotation with HOWARD for minimal VCF annotation (rule howard_minimal)
HOWARD_ANNOTATION_MINIMAL="null"
# Default annotation with HOWARD for report
HOWARD_ANNOTATION_REPORT="null"
# Default annotation with HOWARD for whole analysis
HOWARD_ANNOTATION_ANALYSIS="null" # no more annotation

# CALCULATION
# Default calculation with HOWARD for all VCF/pipelines
HOWARD_CALCULATION="VARTYPE,NOMEN"
# Default minimal calculation with HOWARD for final VCF report
HOWARD_CALCULATION_MINIMAL="VARTYPE,NOMEN"
# Default calculation with HOWARD for final VCF report
HOWARD_CALCULATION_REPORT="FindByPipelines,VAF_STATS,DP_STATS,VARTYPE,NOMEN"

# METRICS SNPEFF (default 0)
# Generate snpEff variant metrics from VCF
METRICS_SNPEFF=0

# Report Sections
REPORT_SECTIONS="results_summary sequencing_mapping depth coverage variant_stats annex_coverage annex_depth annex_genes_coverage"