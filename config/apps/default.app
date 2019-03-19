#!/bin/bash
## STARK application DEFAULT

# DEFAULT ENV
######################
#source_app $CONFIG_HEADER

# APPLICATION INFOS
#####################
APP_NAME="DEFAULT"
APP_RELEASE="1.1"
APP_DESCRIPTION="Default application"
APP_GROUP=""
APP_PROJECT=""


# FOLDERS
###########

# TOOLS FOLDER
# tools: All tools needed for STARK, and more, including STARK
FOLDER_TOOLS=/STARK/tools

# DATABASES FOLDER
# genomes: All references genomes. format: $FOLDER_GENOMES/$ASSEMBLY/$ASSEMBLY.fa (with ASSEMBLY=hg19, hg38, mmu19...)
# db: Folder with mandatory databases: dbSNP database for variant calling, such as "dbsnp_138.hg19.vcf.gz" (mandatory, depending on ASSEMBLY), VCF databases for recalibration such as "dbsnp_137.hg19.vcf" (mandatory, depending on ASSEMBLY)
FOLDER_DATABASES=/STARK/databases

# RUN FOLDER
# Illumina Sequencer repository Folder. Subfolder as runs
FOLDER_RUN=/STARK/runs

# MANIFEST FOLDER
# Illumina Manifests repository.
# Files to provide in the SampleSheet of each run
FOLDER_MANIFEST=/STARK/manifests

# RESULTS FOLDER
# All results will be generated in this folder :
# RES: RUN files such as BAM, VCF, metrics
# DEM: Demultiplexing folder
# LOG: log files
# TMP: temporary files
FOLDER_RESULTS=/STARK/results

# REPOSITORY FOLDER
# Results data can be copy in a repository folder. leave it blank for no copy
FOLDER_REPOSITORY=/STARK/repository


# PARAMETERS
######################


# Rules
# Add specific rules to load.
# These files will be added to the list of rules files from APPS folder
# example: RULES_APP="MYGROUP/*.rules.mk" "$APP_FOLDER/*.rules.mk"
# for inheritance, use RULES_APP="$RULES_APP MYGROUP/*.rules.mk"
#RULES_APP=""


# ASSEMBLY (default hg19)
# Assembly is automatically detected in manifest file (configured in the SampleSheet of each run), if any
ASSEMBLY=hg19


# PIPELINES (default "bwamem.gatkHC.howard")
# pipelines to use for the analysis
# Format: "ALIGNER1.CALLER1.ANNOTATOR ALIGNER1.CALLER2.ANNOTATOR1 ALIGNER2.CALLER1.ANNOTATOR1"
# If variables ALIGNERS, CALLERS, and ANNOTATORS are defined (see below), PIPELINES variable will be automatically generated
# This variable can be a additionnal pipeline that those defined by the combinasion of ALIGNERS, CALLERS, and ANNOTATORS
# If no pipeline is finally defined, the default pipeline will be applied
PIPELINES="bwamem.gatkHC.howard"

# ALIGNERS (default "")
# Aligners to use for the analysis
# Example of available aligners: bwamem  bwasw bwaaln
#ALIGNERS="bwamem"

# CALLERS (default "")
# Callers to use for the analysis
# Example of available callers: gatkHC gatkUG VarScan samtools
#CALLERS="gatkHC"

# ANNOTATORS (default "")
# Annotator to use for the analysis
# Example of available annotators: howard snpeff
#ANNOTATORS="howard"

# BLANK samples
# Used to reject for CNV analysis
BLANK="BlcADN,blanc,BlcPCR,blcPCR,T_NTC,Z_NTC"

# BARCODE_MISMATCHES (default 1)
# Used to demultiplex
BARCODE_MISMATCHES=1

# BAM METRICS (default 1/TRUE/YES/Y)
# Performs BAM METRICS (1/TRUE/YES/Y or 0/FALSE/NO/N). Time and space consuming. Switch off for exome/genome for better perfomances
BAM_METRICS=0

# VARIANT RECALIBRATION (default 0/FALSE/NO/N)
# Performs Variant recalibration after Calling (1/TRUE/YES/Y or 0/FALSE/NO/N).
# If the recalibration fail (usually due to lack of data for statistic calculation), nothing will be done
VARIANT_RECALIBRATION=0

# INTERVAL_PADDING (default 0)
# Add some “padding” to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=0

# COVERAGE CRITERIA (default "1,30")
# For gene coverage metrics
# the criteria to calculate the percent of bases over $COVERAGE_CRITERIA X (eg 30 for 30X)
#COVERAGE_CRITERIA="1,30"
COVERAGE_CRITERIA="1,30,100"

# COVERAGE DP THRESHOLD (default "30" "100" "1")
# For gene coverage metrics
# the criteria test if genes failed (or just warning) the coverage threshold
DP_FAIL="30"	# fail DP threshold (default 30X)
DP_WARN="100"	# warn DP threshold (default 100X)
DP_THRESHOLD="1" # threshold percentage of bases over the DP threshold 

# NB_BASES_AROUND (default 0)
# For gene coverage metrics
# the number of bases to look around the exons from the given bed file
NB_BASES_AROUND=0

# BEDFILE_GENES (default "")
# For gene coverage metrics
# the bed file containing the 5'UTR, 3'UTR and genomic coding coordinates.
BEDFILE_GENES=""

# VARANK ANALYSIS (default 0)
# Performs VARANK ANALYSIS with Alamut (1 (for true) or 0 (for false))
VARANK_ANALYSIS=0
# if yes, define a folder to store the results
# VARANK_FOLDER=$FOLDER_RESULTS/VARANK

# BAM CHECK (default 0)
# Check BAM for each manipulation step (clipping, realignment...)
# Time consuming, but will stop the analysis in case of missing reads.
# A Metrics on each BAM check the BAM decrependy in any case
BAM_CHECK_STEPS=0

# METRICS SNPEFF (default 0)
# Generate snpEff variant metrics from VCF
METRICS_SNPEFF=0

# PIPELINES PRIORITIZATION
# List of pipelines to prioritize for the report (final.vcf)
PRIORITIZE_PIPELINES_LIST=""


# POST ALIGNEMENT STEPS (default "sorting realignment clipping compress")
# All steps after alignement and before calling
# This sequence correspond to the BAM file generated jsut after the alignemnt
# Format: "step1 step2 step3"
# Example: "sorting realignment clipping compress"
#    This sequence will generate the file $ALIGNER.compress.clipping.realigned.sorting.bam whose will be processed
#    Then, this BAM file will be 1/ sorted, 2/ realigned, 3/ clipped and 4/ compressed
# The steps are defined as makefiles rules
# Check available steps by using the command: STARK --pipelines_infos
# Available steps (not up-to-date:
#    sorting: BAM sorting
#    compress: BAM compression (see $BAM_COMPRESSION variable)
#    realignment: local realignment
#    markduplicates: BAM Mark Duplicates
#    clipping: BAM Clipping according to primer definition in manifest file, if any
# Usually:
#    "sorting realignment clipping compress" for Amplicon technology
#    "sorting markduplicates realignment compress" for Capture technology
#POST_ALIGNMENT_STEPS="sorting realignment recalibration clipping compress"
POST_ALIGNMENT_STEPS="sorting markduplicates realignment recalibration compress"

# BAM COMPRESSION
# Final BAM copression level (unaligned.bam, ALIGNER.bam)
BAM_COMPRESSION=5


# THREADS (default AUTO)
# Number of threads to use for the analysis
# AUTO will considere CORE-1 threads to use
# The number of threads need to be between 1 and the total number of cores available (autoadjusting if bad value)
THREADS=AUTO


# HOWARD ANNOTATION/PRIOTITIZATION/TRANSLATION CONFIGURATION

# HOWARD Configuration files for Annotation and Prioritization
# Use APP_FOLDER if necessary
# Example: HOWARD_CONFIG_ANNOTATION=$APP_FOLDER/config.annotation.stark.ini
# Example: HOWARD_CONFIG_ANNOTATION=$STARK_FOLDER_APPS/MY_APP_GROUP/config.annotation.stark.ini
#HOWARD_CONFIG_ANNOTATION=$APP_FOLDER/config.annotation.stark.ini
#HOWARD_CONFIG_PRIORITIZATION=$APP_FOLDER/config.prioritization.stark.ini


# ANNOTATION
# Default annotation with HOWARD for intermediate VCF (for each caller) used by default with annotation rule "howard"
#ANNOTATION_TYPE="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs" "core,symbol,location,outcome,hgvs,snpeff,snpeff_hgvs,snpeff_split"
HOWARD_ANNOTATION="core,snpeff_split"
# Default annotation with HOWARD for minimal VCF annotation (rule howard_minimal)
HOWARD_ANNOTATION_MINIMAL="core,snpeff_split"
# Default annotation with HOWARD for report
HOWARD_ANNOTATION_REPORT="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs,snpeff_split"
#HOWARD_ANNOTATION_REPORT=""


# CALCULATION
# Default calculation with HOWARD for all VCF/pipelines
HOWARD_CALCULATION="VARTYPE,NOMEN"
# Default minimal calculation with HOWARD for final VCF report
HOWARD_CALCULATION_MINIMAL="VARTYPE,NOMEN"
# Default calculation with HOWARD for final VCF report
HOWARD_CALCULATION_REPORT="FindByPipelines,GenotypeConcordance,VAF,VAF_STATS,CALLING_QUALITY,CALLING_QUALITY_EXPLODE,VARTYPE,NOMEN"


# PRIORITIZATION
# Default filter to prioritize/rank variant.
# This option create ranking scores in VCF and comment in TXT (after translation).
# Scores can be used to sort variant in the TXT
# HOWARD_FILTER_DEFAULT="default" # in env_header.sh
# Default calculation with HOWARD
HOWARD_PRIORITIZATION=$HOWARD_PRIORITIZATION_DEFAULT # "default"
# Minimal calculation with HOWARD
HOWARD_PRIORITIZATION_MINIMAL=$HOWARD_PRIORITIZATION_DEFAULT # "default"
# Default calculation with HOWARD for Report (full/final VCF)
HOWARD_PRIORITIZATION_REPORT=$HOWARD_PRIORITIZATION_DEFAULT # "default"


# TRANSLATION
# List of fields to show in the TSV file
# use ALL to show ALL "other" annotations
# Default filter to prioritize/rank variant, Sort variant in the TXT using 2 fields, Order fields in variant ranking
HOWARD_FIELDS="NOMEN,PZFlag,PZScore,PZComment,CNOMEN,PNOMEN,location,outcome,VAF_average,dbSNP,dbSNPNonFlagged,popfreq,snpeff_impact,ALL"
HOWARD_SORT_BY="PZFlag,PZScore"
HOWARD_ORDER_BY="DESC,DESC"
