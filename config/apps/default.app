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

# MAIN STARK FOLDER
STARK_FOLDER_MAIN="/STARK"

# TOOLS FOLDER
# tools: All tools needed for STARK, and more, including STARK
FOLDER_TOOLS=$STARK_FOLDER_MAIN/tools

# DATABASES FOLDER
# genomes: All references genomes. format: $FOLDER_GENOMES/$ASSEMBLY/$ASSEMBLY.fa (with ASSEMBLY=hg19, hg38, mmu19...)
# db: Folder with mandatory databases: dbSNP database for variant calling, such as "dbsnp_138.hg19.vcf.gz" (mandatory, depending on ASSEMBLY), VCF databases for recalibration such as "dbsnp_137.hg19.vcf" (mandatory, depending on ASSEMBLY)
FOLDER_DATABASES=$STARK_FOLDER_MAIN/databases


# INPUT FOLDER
FOLDER_INPUT=$STARK_FOLDER_MAIN/input
# Illumina Sequencer repository Folder. Subfolder as runs
#FOLDER_RUN=$FOLDER_INPUT/runs
# Illumina Manifests repository.
# Files to provide in the SampleSheet of each run
#FOLDER_MANIFEST=$FOLDER_INPUT/manifests
# Pedigree repository.
#FOLDER_PEDIGREE=$FOLDER_INPUT/pedigree


# OUTPUT FOLDER
# All results will be generated in this folder :
FOLDER_OUTPUT=$STARK_FOLDER_MAIN/output
# RES: RUN files such as BAM, VCF, metrics
#FOLDER_RESULTS=$FOLDER_OUTPUT/results
# DEM: Demultiplexing folder
#FOLDER_DEMULTIPLEXING=$FOLDER_OUTPUT/demulitplexing
# LOG: log files
#FOLDER_LOG=$FOLDER_OUTPUT/log
# TMP: temporary files
#FOLDER_TMP=$FOLDER_OUTPUT/tmp

# REPOSITORY and ARCHIVES folder
# Results data can be copy in a repository folder. leave it blank for no copy
FOLDER_REPOSITORY=$FOLDER_OUTPUT/repository
# Results data can be copy in a archives folder. leave it blank for no copy
FOLDER_ARCHIVES=$FOLDER_OUTPUT/archives
# Results data can be copy in a favorites folder. leave it blank for no copy
# Configurations:
#    - to NOT use favorites folder: FOLDER_FAVORITES=
#    - to copy favorites within repository folder: FOLDER_FAVORITES=$FOLDER_REPOSITORY
#    - to copy favorites within default favorites folder: FOLDER_FAVORITES=$FOLDER_OUTPUT/favorites
FOLDER_FAVORITES=


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

# GLOBAL METRICS VARIABLES
# Values for BAM metrics
# Minimum mapping quality to consider in the metrics BAM
METRICS_MINIMUM_MAPPING_QUALITY=10

# Minimum bases quality to consider in the metrics BAM
METRICS_MINIMUM_BASE_QUALITY=10

# Clipping overlapping reads in the metrics BAM
CLIP_OVERLAPPING_READS=1

# Flagged reads in the metrics BAM (mpileup format)
# default UNMAP,SECONDARY,QCFAIL,DUP
METRICS_FLAGS="UNMAP,SECONDARY,QCFAIL,DUP"

# Flagged reads in the metrics BAM (samtools format)
# Generated from mpileup format (if empty)
#SAMTOOLS_METRICS_FLAG_PARAM=" -F 0x4 -F 0x100 -F 0x200 -F 0x400"


# INTERVAL_PADDING (default 0)
# Add some “padding” to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=0


# COVERAGE CRITERIA (default "1,30")
# For gene coverage metrics
# the criteria to calculate the percent of bases over $COVERAGE_CRITERIA X (eg 30 for 30X)
#COVERAGE_CRITERIA="1,30"
COVERAGE_CRITERIA="1,5,10,20,30,50,100,200,300"

# COVERAGE DP THRESHOLD (default "30" "100" "1")
# For gene coverage metrics
# the criteria test if genes failed (or just warning) the coverage threshold
SEQUENCING_DEPTH="1" # Sequencing depth threshold
SEQUENCING_COVERAGE_THRESHOLD="1" # Sequencing coverage threshold
MINIMUM_DEPTH="30" # fail DP threshold (default 30X)
EXPECTED_DEPTH="100" # warn DP threshold (default 100X)
DEPTH_COVERAGE_THRESHOLD="0.95" # threshold percentage of bases over the DP threshold

# CLIP_OVERLAPPING_READS (default 1)
# From PICARD: For paired reads, soft clip the 3' end of each read if necessary so that it does not extend past the 5' end of its mate
CLIP_OVERLAPPING_READS=1

# NB_BASES_AROUND (default 0)
# For gene coverage metrics
# the number of bases to look around the exons from the given bed file
NB_BASES_AROUND=0

# GENESCOVERAGE_PRECISION (default 2)
# Genes Coverage calculation precision
GENESCOVERAGE_PRECISION=2

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
# Only for Report final VCF
METRICS_SNPEFF=1

# PIPELINES PRIORITIZATION
# List of pipelines to prioritize for the report (final.vcf)
PRIORITIZE_PIPELINES_LIST=""



# FASTQ Processing

# Set mask for demultiplexing
# e.g. "" (auto from SampleSheet), "Y150,I10,Y10,Y150", "Y150,I8,Y10,Y150" (UMI index2)
# default ""
STARK_DEMULTIPLEXING_BASES_MASK=""

# Set short read size for demultiplexing
# If demultiplexing UMI within a read must be set to 0 for Agilent XTHS kits
# default "" (auto)
STARK_DEMULTIPLEXING_MASK_SHORT_ADAPTATER_READ=""

# Set read mapping
# Redefine FASTQ files in order to identify R1, R2, I1 and I2
# Order: R1 I1 I2 R2
# e.g. "" (default, corresponding to "R1 I1 I2 R2"), "R1 I1 R2 R3" (UMI index2)
# default "R1 I1 I2 R2"
STARK_DEMULTIPLEXING_READS_MAPPING=""

# Demultiplexing adaptated stringency
# For BCL2FASTQ demultiplexing (see doc)
ADAPTER_STRINGENCY=0.9

# Demultiplexing options
# For BCL2FASTQ demultiplexing (see doc)
# Usually: "--no-lane-splitting --create-fastq-for-index-reads"
# Default: ""
STARK_DEMULTIPLEXING_BCL2FASTQ_OPTIONS="--no-lane-splitting --create-fastq-for-index-reads"

# FASTQ compression level for demultiplexing FASTQ files
# zlib compression level (1-9) used for FASTQ files during demultiplexing
# Used by BCL2FASTQ
# If FASTQ_DEMULTIPLEXING_KEEP=1, we suggest a high level of compression (at least 5)
FASTQ_DEMULTIPLEXING_COMPRESSION_LEVEL=1

# FASTQ compression level for main FASTQ files
# zlib compression level (1-9) used for FASTQ files
# Used by FASTP
FASTQ_COMPRESSION_LEVEL=1

# ENABLE_ADAPTER_TRIMMING
# Trim adapter and autodetect adapter for paired end
# Either 0 or 1
# Default: 0 (i.e. adapter trimming is disabled)
ENABLE_ADAPTER_TRIMMING=0

# FASTQ Read quality filtering
# Read Quality threshold. Read quality below will be removed
# Default: 0 (i.e. disable)
FASTQ_QUALITY_FILTERING=0

# polyG tail trimming
# force polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
# the minimum length to detect polyG in the read tail.
# Default: 0 (i.e. disable)
POLY_G_MIN_LEN=0

# Read length filtering
# reads shorter than length_required will be discarded
# Default: 0 (i.e. disable)
READ_LENGTH_REQUIRED=0

# UMI extract location
# Set the UMI location
# If not null, NO UMI extraction and analysis
# Available locations: 
#    index1: the first index is used as UMI. If the data is PE, this UMI will be used for both read1/read2.
#    index2 the second index is used as UMI. PE data only, this UMI will be used for both read1/read2.
#    read1 the head of read1 is used as UMI. If the data is PE, this UMI will be used for both read1/read2.
#    read2 the head of read2 is used as UMI. PE data only, this UMI will be used for both read1/read2.
#    per_index read1 will use UMI extracted from index1, read2 will use UMI extracted from index2.
#    per_read read1 will use UMI extracted from the head of read1, read2 will use UMI extracted from the head of read2.
# e.g.: UMI_LOC="index2"
# See FASTP/UMI TOOLS documentation for more information
UMI_LOC=""

# UMI extract tag
# Set the UMI Barcode pattern
# If not null, STARK will prepare fastq containg UMIs +/- cell barcodes for alignment
# e.g.: UMI_BARCODE_PATTERN="NNNNNNNNNN" for simplex
# e.g.: UMI_BARCODE_PATTERN="NNNNN-NNNNN" for duplex
# if UMI_LOC is "per_index" or "per_read", and UMI_BARCODE_PATTERN is defined as "NNNNN", it with be redefined as "NNNNN-NNNNN"
# Only length of the first part of the duplex barcode will be used with FASTP
# See FASTP/UMI TOOLS documentation for more information
UMI_BARCODE_PATTERN=""

# Barcode tag
# Barcode to use for Mark Duplicates
# If not null, Mark Duplicates will consider this tag (default null)
# e.g.: BARCODE_TAG="BC" for 10X Genomics, BARCODE_TAG="BX" for UMI
# See PICARD documentation for more information
BARCODE_TAG=""

# set variable to "READ_NAME_REGEX=null" to disable optical deduplication
# PICARD_MARKDUP_OPTICAL_DEDUP="-READ_NAME_REGEX null"
PICARD_MARKDUP_OPTICAL_DEDUP=""

# Keep demultiplexing FASTQ
# Keep fastq demultiplexed or from input reads/reads2
FASTQ_DEMULTIPLEXING_KEEP=0

# Folder for demultiplexing FASTQ
# within $SAMPLE.sequencing folder
SEQUENCING_DEMULTIPLEXING_FOLDER=demultiplexing


# FASTQ_PROCESSING_STEPS
# All steps to process input FASTQ files, after sequencing and demultiplexing (if any)
# Format: "step1 step2 step3"
# Example (default): fastq_reheader sort fastp fastq_clean_header compress
# Example (UMItools): fastq_reheader sort umi_tools fastp fastq_clean_header compress
# Available steps:
#    fastq_reheader: FASTQ reheader to integreate index within FASTQ comment Illumina tag (e.g. 1:N:0:xxx). Nothing done if already integrated (same header or tag BC or RX exists)
#    fastq_clean_header: FASTQ read head formatting, especially SAMTOOLS tags. Nothing done if no needs
#    compress: FASTQ files compression (see FASTQ_COMPRESSION_LEVEL)
#    sort: sort FASTQ using read name
#    fastp: process FASTP algorithm and report, UMI extraction (if any, see UMI_LOC and UMI_BARCODE_PATTERN), quality filtration...
#    umi_tools: process UMITools algorithm for UMI extraction (if any, see UMI_LOC and UMI_BARCODE_PATTERN)
# Usually:
#    "fastq_reheader sort fastp fastq_clean_header compress" for UMI technology
# dafault:
#    "sort compress" for sorting and compression

#FASTQ_PROCESSING_STEPS="fastq_reheader sort umi_tools fastp fastq_clean_header compress"
FASTQ_PROCESSING_STEPS="fastq_reheader sort fastp fastq_clean_header compress"



# POST SEQUENCING STEPS (default '')
# All steps and before alignment
# This sequence correspond to the FASTQ file processing before the alignemnt (trimming...)
# Format: "step1 step2 step3"
# The steps are defined as makefiles rules
# Check available steps by using the command: STARK --pipelines_infos
# Available steps (not up-to-date):
# Usually:
#    "" nothing to do
POST_SEQUENCING_STEPS=""



# POST ALIGNEMENT STEPS (default "sorting realignment clipping compress")
# All steps after alignement and before calling
# This sequence correspond to the BAM file generated jsut after the alignemnt
# Format: "step1 step2 step3"
# Example: "sorting realignment clipping compress"
#    This sequence will generate the file $ALIGNER.compress.clipping.realigned.sorting.bam whose will be processed
#    Then, this BAM file will be 1/ sorted, 2/ realigned, 3/ clipped and 4/ compressed
# The steps are defined as makefiles rules
# Check available steps by using the command: STARK --pipelines_infos
# Available steps (not up-to-date):
#    sorting: BAM sorting
#    compress: BAM compression (see $BAM_COMPRESSION variable)
#    realignment: local realignment
#    recalibration: reads recalibration
#    gencore: gencore is a tool for fast and powerful deduplication for paired-end next-generation sequencing
#    markduplicates: Mark duplicated reads in BAM with PICARD MarkDuplicates. Use BARCODE_TAG to specify tag
#    clipping: BAM Clipping according to primer definition in manifest file, if any
# Usually:
#    "sorting realignment clipping compress" for Amplicon technology
#    "sorting markduplicates realignment compress" for Capture technology
#    "sorting gencore markduplicates realignment compress" for UMI technology
#POST_ALIGNMENT_STEPS="sorting realignment recalibration clipping compress"
POST_ALIGNMENT_STEPS="sorting markduplicates realignment recalibration compress"



# POST CALLING STEPS (default " ")
# All steps after calling
# This sequence correspond to the VCF file generated just after the calling
# Format: "step1 step2 step3"
# Example: "recalibration filtration"
#    This sequence will generate the file $CALLER.filtration.recalibration.vcf whose will be processed
#    Then, this VCF file will be 1/ recalibrated, 2/ filtrered
# The steps are defined as makefiles rules
# Check available steps by using the command: STARK --pipelines_infos
# Available steps (not up-to-date):
#    sorting: VCF sort
#    normalization: VCF normalization
#    variantrecalibration: VCF recalibration (using GATK4). Include variantfiltration if no recalibration possible
#    variantfiltration: VCF filtration (using GATK4)
# Usually:
#    " " to avoid at this pipeline step (see POST_CALLING_MERGING_STEPS)
#    "normalization variantfiltration" for gene panel
#    "normalization variantrecalibration" for exome or genome
#POST_CALLING_STEPS="normalization variantfiltration"
POST_CALLING_STEPS=" "



# POST CALLING MERGING STEPS (default "sorting normalization variantrecalibration")
# All steps after merging calling
# This sequence correspond to the VCF file generated after the merge of VCF calling
# Format: "step1 step2 step3"
# Example: "recalibration filtration"
#    This sequence will generate the file $CALLER.filtration.recalibration.vcf whose will be processed
#    Then, this VCF file will be 1/ recalibrated, 2/ filtrered
# The steps are defined as makefiles rules
# Check available steps by using the command: STARK --pipelines_infos
# Available steps (not up-to-date):
#    sorting: VCF sort
#    normalization: VCF normalization
#    variantrecalibration: VCF recalibration (using GATK4). Include variantfiltration if no recalibration possible
#    variantfiltration: VCF filtration (using GATK4)
# Usually:
#    "sorting normalization variantrecalibration" for exome or genome
#    "sorting normalization variantfiltration" for gene panel
POST_CALLING_MERGING_STEPS="sorting normalization variantrecalibration"



# POST ANNOTATION STEPS (default " ")
# All steps after annotation
# This sequence correspond to the VCF file generated just after the annotation
# Format: "step1 step2 step3"
# Example: "sorting normalization"
#    This sequence will generate the file $ANNOTATION.sorting.vcf whose will be processed
#    Then, this VCF file will be 1/ sorted, 2/ normalized
# The steps are defined as makefiles rules
# Check available steps by using the command: STARK --pipelines_infos
# Available steps (not up-to-date):
#    sorting: VCF sorting
# Usually:
#    " " to avoid at this pipeline step (see POST_CALLING_MERGING_STEPS)
#    "sorting"
POST_ANNOTATION_STEPS=" "



# BAM COMPRESSION
# Final BAM compression level (ALIGNER.bam)
BAM_COMPRESSION=9


# BAM VALIDATION COMPRESSION
# Validation BAM compression level (ALIGNER.validation.bam)
BAM_VALIDATION_COMPRESSION=5


### GENCORE

# supporting_reads ; set to 2 to keep the clusters with 2 or more supporting reads / ultrasensible filter
# set to 1 to replace Picard Markduplicates dedup
# default "" (corresponding to "--supporting_reads 1", see gencore doc)
GENCORE_SUP_READS=""

# --score_threshold
# set to 8 recommanded for dup-rate < 50% if you want to keep all the DNA fragments, and for each output read you want to discard all the low quality unoverlapped mutations to obtain a relative clean data
# default "" (corresponfing to "--score_threshold 6", see gencore doc)
GENCORE_SCORE_THREESHOLD=""

# --ratio_threshold
# if the ratio of the major base in a cluster is less than <ratio_threshold>, it will be further compared to the reference.
# The valud should be 0.5~1.0, and the default value is 0.8
# default "" (corresponding to "--ratio_threshold 0.8", see gencore doc)
GENCORE_RATIO_THREESHOLD=""

# --umi_diff_threshold
# if two reads with identical mapping position have UMI difference <= <umi_diff_threshold>, then they will be merged to generate a consensus read
# Default "" (corresponding to "--umi_diff_threshold 2", see gencore doc)
GENCORE_DIFF_THREESHOLD=""

# Quality Threeshold
# --high_qual (Q30) ; --moderate_qual (Q20) ; --low_qual (Q15)
# --high_qual : the threshold for a quality score to be considered as high quality. Default 30 means Q30. (int [=30])
# --moderate_qual : the threshold for a quality score to be considered as moderate quality. Default 20 means Q20. (int [=20])
# --low_qual : the threshold for a quality score to be considered as low quality. Default 15 means Q15. (int [=15])
# e.g. "--moderate_qual 20", "--high_qual 20 --moderate_qual 15 --low_qual 10"
# default "" (corresponding to default quality, see gencore doc)
GENCORE_QUAL_THREESHOLD=""

# --coverage_sampling
# the sampling rate for genome scale coverage statistics. Default 10000 means 1/10000
# for statistics purpose, can be reduce to gain performance
# not included in the mk file
# default "" (corresponding to "--coverage_sampling=10000", see gencore doc)
GENCORE_COVERAGE_SAMPLING=""


### CRAM

# CRAM OPTIONS
# Final CRAM options for compression (archive.cram)
# example: CRAM_OPTIONS="version=3.0,level=9,no_ref,use_lzma,use_bzip2,use_fqz,seqs_per_slice=100000"
CRAM_OPTIONS="version=3.0,level=9,no_ref"


# CRAM REMOVE TAGS
# Final CRAM options for tags (archive.cram)
# example: CRAM_REMOVE_TAGS="BD,BI,OQ"
CRAM_REMOVE_TAGS="BD,BI"



# THREADS (default AUTO)
# Number of threads to use for the analysis
# AUTO will considere CORE-1 threads to use
# The number of threads need to be between 1 and the total number of cores available (autoadjusting if bad value)
THREADS=AUTO

# THREADS_LOADING (default THREADS)
# Number of threads used for loading demultiplexed data
# AUTO will considere THREADS threads to use
# The number of threads need to be between 1 and the total number of cores available (autoadjusting if bad value)
#THREADS_LOADING=

# THREADS_WRITING (default THREADS)
# Number of threads used for writing demultiplexed data
# AUTO will considere THREADS threads to use
# The number of threads need to be between 1 and the total number of cores available (autoadjusting if bad value)
#THREADS_WRITING=

# THREADS_COPY (default 1)
# Number of threads used for copy files in repositories
# The number of threads need to be between 1 and the total number of cores available (autoadjusting if bad value)
#THREADS_COPY=


# MEMORY (default AUTO)

# Total memory to use by thread
# default: MemTotal from /proc/meminfo divided by number of threads
#MEMORY=

# Maximum memory to use for Java
# default: MEMORY
#JAVA_MEMORY=



# HOWARD ANNOTATION/PRIOTITIZATION/TRANSLATION CONFIGURATION


# DATABASES FOLDER for HOWARD/ANNOVAR/SNPEFF/DEJAVU
# Define ANNOVAR and SNPEFF and DEJAVU folders by default
#FOLDER_DATABASES_ANNOVAR=$FOLDER_DATABASES/annovar
#FOLDER_DATABASES_SNPEFF=$FOLDER_DATABASES/snpeff/current
#FOLDER_DATABASES_DEJAVU_ANNOVAR=$FOLDER_DATABASES/annovar
#FOLDER_DATABASES_DEJAVU_ANNOVAR=$FOLDER_DATABASES/dejavu/latest/annovar


# HOWARD Configuration files for Annotation and Prioritization
# Use APP_FOLDER if necessary
# default: $HOWARD_FOLDER_CONFIG/config.annotation.ini and $HOWARD_FOLDER_CONFIG/config.prioritization.ini
# Example: HOWARD_CONFIG_ANNOTATION=$APP_FOLDER/config.annotation.ini
# Example: HOWARD_CONFIG_PRIORITIZATION=$STARK_FOLDER_APPS/MY_APP_GROUP/config.prioritization.ini
#HOWARD_CONFIG_ANNOTATION=$HOWARD_FOLDER_CONFIG/config.annotation.ini
#HOWARD_CONFIG_PRIORITIZATION=$HOWARD_FOLDER_CONFIG/config.prioritization.ini

# HOWARD DEJAVU Configuration files for Annotation
# Use APP_FOLDER if necessary
# default: $HOWARD_FOLDER_CONFIG/config.annotation.ini
# Example: HOWARD_CONFIG_DEJAVU_ANNOTATION=$APP_FOLDER/config.annotation.ini
# Example: HOWARD_CONFIG_DEJAVU_ANNOTATION=$STARK_FOLDER_APPS/MY_APP_GROUP/config.annotation.ini
#HOWARD_CONFIG_DEJAVU_ANNOTATION=$HOWARD_FOLDER_CONFIG/config.annotation.ini


# ANNOTATION
# Default annotation with HOWARD for intermediate VCF (for each caller) used by default with annotation rule "howard"
#ANNOTATION_TYPE="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs" "core,symbol,location,outcome,hgvs,snpeff,snpeff_hgvs,snpeff_split"
HOWARD_ANNOTATION="symbol,location,outcome,hgvs"
# Default annotation with HOWARD for minimal VCF annotation (rule howard_minimal)
HOWARD_ANNOTATION_MINIMAL="core,snpeff_split"
# Default annotation with HOWARD for report
HOWARD_ANNOTATION_REPORT="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs,snpeff_split"
# Default annotation with HOWARD for whole analysis
HOWARD_ANNOTATION_ANALYSIS="null" # no more annotation


# CALCULATION
# Default calculation with HOWARD for all VCF/pipelines
HOWARD_CALCULATION="VAF_STATS,DP_STATS,VARTYPE,NOMEN,BARCODE"
# Default minimal calculation with HOWARD for final VCF report
HOWARD_CALCULATION_MINIMAL="VAF_STATS,DP_STATS,VARTYPE,NOMEN,BARCODE"
# Default calculation with HOWARD for final VCF report
HOWARD_CALCULATION_REPORT="FindByPipelines,GenotypeConcordance,VAF,VAF_STATS,DP_STATS,VARTYPE,NOMEN,BARCODE"
# Default calculation with HOWARD for whole analysis (calculation forced, no transcripts list available)
HOWARD_CALCULATION_ANALYSIS="VAF_STATS,DP_STATS,BARCODE"
# List of annotation fields to extract NOMEN annotation (default 'hgvs', see HOWARD docs)
HOWARD_NOMEN_FIELDS="hgvs"


# PRIORITIZATION
# Default filter to prioritize/rank variant.
# This option create ranking scores in VCF and comment in TXT (after translation).
# Scores can be used to sort variant in the TXT
# HOWARD_FILTER_DEFAULT="default" # in env_header.sh
# Default prioritization with HOWARD
HOWARD_PRIORITIZATION=$HOWARD_PRIORITIZATION_DEFAULT # "default"
# Minimal prioritization with HOWARD
HOWARD_PRIORITIZATION_MINIMAL=$HOWARD_PRIORITIZATION_DEFAULT # "default"
# Default prioritization with HOWARD for Report (full/final VCF)
HOWARD_PRIORITIZATION_REPORT=$HOWARD_PRIORITIZATION_DEFAULT # "default"
# Default prioritization with HOWARD for whole analysis (prioritization forced)
HOWARD_PRIORITIZATION_ANALYSIS="" # "none"
# Default prioritization with HOWARD for VaRank score mode
HOWARD_PRIORITIZATION_VARANK=VaRank # "default"


# TRANSLATION
# List of fields to show in the TSV file
# use ALL to show ALL "other" annotations
# Default filter to prioritize/rank variant, Sort variant in the TXT using 2 fields, Order fields in variant ranking
HOWARD_FIELDS="NOMEN,PZFlag,PZScore,PZComment,CNOMEN,PNOMEN,location,outcome,snpeff_impact,VAF_average,dbSNP,dbSNPNonFlagged,popfreq,ALL"
HOWARD_SORT="PZFlag::DESC,PZScore:n:DESC"
HOWARD_SORT_BY="PZFlag,PZScore"
HOWARD_ORDER_BY="DESC,DESC"
# Minimal
HOWARD_FIELDS_MINIMAL="NOMEN,PZFlag,PZScore,PZComment,CNOMEN,PNOMEN,location,outcome,snpeff_impact,VAF_average,dbSNP,dbSNPNonFlagged,popfreq"
HOWARD_SORT_MINIMAL="PZFlag::DESC,PZScore:n:DESC"
HOWARD_SORT_BY_MINIMAL="PZFlag,PZScore"
HOWARD_ORDER_BY_MINIMAL="DESC,DESC"
# REPORT
HOWARD_FIELDS_REPORT="NOMEN,PZFlag,PZScore,PZComment,CNOMEN,PNOMEN,location,outcome,snpeff_impact,VAF_average,dbSNP,dbSNPNonFlagged,popfreq"
HOWARD_SORT_REPORT="PZFlag::DESC,PZScore:n:DESC"
HOWARD_SORT_BY_REPORT="PZFlag,PZScore"
HOWARD_ORDER_BY_REPORT="DESC,DESC"




# GATK4 Recalibrator and Filtration

# Variant Filtration
# Filter variant calls based on INFO and/or FORMAT annotations

# Variant Filtration main option (see documentation guide for more info)
# default: VARIANTFILTRATION_OPTIONS=
#VARIANTFILTRATION_OPTIONS=

# One or more expression used with INFO fields to filter SNP (see documentation guide for more info)
# default: VARIANTFILTRATION_SNP_FILTER_OPTION=''
# example: VARIANTFILTRATION_SNP_FILTER_OPTION='--filter-name "SNP_filter_QD" --filter-expression "QD < 2.0" --filter-name "SNP_filter_FS" --filter-expression "FS > 60.0" --filter-name "SNP_filter_MQ" --filter-expression "MQ < 40.0" --filter-name "SNP_filter_MQRankSum" --filter-expression "MQRankSum < -12.5" --filter-name "SNP_filter_ReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0"'
# example: VARIANTFILTRATION_SNP_FILTER_OPTION='--filter-name "HARD_TO_VALIDATE" --filter-expression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)"    --filter-name "VeryVeryLowDepth" --filter-expression "DP == 0"    --filter-name "VeryLowDepth" --filter-expression "DP > 0 && DP < 10"    --filter-name "LowDepth" --filter-expression "DP >= 10 && DP < 30"    --filter-name "VeryVeryLowQual" --filter-expression "QUAL == 0"    --filter-name "VeryLowQual" --filter-expression "QUAL > 0 && QUAL < 30.0"    --filter-name "LowQual" --filter-expression "QUAL >= 30.0 && QUAL < 50.0"    --filter-name "LowQD" --filter-expression "QD >= 0.0 && QD < 1.5"'
VARIANTFILTRATION_SNP_FILTER_OPTION='--filter-name "SNP_filter_QD" --filter-expression "QD < 2.0" --filter-name "SNP_filter_FS" --filter-expression "FS > 60.0" --filter-name "SNP_filter_MQ" --filter-expression "MQ < 40.0" --filter-name "SNP_filter_MQRankSum" --filter-expression "MQRankSum < -12.5" --filter-name "SNP_filter_ReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0"'

# One or more expression used with FORMAT (sample/genotype-level) fields to filter SNP (see documentation guide for more info)
# default: VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION=''
# example: VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION='--genotype-filter-expression "GQ == 0" --genotype-filter-name "genotype_GQ_filter_VeryVeryLow" --genotype-filter-expression "GQ > 0 && GQ < 50.0" --genotype-filter-name "genotype_GQ_filter_VeryLow"  --genotype-filter-expression "GQ >= 50.0 && GQ < 90.0" --genotype-filter-name "genotype_GQ_filter_Low" --genotype-filter-expression "DP == 0" --genotype-filter-name "genotype_DP_filter_VeryVeryLow" --genotype-filter-expression "DP >= 0 && DP < 10" --genotype-filter-name "genotype_DP_filter_VeryLow"  --genotype-filter-expression "DP >= 10 && DP < 30" --genotype-filter-name "genotype_GQ_filter_Low"'
# example: VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION='--genotype-filter-name "VeryVeryLowGQ" --genotype-filter-expression "GQ == 0"    --genotype-filter-name "VeryLowGQ" --genotype-filter-expression "GQ > 0 && GQ < 50.0"    --genotype-filter-name "LowGQ" --genotype-filter-expression "GQ >= 50.0 && GQ < 90.0"    --genotype-filter-name "VeryVeryLowDP" --genotype-filter-expression "DP == 0"    --genotype-filter-name "VeryLowDP" --genotype-filter-expression "DP >= 0 && DP < 10"    --genotype-filter-name "LowDP" --genotype-filter-expression "DP >= 10 && DP < 30" '
VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION='--genotype-filter-expression "GQ == 0" --genotype-filter-name "genotype_GQ_filter_VeryVeryLow" --genotype-filter-expression "GQ > 0 && GQ < 50.0" --genotype-filter-name "genotype_GQ_filter_VeryLow"  --genotype-filter-expression "GQ >= 50.0 && GQ < 90.0" --genotype-filter-name "genotype_GQ_filter_Low" --genotype-filter-expression "DP == 0" --genotype-filter-name "genotype_DP_filter_VeryVeryLow" --genotype-filter-expression "DP >= 0 && DP < 10" --genotype-filter-name "genotype_DP_filter_VeryLow"  --genotype-filter-expression "DP >= 10 && DP < 30" --genotype-filter-name "genotype_GQ_filter_Low"'

# One or more expression used with INFO fields to filter INDEL (see documentation guide for more info)
# default: VARIANTFILTRATION_INDEL_FILTER_OPTION=''
# example: VARIANTFILTRATION_INDEL_FILTER_OPTION='--filter-name "INDEL_filter_QD" --filter-expression "QD < 2.0" --filter-name "INDEL_filter_FS" --filter-expression "FS > 200.0" --filter-name "INDEL_filter_ReadPosRankSum" --filter-expression "ReadPosRankSum < -20.0"'
# example: VARIANTFILTRATION_INDEL_FILTER_OPTION='--filter-name "HARD_TO_VALIDATE" --filter-expression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)"    --filter-name "VeryVeryLowDepth" --filter-expression "DP == 0"    --filter-name "VeryLowDepth" --filter-expression "DP > 0 && DP < 10"    --filter-name "LowDepth" --filter-expression "DP >= 10 && DP < 30"    --filter-name "VeryVeryLowQual" --filter-expression "QUAL == 0"    --filter-name "VeryLowQual" --filter-expression "QUAL > 0 && QUAL < 30.0"    --filter-name "LowQual" --filter-expression "QUAL >= 30.0 && QUAL < 50.0"    --filter-name "LowQD" --filter-expression "QD >= 0.0 && QD < 1.5"'
VARIANTFILTRATION_INDEL_FILTER_OPTION='--filter-name "INDEL_filter_QD" --filter-expression "QD < 2.0" --filter-name "INDEL_filter_FS" --filter-expression "FS > 200.0" --filter-name "INDEL_filter_ReadPosRankSum" --filter-expression "ReadPosRankSum < -20.0"'

# One or more expression used with FORMAT (sample/genotype-level) fields to filter INDEL (see documentation guide for more info)
# default: VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION=''
# example: VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION='--genotype-filter-expression "GQ == 0" --genotype-filter-name "genotype_GQ_filter_VeryVeryLow" --genotype-filter-expression "GQ > 0 && GQ < 50.0" --genotype-filter-name "genotype_GQ_filter_VeryLow"  --genotype-filter-expression "GQ >= 50.0 && GQ < 90.0" --genotype-filter-name "genotype_GQ_filter_Low" --genotype-filter-expression "DP == 0" --genotype-filter-name "genotype_DP_filter_VeryVeryLow" --genotype-filter-expression "DP >= 0 && DP < 10" --genotype-filter-name "genotype_DP_filter_VeryLow"  --genotype-filter-expression "DP >= 10 && DP < 30" --genotype-filter-name "genotype_GQ_filter_Low"'
# example: VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION='--genotype-filter-name "VeryVeryLowGQ" --genotype-filter-expression "GQ == 0"    --genotype-filter-name "VeryLowGQ" --genotype-filter-expression "GQ > 0 && GQ < 50.0"    --genotype-filter-name "LowGQ" --genotype-filter-expression "GQ >= 50.0 && GQ < 90.0"    --genotype-filter-name "VeryVeryLowDP" --genotype-filter-expression "DP == 0"    --genotype-filter-name "VeryLowDP" --genotype-filter-expression "DP >= 0 && DP < 10"    --genotype-filter-name "LowDP" --genotype-filter-expression "DP >= 10 && DP < 30" '
VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION='--genotype-filter-expression "GQ == 0" --genotype-filter-name "genotype_GQ_filter_VeryVeryLow" --genotype-filter-expression "GQ > 0 && GQ < 50.0" --genotype-filter-name "genotype_GQ_filter_VeryLow"  --genotype-filter-expression "GQ >= 50.0 && GQ < 90.0" --genotype-filter-name "genotype_GQ_filter_Low" --genotype-filter-expression "DP == 0" --genotype-filter-name "genotype_DP_filter_VeryVeryLow" --genotype-filter-expression "DP >= 0 && DP < 10" --genotype-filter-name "genotype_DP_filter_VeryLow"  --genotype-filter-expression "DP >= 10 && DP < 30" --genotype-filter-name "genotype_GQ_filter_Low"'

# Remove previous filters applied to the VCF
# within makefile rule: VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS_OPTION?=$(shell if (( $(VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS) )); then echo " --invalidate-previous-filters "; fi )
# default: VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS=0
VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS=1


# Variant Recalibrator

# Variant Recalibrator main option (see documentation guide for more info)
# default: VARIANTRECALIBRATOR_OPTIONS=
#VARIANTRECALIBRATOR_OPTIONS=

# Variant Recalibrator SNP resources option (see documentation guide for more info)
# These resources need to be available on STARK Databases folder for GATK
# default:
# VARIANTRECALIBRATION_SNP_RESOURCES="
#   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.vcf.gz
#   -resource:omni,known=false,training=true,truth=true,prior=12.0 1000G_omni2.5.b37.vcf.gz
#   -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.b37.vcf.gz 
#   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.b37.vcf.gz
# "
VARIANTRECALIBRATION_SNP_RESOURCES="
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.vcf.gz
    -resource:omni,known=false,training=true,truth=true,prior=12.0 1000G_omni2.5.b37.vcf.gz
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.b37.vcf.gz 
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.b37.vcf.gz
"

# Variant Recalibrator INDEL resources option (see documentation guide for more info)
# These resources need to be available on STARK Databases folder for GATK
# default:
# VARIANTRECALIBRATION_INDEL_RESOURCES="
#   -resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.b37.vcf.gz
#   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.b37.vcf.gz
# "
VARIANTRECALIBRATION_INDEL_RESOURCES="
    -resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.b37.vcf.gz
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.b37.vcf.gz
"

# Variant Recalibrator SNP annotations option (see documentation guide for more info)
# default: VARIANTRECALIBRATION_SNP_ANNOTATIONS="-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP"
VARIANTRECALIBRATION_SNP_ANNOTATIONS="-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP"

# Variant Recalibrator INDEL annotations option (see documentation guide for more info)
# default: VARIANTRECALIBRATION_INDEL_ANNOTATIONS="-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum"
VARIANTRECALIBRATION_INDEL_ANNOTATIONS="-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum"

# Variant Recalibrator SNP tranches option (see documentation guide for more info)
# default: VARIANTRECALIBRATION_SNP_TRANCHES="-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.1 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0"
VARIANTRECALIBRATION_SNP_TRANCHES="-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.1 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0"

# Variant Recalibrator INDEL tranches option (see documentation guide for more info)
# default: VARIANTRECALIBRATION_INDEL_TRANCHES="-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.1 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0"
VARIANTRECALIBRATION_INDEL_TRANCHES="-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.1 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0"


# Variant Recalibrator Optional Variant Filtration 
# If Variant Recalibrator failed, usually due to lack of variant in the input callset, a Variant Filtration is optional
# default: no filtration
# example: Use empty value to switch off this option (no filtration)
# VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_OPTION=
# VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION=
# VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_OPTION=
# VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION=
# example: Use same options than Variant Filtration in stand alone (see above)
# VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_OPTION=$VARIANTFILTRATION_SNP_FILTER_OPTION
# VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION=$VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION
# VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_OPTION=$VARIANTFILTRATION_INDEL_FILTER_OPTION
# VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION=$VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION
# Filter variant calls based on INFO and/or FORMAT annotations for SNP and INDEL 
VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_OPTION=$VARIANTFILTRATION_SNP_FILTER_OPTION
VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION=$VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION
VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_OPTION=$VARIANTFILTRATION_INDEL_FILTER_OPTION
VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION=$VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION



# REPORT
# Report variables

# Report Sections
# List of sections to show in the report (default "ALL")
# Sections :
#   Report Sections: results_summary sequencing_mapping depth coverage variant_calling variant_stats
#   Report Annex Sections: annex_coverage annex_depth annex_genes_coverage annex_variants annex_annotations
REPORT_SECTIONS="ALL"


### REPORT for run Files
# Generate variants files from run with full VCF (include all calling information)
REPORT_VARIANTS_FULL=0


# REPOSITORY and ARCHIVES
# Copy files into folders for repository and archives
# Repository files patterns to add
REPOSITORY_FILE_PATTERNS=' $SAMPLE.*.validation.bam $SAMPLE.*.validation.bam.bai $SAMPLE.*.bam.metrics/$SAMPLE.*.validation.flags.Design.bed $SAMPLE.reports/$SAMPLE.full.Design.vcf.gz $SAMPLE.reports/$SAMPLE.full.Design.vcf.gz.tbi $SAMPLE.reports/$SAMPLE.full.Design.tsv '
# Repository files patterns to exclude on results subfolder (usually 'STARK'). Useful to reduce storage and exclude repeated files (in others files patterns varaibles)
REPOSITORY_FILE_SUBFOLDER_PATTERNS=' $SAMPLE.*fastq.gz $SAMPLE.*.validation.bam $SAMPLE.*.validation.bam.bai '
# Archives files patterns to add
ARCHIVES_FILE_PATTERNS=' $SAMPLE.reports/$SAMPLE.full.vcf.gz $SAMPLE.reports/$SAMPLE.full.vcf.gz.tbi $SAMPLE.reports/$SAMPLE.final.tsv $SAMPLE.*.bam.metrics/$SAMPLE.*.validation.flags.*.bed '
# Favorites files patterns to add
FAVORITES_FILE_PATTERNS=''


# IGV SESSION
# IGV session xml file (*igv_session.xml) is generated with Samples and Runs files (BAM, VCF, BED)
# Parameters select files through patterns and folder depth
# Be aware to coordinate these parameters with availabled files within repositories folders

# SAMPLE
IGV_SESSION_PATTERNS_SAMPLE_DEFAULT_APP=' $SAMPLE.final.vcf.gz $SAMPLE*Design*bed $SAMPLE*Panel*bed *Design*vcf.gz *Panel*vcf.gz *bam *cram '
IGV_SESSION_MINDEPTH_SAMPLE_DEFAULT_APP=0
IGV_SESSION_MAXDEPTH_SAMPLE_DEFAULT_APP=1

REPOSITORY_SAMPLE_IGV_SESSION_PATTERNS=$IGV_SESSION_PATTERNS_SAMPLE_DEFAULT_APP
REPOSITORY_SAMPLE_IGV_SESSION_MINDEPTH=$IGV_SESSION_MINDEPTH_SAMPLE_DEFAULT_APP
REPOSITORY_SAMPLE_IGV_SESSION_MAXDEPTH=$IGV_SESSION_MAXDEPTH_SAMPLE_DEFAULT_APP


ARCHIVES_SAMPLE_IGV_SESSION_PATTERNS=$IGV_SESSION_PATTERNS_SAMPLE_DEFAULT_APP
ARCHIVES_SAMPLE_IGV_SESSION_MINDEPTH=$IGV_SESSION_MINDEPTH_SAMPLE_DEFAULT_APP
ARCHIVES_SAMPLE_IGV_SESSION_MAXDEPTH=$IGV_SESSION_MAXDEPTH_SAMPLE_DEFAULT_APP

FAVORITES_SAMPLE_IGV_SESSION_PATTERNS=$IGV_SESSION_PATTERNS_SAMPLE_DEFAULT_APP
FAVORITES_SAMPLE_IGV_SESSION_MINDEPTH=$IGV_SESSION_MINDEPTH_SAMPLE_DEFAULT_APP
FAVORITES_SAMPLE_IGV_SESSION_MAXDEPTH=$IGV_SESSION_MAXDEPTH_SAMPLE_DEFAULT_APP

# RUN
IGV_SESSION_PATTERNS_RUN_DEFAULT_APP=' *.final.vcf.gz *Design*bed *Panel*bed *Design*vcf.gz *Panel*vcf.gz *bam *cram '
IGV_SESSION_MINDEPTH_RUN_DEFAULT_APP=0
IGV_SESSION_MAXDEPTH_RUN_DEFAULT_APP=2

REPOSITORY_RUN_IGV_SESSION_PATTERNS=$IGV_SESSION_PATTERNS_RUN_DEFAULT_APP
REPOSITORY_RUN_IGV_SESSION_MINDEPTH=$IGV_SESSION_MINDEPTH_RUN_DEFAULT_APP
REPOSITORY_RUN_IGV_SESSION_MAXDEPTH=$IGV_SESSION_MAXDEPTH_RUN_DEFAULT_APP

ARCHIVES_RUN_IGV_SESSION_PATTERNS=$IGV_SESSION_PATTERNS_RUN_DEFAULT_APP
ARCHIVES_RUN_IGV_SESSION_MINDEPTH=$IGV_SESSION_MINDEPTH_RUN_DEFAULT_APP
ARCHIVES_RUN_IGV_SESSION_MAXDEPTH=$IGV_SESSION_MAXDEPTH_RUN_DEFAULT_APP

FAVORITES_RUN_IGV_SESSION_PATTERNS=$IGV_SESSION_PATTERNS_RUN_DEFAULT_APP
FAVORITES_RUN_IGV_SESSION_MINDEPTH=$IGV_SESSION_MINDEPTH_RUN_DEFAULT_APP
FAVORITES_RUN_IGV_SESSION_MAXDEPTH=$IGV_SESSION_MAXDEPTH_RUN_DEFAULT_APP


# IGV SESSION DB
# Additionnal databases for IGV session (see IGV doc)
# Example :
# IGV_SESSION_RESSOURCES='
# 	<Resource index="https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz.tbi" path="https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz" type="vcf"/>
# 	<Resource index="http://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi" path="http://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz" type="vcf"/>
# 	<Resource name="GC Percentage" path="http://www.broadinstitute.org/igvdata/annotations/hg19/hg19.gc5base.tdf" type="tdf"/>
# 	<Resource index="https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz.tbi" name="Refseq Genes" path="https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz" type="refgene"/>
# 	<Resource hyperlink="http://www.gencodegenes.org/" index="https://s3.amazonaws.com/igv.broadinstitute.org/annotations/hg19/genes/gencode.v18.annotation.sorted.gtf.gz.tbi" name="Gencode Genes (v18)" path="https://s3.amazonaws.com/igv.broadinstitute.org/annotations/hg19/genes/gencode.v18.annotation.sorted.gtf.gz" trackLine="visibilitywindow:10000000" type="gtf"/>
# '
# IGV_SESSION_DATAPANEL='
# 	<Track attributeKey="00-All.vcf.gz" clazz="org.broad.igv.variant.VariantTrack" featureVisibilityWindow="100100" fontSize="10" groupByStrand="false" id="https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz" name="dbSNP" siteColorMode="ALLELE_FREQUENCY" displayMode="COLLAPSE" visible="true"/>
# 	<Track attributeKey="clinvar.vcf.gz" clazz="org.broad.igv.variant.VariantTrack" featureVisibilityWindow="100100" fontSize="10" groupByStrand="false" id="http://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz" name="ClinVar" siteColorMode="ALLELE_FREQUENCY" displayMode="COLLAPSE" visible="true"/>
# '
# IGV_SESSION_FEATUREPANEL='
# 	<Track attributeKey="GC Percentage" autoScale="false" clazz="org.broad.igv.track.DataSourceTrack" colorScale="ContinuousColorScale;0.0;80.0;255,255,255;0,0,178" fontSize="10" height="20" id="http://www.broadinstitute.org/igvdata/annotations/hg19/hg19.gc5base.tdf" name="GC Percentage" renderer="HEATMAP" visible="true" windowFunction="mean">
# 		<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="80.0" minimum="0.0" type="LINEAR"/>
# 	</Track>
# 	<Track attributeKey="Refseq Genes" clazz="org.broad.igv.track.FeatureTrack" featureVisibilityWindow="63761233" fontSize="10" groupByStrand="false" id="https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz" name="Refseq Genes" visible="true"/>
# 	<Track attributeKey="Gencode Genes (v18)" clazz="org.broad.igv.track.FeatureTrack" featureVisibilityWindow="1000000" fontSize="10" groupByStrand="false" id="https://s3.amazonaws.com/igv.broadinstitute.org/annotations/hg19/genes/gencode.v18.annotation.sorted.gtf.gz" name="Gencode Genes (v18)" visible="true"/>
# '
# Default Refseq Gene
IGV_SESSION_RESSOURCES='<Resource index="https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz.tbi" name="Refseq Genes" path="https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz" type="refgene"/>'
IGV_SESSION_DATAPANEL=''
IGV_SESSION_FEATUREPANEL='<Track attributeKey="Refseq Genes" clazz="org.broad.igv.track.FeatureTrack" featureVisibilityWindow="63761233" fontSize="10" groupByStrand="false" id="https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz" name="Refseq Genes" visible="true"/>'


# IGV SESSION JSON
# Example:
# IGV_SESSION_JSON_TRACKS_ADDITIONAL='
# {
# 	"type": "bed",
# 	"url": "https://s3.amazonaws.com/igv.org.test/data/gencode.v18.collapsed.bed",
# 	"indexURL": "https://s3.amazonaws.com/igv.org.test/data/gencode.v18.collapsed.bed.idx",
# 	"name": "Gencode V18",
#   "type": "annotation",
#   "displayMode": "COLLAPSED"
# }'
IGV_SESSION_DAS="http://localhost:4201/static/data/public/repositories"
IGV_SESSION_DAS_REPOSITORIES="$IGV_SESSION_DAS/Repository $IGV_SESSION_DAS/Archives $IGV_SESSION_DAS/Favorites"
IGV_SESSION_JSON_TRACKS_ADDITIONAL=""

