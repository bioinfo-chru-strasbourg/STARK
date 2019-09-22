#!/bin/bash
#################################
## STARK environment
#################################



# APPLICATION INFOS
############

# Name of the APP. usefull to include specific rules if exists (in $STARK/$APP_NAME.rules.mk/*rules.mk)
# AUTO detect: $(basename ${BASH_SOURCE[0]} | sed 's/^env//' | sed 's/\.sh$//gi' | cut -d. -f2-)
if [ "$APP_NAME" == "" ]; then
	APP_NAME="DEFAULT"
fi;
export APP_NAME

# Release of the APP
# AUTO detect: $(basename ${BASH_SOURCE[0]} | sed 's/^env//' | sed 's/\.sh$//gi' | cut -d. -f2-)
if [ "$APP_RELEASE" == "" ]; then
	APP_RELEASE=""
fi;
export APP_RELEASE

# GROUP and PROJECT Associated with the APP
# Use to structure data in the repository folder
if [ "$GROUP" == "" ]; then
	GROUP=""
fi;
export GROUP
if [ "$PROJECT" == "" ]; then
	PROJECT=""
fi;
export PROJECT


# SCRIPT DIR
##############

# Folder of STARK ENV
#export STARK_FOLDER_CONFIG="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )"



# FOLDERS
#################
# Checking and adaptation of variables for the workflow


# MAIN FOLDERS
[ -z $STARK_FOLDER_CONFIG ] && export STARK_FOLDER_CONFIG=$STARK_FOLDER_ROOT/config
export STARK_FOLDER_ROOT=$(dirname $STARK_FOLDER_CONFIG)
#export STARK_FOLDER_CONFIG=$STARK_FOLDER_ROOT/config
export STARK_FOLDER_BIN=$STARK_FOLDER_ROOT/bin
export STARK_FOLDER_DOCS=$STARK_FOLDER_ROOT/docs
export STARK_FOLDER_TOOLBOX=$STARK_FOLDER_ROOT/toolbox
export STARK_FOLDER_APPS=$STARK_FOLDER_CONFIG/apps
export STARK_FOLDER_RULES=$STARK_FOLDER_CONFIG/rules


# FOLDER INFRASTRUCTURE

# MAIN STARK FOLDER
STARK_FOLDER_MAIN="/STARK"

# INPUT FOLDER
[ "$INPUT" != "" ] && FOLDER_INPUT=$INPUT && FOLDER_RUN="$FOLDER_INPUT/runs" && FOLDER_MANIFEST="$FOLDER_INPUT/manifests"
if [ "$FOLDER_INPUT" == "" ]; then
	FOLDER_INPUT="$STARK_FOLDER_MAIN/input"
fi;

# RUN FOLDER
if [ "$FOLDER_RUN" == "" ]; then
	FOLDER_RUN="$FOLDER_INPUT/runs"
fi;
if [ "$FOLDER_MANIFEST" == "" ]; then
	FOLDER_MANIFEST="$FOLDER_INPUT/manifests"
fi;

# OUTPUT FOLDER
[ "$OUTPUT" != "" ] && FOLDER_OUTPUT=$OUTPUT
[ "$OUTPUT" != "" ] && FOLDER_OUTPUT=$OUTPUT && FOLDER_RESULTS="$FOLDER_OUTPUT/results" && FOLDER_DEMULTIPLEXING="$FOLDER_OUTPUT/demultiplexing"  && FOLDER_LOG="$FOLDER_OUTPUT/log"  && FOLDER_TMP="$FOLDER_OUTPUT/tmp"

if [ "$FOLDER_OUTPUT" == "" ]; then
	FOLDER_OUTPUT="$STARK_FOLDER_MAIN/output"
fi;

# RESULTS FOLDER
[ "$RESULTS" != "" ] && FOLDER_RESULTS=$RESULTS
if [ "$FOLDER_RESULTS" == "" ]; then
	FOLDER_RESULTS="$FOLDER_OUTPUT/results"
fi;

# DEMULTIPLEXING FOLDER
[ "$DEMULTIPLEXING" != "" ] && FOLDER_DEMULTIPLEXING=$DEMULTIPLEXING
if [ "$FOLDER_DEMULTIPLEXING" == "" ]; then
	FOLDER_DEMULTIPLEXING="$FOLDER_OUTPUT/demultiplexing"
fi;

# ANALYSIS/LOG FOLDER
[ "$LOG" != "" ] && FOLDER_LOG=$LOG
if [ "$FOLDER_LOG" == "" ]; then
	FOLDER_LOG="$FOLDER_OUTPUT/log"
fi;

# TMP FOLDER
[ "$TMP" != "" ] && FOLDER_TMP=$TMP
if [ "$FOLDER_TMP" == "" ]; then
	FOLDER_TMP="$FOLDER_OUTPUT/tmp"
fi;

# EXPORT FOLDER REPOSITORY
export FOLDER_REPOSITORY

# EXPORT FOLDER ARCHIVE
export FOLDER_ARCHIVE


# CONFIG FOLDERS

if [ "$FOLDER_TOOLS" == "" ]; then
	FOLDER_TOOLS="$STARK_FOLDER_MAIN/tools"
fi;

[ "$DATABASES" != "" ] && [ -d $DATABASES ] && FOLDER_DATABASES=$DATABASES
if [ "$FOLDER_DATABASES" == "" ]; then
	FOLDER_DATABASES="$STARK_FOLDER_MAIN/databases"
fi;


export MISEQ_FOLDER=$FOLDER_RUN				# MISEQ RAW folder
MSR_SUBFOLDER=Data/Intensities/BaseCalls		# SUBFOLDER Illumina Sequencer Data subfolder


# MANIFEST FOLDER
export MANIFEST_FOLDER=$FOLDER_MANIFEST			# MISEQ RAW folder


# RESULTS FOLDER
#export DEMULTIPLEXING_FOLDER=$FOLDER_RESULTS/DEM	# DEMULTIPLEXING folder
export DEMULTIPLEXING_FOLDER=$FOLDER_DEMULTIPLEXING	# DEMULTIPLEXING folder
export RESULTS_FOLDER=$FOLDER_RESULTS		# RESULTS folder ALL
#export ANALYSIS_FOLDER=$FOLDER_RESULTS/ANA		# ANALYSIS folder
export ANALYSIS_FOLDER=$FOLDER_LOG		# ANALYSIS folder
#export TMP_FOLDER_TMP=$FOLDER_RESULTS/TMP		# TEMPORARY folder tmp
export TMP_FOLDER_TMP=$FOLDER_TMP		# TEMPORARY folder tmp
#export VALIDATION_FOLDER=$FOLDER_RESULTS/VAL		# VALIDATION folder
#export VALIDATION_FOLDER=$FOLDER_VALIDATION		# VALIDATION folder
export TMP_FOLDER=$FOLDER_OUTPUT			# TEMPORARY folder = MAIN V2 folder
export MAIN_FOLDER=$FOLDER_OUTPUT			# MAIN V2 Folder (RES, DEM...)

# Create directories if does not exist
mkdir -p $DEMULTIPLEXING_FOLDER $RESULTS_FOLDER $ANALYSIS_FOLDER $TMP_FOLDER_TMP $TMP_FOLDER $MAIN_FOLDER 2>/dev/null
# TOOLS FOLDER
export NGS_TOOLS=$FOLDER_TOOLS				# TOOLS
export DBFOLDER=$FOLDER_DATABASES			# DB
export NGS_GENOMES=$DBFOLDER/genomes			# GENOMES folder
export RULES=$STARK_FOLDER_RULES/*.rules.mk		# RULES
export APPS=$STARK_FOLDER_APPS				# APPS
export GENOMES=$NGS_GENOMES				# Genomes folder in NGS folder
export TMP_SYS_FOLDER=/tmp				# TEMPORARY SYSTEM folder (tmpfs)
export NGS_FOLDER=$FOLDER_TOOLS
export NGS_SCRIPTS=$STARK_FOLDER_BIN			#"$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )"


# REPOSITORY FOLDER
export RESULTS_FOLDER_BY_GROUP_PROJECT_COPY=$FOLDER_REPOSITORY	# Copy data into group and project (if defined)
export RESULTS_SUBFOLDER_DATA="DATA";				# Copy sample results in a SUBFOLDER
# Copy some sample results files in the root sample folder, if any SUBFOLDER defined
#export RESULTS_SUBFOLDER_ROOT_FILE_PATTERNS=$RESULTS_SUBFOLDER_ROOT_FILE_PATTERNS' $SAMPLE.reports/$SAMPLE.final.vcf $SAMPLE.reports/$SAMPLE.full.vcf $SAMPLE.reports/$SAMPLE.final.txt $SAMPLE.reports/$SAMPLE.full.txt *.reports/*.report.pdf *.reports/latex*/*pdf';
#export RESULTS_SUBFOLDER_ROOT_FILE_PATTERNS=$RESULTS_SUBFOLDER_ROOT_FILE_PATTERNS' $SAMPLE.reports/$SAMPLE.final.vcf $SAMPLE.reports/$SAMPLE.final.vcf.idx $SAMPLE.reports/$SAMPLE.final.vcf.gz $SAMPLE.reports/$SAMPLE.final.vcf.gz.tbi $SAMPLE.reports/$SAMPLE.full.vcf $SAMPLE.reports/$SAMPLE.full.vcf.idx $SAMPLE.reports/$SAMPLE.full.vcf.gz $SAMPLE.reports/$SAMPLE.full.vcf.gz.tbi $SAMPLE.reports/$SAMPLE.final.tsv $SAMPLE.reports/$SAMPLE.full.tsv *.reports/*.report.pdf $SAMPLE.bed';
export REPOSITORY_FILE_PATTERNS=$(echo $REPOSITORY_FILE_PATTERNS' $SAMPLE.reports/$SAMPLE.final.vcf.gz $SAMPLE.reports/$SAMPLE.final.vcf.gz.tbi $SAMPLE.reports/$SAMPLE.final.tsv $SAMPLE.reports/*.report*.html $SAMPLE.reports/*.report.html.folder:FOLDER $SAMPLE.bed' | tr "," " " | tr " " "\n" | sort -u | tr "\n" " ");
export ARCHIVE_FILE_PATTERNS=$(echo $ARCHIVE_FILE_PATTERNS' $SAMPLE.reports/$SAMPLE.final.vcf.gz $SAMPLE.reports/$SAMPLE.final.tsv $SAMPLE.reports/*.report*.html $SAMPLE.reports/*.report.html.folder:FOLDER $SAMPLE.bed $SAMPLE.manifest $SAMPLE*.genes $SAMPLE*.transcripts $SAMPLE*.tag $SAMPLE.launch.json $SAMPLE.archive.cram' | tr "," " " | tr " " "\n" | sort -u | tr "\n" " ");


#Rules defined in the app
for RA in $(echo "$RULES_APP" | tr "," " " | tr " " "\n"); do
	RULES_APP=$RULES_APP" "$(find $STARK_FOLDER_APPS/$RA -maxdepth 0 -type f 2>/dev/null)
done;

# RULES for the application
[ "$APP_NAME" != "" ] && RULES_APP="$RULES_APP $STARK_FOLDER_RULES/$APP_NAME.app/*rules.mk";

# RULES for the group
[ "$APP_GROUP" != "" ] && RULES_APP="$RULES_APP $STARK_FOLDER_RULES/$APP_GROUP/*rules.mk";

# RULES for the project
[ "$APP_PROJECT" != "" ] && RULES_APP="$RULES_APP $STARK_FOLDER_RULES/$APP_PROJECT/*rules.mk";


#export RULES="$RULES $RULES_APP"
#export RULES=$(ls -N $RULES $RULES_APP)
export RULES=$(find $RULES $RULES_APP -maxdepth 0 -type f 2>/dev/null | sort -u | tr "\n" " ")
#export RULES=$(find $RULES $RULES_APP)
#echo "RULLES: $RULES"

#SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )"
#echo $SCRIPT_DIR
#echo $APP_NAME
#echo $APP_GROUP
#echo $APP_PROJECT


# ASSEMBLY
##########

# default Assembly
if [ -z $ASSEMBLY ] || [ "$ASSEMBLY" == "" ]; then
	ASSEMBLY=hg19
fi;
export ASSEMBLY
export REF=$GENOMES/$ASSEMBLY/$ASSEMBLY.fa	# Defautl REF



# TOOLS
#########

#source $STARK_FOLDER_CONFIG/tools.app
source $CONFIG_TOOLS

# DATABASES
#############

#source $STARK_FOLDER_CONFIG/databases.app
source $CONFIG_DATABASES

# databases release
DATABASES_RELEASE_FILE=$FOLDER_DATABASES/RELEASE
[ -s $DATABASES_RELEASE_FILE ] && DATABASES_RELEASE=$(grep "^RELEASE" $DATABASES_RELEASE_FILE | cut -f2) || DATABASES_RELEASE="UNKNOWN";
export DATABASES_RELEASE


# PIPELINES
##############


# PIPELINES (default "bwamem.gatkHC.howard")
# pipelines to use for the analysis
# Format: "ALIGNER1.CALLER1.ANNOTATOR ALIGNER1.CALLER2.ANNOTATOR1 ALIGNER2.CALLER1.ANNOTATOR1"
# If variables ALIGNERS, CALLERS, and ANNOTATORS are defined (see below), PIPELINES variable will be automatically generated
# This variable can be a additionnal pipeline that those defined by the combinasion of ALIGNERS, CALLERS, and ANNOTATORS
# If no pipeline is finally defined, the default pipeline will be applied
#PIPELINES=""

# ALIGNERS (default "")
# Aligners to use for the analysis
# Example of available aligners: bwamem  bwasw bwaaln
#ALIGNERS=""
export ALIGNERS

# CALLERS (default "")
# Callers to use for the analysis
# Example of available callers: gatkHC gatkUG VarScan samtools
#CALLERS=""
export CALLERS

# ANNOTATORS (default "")
# Annotator to use for the analysis
# Example of available annotators: howard snpeff
#ANNOTATORS=""
export ANNOTATORS

# CREATE PIPELINES from ALIGNERS/CALLERS/ANNOTATORS
for ALIGNER in $ALIGNERS; do
	for CALLER in $CALLERS; do
		for ANNOTATOR in $ANNOTATORS; do
			PIPELINES="$PIPELINES $ALIGNER.$CALLER.$ANNOTATOR"
		done;
	done;
done;
#PIPELINES=$(echo $PIPELINES | tr "," "\n" | tr " " "\n" | sort | uniq | tr "\n" " " | sed -e 's/^[[:space:]]*//'  | sed -e 's/[[:space:]]*$//')
PIPELINES=$(echo $PIPELINES | tr "," "\n" | tr " " "\n" | awk '!x[$0]++' | tr "\n" " " | sed -e 's/^[[:space:]]*//'  | sed -e 's/[[:space:]]*$//')
if [ -z "$PIPELINES" ]; then
	PIPELINES="bwamem.gatkHC.howard"
fi;
export PIPELINES

# MANIFEST
export MANIFEST

# BLANK samples
export BLANK

# BARCODE_MISMATCHES (default 1)
# Used to demultiplex
if [ -z $BARCODE_MISMATCHES ]; then
	BARCODE_MISMATCHES=1
fi;
export BARCODE_MISMATCHES

# BAM METRICS (default 1/TRUE/YES/Y)
# Performs BAM METRICS (1/TRUE/YES/Y or 0/FALSE/NO/N). Time and space consuming. Switch off for exome/genome for better perfomances
if [ "${BAM_METRICS^^}" == "FALSE" ] || [ "${BAM_METRICS^^}" == "NO" ] || [ "${BAM_METRICS^^}" == "N" ] || [ "$BAM_METRICS" == "0" ]; then
	BAM_METRICS=0
elif [ -z $BAM_METRICS ] || [ "${BAM_METRICS^^}" == "TRUE" ] || [ "${BAM_METRICS^^}" == "YES" ] || [ "${BAM_METRICS^^}" == "Y" ] || [ "$BAM_METRICS" == "1" ]; then
	BAM_METRICS=1
else
	BAM_METRICS=1
fi;
export BAM_METRICS


# GLOBAL METRICS VARIABLES
# Values for BAM metrics
# Minimum mapping quality to consider in the metrics BAM
if [ -z "$METRICS_MINIMUM_MAPPING_QUALITY" ] || ! [[ $METRICS_MINIMUM_MAPPING_QUALITY =~ ^[0-9]+$ ]]; then
	METRICS_MINIMUM_MAPPING_QUALITY=20
fi;
export METRICS_MINIMUM_MAPPING_QUALITY

# Minimum bases quality to consider in the metrics BAM
if [ -z "$METRICS_MINIMUM_BASE_QUALITY" ] || ! [[ $METRICS_MINIMUM_BASE_QUALITY =~ ^[0-9]+$ ]]; then
	METRICS_MINIMUM_BASE_QUALITY=20
fi;
export METRICS_MINIMUM_BASE_QUALITY

# Clipping overlapping reads in the metrics BAM
if [ -z "$CLIP_OVERLAPPING_READS" ]; then
	CLIP_OVERLAPPING_READS=1
fi;
export CLIP_OVERLAPPING_READS

# Flagged reads in the metrics BAM (mpileup format)
# default UNMAP,SECONDARY,QCFAIL,DUP
if [ -z "$METRICS_FLAGS" ]; then
	METRICS_FLAGS="UNMAP,SECONDARY,QCFAIL,DUP"
fi;
METRICS_FLAGS=$(echo $METRICS_FLAGS | tr -d " ")
export METRICS_FLAGS

# Flagged reads in the metrics BAM (samtools format)
# Generated from mpileup format (if empty)
if [ -z "$SAMTOOLS_METRICS_FLAG_PARAM" ]; then
	#SAMTOOLS_METRICS_FLAG_PARAM=" -F 0x4 -F 0x100 -F 0x200 -F 0x400"
	SAMTOOLS_METRICS_FLAG_PARAM=""
	for F in $(echo $METRICS_FLAGS | tr "," " "); do
		[ "$F" == "UNMAP" ] && SAMTOOLS_METRICS_FLAG_PARAM=$SAMTOOLS_METRICS_FLAG_PARAM" -F 0x4 "
		[ "$F" == "SECONDARY" ] && SAMTOOLS_METRICS_FLAG_PARAM=$SAMTOOLS_METRICS_FLAG_PARAM" -F 0x100 "
		[ "$F" == "QCFAIL" ] && SAMTOOLS_METRICS_FLAG_PARAM=$SAMTOOLS_METRICS_FLAG_PARAM" -F 0x200 "
		[ "$F" == "DUP" ] && SAMTOOLS_METRICS_FLAG_PARAM=$SAMTOOLS_METRICS_FLAG_PARAM" -F 0x400 "
	done;
fi;
export SAMTOOLS_METRICS_FLAG_PARAM


# VARIANT RECALIBRATION (default 0/FALSE/NO/N)
# Performs Variant recalibration after Calling (1/TRUE/YES/Y or 0/FALSE/NO/N).
# If the recalibration fail (usually due to lack of data for statistic calculation), nothing will be done
if [ -z "$VARIANT_RECALIBRATION" ] || [ "${VARIANT_RECALIBRATION^^}" == "FALSE" ] || [ "${VARIANT_RECALIBRATION^^}" == "NO" ] || [ "${VARIANT_RECALIBRATION^^}" == "N" ]  || [ "$VARIANT_RECALIBRATION" == "0" ]; then
	VARIANT_RECALIBRATION=0
elif [ "${VARIANT_RECALIBRATION^^}" == "TRUE" ] || [ "${VARIANT_RECALIBRATION^^}" == "YES" ] || [ "${VARIANT_RECALIBRATION^^}" == "Y" ]  || [ "$VARIANT_RECALIBRATION" == "1" ]; then
	VARIANT_RECALIBRATION=1

else
	VARIANT_RECALIBRATION=0
fi;
export VARIANT_RECALIBRATION

# INTERVAL_PADDING (default 0)
# Add some “padding” to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
if [ -z "$INTERVAL_PADDING" ] || ! [[ $INTERVAL_PADDING =~ ^[0-9]+$ ]]; then
	INTERVAL_PADDING=0
fi;
export INTERVAL_PADDING

# COVERAGE CRITERIA (default "1,30")
# For gene coverage metrics
# the criteria to calculate the percent of bases over $COVERAGE_CRITERIA X (eg 30 for 30X)
if [ -z "$COVERAGE_CRITERIA" ]; then
	COVERAGE_CRITERIA="1,5,10,20,30,50,100,200,300"
fi;
export COVERAGE_CRITERIA


# COVERAGE DP THRESHOLD (default "30" "100" "1")
# For gene coverage metrics
# the criteria test if genes failed (or just warning) the coverage threshold
# Sequencing depth threshold
if [ -z "$SEQUENCING_DEPTH" ]; then
	SEQUENCING_DEPTH="1"
fi;
export SEQUENCING_DEPTH
# Sequencing coverage threshold
if [ -z "$SEQUENCING_COVERAGE_THRESHOLD" ]; then
	SEQUENCING_COVERAGE_THRESHOLD="1"
fi;
export SEQUENCING_COVERAGE_THRESHOLD
# fail DP threshold (default 30X)
if [ -z "$MINIMUM_DEPTH" ]; then
	MINIMUM_DEPTH="30"
fi;
export MINIMUM_DEPTH
# warn DP threshold (default 100X)
if [ -z "$EXPECTED_DEPTH" ]; then
	EXPECTED_DEPTH="100"
fi;
export EXPECTED_DEPTH
# threshold percentage of bases over the DP threshold
if [ -z "$DEPTH_COVERAGE_THRESHOLD" ]; then
	DEPTH_COVERAGE_THRESHOLD="0.95"
fi;
export DEPTH_COVERAGE_THRESHOLD


# COVERAGE HARMONIZATION
COVERAGE_CRITERIA=$(echo "$SEQUENCING_DEPTH,$MINIMUM_DEPTH,$EXPECTED_DEPTH,$COVERAGE_CRITERIA" | tr "," "\n" | tr " " "\n" | grep -v "^$" | sort -u -n | tr "\n" ","  | sed "s/,$//")


# CLIP_OVERLAPPING_READS (default 1)
# From PICARD: For paired reads, soft clip the 3' end of each read if necessary so that it does not extend past the 5' end of its mate
if [ -z "$CLIP_OVERLAPPING_READS" ]; then
	CLIP_OVERLAPPING_READS=1
fi;
export CLIP_OVERLAPPING_READS


# NB_BASES_AROUND (default 0)
# For gene coverage metrics
# the number of bases to look around the exons from the given bed file
if [ -z $NB_BASES_AROUND ] || ! [[ $NB_BASES_AROUND =~ ^[0-9]+$ ]]; then
	NB_BASES_AROUND=0
fi;
export NB_BASES_AROUND

# BEDFILE_GENES (default "")
# For gene coverage metrics
# the bed file containing the 5'UTR, 3'UTR and genomic coding coordinates.
#BEDFILE_GENES=""
#if [ "$BEDFILE_GENES" != "" ] && [ ! -s $BEDFILE_GENES ]; then
	export BEDFILE_GENES
#fi;

# VARANK ANALYSIS (default 0/FALSE/NO/N)
# Performs VARANK ANALYSIS with Alamut (1/TRUE/YES/Y or 0/FALSE/NO/N).
if [ -z $VARANK_ANALYSIS ] || [ "${VARANK_ANALYSIS^^}" == "FALSE" ] || [ "${VARANK_ANALYSIS^^}" == "NO" ] || [ "${VARANK_ANALYSIS^^}" == "N" ]  || [ "$VARANK_ANALYSIS" == "0" ]; then
	VARANK_ANALYSIS=0
elif [ "${VARANK_ANALYSIS^^}" == "TRUE" ] || [ "${VARANK_ANALYSIS^^}" == "YES" ] || [ "${VARANK_ANALYSIS^^}" == "Y" ]  || [ "$VARANK_ANALYSIS" == "1" ]; then
	VARANK_ANALYSIS=1
else
	VARANK_ANALYSIS=0
fi;
export VARANK_ANALYSIS
# if yes, define a folder to store the results
#VARANK_FOLDER=$FOLDER_RESULTS/VARANK
export VARANK_FOLDER

# BAM CHECK (default 0)
# Check BAM for each manipulation step (clipping, realignment...)
# Time consuming, but will stop the analysis in case of missing reads.
# A Metrics on each BAM check the BAM decrependy in any case
if [ -z $BAM_CHECK_STEPS ] || [ "${BAM_CHECK_STEPS^^}" == "FALSE" ] || [ "${BAM_CHECK_STEPS^^}" == "NO" ] || [ "${BAM_CHECK_STEPS^^}" == "N" ]  || [ "$BAM_CHECK_STEPS" == "0" ]; then
	BAM_CHECK_STEPS=0
elif [ "${BAM_CHECK_STEPS^^}" == "TRUE" ] || [ "${BAM_CHECK_STEPS^^}" == "YES" ] || [ "${BAM_CHECK_STEPS^^}" == "Y" ]  || [ "$BAM_CHECK_STEPS" == "1" ]; then
	BAM_CHECK_STEPS=1
else
	BAM_CHECK_STEPS=0
fi;
export BAM_CHECK_STEPS

# METRICS SNPEFF (default 0)
# Generate snpEff variant metrics from VCF
if [ -z $METRICS_SNPEFF ] || [ "${METRICS_SNPEFF^^}" == "FALSE" ] || [ "${METRICS_SNPEFF^^}" == "NO" ] || [ "${METRICS_SNPEFF^^}" == "N" ]  || [ "$METRICS_SNPEFF" == "0" ]; then
	METRICS_SNPEFF=0
elif [ "${METRICS_SNPEFF^^}" == "TRUE" ] || [ "${METRICS_SNPEFF^^}" == "YES" ] || [ "${METRICS_SNPEFF^^}" == "Y" ]  || [ "$METRICS_SNPEFF" == "1" ]; then
	METRICS_SNPEFF=1
else
	METRICS_SNPEFF=0
fi;
export METRICS_SNPEFF


# PIPELINES PRIORITIZATION
# List of pipelines to prioritize for the report (final.vcf)
#PRIORITIZE_PIPELINES_LIST=""
export PRIORITIZE_PIPELINES_LIST

# POST ALIGNEMENT STEPS (default .uncompressed.unclipped.unrealigned.unsorted)
# All steps after alignement and before calling
# This sequence correspond to the BAM file generated jsut after the alignemnt
# Format: ".step3.step2.step1"
# Example: .uncompressed.unrealigned.unsorted
#    This sequence will generate the file $ALIGNER.uncompressed.unrealigned.unsorted.bam
#    Then, this BAM file will be 1/ sorted, 2/ realigned and 3/ compressed
# The steps are defined as makefiles rules
# Available steps:
#    sorting: BAM sorting
#    compress: BAM compression (see $BAM_COMPRESSION variable)
#    realignment: local realignment
#    markduplicates: BAM Mark Duplicates
#    clipping: BAM Clipping according to primer definition in manifest file, if any
if [ -z "$POST_ALIGNMENT_STEPS" ]; then
	POST_ALIGNMENT_STEPS="sorting realignment clipping compress" #.uncompressed.unclipped.unrealigned.unsorted
fi;
#echo $POST_ALIGNMENT;
# Create POST_ALIGNEMENT variable
#echo "."$(echo $(echo $test | tr "," " " | tr "." " " | tr " " "\n" | tac | tr "\n" " ") | tr " " ".")
POST_ALIGNMENT="."$(echo $POST_ALIGNMENT_STEPS | tr "," " " | tr "." " " | tr " " "\n" | sed '/^$/d' | tac  | tr "\n" " " | sed s/\.$//g | tr " " ".")
#echo $POST_ALIGNMENT;
#exit 1;
export POST_ALIGNMENT

# BAM COMPRESSION
# Final BAM copression level (unaligned.bam, ALIGNER.bam)
if [ -z $BAM_COMPRESSION ] || ! [[ $BAM_COMPRESSION =~ ^[0-9]$ ]]; then
	BAM_COMPRESSION=9
fi;
export BAM_COMPRESSION




# ANNOTATION
#############

# DEFAULT

# Default annotation with HOWARD for HOWARD for intermediate VCF (for each caller) used by default with annotation rule "howard"
if [ -z "$HOWARD_ANNOTATION" ]; then
	#HOWARD_ANNOTATION="core,symbol,location,outcome,hgvs,snpeff,snpeff_hgvs"
	#HOWARD_ANNOTATION="core,snpeff_hgvs"
	HOWARD_ANNOTATION="snpeff_split"
	# ANNOTATINO FULL
	#HOWARD_ANNOTATION="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs"
fi;
# DEJAVU
if [ -s $ANNOVAR_DATABASES/$ASSEMBLY"_dejavu."$GROUP.$PROJECT.txt ]; then
	HOWARD_ANNOTATION="$HOWARD_ANNOTATION,dejavu.$GROUP.$PROJECT"
fi
# DEJAVU for all the GROUP
for DEJAVU_DATABASE_ONE in $(find $ANNOVAR_DATABASES -name $ASSEMBLY"_dejavu."$GROUP"*.txt" 2>/dev/null); do
	HOWARD_ANNOTATION=$HOWARD_ANNOTATION","$(basename $DEJAVU_DATABASE_ONE | sed s/^$ASSEMBLY"_"//g | sed s/.txt$//g );
done;
export HOWARD_ANNOTATION

# MINIMAL

# Minimal annotation with HOWARD for minimal VCF annotation  (rule howard_minimal)
if [ -z "$HOWARD_ANNOTATION_MINIMAL" ]; then
	HOWARD_ANNOTATION_MINIMAL=$HOWARD_ANNOTATION
fi;
export HOWARD_ANNOTATION_MINIMAL

# REPORT
# Default annotation with HOWARD for Report (/fullfinal VCF)
if [ -z "$HOWARD_ANNOTATION_REPORT" ]; then
	HOWARD_ANNOTATION_REPORT=$HOWARD_ANNOTATION
fi;
export HOWARD_ANNOTATION_REPORT

# ANALYSIS
# Default annotation with HOWARD for whole analysis (annotation NOT forced)
if [ -z "$HOWARD_ANNOTATION_ANALYSIS" ]; then
	HOWARD_ANNOTATION_ANALYSIS=$HOWARD_ANNOTATION
fi;
export HOWARD_ANNOTATION_ANALYSIS



# CALCULATION
###############

# DEFAULT
# Default calculation with HOWARD
if [ -z "$HOWARD_CALCULATION" ]; then
	HOWARD_CALCULATION="VARTYPE,NOMEN"
fi;
export HOWARD_CALCULATION

# MINIMAL
# Minimal calculation with HOWARD
if [ -z "$HOWARD_CALCULATION_MINIMAL" ]; then
	HOWARD_CALCULATION_MINIMAL=$HOWARD_CALCULATION
fi;
export HOWARD_CALCULATION_MINIMAL

# REPORT
# Default calculation with HOWARD for Report (full/final VCF)
if [ -z "$HOWARD_CALCULATION_REPORT" ]; then
	HOWARD_CALCULATION_REPORT=$HOWARD_CALCULATION
fi;
export HOWARD_CALCULATION_REPORT

# ANALYSIS
# Default calculation with HOWARD for whole analysis (prioritization forced)
if [ -z "$HOWARD_CALCULATION_ANALYSIS" ]; then
	HOWARD_CALCULATION_ANALYSIS=$HOWARD_CALCULATION
fi;
export HOWARD_CALCULATION_ANALYSIS


# List of annotation fields to extract NOMEN annotation (default 'hgvs', see HOWARD docs)
if [ -z "$HOWARD_NOMEN_FIELDS" ]; then
	HOWARD_NOMEN_FIELDS="hgvs"
fi;
export HOWARD_NOMEN_FIELDS



# PRIORITIZATION
###################
# This option create ranking scores in VCF and comment in TXT (after translation).
# Scores can be used to sort variant in the TXT

# DEFAULT
# Default filter to prioritize/rank variant.
if [ -z "$HOWARD_PRIORITIZATION" ]; then
	HOWARD_PRIORITIZATION=$HOWARD_PRIORITIZATION_DEFAULT
fi;
if [ ! -z "$APP_NAME" ]; then
	HOWARD_PRIORITIZATION="$HOWARD_PRIORITIZATION $APP_NAME"
fi;
HOWARD_PRIORITIZATION=$(echo $HOWARD_PRIORITIZATION | tr "," " " | tr " " "\n" | sort -u | tr "\n" "," | sed s/,$//)
export HOWARD_PRIORITIZATION

# MINIMAL
# Minimal filter to prioritize/rank variant.
if [ -z "$HOWARD_PRIORITIZATION_MINIMAL" ]; then
	HOWARD_PRIORITIZATION_MINIMAL=$HOWARD_PRIORITIZATION
fi;
if [ ! -z "$APP_NAME" ]; then
	HOWARD_PRIORITIZATION_MINIMAL="$HOWARD_PRIORITIZATION_MINIMAL $APP_NAME"
fi;
HOWARD_PRIORITIZATION_MINIMAL=$(echo $HOWARD_PRIORITIZATION_MINIMAL | tr "," " " | tr " " "\n" | sort -u | tr "\n" "," | sed s/,$//)
export HOWARD_PRIORITIZATION_MINIMAL

#if [ -z "$HOWARD_PRIORITIZATION_MINIMAL" ]; then
#	HOWARD_PRIORITIZATION_MINIMAL=$HOWARD_PRIORITIZATION
#fi;
#export HOWARD_PRIORITIZATION_MINIMAL

# Report
# Default filter to prioritize/rank variant for Report
if [ -z "$HOWARD_PRIORITIZATION_REPORT" ]; then
	HOWARD_PRIORITIZATION_REPORT=$HOWARD_PRIORITIZATION
fi;
if [ ! -z "$APP_NAME" ]; then
	HOWARD_PRIORITIZATION_REPORT="$HOWARD_PRIORITIZATION_REPORT $APP_NAME"
fi;
HOWARD_PRIORITIZATION_REPORT=$(echo $HOWARD_PRIORITIZATION_REPORT | tr "," " " | tr " " "\n" | sort -u | tr "\n" "," | sed s/,$//)
export HOWARD_PRIORITIZATION_REPORT

#if [ -z "$HOWARD_PRIORITIZATION_REPORT" ]; then
#	HOWARD_PRIORITIZATION_REPORT=$HOWARD_PRIORITIZATION
#fi;
#export HOWARD_PRIORITIZATION_REPORT


# ANALYSIS
# Default filter to prioritize/rank variant for whole analysis (calculation forced)
if [ -z "$HOWARD_PRIORITIZATION_ANALYSIS" ]; then
	HOWARD_PRIORITIZATION_ANALYSIS=$HOWARD_PRIORITIZATION
fi;
if [ ! -z "$APP_NAME" ]; then
	HOWARD_PRIORITIZATION_ANALYSIS="$HOWARD_PRIORITIZATION_ANALYSIS $APP_NAME"
fi;
HOWARD_PRIORITIZATION_ANALYSIS=$(echo $HOWARD_PRIORITIZATION_ANALYSIS | tr "," " " | tr " " "\n" | sort -u | tr "\n" "," | sed s/,$//)
export HOWARD_PRIORITIZATION_ANALYSIS


# TRANSLATION
################
# List of fields to show in the TXT file
# use ALL to show ALL "other" annotations

# DEFAULT
# Fields to show after translation
if [ -z $HOWARD_FIELDS ]; then
	HOWARD_FIELDS="NOMEN,PZFlag,PZScore,PZComment,CNOMEN,PNOMEN,location,outcome,snpeff_impact,VAF_average,dbSNP,dbSNPNonFlagged,popfreq,ALL"
fi;
export HOWARD_FIELDS
# Sort variant in the TXT using 2 fields
if [ -z $HOWARD_SORT_BY ]; then
	HOWARD_SORT_BY="PZFlag,PZScore"
fi;
export HOWARD_SORT_BY
# Order fields in variant ranking
if [ -z $HOWARD_ORDER_BY ]; then
	HOWARD_ORDER_BY="DESC,DESC"
fi;
export HOWARD_ORDER_BY

# MINIMAL
# Fields to show after minimal translation
if [ -z $HOWARD_FIELDS_MINIMAL ]; then
	#HOWARD_FIELDS_MINIMAL=$HOWARD_FIELDS
	HOWARD_FIELDS_MINIMAL=$(echo $HOWARD_FIELDS | sed "s/ALL//")
fi;
export HOWARD_FIELDS_MINIMAL
# Sort variant in the TXT using 2 fields
if [ -z $HOWARD_SORT_BY_MINIMAL ]; then
	HOWARD_SORT_BY_MINIMAL=$HOWARD_SORT_BY
fi;
export HOWARD_SORT_BY_MINIMAL
# Order fields in variant ranking
if [ -z $HOWARD_ORDER_BY_MINIMAL ]; then
	HOWARD_ORDER_BY_MINIMAL=$HOWARD_ORDER_BY
fi;
export HOWARD_ORDER_BY_MINIMAL

# REPORT
# Fields to show after minimal translation
if [ -z $HOWARD_FIELDS_REPORT ]; then
	#HOWARD_FIELDS_REPORT=$HOWARD_FIELDS
	HOWARD_FIELDS_REPORT=$(echo $HOWARD_FIELDS | sed "s/ALL//")
fi;
export HOWARD_FIELDS_REPORT
# Sort variant in the TXT using 2 fields
if [ -z $HOWARD_SORT_BY_REPORT ]; then
	HOWARD_SORT_BY_REPORT=$HOWARD_SORT_BY
fi;
export HOWARD_SORT_BY_REPORT
# Order fields in variant ranking
if [ -z $HOWARD_ORDER_BY_REPORT ]; then
	HOWARD_ORDER_BY_REPORT=$HOWARD_ORDER_BY
fi;
export HOWARD_ORDER_BY_REPORT


# DAEMON CONFIG
##########

# Config
export TMP_REMOVE_TIME=10
export LOG_FILE=analysis.log
export RUN_ANALYSIS_CONDITIONS="RTAComplete.txt SampleSheet.csv"


# DATABASE PEPPER CONFIG (depreciated)
########################################

CONFIG=$PEPPER_CONFIG

# RELEASE
###########

README_FILE=$STARK_FOLDER_CONFIG/README
RELEASE_FILE=$STARK_FOLDER_CONFIG/RELEASE

# INFRASTRUCTURE
##################

# THREADS is CORES - 1

if [ -z $THREADS ] || [ "${THREADS^^}" == "AUTO" ] || ! [[ $THREADS =~ ^[0-9]+$ ]]; then
	export CORES=$(ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w)	# NB of cores in the server
	export CORES_FREE=1							# Number of threads free for other command
	export CORES_TO_USE=$(($CORES-$CORES_FREE))				# Nb of cores to use
	THREADS=$CORES_TO_USE
	#THREADS=2
	#if ! [[ $THREADS =~ "^\d+" ]]; then # || [ "$THREADS" == "" ]; then
	#if [[ -z $THREADS ]]; then # || [ "$THREADS" == "" ]; then
	#	THREADS=1
	#fi;
	#echo "|$THREADS|"
	if [ $THREADS -gt $CORES ]; then
		THREADS=$CORES
	fi;

fi;
export THREADS
export CORES_TO_USE=$THREADS

#export THREADS_BY_SAMPLE=1				# Number of threads for a sample's command
#export THREADS_BY_PIPELINE=1				# Number of threads for a pipeline's command
#export THREADS_BY_ALIGNER=1				# Number of threads for a aligner's command
#export THREADS_BY_CALLER=1				# Number of threads for a caller's command
#export THREADS=$(($CORES_TO_USE/$THREADS_BY_SAMPLE))	# Total nb of threads to use for Makefile
#export THREADS_BWA=$THREADS_BY_SAMPLE			# NB of threads for BWA command
#export THREADS_RTC=$THREADS_BY_SAMPLE			# NB of threads for GATK realigment RealignerTargetCreator command

export THREADS_AUTO=1					# Automatize nb of threads depending on the input and number of cores to use (TRUE if = 1)

# MEMORY
MEMTOTAL=$(cat /proc/meminfo | grep MemTotal | awk '{print $2}')	# MEMORY in octet
#export MEMORY=$(($MEMTOTAL/$THREADS/1024/1024))			# MEMORY in Go
MEMORY=$(($MEMTOTAL/$CORES_TO_USE/1024/1024))			# MEMORY in Go
if [ "$MEMORY" == "" ] || [ $MEMORY -lt 1 ]; then
	MEMORY=1;
fi;
export MEMORY				# MEMORY in Go
export MEMORY_JAVA=$MEMORY		# 4Go
export JAVA_MEMORY=$MEMORY_JAVA		# for naming...


# JAVA_MEMORY by sample and by caller
###############

if ((1)); then

	NB_SAMPLE=1
	NB_PIPELINES=1
	NB_ALIGNERS=1
	NB_CALLERS=1

	# Default
	JAVA_MEMORY_BY_SAMPLE=$JAVA_MEMORY
	JAVA_MEMORY_BY_CALLER=$JAVA_MEMORY

	# Calcul
	JAVA_MEMORY_BY_SAMPLE=$(($MEMTOTAL/$NB_SAMPLE/1024/1024))				# Number of threads by sample
	JAVA_MEMORY_BY_PIPELINE=$(($JAVA_MEMORY_BY_SAMPLE/$NB_PIPELINES))	# Number of threads for a pipeline's command
	JAVA_MEMORY_BY_ALIGNER=$(($JAVA_MEMORY_BY_SAMPLE/$NB_ALIGNERS))	# Number of threads for a aligner's command
	JAVA_MEMORY_BY_CALLER=$(($JAVA_MEMORY_BY_SAMPLE/$NB_CALLERS))		# Number of threads for a caller's command

	# Test
	if [ "$JAVA_MEMORY_BY_SAMPLE" == "" ] || [ $JAVA_MEMORY_BY_SAMPLE -lt 1 ]; then JAVA_MEMORY_BY_SAMPLE=1; fi;
	if [ "$JAVA_MEMORY_BY_PIPELINE" == "" ] || [ $JAVA_MEMORY_BY_PIPELINE -lt 1 ]; then JAVA_MEMORY_BY_PIPELINE=1; fi;
	if [ "$JAVA_MEMORY_BY_ALIGNER" == "" ] || [ $JAVA_MEMORY_BY_ALIGNER -lt 1 ]; then JAVA_MEMORY_BY_ALIGNER=1; fi;
	if [ "$JAVA_MEMORY_BY_CALLER" == "" ] || [ $JAVA_MEMORY_BY_CALLER -lt 1 ]; then JAVA_MEMORY_BY_CALLER=1; fi;

fi;

# JAVA FLAGS
##############
export JAVA_FLAGS_TMP_FOLDER=" -Dorg.xerial.snappy.tempdir=$TMP_FOLDER_TMP -Djava.io.tmpdir=$TMP_FOLDER_TMP";
export JAVA_FLAGS_OTHER_PARAM=" -Dsnappy.disable=true -Dsamjdk.try_use_intel_deflater=false ";
export JAVA_FLAGS_DEFAULT=" -Xmx"$JAVA_MEMORY"g $JAVA_FLAGS_OTHER_PARAM $JAVA_FLAGS_TMP_FOLDER";
#if [ "$JAVA_FLAGS" == "" ]; then
#	JAVA_FLAGS=$JAVA_FLAGS_DEFAULT;
#fi;
JAVA_FLAGS=$JAVA_FLAGS_DEFAULT;
export JAVA_FLAGS;
#if [ "$JAVA_FLAGS_BY_SAMPLE" == "" ]; then
#	JAVA_FLAGS_BY_SAMPLE=" -Xmx"$JAVA_MEMORY_BY_SAMPLE"g $JAVA_FLAGS_OTHER_PARAM $JAVA_FLAGS_TMP_FOLDER";
#fi;
#export JAVA_FLAGS_BY_SAMPLE;


# CHECK
#re='^[0-9]+$'
#if ! [[ $THREADS =~ $re ]] || [ "$THREADS" == "" ]; then
#	THREADS=1
#fi;

# TEST
#THREADS=16
