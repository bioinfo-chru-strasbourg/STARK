#################################
## STARK environment
#################################



# APPLICATION INFOS
############

# Name of the APP. usefull to include specific rules if exists (in $STARK/$APP_NAME.rules.mk/*rules.mk)
# AUTO detect: $(basename ${BASH_SOURCE[0]} | sed 's/^env//' | sed 's/\.sh$//gi' | cut -d. -f2-)
if [ "$APP_NAME" == "" ]; then
	APP_NAME="UNKNOWN"
fi;
export APP_NAME

# GROUP and PROJECT Associated with the APP
# Use to structure data in the repository folder
if [ "$GROUP" == "" ]; then
	GROUP="UNKNOWN"
fi;
export GROUP
if [ "$PROJECT" == "" ]; then
	PROJECT="UNKNOWN"
fi;
export PROJECT


# SCRIPT DIR
##############

# Folder of STARK ENV
STARK_FOLDER="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )"


# FOLDERS
#################
# Checking and adaptation of variables for the workflow


# RUN FOLDER
export MISEQ_FOLDER=$FOLDER_RUN				# MISEQ RAW folder
MSR_SUBFOLDER=Data/Intensities/BaseCalls		# SUBFOLDER Illumina Sequencer Data subfolder


# MANIFEST FOLDER
export MANIFEST_FOLDER=$FOLDER_MANIFEST			# MISEQ RAW folder


# RESULTS FOLDER
export DEMULTIPLEXING_FOLDER=$FOLDER_RESULTS/DEM	# DEMULTIPLEXING folder
export RESULTS_FOLDER=$FOLDER_RESULTS/RES		# RESULTS folder ALL
export ANALYSIS_FOLDER=$FOLDER_RESULTS/ANA		# ANALYSIS folder
export TMP_FOLDER_TMP=$FOLDER_RESULTS/TMP		# TEMPORARY folder tmp
export VALIDATION_FOLDER=$FOLDER_RESULTS/VAL		# VALIDATION folder
export TMP_FOLDER=$FOLDER_RESULTS			# TEMPORARY folder = MAIN V2 folder
export MAIN_FOLDER=$FOLDER_RESULTS			# MAIN V2 Folder (RES, DEM...)

# Create directories if does not exist
mkdir -p $DEMULTIPLEXING_FOLDER $RESULTS_FOLDER $ANALYSIS_FOLDER $TMP_FOLDER_TMP $VALIDATION_FOLDER $TMP_FOLDER $MAIN_FOLDER

# TOOLS FOLDER
export NGS_TOOLS=$FOLDER_ENV/tools			# TOOLS
export NGS_GENOMES=$FOLDER_ENV/genomes			# GENOMES folder
export BDFOLDER=$FOLDER_ENV/db				# DB
export RULES=$STARK_FOLDER/rules.mk/*.rules.mk		# RULES
export GENOMES=$NGS_GENOMES				# Genomes folder in NGS folder
export TMP_SYS_FOLDER=/tmp				# TEMPORARY SYSTEM folder (tmpfs)
export NGS_FOLDER=$FOLDER_ENV
export NGS_SCRIPTS=$STARK_FOLDER			#"$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )"


# REPOSITORY FOLDER
export RESULTS_FOLDER_BY_GROUP_PROJECT_COPY=$FOLDER_REPOSITORY	# Copy data into group and project (if defined)
export RESULTS_SUBFOLDER_DATA="DATA";				# Copy sample results in a SUBFOLDER
# Copy some sample results files in the root sample folder, if any SUBFOLDER defined
export RESULTS_SUBFOLDER_ROOT_FILE_PATTERNS=$RESULTS_SUBFOLDER_ROOT_FILE_PATTERNS' $SAMPLE.reports/$SAMPLE.final.vcf $SAMPLE.reports/$SAMPLE.full.vcf $SAMPLE.reports/$SAMPLE.final.txt $SAMPLE.reports/$SAMPLE.full.txt *.reports/*.report.pdf *.reports/latex*/*pdf';	 
								

# RULES for the application
if [ ! -z "$(ls $STARK_FOLDER/$APP_NAME.rules.mk/*rules.mk 2>/dev/null)" ]; then
	RULES_APP="$RULES_APP $STARK_FOLDER/$APP_NAME.rules.mk/*.rules.mk"
fi;
export RULES="$RULES $RULES_APP"

#echo ${BASH_SOURCE[2]}; return 0;

# TOOLS
#########

source $STARK_FOLDER/env_tools.sh


# ASSEMBLY
##########

# default Assembly
if [ -z $ASSEMBLY ] || [ "$ASSEMBLY" == "" ]; then
	ASSEMBLY=hg19				
fi;
export ASSEMBLY
export REF=$GENOMES/$ASSEMBLY/$ASSEMBLY.fa	# Defautl REF
export REFSEQ_GENES=$BDFOLDER/RefSeq.$ASSEMBLY.bed


# DATABASES
#############

# Mandatory DB (for calling with GATK variant recalibration)

# BDSNP DB (used for calling, espacially with GATK)
export VCFDBSNP=$BDFOLDER/dbsnp_138.$ASSEMBLY.vcf.gz	#snp138.vcf.gz
if [ ! -e $VCFDBSNP ]; then
	echo "#[ERROR] No VCFDBSNP '$VCFDBSNP' in the database. Calling step impossible. Please check '$BDFOLDER' folder or configuration file"
	exit 1;
fi;

# VCF DB for recalibration
#export VCFDBSNP_RECALIBRATION=$BDFOLDER/dbsnp_137.$ASSEMBLY.vcf
export VCFDBSNP_RECALIBRATION=$VCFDBSNP
if (($VARIANT_RECALIBRATION)) && [ ! -e $VCFDBSNP_RECALIBRATION ]; then
	echo "#[ERROR] No VCFDBSNP_RECALIBRATION '$VCFDBSNP_RECALIBRATION' in the database. Calling step impossible. Please check '$BDFOLDER' folder or configuration file"
	exit 1;
fi;

# Other databases
export VCFMILLS1000G=$BDFOLDER/Mills_and_1000G_gold_standard.indels.$ASSEMBLY.sites.vcf
export COSMIC=$DBFOLDER/cosmic.$ASSEMBLY.vcf
export KNOWN_ALLELES=$VCFMILLS1000G
export HAPMAP=$BDFOLDER/hapmap_3.3.$ASSEMBLY.sites.vcf
export OMNI=$BDFOLDER/1000G_omni2.5.$ASSEMBLY.vcf
export PHASE1_1000G=$BDFOLDER/1000G_phase1.snps.high_confidence.$ASSEMBLY.sites.vcf




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
PIPELINES=$(echo $PIPELINES | tr " " "\n" | sort | uniq | tr "\n" " " | sed -e 's/^[[:space:]]*//'  | sed -e 's/[[:space:]]*$//')
if [ -z "$PIPELINES" ]; then
	PIPELINES="bwamem.gatkHC.howard"
fi;
export PIPELINES

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
	COVERAGE_CRITERIA="1,30"
fi;
export COVERAGE_CRITERIA

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
export BEDFILE_GENES

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
	BAM_COMPRESSION=5
fi;
export BAM_COMPRESSION



# HOWARD ANNOTATION/PRIOTITIZATION/TRANSLATION CONFIGURATION
# ANNOTATION
# Default annotation with HOWARD
if [ -z "$ANNOTATION_TYPE" ]; then
	ANNOTATION_TYPE="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs"
fi;
export ANNOTATION_TYPE

# PRIORITIZATION
# Default filter to prioritize/rank variant.
# This option create ranking scores in VCF and comment in TXT (after translation).
# Scores can be used to sort variant in the TXT
if [ -z "$HOWARD_FILTER" ]; then
	if [ ! -z "$APP_NAME" ]; then
		HOWARD_FILTER="$APP_NAME"
	else
		HOWARD_FILTER="default"
	fi;
fi;
export HOWARD_FILTER


# TRANSLATION
# List of fields to show in the TXT file
# use ALL to show ALL "other" annotations
if [ -z "$HOWARD_FIELDS" ]; then
	HOWARD_FIELDS="PZScore,PZFlag,PZComment,Symbol,hgvs,location,outcome,AlleleFrequency,AD,dbSNP,dbSNPNonFlagged,popfreq,ALL"
fi;
export HOWARD_FIELDS

# Sort variant in the TXT using 2 fields
if [ -z "$HOWARD_SORT_BY" ]; then
	HOWARD_SORT_BY="PZFlag,PZScore"
fi;
export HOWARD_SORT_BY

# Order fields in variant ranking
if [ -z "$HOWARD_ORDER_BY" ]; then
	HOWARD_ORDER_BY="DESC,DESC"
fi;
export HOWARD_ORDER_BY




# DAEMON CONFIG
##########

# Config
export TMP_REMOVE_TIME=10
export LOG_FILE=analysis.log
export RUN_ANALYSIS_CONDITIONS="RTAComplete.txt SampleSheet.csv"


# DATABASE CONFIG
###################

CONFIG=$PEPPER_CONFIG

# RELEASE
###########

README_FILE=$STARK_FOLDER/README
RELEASE_FILE=$STARK_FOLDER/RELEASE

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
if [ "$JAVA_FLAGS" == "" ]; then
	JAVA_FLAGS=$JAVA_FLAGS_DEFAULT;
fi;
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


