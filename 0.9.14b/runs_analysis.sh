#!/bin/bash
#################################
##
## NGS environment
## Analysis of multiple RUNS
##
## author: Antony Le Bechec
##
#################################

SCRIPT_NAME="RUNS Analysis"
SCRIPT_DESCRIPTION="Analysis of multiple RUNS"
SCRIPT_RELEASE="0.9.8.3"
SCRIPT_DATE="18/12/2015"
SCRIPT_AUTHOR="Antony Le Bechec / Amandine Velt"
SCRIPT_COPYRIGHT="IRC"

# Realse note
#RELEASE_NOTES="# Release Notes\n"
RELEASE_NOTES=$RELEASE_NOTES"0.9.6.1-07/01/2015: Copying temporary files: Modification of user's permissions. Copy\n"
RELEASE_NOTES=$RELEASE_NOTES"0.9.7-30/03/2015: Temporary space debugged\n"
RELEASE_NOTES=$RELEASE_NOTES"0.9.7.1-02/06/2015: Add permissions on folders when not using temporary folder\n"
RELEASE_NOTES=$RELEASE_NOTES"0.9.7.2-20/08/2015: Add ALIGNERS/CALLERS/ANNOTATORS in ENV\n"
RELEASE_NOTES=$RELEASE_NOTES"0.9.8-08/12/2015: Parameters for analysis from the SampleSheet\n"
RELEASE_NOTES=$RELEASE_NOTES"0.9.8.1-10/12/2015: Debug, param.sh generation\n"
RELEASE_NOTES=$RELEASE_NOTES"0.9.8.2-18/12/2015: Add copy files by samples, option by group/project\n"
RELEASE_NOTES=$RELEASE_NOTES"0.9.8.3-29/02/2016: Bug on Sample Project correction\n"

# Header
function header () {
	echo "#######################################";
	echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
	echo "# $SCRIPT_DESCRIPTION ";
	echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © $SCRIPT_LICENCE";
	echo "#######################################";
}

# Release
function release () {
	echo "# RELEASE NOTES:";
	echo -e $RELEASE_NOTES
}

# Usage
function usage {
	echo "# USAGE: $(basename $0) --fastq=<FASTQ> [options...]";
	echo "# -e/--env              ENV file configuration (default: defined in the RUN SampleSheet, or env.sh if not defined)";
	echo "# -r/--runs             List of RUNs to analyse, from Illumina sequencers (mandatory).";
	echo "# -a/--aligners         List of aligner to use for the analysis (optional, e.g. 'bwamem bwasw', default: in ENV or 'bwamem' if not defined).";
	echo "# -c/--callers          List of caller to use for the analysis (optional, e.g. 'gatkHC gatkUG VarScan samtools', default: in ENV or 'gatkHC' if not defined).";
	echo "# -k/--annotators       List of annotators to use for the analysis (optional, e.g. 'howard snpeff', default: in ENV or 'howard' if not defined).";
	echo "# -p/--pipelines        PIPELINES to launch, in the format ALIGNER.CALLER.ANNOTATOR, separated by a comma (optional, defaul: 'bwamem.gatkHC.howard'). Automatically defined in the environment file --env if any. This options has prior over options --aligners, --callers and --annotators";
	echo "# -f/--filter_samples   Samples to use for the analysis (optional, default: all samples).";
	echo "# -s/--samplesheet      Illumina SampleSheet.csv file (optional, default: found in RUN folder).";
	echo "# -z/--parallelization  Parallelization of RUN analyses (optional, default: FALSE). Nota Bene: not fully tested.";
	echo "# -b/--by_sample        Split analysis by SAMPLE, all threads on each sample, one by one (optional, default FALSE).";
	echo "# -v/--remove           Files in RES are removed (if any) before the RUN analysis of the run (optional).";
	
	echo "# -v/--verbose          VERBOSE option";
	echo "# -d/--debug            DEBUG option";
	echo "# -n/--release          RELEASE option";
	echo "# -h/--help             HELP option";
	echo "#";
}

# header
header;

PARAM=$@
### Script Structure
# 1. Configuration, Input parameters
#    1.1. Configuration
#    1.2. Input Parameters
#    1.3. Input Parameters test and default values
#    1.4. Variables
# 2. Preprocessing for each RUN
#    2.1. Check if RUN exists (Folder, SampleSheet)
#    2.2. Check previous analysis, if TMP space (TODO)
#    2.3. Demultiplexing (if needed)
#         2.3.1. Demultiplexing configuration
#         2.3.2. Demultiplexing processing
#    2.4. RUN analysis configuration
#         2.4.1. Makefile and shell from xml
#         2.4.2. Additional parameters to Makefile
# 3. Mutli Run Analysis Process
#    3.1. Combination of all RUN's Makefile
#    3.2. Main analysis process
# 4. Postprocessing for each RUN
#    4.1. Report process
#    4.2. Validation process
#    4.3. Distribution into Groups process (TODO)
# 5. Mutli Run PostProcessing
#    5.1. Report process
#    5.2. Validation process
# 6. Copying temporary files (TODO)
# 7. END

######################################
# 1. Configuration, Input parameters #
######################################

# 1.1. Configuration
#####################

####################################################################################################################################
# Define the function to print the usage of the script
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    sh runs_analysis.sh -e environment -r runs -a aligners -c callers -n annotators -f filter_samples -i intersect -u use_tmp_folder -d database_integration_option -m multi_run_analysis_process -s samplesheet [-h]

		Description:
		    This script is launch the demultiplexing step and the make steps

		Options:
		 	-e, --env
		 	This option is optionnal. It can be automatically defined by Stark via the SampleSheet.csv file.
		 	-r, --runs
		 	/!\ This option is required /!\ List of RUNs to analyse. It is of type "160430_NB551027_0001_AHC5LTAFXX" for Illumina sequencers.
		 	-a --aligners
		 	This option is optionnal. It can be automatically defined in the environment file. List of aligner to use for the analysis (e.g. "bwamem bwasw"). Default "bwamem".
		 	-c --callers
		 	This option is optionnal. It can be automatically defined in the environment file. List of caller to use for the analysis (e.g. "gatkHC gatkUG mutect samtools"). Default "gatkHC gatkUG".
		 	-n, --annotators
		 	This option is optionnal. It can be automatically defined in the environment file. List of annotators to use for the analysis (e.g. "howard"). (default "howard").
		 	-f --filter_samples
		 	This option is optionnal. Samples to use for the analysis. If empty, all the samples will be analysed.
		 	-i --intersect
		 	This option is optionnal. Number of workflow to consider to calculate the intesection list of variants. Default "2".
		 	-u --use_tmp_folder
		 	This option is optionnal. 1 for use it, 0 for not use it. Default "1", if TMP folder different than RESULTS or RAW folders.
		 	-d --database_integration_option
		 	This option is optionnal.
		 	-m --multi_run_analysis_process
		 	This option is optionnal. 0=Analysis of RUNs one by one, 1=Analysis of all RUNs together, 2=Reporting of all RUNs together. If 1 RUN, value=0. ERROR! force to value=0
		 	-s --samplesheet
		 	This option is optionnal. It can be automatically find by Stark via the environment variables.
		  	-h, --help
		 	Print this message and exit the program.
		__EOF__
}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "e:r:a:c:k:p:f:i:u:g:m:s:vdnh" --long "env:,runs:,aligners:,callers:,annotators:,pipelines:,filter_samples:,intersect:,use_tmp_folder:,database_integration_option:,multi_run_analysis_process:,samplesheet:,help" -- "$@" 2> /dev/null)
[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && \
	exit 1
eval set -- "$ARGS"
while true
do
	case "$1" in
		-e|--env)
			ENV="$2"
			shift 2 
			;;
		-r|--runs)
			RUNS="$2"
			shift 2 
			;;
		-a|--aligners)
			ALIGNERS_INPUT="$2"
			shift 2 
			;;
		-c|--callers)
			CALLERS_INPUT="$2"
			shift 2 
			;;
		-k|--annotators)
			ANNOTATORS_INPUT="$2"
			shift 2 
			;;
		-p|--pipelines)
			PIPELINES_INPUT="$2"
			shift 2 
			;;
		-f|--filter_samples)
			FILTER_SAMPLES="$2"
			shift 2 
			;;
		-i|--intersect)
			INTERSEC="$2"
			shift 2 
			;;
		-u|--use_tmp_folder)
			USE_TMP_FOLDER="$2"
			shift 2 
			;;
		-g|--database_integration_option)
			DATABASE_INTEGRATION_SCRIPT="$2"
			shift 2 
			;;
		-m|--multi_run_analysis_process)
			MULTI_RUN_ANALYSIS_PROCESS="$2"
			shift 2 
			;;
		-s|--samplesheet)
			SAMPLESHEET="$2"
			shift 2 
			;;
		-v|--verbose)
			VERBOSE=1
			shift 1
			;;
		-d|--debug)
			VERBOSE=1
			DEBUG=1
			shift 1
			;;
		-n|--release)
			release;
			exit 0
			;;
		-h|--help)
			usage
			exit 0
			;;
		--) shift
			break 
			;;
		*) 	echo "Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Checking the input parameter
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
[ "$RUNS" == "" ] && \
	echo "Option --runs is required. " "Use -h or --help to display the help." && exit 1;
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Other input parameter
ANALYSIS="requeue"  # "requeue" "report" (TODO?)

# 1.3. Input Parameters test and default values
################################################

# MULTI_RUN_ANALYSIS_PROCESS
#if [ $MULTI_RUN_ANALYSIS_PROCESS -ne 1 ] && [ $MULTI_RUN_ANALYSIS_PROCESS -ne 0 ] && [ $MULTI_RUN_ANALYSIS_PROCESS -ne 2 ]; then
if [ "$MULTI_RUN_ANALYSIS_PROCESS" != "1" ] && [ "$MULTI_RUN_ANALYSIS_PROCESS" != "0" ] && [ "$MULTI_RUN_ANALYSIS_PROCESS" != "2" ]; then
	MULTI_RUN_ANALYSIS_PROCESS=0
fi;

if [ $(echo $RUNS | wc -w) -eq 1 ]; then
	MULTI_RUN_ANALYSIS_PROCESS=0
fi;

#MULTI_RUN_ANALYSIS_PROCESS=0 # Forced value = 0 : ERROR in multi run processing when requeue...

# FILTER_SAMPLES
if [ "$FILTER_SAMPLES" == "" ]; then
	FILTER_SAMPLES=""
fi;
# INTERSECTION $USE_TMP_FOLDER
if [ "$INTERSEC" == "" ]; then
	INTERSEC="2"
fi;
# ENV
#if [ "$ENV" == "" ] || [ ! -s $ENV ]; then
#	ENV="env.sh"
#fi;
if [ -s $ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$ENV;
elif [ -s $SCRIPT_DIR/$ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$SCRIPT_DIR/$ENV;
elif [ "$ENV" == "" ] || [ ! -s $ENV ]; then
	if [ -s $SCRIPT_DIR/"env.sh" ]; then
		ENV=$SCRIPT_DIR/"env.sh";
	else
		ENV="";
	fi;
fi;

source $ENV

# DEFAULT ALIGNERS/CALLERS/ANNOTATORS in ENV!?


if ((1)); then

	# ALIGNERS # defined in ENV
	ALIGNERS_DEFAULT="bwamem";
	if [ "$ALIGNERS_INPUT" != "" ]; then # defined in INPUT
		ALIGNERS=$ALIGNERS_INPUT 
	elif [ "$ALIGNERS" == "" ]; then # default if not defined in ENV
		ALIGNERS=$ALIGNERS_DEFAULT # bwasw bwamemUnclipped
		echo "#[WARNING] NO ALIGNERS defined. Default ALIGNERS '$ALIGNERS' will be used."
	fi;
	ALIGNERS=$(echo $ALIGNERS | tr -d "\"")
	
	# CALLERS # defined in ENV
	CALLERS_DEFAULT="gatkHC";				
	if [ "$CALLERS_INPUT" != "" ]; then # defined in INPUT
		CALLERS=$CALLERS_INPUT
	elif [ "$CALLERS" == "" ]; then	# default if not defined in ENV
		CALLERS=$CALLERS_DEFAULT
		echo "#[WARNING] NO CALLERS defined. Default CALLERS '$CALLERS' will be used."
	fi;
	CALLERS=$(echo $CALLERS | tr -d "\"")
	
	# ANNOTATORS # defined in ENV
	ANNOTATORS_DEFAULT="howard";
	if [ "$ANNOTATORS_INPUT" != "" ]; then # defined in INPUT
		ANNOTATORS=$ANNOTATORS_INPUT
	elif [ "$ANNOTATORS" == "" ]; then # default if not defined in ENV
		ANNOTATORS=$ANNOTATORS_DEFAULT
		echo "#[WARNING] NO ANNOTATORS defined. Default ANNOTATORS '$ANNOTATORS' will be used."
	fi;
	ANNOTATORS=$(echo $ANNOTATORS | tr -d "\"")

	# PIPELINES # defined in ENV
	if [ "$ALIGNERS_INPUT" != "" ] || [ "$CALLERS_INPUT" != "" ] || [ "$ANNOTATORS_INPUT" != "" ] || [ "$PIPELINES_INPUT" != "" ]; then
		PIPELINES=$PIPELINES_INPUT #$ALIGNER_DEFAULT"."$CALLER_DEFAULT"."$ANNOTATOR_DEFAULT
	fi;
	# CREATE PIPELINES from ALIGNERS/CALLERS/ANNOTATORS
	if [ "$PIPELINES_INPUT" == "" ]; then
		for ALIGNER in $ALIGNERS; do
			for CALLER in $CALLERS; do
				for ANNOTATOR in $ANNOTATORS; do
					PIPELINES="$PIPELINES $ALIGNER.$CALLER.$ANNOTATOR"
				done;
			done;
		done;
		PIPELINES=$(echo $PIPELINES | tr " " "\n" | sort | uniq | tr "\n" " " | sed "s/^ //g" | sed "s/ $//g")
	fi;
	#if [ "$PIPELINES_INPUT" != "" ]; then # defined in INPUT
	#	PIPELINES=$PIPELINES_INPUT
	#fi;
	#if [ "$PIPELINES" == "" ]; then # default if not defined in ENV
	#	PIPELINES=$PIPELINES_DEFAULT;
	#	echo "#[WARNING] NO PIPELINES defined. Default PIPELINES '$PIPELINES' will be used."
	#fi;
	PIPELINES=$(echo $PIPELINES | tr -d "\"")
	
fi;

if ((0)); then
	# ALIGNERS
	if [ "$ALIGNERS_INPUT" == "" ]; then
		if [ "$ALIGNERS" == "" ]; then
			ALIGNERS="bwamem bwamemUnclipped" # bwasw bwamemUnclipped
			echo "#[WARNING] NO ALIGNERS defined. Default ALIGNERS '$ALIGNERS' will be used."
		fi;
	else
		ALIGNERS=$ALIGNERS_INPUT # INPUT
	fi;

	# CALLERS
	if [ "$CALLERS_INPUT" == "" ]; then
		if [ "$CALLERS" == "" ]; then
			CALLERS="gatkHC gatkUG"
			echo "#[WARNING] NO CALLERS defined. Default CALLERS '$CALLERS' will be used."
		fi;
	else
		CALLERS=$CALLERS_INPUT # INPUT
	fi;

	# ANNOTATORS
	if [ "$ANNOTATORS_INPUT" == "" ]; then
		if [ "$ANNOTATORS" == "" ]; then
			ANNOTATORS="howard"
			echo "#[WARNING] NO ANNOTATORS defined. Default ANNOTATORS '$ANNOTATORS' will be used."
		fi;
	else
		ANNOTATORS=$ANNOTATORS_INPUT # INPUT
	fi;
fi;

if [ "BAM_METRICS" == "" ]; then
	 BAM_METRICS=1
	echo "#[WARNING] No BAM_METRICS option defined. Default METRICS ('$BAM_METRICS') will be used."
fi;

#if [ "CLIPPING" == "" ]; then
#	 CLIPPING=1
#	echo "#[WARNING] No CLIPPING option defined. Default CLIPPING ('$CLIPPING') will be used."
#fi;

if [ "VARIANT_RECALIBRATION" == "" ]; then
	VARIANT_RECALIBRATION=0
	echo "#[WARNING] No VARIANT_RECALIBRATION option defined. Default VARIANT_RECALIBRATION ('$VARIANT_RECALIBRATION') will be used."
fi;

if [ "VARANK_ANALYSIS" == "" ]; then
	VARANK_ANALYSIS=0
	echo "#[WARNING] No VARANK_ANALYSIS option defined. Default VARANK_ANALYSIS ('$VARANK_ANALYSIS') will be used."
fi;

if [ "INTERVAL_PADDING" == "" ]; then
	INTERVAL_PADDING=0
	echo "#[WARNING] No INTERVAL_PADDING option defined. Default INTERVAL_PADDING ('$INTERVAL_PADDING') will be used."
fi;

if [ "PRIORITIZE_PIPELINES_LIST" == "" ]; then
	 PRIORITIZE_PIPELINES_LIST=""
	echo "#[WARNING] No PRIORITIZE_PIPELINES_LIST option defined. Default PRIORITIZE_PIPELINES_LIST ('$PRIORITIZE_PIPELINES_LIST') will be used."
fi;

#echo "["`date '+%Y%m%d-%H%M%S'`"] ALIGNERS '$ALIGNERS', CALLERS '$CALLERS',  ANNOTATORS '$ANNOTATORS', BAM_METRICS '$BAM_METRICS', CLIPPING '$CLIPPING',VARIANT_RECALIBRATION '$VARIANT_RECALIBRATION', INTERVAL_PADDING '$INTERVAL_PADDING', PRIORITIZE_PIPELINES_LIST '$PRIORITIZE_PIPELINES_LIST''"
echo "["`date '+%Y%m%d-%H%M%S'`"] ALIGNERS '$ALIGNERS', CALLERS '$CALLERS',  ANNOTATORS '$ANNOTATORS',  PIPELINES '$PIPELINES', BAM_METRICS '$BAM_METRICS',VARIANT_RECALIBRATION '$VARIANT_RECALIBRATION', INTERVAL_PADDING '$INTERVAL_PADDING', PRIORITIZE_PIPELINES_LIST '$PRIORITIZE_PIPELINES_LIST''"

# PREVIEW
PREVIEW=0

# USE_TMP_FOLDER
if [ "$USE_TMP_FOLDER" == "" ]; then
	USE_TMP_FOLDER=1
fi;
if [ "$TMP_FOLDER" == "$MAIN_FOLDER" ] || [ "$TMP_FOLDER" == "$RAW_FOLDER" ]; then
	USE_TMP_FOLDER=0
fi;

#echo "DATABASE_INTEGRATION_SCRIPT0=$DATABASE_INTEGRATION_SCRIPT";
# DATABASE_INTEGRATION_SCRIPT
if [ "$DATABASE_INTEGRATION_SCRIPT" == "" ] || [ ! -s "$DATABASE_INTEGRATION_SCRIPT" ]; then
	DATABASE_INTEGRATION_SCRIPT="";
fi;

#echo "DATABASE_INTEGRATION_SCRIPT1=$DATABASE_INTEGRATION_SCRIPT";
#echo "USE_TMP_FOLDER=$USE_TMP_FOLDER"; exit 0;


# 1.4. Variables
##################

# 4.1.0 DATE

DATE_DAY=`date '+%Y%m%d'`
DATE_MIN=`date '+%Y%m%d-%H%M%S'`


# 1.4.1 CORES & Threads

re='^[0-9]+$'

# Cores to let free
if ! [[ $CORES_FREE =~ $re ]] || [ "$CORES_FREE" == "" ]; then
	CORES_FREE=1
fi;

# Cores to use
if ! [[ $CORES_TO_USE =~ $re ]] || [ "$CORES_TO_USE" == "" ]; then

	CORES=$(ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w)	# NB of cores in the server
	CORES_FREE=1												# Number of threads free for other command
	CORES_TO_USE=$(($CORES-$CORES_FREE))					# Nb of cores to use

fi;

# Threads
THREADS=$CORES_TO_USE


# 4.1.2. COMMAND COPY

COMMAND_COPY="rsync -aucqpAXoghi" # "cp -auv" or "rsync -auv" # auvpAXog
# while ! rsync auczqpAXoghi $SRC/* $DEST ; do : ; done;

# 4.1.3. Permissions

PERMS="a+rwx"


# 1.5. Log
###########

echo "["`date '+%Y%m%d-%H%M%S'`"] *** Start Analysis '$DATE_MIN'"
echo "["`date '+%Y%m%d-%H%M%S'`"] * Scripts Configuration "
echo "["`date '+%Y%m%d-%H%M%S'`"] Bin folder		'$NGS_BIN'"
echo "["`date '+%Y%m%d-%H%M%S'`"] * Folder Configuration"
echo "["`date '+%Y%m%d-%H%M%S'`"] Raw data folder 	'$MISEQ_FOLDER'"
echo "["`date '+%Y%m%d-%H%M%S'`"] Demultiplexing folder 	'$DEMULTIPLEXING_FOLDER'"
echo "["`date '+%Y%m%d-%H%M%S'`"] Results folder 		'$RESULTS_FOLDER'"
echo "["`date '+%Y%m%d-%H%M%S'`"] Analysis folder 	'$ANALYSIS_FOLDER'"



# 1.6. Temporary Folder
#######################

if [ $USE_TMP_FOLDER -gt 0 ]; then
	echo "["`date '+%Y%m%d-%H%M%S'`"] Temporary space used	'$TMP_FOLDER'"
	if [ "$TMP_FOLDER" == "" ]; then
		TMP_DIR=/tmp
	else
		TMP_DIR=$TMP_FOLDER
	fi;
	TMP_DEM_DIR=$TMP_DIR/DEM
	TMP_RES_DIR=$TMP_DIR/RES/ALL
	TMP_ANA_DIR=$TMP_DIR/ANA
	mkdir -p $TMP_DIR
	mkdir -p $TMP_DEM_DIR
	mkdir -p $TMP_RES_DIR
	mkdir -p $TMP_ANA_DIR
else
	TMP_DEM_DIR=$DEMULTIPLEXING_FOLDER
	TMP_RES_DIR=$RESULTS_FOLDER
	TMP_ANA_DIR=$ANALYSIS_FOLDER
fi;
mkdir -p $TMP_DEM_DIR
mkdir -p $TMP_RES_DIR
mkdir -p $TMP_ANA_DIR

# 1.7. FOLDERS and FILES
#########################

ANALYSIS_DIR=$TMP_ANA_DIR
MAKEFILE_ANALYSIS=$ANALYSIS_DIR/run_analysis.V$DATE_MIN.param.mk
FINAL_REPORT=$ANALYSIS_DIR/run_analysis.V$DATE_MIN.report
FINAL_REPORT_FULL=$ANALYSIS_DIR/run_analysis.V$DATE_MIN.full.report
RELEASE=$ANALYSIS_DIR/run_analysis.V$DATE_MIN.release
FINAL_REPORT_VALIDATION=$ANALYSIS_DIR/run_analysis.V$DATE_MIN.validation.report
LOGFILE_ANA=$ANALYSIS_DIR/run_analysis.V$DATE_MIN.log
LOGFILE_CLEANING=$ANALYSIS_DIR/run_analysis.V$DATE_MIN.cleaning.log
LOGFILE_VALIDATION=$ANALYSIS_DIR/run_analysis.V$DATE_MIN.validation.log
SHELL_LIST=
MAKEFILE_LIST=
DATABASE_INTEGRATION=0


# 1.5. Temporary folder cleaning
################################
# Use date and/or space left...
#/media/data/NGS/e107_2.0/e107_plugins/NGS2_plugin/tmp
if [ "$TMP_REMOVE_TIME" == "" ]; then
	#TMP_REMOVE_TIME=10 # 10 days by default
	TMP_REMOVE_TIME=2 # 10 days by default
fi;
if [ $USE_TMP_FOLDER -gt 0 ]; then
	echo "["`date '+%Y%m%d-%H%M%S'`"] Temporary space cleaning..."

	echo "["`date '+%Y%m%d-%H%M%S'`"] Temporary space cleaning..." > $LOGFILE_CLEANING

	# Remove DEM folder
	#echo "find $TMP_DIR/DEM -maxdepth 1 -mtime +$TMP_REMOVE_TIME -type d -delete"
	find $TMP_DIR/DEM -maxdepth 1 -mtime +$TMP_REMOVE_TIME -type d | xargs rm -Rf 1>> $LOGFILE_CLEANING 2>> $LOGFILE_CLEANING
	# Remove ANA folder
	#echo "find $TMP_DIR/ANA -maxdepth 1 -mtime +$TMP_REMOVE_TIME -type d -delete"
	find $TMP_DIR/ANA -maxdepth 1 -mtime +$TMP_REMOVE_TIME -type d | xargs rm -Rf 1>> $LOGFILE_CLEANING 2>> $LOGFILE_CLEANING
	# Remove RES folder
	#echo "find $TMP_DIR/RES -mindepth 2 -maxdepth 2 -mtime +$TMP_REMOVE_TIME -type d -delete"
	find $TMP_DIR/RES -mindepth 2 -maxdepth 2 -mtime +$TMP_REMOVE_TIME -type d | xargs rm -Rf 1>> $LOGFILE_CLEANING 2>> $LOGFILE_CLEANING

	## DEBUG
	NOT_NAME="";
	for RUN in $RUNS;
	do
		NOT_NAME=$NOT_NAME" -not -name $RUN "
	done;
	find $TMP_DIR/DEM/ -mindepth 1 -maxdepth 1 -type d $NOT_NAME | xargs rm -Rf 1>> $LOGFILE_CLEANING 2>> $LOGFILE_CLEANING
	find $TMP_DIR/RES/ALL/ -mindepth 1 -maxdepth 1 -type d $NOT_NAME | xargs rm -Rf 1>> $LOGFILE_CLEANING 2>> $LOGFILE_CLEANING

	# ERROR
	if [ "`grep ^ $LOGFILE_CLEANING -c`" != "1" ]; then
		echo "["`date '+%Y%m%d-%H%M%S'`"] [ERROR] Cleaning finished with errors. See '$LOGFILE_CLEANING'";
		#exit 1;
	fi;

fi;




#################################
# 2. Preprocessing for each RUN #
#################################

for RUN in $RUNS;
do

	# Preprocessing
	echo "["`date '+%Y%m%d-%H%M%S'`"] * Preprocessing for RUN '$RUN'"

	# Test demultiplexing parallelized using makefile: #make RUNS=$RUN -f $NGS_SCRIPTS/Demultiplexing.mk

	# FOLDERS and FILES
	INPUT_DIR=$MISEQ_FOLDER/$RUN/$MSR_SUBFOLDER
	OUTPUT_DEM_DIR=$TMP_DEM_DIR/$RUN # Force to use tmp folder for demultiplexing. Original: OUTPUT_DEM_DIR=$DEMULTIPLEXING_FOLDER/$RUN
	OUTPUT_RES_DIR=$TMP_RES_DIR/$RUN # Force to use tmp folder for analysis. Original: RES_DIR=$RESULTS_FOLDER/$RUN
	# Create output folder
	mkdir -p $OUTPUT_DEM_DIR
	chmod $PERMS -R $OUTPUT_DEM_DIR 1>/dev/null 2>/dev/null
	mkdir -p $OUTPUT_RES_DIR
	chmod $PERMS -R $OUTPUT_RES_DIR 1>/dev/null 2>/dev/null

	# TMP DIR

	# 2.1. Check if RUN exists (Folder, SampleSheet)
	#################################################

	# Test if RUN folder exists in RAW data
	if [ ! -d $MISEQ_FOLDER/$RUN ]; then
		echo "["`date '+%Y%m%d-%H%M%S'`"] ERROR: RUN '$RUN' missing! Folder '$MISEQ_FOLDER/$RUN' doesn't exist"
		exit 1;
	fi

	# Test if SampleSheet.txt file exists in RAW data
	#if [ -z $SAMPLESHEET ]; then
	SAMPLE_SHEET_MISEQ_ORIGINAL=$MISEQ_FOLDER/$RUN/SampleSheet.csv
	#fi;
	#Data\Intensities\BaseCalls\Alignment
	SAMPLE_SHEET_MISEQ_EXISTS=0
	if [ -e $SAMPLESHEET ]; then
		SAMPLE_SHEET_MISEQ_ORIGINAL=$SAMPLESHEET
		SAMPLE_SHEET_MISEQ_EXISTS=1
	elif [ -e $MISEQ_FOLDER/$RUN/SampleSheet.csv ]; then
		SAMPLE_SHEET_MISEQ_ORIGINAL=$MISEQ_FOLDER/$RUN/SampleSheet.csv
		SAMPLE_SHEET_MISEQ_EXISTS=1
	elif [ -e $MISEQ_FOLDER/$RUN/Data/Intensities/BaseCalls/Alignment/SampleSheetUsed.csv ]; then
		SAMPLE_SHEET_MISEQ_ORIGINAL=$MISEQ_FOLDER/$RUN/Data/Intensities/BaseCalls/Alignment/SampleSheetUsed.csv
		SAMPLE_SHEET_MISEQ_EXISTS=1
	fi;
	if [ $SAMPLE_SHEET_MISEQ_EXISTS -eq 0 ]; then
		echo "["`date '+%Y%m%d-%H%M%S'`"] ERROR: SampleSheet.csv for RUN '$RUN' missing! File '$SAMPLE_SHEET_MISEQ_ORIGINAL' doesn't exist"
		exit 1;
	fi

	# Test if RUN is completed in RAW data (TODO)

	# 2.2. Check previous analysis, if TMP folder (TODO-TEST)
	##############################################

	#DEM_FINAL_FILE="$DEMULTIPLEXING_FOLDER/$RUN/Basecall_Stats_$RUN/Demultiplex_Stats.htm"
	DEM_FINAL_FILE=$OUTPUT_DEM_DIR/STARKComplete.txt

	if [ -e $DEM_FINAL_FILE ]; then
		echo "["`date '+%Y%m%d-%H%M%S'`"] RUN '$RUN' already demultiplexed"
		OUTPUT_DEM_DIR=$DEMULTIPLEXING_FOLDER/$RUN
		SAMPLE_SHEET_MISEQ=$OUTPUT_DEM_DIR/SampleSheet.csv
	else

		echo "["`date '+%Y%m%d-%H%M%S'`"] Start Demultiplexing, File '$DEM_FINAL_FILE' doesn't exist" 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA

		if [ $USE_TMP_FOLDER -gt 0 ]; then

			# 2.2.1. Send Demultiplexing data to DEM

			if [ -d $OUTPUT_DEM_DIR/$RUN ]; then #if [ ! -d $TMP_DEM_DIR ]; then # Synchronize if folder exists
				echo "["`date '+%Y%m%d-%H%M%S'`"] RUN '$RUN' files for Demultiplexing available..."
			else
				if [ -d $DEMULTIPLEXING_FOLDER/$RUN ]; then #if [ ! -d $TMP_DEM_DIR ]; then # Synchronize if folder exists
					echo "["`date '+%Y%m%d-%H%M%S'`"] Copying RUN '$RUN' files from Demultiplexing folder to Temporary folder..."
				 	#$COMMAND_COPY $DEMULTIPLEXING_FOLDER/$RUN $TMP_DEM_DIR 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA #&
				 	while ! $COMMAND_COPY $DEMULTIPLEXING_FOLDER/$RUN $TMP_DEM_DIR 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA ; do : ; done;
				fi;
			fi;

		fi;



		# FOLDERS and FILES
		SAMPLE_SHEET_MISEQ_ORIGINAL=$MISEQ_FOLDER/$RUN/SampleSheet.csv
		LOGFILE_DEM=$OUTPUT_DEM_DIR/demultiplexing.V$DATE_MIN.log
		SAMPLE_SHEET_MISEQ=$OUTPUT_DEM_DIR/SampleSheet.csv
		SAMPLE_SHEET_CASAVA=$OUTPUT_DEM_DIR/SampleSheet.casava.csv
		READS_LENGTH=$OUTPUT_DEM_DIR/readsLength.txt
		MASK=$OUTPUT_DEM_DIR/mask.txt
		ADAPTERS=$OUTPUT_DEM_DIR/adapters.txt
		MAKEFILE=$OUTPUT_DEM_DIR/Makefile
		DEM_STOP=$OUTPUT_DEM_DIR/STARKComplete.txt

		# 2.3. Demultiplexing (if needed)
		##################################

		# IF DEM/RUN Doesn't exists THEN DO Demultiplexing
		# IF DEM/RUN InComplete THEN Copy Files to TMP and DO Demultiplexing
		# IF DEM/RUN Complete THEN DO NOT Demultiplex and Switch DEM_SOURCE to DEM/RUN (for the Analysis Step)
		# Need RunInfo.xml in the root folder of the RUN !!! TODO

		

		# Demultiplexing Step
		#if [ ! -e $MAKEFILE ] || [ ! -e $OUTPUT_DEM_DIR/DemultiplexedBustardSummary.xml ];
		if [ ! -e $MAKEFILE ] || [ ! -e $DEM_FINAL_FILE ];
		then

			echo "["`date '+%Y%m%d-%H%M%S'`"] Demultiplexing RUN '$RUN'"

			# 2.3.1. Demultiplexing configuration
			if [ ! -e $MAKEFILE ];
			then
				echo "["`date '+%Y%m%d-%H%M%S'`"] Demultiplexing RUN '$RUN' Configuration (see '$LOGFILE_DEM')..." # Create $MAKEFILE
				$NGS_SCRIPTS/demultiplexing.sh "$ENV" "$RUN" "$TMP_DEM_DIR" "" "$SAMPLESHEET" 1>>$LOGFILE_DEM 2>>$LOGFILE_DEM
			fi

			# 2.3.2. Demultiplexing processing
			#if  [ ! -e $OUTPUT_DEM_DIR/DemultiplexedBustardSummary.xml ]; #DemultiplexedBustardSummary.xml
			#then
			#	echo "["`date '+%Y%m%d-%H%M%S'`"] Demultiplexing RUN '$RUN' Process (see '$LOGFILE_DEM')..."
			#	make -j $THREADS --directory=$OUTPUT_DEM_DIR -f $MAKEFILE 1>>$LOGFILE_DEM 2>>$LOGFILE_DEM  # Create $OUTPUT_DEM_DIR/DemultiplexedBustardSummary.xml
			#fi

		fi

		# 2.4. Send Demultiplexing data to DEM
		#######################################
		if [ -d $TMP_DEM_DIR ] && [ "$TMP_DEM_DIR" != "$DEMULTIPLEXING_FOLDER" ]; then
			echo "["`date '+%Y%m%d-%H%M%S'`"] Copying RUN '$RUN' files from Temporary folder '$TMP_DEM_DIR/*' to Demultiplexing folder '$DEMULTIPLEXING_FOLDER'..."
			#$COMMAND_COPY $TMP_DEM_DIR/$RUN $DEMULTIPLEXING_FOLDER 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA #&
			while ! $COMMAND_COPY $TMP_DEM_DIR/$RUN $DEMULTIPLEXING_FOLDER 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA ; do : ; done;
		fi;


	fi; # DEM_FINAL_FILE exists!

	echo "["`date '+%Y%m%d-%H%M%S'`"] Demultiplexing RUN folder '$OUTPUT_DEM_DIR'"


	# 2.5. RUN analysis configuration
	##################################

	# Config files
	MAKEFILE_ANALYSIS_RUN=$OUTPUT_RES_DIR/run_analysis.V$DATE_MIN.param.mk
	SHELL_ANALYSIS_RUN=$OUTPUT_RES_DIR/run_analysis.V$DATE_MIN.param.sh

	# Analysis Configuration Step
	echo "["`date '+%Y%m%d-%H%M%S'`"] Analysis Parameters Configuration (see $MAKEFILE_ANALYSIS_RUN)..."

	# 2.5.1. Parameters created from SampleSeet TODO

	RUNS_SAMPLES=""
	NB_SAMPLE=0
	RUN_Investigator_Name=$(grep "Investigator Name" $SAMPLE_SHEET_MISEQ | tr -d "\r\n" | awk -F, '{print $2}')
	DATA_SECTION_LINE=$(grep "^\[Data\]" $SAMPLE_SHEET_MISEQ -n | awk -F: '{print $1}')
	SAMPLE_SECTION_FIRST_LINE=$(($DATA_SECTION_LINE+1))
	SAMPLE_LINES=$(tail -n $(($SAMPLE_SECTION_FIRST_LINE-$(grep ^ -c $SAMPLE_SHEET_MISEQ)))  $SAMPLE_SHEET_MISEQ)
 	SAMPLE_PROJECT_COL=$(grep -i ^Sample_ID $SAMPLE_SHEET_MISEQ | tr -d '\r\n' | sed -e 's/,/\n/g' | grep -i -n Sample_Project | cut -d \: -f 1)
	for L in $(tail -n $(($SAMPLE_SECTION_FIRST_LINE-$(grep ^ -c $SAMPLE_SHEET_MISEQ)))  $SAMPLE_SHEET_MISEQ);
	do
		if [ "$L" == "" ] || [[ $L =~ ^\[.*\] ]]; then
			break;
		else
			((NB_SAMPLE++))
			SAMPLE_ID=$(echo $L | awk -F, '{print $1}')
			SAMPLE_PROJECT=$(echo $L | tr -d '\r\n' | cut -d \, -f $SAMPLE_PROJECT_COL )
			if [ "$SAMPLE_PROJECT" == "" ]; then SAMPLE_PROJECT=$RUN_Investigator_Name; fi;
			RUNS_SAMPLES="$RUNS_SAMPLES$RUN:$SAMPLE_ID:$SAMPLE_PROJECT "
		fi;
	done;

	# Additional parameters to Makefile Header
	echo "" >> $MAKEFILE_ANALYSIS_RUN
	echo "#############" >> $MAKEFILE_ANALYSIS_RUN
	echo "# Parameters #" >> $MAKEFILE_ANALYSIS_RUN
	echo "#############" >> $MAKEFILE_ANALYSIS_RUN
	echo "" >> $MAKEFILE_ANALYSIS_RUN
	echo "RUNS=$RUN" >> $MAKEFILE_ANALYSIS_RUN
	echo "RUNS_SAMPLES=$RUNS_SAMPLES" >> $MAKEFILE_ANALYSIS_RUN
	echo "" >> $MAKEFILE_ANALYSIS_RUN

	# Additional parameters to Shell file Header
	echo "" >> $SHELL_ANALYSIS_RUN
	echo "#########################" >> $SHELL_ANALYSIS_RUN
	echo "# Additional parameters #" >> $SHELL_ANALYSIS_RUN
	echo "#########################" >> $SHELL_ANALYSIS_RUN
	echo "" >> $SHELL_ANALYSIS_RUN
	echo "RUNS=\"$RUN\"" >> $SHELL_ANALYSIS_RUN
	echo "RUNS_SAMPLES=\"$RUNS_SAMPLES\"" >> $SHELL_ANALYSIS_RUN
	echo "" >> $SHELL_ANALYSIS_RUN



	# TODO!!!


	# 2.4.1. Makefile and shell from xml
	#xsltproc $NGS_SCRIPTS/DemultiplexConfig.xml2sh.param.xsl $OUTPUT_DEM_DIR/DemultiplexConfig.xml > $SHELL_ANALYSIS_RUN # Create $MAKEFILE_ANALYSIS
	#xsltproc $NGS_SCRIPTS/DemultiplexConfig.xml2mk.param.xsl $OUTPUT_DEM_DIR/DemultiplexConfig.xml > $MAKEFILE_ANALYSIS_RUN # Create $MAKEFILE_ANALYSIS

	SHELL_LIST="$SHELL_LIST $SHELL_ANALYSIS_RUN"
	MAKEFILE_LIST="$MAKEFILE_LIST $MAKEFILE_ANALYSIS_RUN"

	# ERROR
	if [ ! -e $MAKEFILE_ANALYSIS_RUN ] || [ ! -e $SHELL_ANALYSIS_RUN ]; then
		echo "[ERROR] Files '$MAKEFILE_ANALYSIS_RUN' or '$SHELL_ANALYSIS_RUN' do NOT exist!!!";
		exit 1;
	fi;

	# 2.4.2. Filter by sample (TODO)
	if [ "$FILTER_SAMPLES" != "" ]; then
		MAKEFILE_ANALYSIS_RUN_HEADER=`cat $MAKEFILE_ANALYSIS_RUN | grep ^RUNS_SAMPLES -v`
		SHELL_ANALYSIS_RUN_HEADER=`cat $SHELL_ANALYSIS_RUN | grep ^RUNS_SAMPLES -v`
		for FILTER_SAMPLE in $FILTER_SAMPLES;
		do
			RUNS_SAMPLES_MAKEFILE_VARIABLE="$RUNS_SAMPLES_MAKEFILE_VARIABLE "`cat $MAKEFILE_ANALYSIS_RUN | grep ^RUNS_SAMPLES  | tr -d \" | awk -F"=" '{print $2}' | tr -s " " "\n" | grep :$FILTER_SAMPLE:`
			RUNS_SAMPLES_SHELL_VARIABLE="$RUNS_SAMPLES_SHELL_VARIABLE "`cat $SHELL_ANALYSIS_RUN | grep ^RUNS_SAMPLES  | tr -d \" | awk -F"=" '{print $2}' | tr -s " " "\n" | grep :$FILTER_SAMPLE:`
		done;
		echo -e "$MAKEFILE_ANALYSIS_RUN_HEADER" "\n" "RUNS_SAMPLES=$RUNS_SAMPLES_MAKEFILE_VARIABLE" | sed 's/ RUNS_SAMPLES= /RUNS_SAMPLES=/' > $MAKEFILE_ANALYSIS_RUN
		echo -e "$SHELL_ANALYSIS_RUN_HEADER" "\n" "RUNS_SAMPLES=\"$RUNS_SAMPLES_SHELL_VARIABLE\"" |  sed 's/ RUNS_SAMPLES=\" /RUNS_SAMPLES="/' > $SHELL_ANALYSIS_RUN
	fi;

	# 2.4.3. Additional parameters to Makefile
	echo "" >> $MAKEFILE_ANALYSIS_RUN
	echo "#########################" >> $MAKEFILE_ANALYSIS_RUN
	echo "# Additional parameters #" >> $MAKEFILE_ANALYSIS_RUN
	echo "#########################" >> $MAKEFILE_ANALYSIS_RUN
	echo "" >> $MAKEFILE_ANALYSIS_RUN
	echo "RUNS=$RUN" >> $MAKEFILE_ANALYSIS_RUN
	echo "ALIGNERS=$ALIGNERS" >> $MAKEFILE_ANALYSIS_RUN
	echo "CALLERS=$CALLERS" >> $MAKEFILE_ANALYSIS_RUN
	echo "ANNOTATORS=$ANNOTATORS" >> $MAKEFILE_ANALYSIS_RUN
	echo "INTERSEC=$INTERSEC" >> $MAKEFILE_ANALYSIS_RUN
	echo "BAM_METRICS=\"$BAM_METRICS\"" >> $MAKEFILE_ANALYSIS_RUN
	echo "CLIPPING=\"$CLIPPING\"" >> $MAKEFILE_ANALYSIS_RUN
	echo "VARIANT_RECALIBRATION=\"$VARIANT_RECALIBRATION\"" >> $MAKEFILE_ANALYSIS_RUN
	echo "INTERVAL_PADDING=\"$INTERVAL_PADDING\"" >> $MAKEFILE_ANALYSIS_RUN
	echo "COVERAGE_CRITERIA=\"$COVERAGE_CRITERIA\"" >> $MAKEFILE_ANALYSIS_RUN
	echo "NB_BASES_AROUND=\"$NB_BASES_AROUND\"" >> $MAKEFILE_ANALYSIS_RUN
	echo "BEDFILE_GENES=\"$BEDFILE_GENES\"" >> $MAKEFILE_ANALYSIS_RUN
	echo "PRIORITIZE_PIPELINES_LIST=\"$PRIORITIZE_PIPELINES_LIST\"" >> $MAKEFILE_ANALYSIS_RUN
	echo "" >> $MAKEFILE_ANALYSIS_RUN

	# 2.4.3. Additional parameters to Shell file
	echo "" >> $SHELL_ANALYSIS_RUN
	echo "#########################" >> $SHELL_ANALYSIS_RUN
	echo "# Additional parameters #" >> $SHELL_ANALYSIS_RUN
	echo "#########################" >> $SHELL_ANALYSIS_RUN
	echo "" >> $SHELL_ANALYSIS_RUN
	echo "RUNS=\"$RUN\"" >> $SHELL_ANALYSIS_RUN
	echo "ALIGNERS=\"$ALIGNERS\"" >> $SHELL_ANALYSIS_RUN
	echo "CALLERS=\"$CALLERS\"" >> $SHELL_ANALYSIS_RUN
	echo "ANNOTATORS=\"$ANNOTATORS\"" >> $SHELL_ANALYSIS_RUN
	echo "INTERSEC=\"$INTERSEC\"" >> $SHELL_ANALYSIS_RUN
	echo "BAM_METRICS=\"$BAM_METRICS\"" >> $SHELL_ANALYSIS_RUN
	echo "CLIPPING=\"$CLIPPING\"" >> $SHELL_ANALYSIS_RUN
	echo "VARIANT_RECALIBRATION=\"$VARIANT_RECALIBRATION\"" >> $SHELL_ANALYSIS_RUN
	echo "INTERVAL_PADDING=\"$INTERVAL_PADDING\"" >>  $SHELL_ANALYSIS_RUN
	echo "COVERAGE_CRITERIA=\"$COVERAGE_CRITERIA\"" >> $SHELL_ANALYSIS_RUN
	echo "NB_BASES_AROUND=\"$NB_BASES_AROUND\"" >> $SHELL_ANALYSIS_RUN
	echo "BEDFILE_GENES=\"$BEDFILE_GENES\"" >> $SHELL_ANALYSIS_RUN
	echo "PRIORITIZE_PIPELINES_LIST=\"$PRIORITIZE_PIPELINES_LIST\"" >> $SHELL_ANALYSIS_RUN
	echo "" >> $SHELL_ANALYSIS_RUN

done;

#################################
# 3. Mutli Run Analysis Process #
#################################

if [ $MULTI_RUN_ANALYSIS_PROCESS -gt 0 ]; # Multi run process, analysis and/or report
then

	echo "["`date '+%Y%m%d-%H%M%S'`"] * Multi Run Process"

	# 3.1. Combination of all RUN's Makefile
	#########################################

	echo "["`date '+%Y%m%d-%H%M%S'`"] Multi Run Parameters Combination Process (see '$MAKEFILE_ANALYSIS') ..."

	# Combination of all RUNs and SAMPLEs
	RUNS_SAMPLES_ALL=
	RUNS_ALL=
	for RUN_PARAM in $SHELL_LIST;
	do
		source $RUN_PARAM
		rm $RUN_PARAM
		RUNS_SAMPLES_ALL="$RUNS_SAMPLES_ALL $RUNS_SAMPLES"
		RUNS_ALL="$RUNS_ALL $RUNS"
	done;

	# Creation of a combined Makefile
	cat $MAKEFILE_LIST >> $MAKEFILE_ANALYSIS
	echo "" >> $MAKEFILE_ANALYSIS
	echo "##################################" >> $MAKEFILE_ANALYSIS
	echo "# Combination of all RUNS_SAMPLES #" >> $MAKEFILE_ANALYSIS
	echo "##################################" >> $MAKEFILE_ANALYSIS
	echo "" >> $MAKEFILE_ANALYSIS
	echo "RUNS_SAMPLES=$RUNS_SAMPLES_ALL" >> $MAKEFILE_ANALYSIS

	# Additional configuration to Makefile
	echo "" >> $MAKEFILE_ANALYSIS
	echo "#########################" >> $MAKEFILE_ANALYSIS
	echo "# Additional parameters #" >> $MAKEFILE_ANALYSIS
	echo "#########################" >> $MAKEFILE_ANALYSIS
	echo "" >> $MAKEFILE_ANALYSIS
	echo "RUNS=$RUNS" >> $MAKEFILE_ANALYSIS
	echo "ALIGNERS=$ALIGNERS" >> $MAKEFILE_ANALYSIS
	echo "CALLERS=$CALLERS" >> $MAKEFILE_ANALYSIS
	echo "ANNOTATORS=$ANNOTATORS" >> $MAKEFILE_ANALYSIS
	echo "INTERSEC=$INTERSEC" >> $MAKEFILE_ANALYSIS
	echo "" >> $MAKEFILE_ANALYSIS

	# 3.2. Multi Run Main Analysis Process
	#######################################

	if [ $MULTI_RUN_ANALYSIS_PROCESS -eq 1 ]; # Multi Run analysis 'full'
	then

		# 3.2.1. Define source of data

		RES_SOURCE=$RESULTS_FOLDER 	# Source of data by default
		NOTHINGTOBEDONE=0		# Something to be done by default
		# IF Folder RES exists
		if [ -d $RESULTS_FOLDER ];
		then
			
			#NOTHINGTOBEDONE=$(make -j $THREADS -e NGSEnv=$NGS_FOLDER ENV=$ENV PARAM=$MAKEFILE_ANALYSIS VALIDATION=1 OUTDIR=$RESULTS_FOLDER -f $NGS_SCRIPTS/NGSWorkflow.mk -n | grep 'Nothing to be done for `all' -c)
			#
			## IF Nothing to be done
			#if [ $NOTHINGTOBEDONE -eq 1 ];
			#then
			#	# Nothing to be done for this Analysis Process: Source of data is RES
			#	echo "["`date '+%Y%m%d-%H%M%S'`"] Nothing to be Done RUN for '$RUN' Analysis Process"
			#	RES_SOURCE=$RESULTS_FOLDER
			#else
				# Something need to be done: Source of data is TMP. Copy Files of ALL RUNs into TMP
				# Improvment: test if analysis need to be done for each RUN. TODO
				for RUN in $RUNS;
				do
					if [ $USE_TMP_FOLDER -gt 0 ]; then
						echo "["`date '+%Y%m%d-%H%M%S'`"] Copying RUN '$RUN' files from Results folder to Temporary folder..."
					 	#$COMMAND_COPY $RESULTS_FOLDER/$RUN $TMP_RES_DIR 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA #&
					 	while ! $COMMAND_COPY $RESULTS_FOLDER/$RUN $TMP_RES_DIR 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA ; do : ; done;
					fi;
				done;
				RES_SOURCE=$TMP_RES_DIR
			#fi;
		else
			# Folder RES/RUN DOESN't exists. Run Analysis on TMP folder
			RES_SOURCE=$TMP_RES_DIR
		fi;
		echo "["`date '+%Y%m%d-%H%M%S'`"] Results folder for Analysis is '$RES_SOURCE'"

		# 3.2.2. Main Analysis process, all RUNs together

		if [ $NOTHINGTOBEDONE -ne 1 ]; # If Something to be done
		then
			echo "["`date '+%Y%m%d-%H%M%S'`"] Multi Run Main Analysis Process (see '$LOGFILE_ANA') ..."
			make -j $THREADS -e NGSEnv=$NGS_FOLDER ENV=$ENV PARAM=$MAKEFILE_ANALYSIS VALIDATION=1  OUTDIR=$RES_SOURCE -f $NGS_SCRIPTS/NGSWorkflow.mk 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA
		fi;

	fi;

fi;


##################################
# 4. Postprocessing for each RUN #
##################################

# FOREACH RUN's Makefile
for MAKEFILE_ANALYSIS_RUN in $MAKEFILE_LIST;
do

	# 4.0 Postprocessing Step
	########################

	RUN=$(basename $(dirname $MAKEFILE_ANALYSIS_RUN))
	echo "["`date '+%Y%m%d-%H%M%S'`"] * Analysis of RUN '$RUN'"
	#echo "["`date '+%Y%m%d-%H%M%S'`"] ALIGNERS '$ALIGNERS', CALLERS '$CALLERS',  ANNOTATORS '$ANNOTATORS'"

	# 4.0.1. FOLDERS and FILES

	TMP_RES_DIR_RUN=$(dirname $MAKEFILE_ANALYSIS_RUN)
	FINAL_REPORT_RUN=$TMP_RES_DIR_RUN"/"$(basename $FINAL_REPORT)
	FINAL_REPORT_FULL_RUN=$TMP_RES_DIR_RUN"/"$(basename $FINAL_REPORT_FULL)
	RELEASE_RUN=$TMP_RES_DIR_RUN"/"$(basename $RELEASE)
	FINAL_REPORT_VALIDATION_RUN=$TMP_RES_DIR_RUN"/"$(basename $FINAL_REPORT_VALIDATION)
	LOGFILE_RES_RUN=$TMP_RES_DIR_RUN"/"$(basename $LOGFILE_ANA)
	LOGFILE_REP_RUN=$TMP_RES_DIR_RUN"/"$(basename $LOGFILE_ANA | sed "s/.log$/.report.log/gi")
	LOGFILE_VALIDATION_RUN=$TMP_RES_DIR_RUN"/"$(basename $LOGFILE_VALIDATION)
	#FINAL_REPORT_FULL_VCF=$TMP_RES_DIR_RUN/$RUN.vcf

	# 4.0.2. Start Log File

	echo "["`date '+%Y%m%d-%H%M%S'`"] Start Main Analysis" >$LOGFILE_RES_RUN

	# 4.0.3. ERROR

	if [ ! -e $LOGFILE_RES_RUN ]; then
		echo "[ERROR] File '$LOGFILE_RES_RUN' does NOT exist!!!";
		exit 1;
	fi;

	# 4.0.4. Thread calculation

	THREAD_PARAMETERS="" # Default in $ENV

	if [ "$THREADS_AUTO" == "1" ]; then

		# NB SAMPLE
		NB_SAMPLE=$(grep "^RUNS_SAMPLES=" $MAKEFILE_ANALYSIS_RUN | wc -w)
		#echo "NB_SAMPLE=$NB_SAMPLE"

		# NB ALIGNERS
		NB_ALIGNERS=$(grep "^ALIGNERS=" $MAKEFILE_ANALYSIS_RUN | wc -w)
		#echo "NB_ALIGNERS=$NB_ALIGNERS"

		# NB CALLERS
		NB_CALLERS=$(grep "^CALLERS=" $MAKEFILE_ANALYSIS_RUN | wc -w)
		#echo "NB_CALLERS=$NB_CALLERS"

		# NB ANNOTATORS
		NB_ANNOTATORS=$(grep "^ANNOTATORS=" $MAKEFILE_ANALYSIS_RUN | wc -w)
		#echo "NB_ANNOTATORS=$NB_ANNOTATORS"

		# NB SAMPLE
		# CREATE PIPELINES from ALIGNERS/CALLERS/ANNOTATORS
		PIPELINES=$(grep "^PIPELINES=" $MAKEFILE_ANALYSIS_RUN | cut -d= -f2)
		for ALIGNER in $ALIGNERS; do
			for CALLER in $CALLERS; do
				for ANNOTATOR in $ANNOTATORS; do
					PIPELINES="$PIPELINES $ALIGNER.$CALLER.$ANNOTATOR"
				done;
			done;
		done;
		PIPELINES=$(echo $PIPELINES | tr " " "\n" | sort | uniq | tr "\n" " ")
		NB_PIPELINES=$(echo $PIPELINES | wc -w)
		#echo "NB_PIPELINES=$NB_PIPELINES"

		# ADJUSTING threads parameters
		THREADS_BY_SAMPLE=$( echo "scale=0;$CORES_TO_USE/$NB_SAMPLE" | bc )              # Number of threads by sample
             	THREADS_BY_PIPELINE=$( echo "scale=0;$THREADS_BY_SAMPLE/$NB_PIPELINES" | bc )    # Number of threads for a pipeline's command
            	THREADS_BY_ALIGNER=$( echo "scale=0;$THREADS_BY_SAMPLE/$NB_ALIGNERS" | bc )      # Number of threads for a aligner's command
             	THREADS_BY_CALLER=$( echo "scale=0;$THREADS_BY_SAMPLE/$NB_CALLERS" | bc )        # Number of threads for a caller's command
             	THREADS_BWA=$THREADS_BY_ALIGNER                                                  # NB of threads for BWA command
             	THREADS_RTC=$THREADS_BY_ALIGNER                                                  # NB of threads for GATK RealignerTargetCreator

		if [ "$THREADS_BY_SAMPLE" == "0" ]; then THREADS_BY_SAMPLE=1; fi;
		if [ "$THREADS_BY_ALIGNER" == "0" ]; then THREADS_BY_ALIGNER=1; fi;
		if [ "$THREADS_BY_CALLER" == "0" ]; then THREADS_BY_CALLER=1; fi;
		if [ "$THREADS_BY_PIPELINE" == "0" ]; then THREADS_BY_PIPELINE=1; fi;
		if [ "$THREADS_BWA" == "0" ]; then THREADS_BWA=1; fi;
		if [ "$THREADS_RTC" == "0" ]; then THREADS_RTC=1; fi;

		# THREAD_PARAMETERS
		THREAD_PARAMETERS=$THREAD_PARAMETERS" THREADS_BY_SAMPLE=$THREADS_BY_SAMPLE "
		THREAD_PARAMETERS=$THREAD_PARAMETERS" THREADS_BY_PIPELINE=$THREADS_BY_PIPELINE "
		THREAD_PARAMETERS=$THREAD_PARAMETERS" THREADS_BY_ALIGNER=$THREADS_BY_ALIGNER "
		THREAD_PARAMETERS=$THREAD_PARAMETERS" THREADS_BY_CALLER=$THREADS_BY_CALLER "
		THREAD_PARAMETERS=$THREAD_PARAMETERS" THREADS_BWA=$THREADS_BWA "
		THREAD_PARAMETERS=$THREAD_PARAMETERS" THREADS_RTC=$THREADS_RTC "

	fi;

	#echo $THREAD_PARAMETERS; exit 0;

	# 4.1. Main Analysis Process
	#############################

	# 4.1.0. Parameters

	VALIDATION=0
	SNAPSHOT=0
	#MAKE_DEBUG=" --debug=b "
	MAKE_DEBUG=" "

	# 4.1.1. Define source of data

	RES_SOURCE=$RESULTS_FOLDER # By default
	NOTHINGTOBEDONE=0
	# Folder RES/RUN exists
	#echo "RESULTS_FOLDER/RUN=$RESULTS_FOLDER/$RUN" >>$LOGFILE_RES_RUN
	if [ -d $RESULTS_FOLDER/$RUN ];
	then
		#TOBEDONE=$(make -j $THREADS -e NGSEnv=$NGS_FOLDER ENV=$ENV PARAM=$MAKEFILE_ANALYSIS_RUN OUTDIR=$RESULTS_FOLDER -f $NGS_SCRIPTS/NGSWorkflow.mk -n) #
		#echo -e $TOBEDONE > $TMP_RES_DIR_RUN"/"$(basename $LOGFILE_ANA).tobedone
		#NOTHINGTOBEDONE=$(echo -e $TOBEDONE | grep 'Nothing to be done for `all' -c)
		
		# IF Nothing to be done
		#if [ $NOTHINGTOBEDONE -eq 1 ];
		#then
		#	# Nothing to be done for this Analysis RUN Process: Source of data is RES
		#	echo "["`date '+%Y%m%d-%H%M%S'`"] Nothing to be Done RUN for '$RUN' Analysis Process"
		#	RES_SOURCE=$RESULTS_FOLDER
		#else
			# Something need to be done: Source of data is TMP. Copy Files of RUN into TMP
			if [ $USE_TMP_FOLDER -gt 0 ]; then
				echo "["`date '+%Y%m%d-%H%M%S'`"] Copying RUN '$RUN' files from Results folder to Temporary folder..."
			 	#$COMMAND_COPY $RESULTS_FOLDER/$RUN $TMP_RES_DIR 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA #&
			 	while ! $COMMAND_COPY $RESULTS_FOLDER/$RUN $TMP_RES_DIR 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA ; do : ; done;
				RES_SOURCE=$TMP_RES_DIR
			fi;
		#fi;
	else
		# Folder RES/RUN DOESN't exists. Run Analysis on TMP folder
		RES_SOURCE=$TMP_RES_DIR
	fi;
	echo "["`date '+%Y%m%d-%H%M%S'`"] Results folder for Analysis is '$RES_SOURCE'"

	# 4.1.2. Main Analysis Process, for the RUN

	if [ $NOTHINGTOBEDONE -ne 1 ]; # If Something to be done
	then

		echo "["`date '+%Y%m%d-%H%M%S'`"] Main Analysis Process for RUN '$(basename $(dirname $MAKEFILE_ANALYSIS_RUN))' (see '$LOGFILE_RES_RUN')..."
		#make -j $THREADS NGSEnv=$NGS_FOLDER PARAM=$MAKEFILE_ANALYSIS_RUN VALIDATION=$VALIDATION OUTDIR=$RES_SOURCE -f $NGS_SCRIPTS/NGSWorkflow.mk summary 1>>$LOGFILE_RES_RUN 2>>$LOGFILE_RES_RUN
		echo "["`date '+%Y%m%d-%H%M%S'`"] Main Analysis Process for RUN '$(basename $(dirname $MAKEFILE_ANALYSIS_RUN))' START" >>$LOGFILE_RES_RUN
		make -j $THREADS -e NGSEnv=$NGS_FOLDER ENV=$ENV PARAM=$MAKEFILE_ANALYSIS_RUN $THREAD_PARAMETERS OUTDIR=$RES_SOURCE RELEASE=$RELEASE_RUN SNAPSHOT=$SNAPSHOT -f $NGS_SCRIPTS/NGSWorkflow.mk $MAKE_DEBUG 1>>$LOGFILE_RES_RUN 2>>$LOGFILE_RES_RUN
		echo "["`date '+%Y%m%d-%H%M%S'`"] Main Analysis Process for RUN '$(basename $(dirname $MAKEFILE_ANALYSIS_RUN))' END" >>$LOGFILE_RES_RUN
		if [ $(grep "\*\*\*" $LOGFILE_RES_RUN -c) -gt 0 ];
		then
			echo "["`date '+%Y%m%d-%H%M%S'`"] ERROR: Main Analysis Process for RUN '$(basename $(dirname $MAKEFILE_ANALYSIS_RUN))'. (see '$LOGFILE_RES_RUN')";
			echo "["`date '+%Y%m%d-%H%M%S'`"] '`grep "\*\*\*" $LOGFILE_RES_RUN -m 1`'";
			if [ "$(dirname $LOGFILE_RES_RUN)" != "$RESULTS_FOLDER/$RUN" ]; then
				echo "["`date '+%Y%m%d-%H%M%S'`"] Copy Log files '$LOGFILE_RES_RUN' into '$RESULTS_FOLDER/$RUN' ";
				mkdir -p $RESULTS_FOLDER/$RUN
				#$COMMAND_COPY $LOGFILE_RES_RUN $RESULTS_FOLDER/$RUN 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA #&
				while ! $COMMAND_COPY $LOGFILE_RES_RUN $RESULTS_FOLDER/$RUN 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA ; do : ; done;
			fi;
			#exit 1;
			continue;
		else
			# RUN Analysis OK. Copy RUN
			RUN_ANALYZED=$RUN_ANALYZED" "$RUN
		fi
		# Database integration
		DATABASE_INTEGRATION=1
	else
		echo "["`date '+%Y%m%d-%H%M%S'`"] Main Analysis 'Nothing to be done' for RUN '$(basename $(dirname $MAKEFILE_ANALYSIS_RUN))' (see '$LOGFILE_RES_RUN')..."
		echo "["`date '+%Y%m%d-%H%M%S'`"] Main Analysis 'Nothing to be done' for RUN '$(basename $(dirname $MAKEFILE_ANALYSIS_RUN))' (see '$LOGFILE_RES_RUN')..." >>$LOGFILE_RES_RUN
		# RUN Analysis OK. Copy RUN
		RUN_ANALYZED=$RUN_ANALYZED" "$RUN
	fi;

	#echo "["`date '+%Y%m%d-%H%M%S'`"] * Postprocessing for RUN '$RUN'"

	# 4.2. Report process
	######################

	# Generate a repport and the release file
	#NOTHINGTOBEDONE=0
	#if [ $NOTHINGTOBEDONE -ne 1 ] || [ 1 -eq 1 ]; # If NOT Something to be done
	#then
		echo "["`date '+%Y%m%d-%H%M%S'`"] Report Process for RUN '$(basename $(dirname $MAKEFILE_ANALYSIS_RUN))' (see '$LOGFILE_REP_RUN')..."
		echo "["`date '+%Y%m%d-%H%M%S'`"] Report Process for RUN '$(basename $(dirname $MAKEFILE_ANALYSIS_RUN))' START" >>$LOGFILE_REP_RUN
		#echo "make -j $THREADS -e NGSEnv=$NGS_FOLDER ENV=$ENV PARAM=$MAKEFILE_ANALYSIS_RUN VALIDATION=$VALIDATION SNAPSHOT=$SNAPSHOT OUTDIR=$RES_SOURCE FINAL_REPORT=$FINAL_REPORT_RUN FINAL_REPORT_FULL=$FINAL_REPORT_FULL_RUN RELEASE=$RELEASE_RUN -f $NGS_SCRIPTS/NGSWorkflow.mk 1>>$LOGFILE_REP_RUN 2>>$LOGFILE_REP_RUN" >>$LOGFILE_REP_RUN
		make -j $THREADS -e NGSEnv=$NGS_FOLDER ENV=$ENV PARAM=$MAKEFILE_ANALYSIS_RUN VALIDATION=$VALIDATION SNAPSHOT=$SNAPSHOT OUTDIR=$RES_SOURCE FINAL_REPORT=$FINAL_REPORT_RUN FINAL_REPORT_FULL=$FINAL_REPORT_FULL_RUN RELEASE=$RELEASE_RUN ANALYSIS_DATE=$DATE_MIN -f $NGS_SCRIPTS/NGSWorkflow.mk $MAKE_DEBUG 1>>$LOGFILE_REP_RUN 2>>$LOGFILE_REP_RUN
		#echo "make -j $THREADS -e NGSEnv=$NGS_FOLDER ENV=$ENV PARAM=$MAKEFILE_ANALYSIS_RUN VALIDATION=$VALIDATION SNAPSHOT=$SNAPSHOT OUTDIR=$RES_SOURCE FINAL_REPORT=$FINAL_REPORT_RUN FINAL_REPORT_FULL=$FINAL_REPORT_FULL_RUN RELEASE=$RELEASE_RUN -f $NGS_SCRIPTS/NGSWorkflow.mk "
		echo "["`date '+%Y%m%d-%H%M%S'`"] Report Process for RUN '$(basename $(dirname $MAKEFILE_ANALYSIS_RUN))' END" >>$LOGFILE_REP_RUN
		
		
		if [ $(grep "\*\*\*" $LOGFILE_REP_RUN -c) -gt 0 ];
		then
			echo "["`date '+%Y%m%d-%H%M%S'`"] ERROR: Report Process for RUN '$(basename $(dirname $MAKEFILE_ANALYSIS_RUN))'. (see '$LOGFILE_REP_RUN')";
			echo "["`date '+%Y%m%d-%H%M%S'`"] '`grep "\*\*\*" $LOGFILE_REP_RUN -m 1`'";
			if [ "$(dirname $LOGFILE_REP_RUN)" != "$RESULTS_FOLDER/$RUN" ]; then
				echo "["`date '+%Y%m%d-%H%M%S'`"] Copy Log files '$LOGFILE_REP_RUN' into '$RESULTS_FOLDER/$RUN' ";
				mkdir -p $RESULTS_FOLDER/$RUN
				#$COMMAND_COPY $LOGFILE_RES_RUN $RESULTS_FOLDER/$RUN 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA #&
				while ! $COMMAND_COPY $LOGFILE_REP_RUN $RESULTS_FOLDER/$RUN 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA ; do : ; done;
			fi;
			#exit 1;
			# remove RUN from RUN_ANALYZED
			RUN_ANALYZED=$(echo $RUN_ANALYZED | sed s/$RUN//gi)
			continue;
		fi
		
		
	#else
	#	echo "["`date '+%Y%m%d-%H%M%S'`"] Report Process 'Nothing to be done' for RUN '$(basename $(dirname $MAKEFILE_ANALYSIS_RUN))' (see '$LOGFILE_REP_RUN')..."
	#	echo "["`date '+%Y%m%d-%H%M%S'`"] Report Process 'Nothing to be done' for RUN '$(basename $(dirname $MAKEFILE_ANALYSIS_RUN))' (see '$LOGFILE_REP_RUN')..." >>$LOGFILE_REP_RUN
	#fi;

	# 4.3. Validation process
	##########################

	#echo "["`date '+%Y%m%d-%H%M%S'`"] Validation Process (see '$LOGFILE_VALIDATION_RUN')..."
	#echo "["`date '+%Y%m%d-%H%M%S'`"] Validation Process for RUN '$(basename $(dirname $MAKEFILE_ANALYSIS_RUN))' START" >>$LOGFILE_RES_RUN
	#make -j $THREADS NGSEnv=$NGS_FOLDER PARAM=$MAKEFILE_ANALYSIS_RUN OUTDIR=$RES_SOURCE FINAL_REPORT_VALIDATION=$FINAL_REPORT_VALIDATION_RUN RELEASE=$RELEASE_RUN -f $NGS_SCRIPTS/NGSWorkflow.validation.mk 1>>$LOGFILE_VALIDATION_RUN 2>>$LOGFILE_VALIDATION_RUN #--debug --
	#echo "["`date '+%Y%m%d-%H%M%S'`"] Validation Process for RUN '$(basename $(dirname $MAKEFILE_ANALYSIS_RUN))' STOP" >>$LOGFILE_RES_RUN

	# 4.4. Distribution into Groups process (TODO)
	#######################################

done;



##############################
# 5. Permissions and Copying temporary files #
##############################

echo "["`date '+%Y%m%d-%H%M%S'`"] * Copying files and Permission"

for RUN in $RUN_ANALYZED;
do

	SAMPLE_INFO_DEFAULT="UNKNOWN";
	if [ -z "$GROUP" ] || [ "$GROUP" == "$SAMPLE_INFO_DEFAULT" ]; then
		SAMPLE_GROUP_DEFAULT=$SAMPLE_INFO_DEFAULT;
	else
		SAMPLE_GROUP_DEFAULT=$GROUP
	fi;
	SAMPLE_PROJECT_DEFAULT=$SAMPLE_INFO_DEFAULT;
	if [ -z "$PROJECT" ] || [ "$PROJECT" == "$SAMPLE_INFO_DEFAULT" ]; then
		SAMPLE_PROJECT_DEFAULT=$SAMPLE_INFO_DEFAULT;
	else
		SAMPLE_PROJECT_DEFAULT=$PROJECT
	fi;
	SAMPLE_USER_DEFAULT=$SAMPLE_INFO_DEFAULT;

	# SAMPLES
	NB_SAMPLE=$(grep "^RUNS_SAMPLES=" $MAKEFILE_ANALYSIS_RUN | wc -w)
	if [ "$VARANK_ANALYSIS" == 1 ] ; then
		i=0
	fi
	for s in $(grep "RUNS_SAMPLES=" $MAKEFILE_ANALYSIS_RUN | cut -d= -f2);
	do
		RUN=$(echo $s | awk -F: '{print $1}');
		SAMPLE=$(echo $s | awk -F: '{print $2}');
		SAMPLE_INFOS=$(echo $s | awk -F: '{print $3}');
		if [ -z $GROUP ]; then
			SAMPLE_GROUP=$(echo $SAMPLE_INFOS | awk -F- '{print $1}');
		else
			SAMPLE_GROUP=$GROUP
		fi;
		if [ -z $PROJECT ]; then
			SAMPLE_PROJECT=$(echo $SAMPLE_INFOS | awk -F- '{print $2}');
		else
			SAMPLE_PROJECT=$PROJECT
		fi;
		#SAMPLE_PROJECT=$(echo $SAMPLE_INFOS | awk -F- '{print $2}');
		SAMPLE_USER=$(echo $SAMPLE_INFOS | awk -F- '{print $3}');
		#echo -e "$RUN\t$SAMPLE\t$SAMPLE_INFOS\t$SAMPLE_GROUP\t$SAMPLE_PROJECT\t$SAMPLE_USER";

		if [ "$SAMPLE_GROUP" == "" ]; then SAMPLE_GROUP=$SAMPLE_GROUP_DEFAULT ; fi;
		if [ "$SAMPLE_PROJECT" == "" ]; then SAMPLE_PROJECT=$SAMPLE_PROJECT_DEFAULT ; fi;
		if [ "$SAMPLE_USER" == "" ]; then SAMPLE_USER=$SAMPLE_USER_DEFAULT ; fi;

		# Temporary Copy
		if [ $USE_TMP_FOLDER -gt 0 ] && [ "$TMP_RES_DIR" != "$RESULTS_FOLDER" ]; then
			echo "["`date '+%Y%m%d-%H%M%S'`"] Copying RUN/SAMPLE '$RUN/$SAMPLE' files from Temporary folder '$TMP_RES_DIR/$RUN/$SAMPLE' to Results folder '$RESULTS_FOLDER/$RUN/$SAMPLE'..."
			chmod $PERMS -R $TMP_RES_DIR/$RUN/$SAMPLE 1>/dev/null 2>/dev/null
			mkdir -p $RESULTS_FOLDER/$RUN/$SAMPLE
		 	chmod $PERMS -R $RESULTS_FOLDER/$RUN/$SAMPLE 1>/dev/null 2>/dev/null
		 	#$COMMAND_COPY $TMP_RES_DIR/$RUN/$SAMPLE/* $RESULTS_FOLDER/$RUN/$SAMPLE 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA #&
		 	while ! $COMMAND_COPY $TMP_RES_DIR/$RUN/$SAMPLE/* $RESULTS_FOLDER/$RUN/$SAMPLE 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA ; do : ; done;
		fi;

		#echo "["`date '+%Y%m%d-%H%M%S'`"] Reporting RUN/SAMPLE '$RUN/$SAMPLE'..."
		#$NGS_SCRIPTS/stark_report.sh -f $RUN -p $SAMPLE_PROJECT -g $SAMPLE_GROUP -u $SAMPLE_USER -s $SAMPLE -e $ENV -d $DATE_MIN
		#$NGS_SCRIPTS/stark_report.sh -f $RUN -p $SAMPLE_PROJECT -g $SAMPLE_GROUP -u $SAMPLE_USER -s $SAMPLE -e $ENV -i $(echo $PIPELINES | tr " " ",") -d $DATE_MIN
		#  -i $(echo $PIPELINES | tr " " ",") 

		# COPY of run/sample
		# List of folders
		RESULTS_FOLDER_COPY_ALL="";
		if [ "$RESULTS_FOLDER_COPY" != "" ] ; then
			#RESULTS_FOLDER_COPY_ALL=$RESULTS_FOLDER_COPY_ALL" $RESULTS_FOLDER_COPY"
			if [ ! -d $RESULTS_FOLDER_COPY ]; then
				mkdir -p $RESULTS_FOLDER_COPY;
			fi;
			if [ -d $RESULTS_FOLDER_COPY ]; then
				RESULTS_FOLDER_COPY_ALL=$RESULTS_FOLDER_COPY_ALL" $RESULTS_FOLDER_COPY"
			fi;
		fi;
		if [ "$RESULTS_FOLDER_BY_GROUP_PROJECT_COPY" != "" ] ; then
			for RESULTS_FOLDER_COPY_FOLDER in $RESULTS_FOLDER_BY_GROUP_PROJECT_COPY; do
				#RESULTS_FOLDER_COPY_ALL=$RESULTS_FOLDER_COPY_ALL" $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT"
				if [ ! -d $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT ]; then
					mkdir -p $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT;
				fi;
				if [ -d $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT ]; then
					RESULTS_FOLDER_COPY_ALL=$RESULTS_FOLDER_COPY_ALL" $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT";
				fi;
			done;
		fi;
		# Copy
		if [ "$RESULTS_FOLDER_COPY_ALL" != "$RESULTS_FOLDER" ] && [ "$RESULTS_FOLDER_COPY_ALL" != "" ] ; then
			for RESULTS_FOLDER_COPY_FOLDER in $RESULTS_FOLDER_COPY_ALL;
			do
				echo "["`date '+%Y%m%d-%H%M%S'`"] Copying RUN/SAMPLE '$RUN/$SAMPLE' files from main result folder '$RESULTS_FOLDER/$RUN/$SAMPLE' to Copy folder '$RESULTS_FOLDER_COPY_FOLDER/$RUN/$SAMPLE'..."

				# Copy SAMPLE files
				chmod $PERMS -R $RESULTS_FOLDER/$RUN/$SAMPLE 1>/dev/null 2>/dev/null
				mkdir -p $RESULTS_FOLDER_COPY_FOLDER/$RUN/$SAMPLE/$RESULTS_SUBFOLDER_DATA
				chmod $PERMS -R $RESULTS_FOLDER_COPY_FOLDER/$RUN/$SAMPLE/$RESULTS_SUBFOLDER_DATA 1>/dev/null 2>/dev/null
				#$COMMAND_COPY $RESULTS_FOLDER/$RUN/$SAMPLE/* $RESULTS_FOLDER_COPY_FOLDER/$RUN/$SAMPLE/$RESULTS_SUBFOLDER_DATA 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA #&
				while ! nohup $COMMAND_COPY $RESULTS_FOLDER/$RUN/$SAMPLE/* $RESULTS_FOLDER_COPY_FOLDER/$RUN/$SAMPLE/$RESULTS_SUBFOLDER_DATA 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA ; do : ; done;

				# Copy RUN files
				#$COMMAND_COPY $(find $RESULTS_FOLDER/$RUN -mindepth 1 -maxdepth 1 -type f) $RESULTS_FOLDER_COPY_FOLDER/$RUN 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA #&
				while ! nohup $COMMAND_COPY $(find $RESULTS_FOLDER/$RUN -mindepth 1 -maxdepth 1 -type f) $RESULTS_FOLDER_COPY_FOLDER/$RUN 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA ; do : ; done;
				chmod $PERMS $RESULTS_FOLDER_COPY_FOLDER/$RUN/* 1>/dev/null 2>/dev/null

				# Copy ROOT FILE PATTERNS
				if [ $RESULTS_SUBFOLDER_DATA != "" ]; then
					for ROOT_FILE_PATTERN in $RESULTS_SUBFOLDER_ROOT_FILE_PATTERNS; do
						if [[ $ROOT_FILE_PATTERN =~ "SAMPLE" ]]
						then
							eval ROOT_FILE_PATTERN_VAR=$ROOT_FILE_PATTERN
							#$COMMAND_COPY "$RESULTS_FOLDER/$RUN/$SAMPLE/$ROOT_FILE_PATTERN_VAR" "$RESULTS_FOLDER_COPY_FOLDER/$RUN/$SAMPLE"
							while ! $COMMAND_COPY "$RESULTS_FOLDER/$RUN/$SAMPLE/$ROOT_FILE_PATTERN_VAR" "$RESULTS_FOLDER_COPY_FOLDER/$RUN/$SAMPLE" 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA ; do : ; done;
						else
							#$COMMAND_COPY $RESULTS_FOLDER/$RUN/$SAMPLE/$ROOT_FILE_PATTERN $RESULTS_FOLDER_COPY_FOLDER/$RUN/$SAMPLE
							while ! $COMMAND_COPY $RESULTS_FOLDER/$RUN/$SAMPLE/$ROOT_FILE_PATTERN $RESULTS_FOLDER_COPY_FOLDER/$RUN/$SAMPLE 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA ; do : ; done;
						fi
					done;
				fi;
			done;
		fi;

		# LINK of run
		if [ "$RESULTS_FOLDER_LINK" != "$RESULTS_FOLDER" ] && [ "$RESULTS_FOLDER_LINK" != "" ] ; then
			for RESULTS_FOLDER_LINK_FOLDER in $RESULTS_FOLDER_LINK;
			do
				echo "["`date '+%Y%m%d-%H%M%S'`"] Link RUN/SAMPLE '$RUN/$SAMPLE' from main result folder '$RESULTS_FOLDER/$RUN/$SAMPLE' to link folder '$RESULTS_FOLDER_LINK_FOLDER/$RUN/$SAMPLE'..."
				chmod $PERMS -R $RESULTS_FOLDER/$RUN/$SAMPLE 1>/dev/null 2>/dev/null
				mkdir -p $RESULTS_FOLDER_LINK_FOLDER/$RUN
				chmod $PERMS -R $RESULTS_FOLDER_LINK_FOLDER/$RUN 1>/dev/null 2>/dev/null
				if [ ! -s $RESULTS_FOLDER_LINK_FOLDER/$RUN/$SAMPLE ]; then
					ln -s $RESULTS_FOLDER/$RUN/$SAMPLE/ $RESULTS_FOLDER_LINK_FOLDER/$RUN/$SAMPLE 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA;
				fi;
				#echo "ln -s $RESULTS_FOLDER/$RUN/$SAMPLE/ $RESULTS_FOLDER_LINK_FOLDER/$RUN/$SAMPLE"
			done;
		fi;

		### VARANK ANALYSIS ###
		# need alamut-batch license #
		#######################
		# copy the final vcf for varank analysis
		if [ "$VARANK_ANALYSIS" == 1 ] ; then
			if [ -s "$RESULTS_FOLDER/$RUN/$SAMPLE/${SAMPLE}.reports/${SAMPLE}.final.vcf" ] ; then
				config_cpt=1;
				if [ ! -d "${VARANK_FOLDER}/${SAMPLE_GROUP}/${SAMPLE_PROJECT}" ]
				then
					mkdir ${VARANK_FOLDER}/${SAMPLE_GROUP}/${SAMPLE_PROJECT}
					### created a new varank config file
					cp $VARANK/configfile ${VARANK_FOLDER}/${SAMPLE_GROUP}/${SAMPLE_PROJECT}/configfile
					chmod 0777 ${VARANK_FOLDER}/${SAMPLE_GROUP}/${SAMPLE_PROJECT}/configfile
					config_cpt=0;
				fi;
				cp $RESULTS_FOLDER/$RUN/$SAMPLE/${SAMPLE}.reports/${SAMPLE}.final.vcf ${VARANK_FOLDER}/${SAMPLE_GROUP}/${SAMPLE_PROJECT}/
				if [ -f "${VARANK_FOLDER}/${SAMPLE_GROUP}/${SAMPLE_PROJECT}/${SAMPLE}.final.vcf.gz" ]
					then
					rm ${VARANK_FOLDER}/${SAMPLE_GROUP}/${SAMPLE_PROJECT}/${SAMPLE}.final.vcf.gz
				fi;
				unset GZIP; gzip ${VARANK_FOLDER}/${SAMPLE_GROUP}/${SAMPLE_PROJECT}/${SAMPLE}.final.vcf
				array[$i]="${SAMPLE_GROUP}/${SAMPLE_PROJECT}"
				else
				echo "$RESULTS_FOLDER/$RUN/$SAMPLE/${SAMPLE}.reports/${SAMPLE}.final.vcf => this sample is not include in the varank analysis"
			fi;
			i=$(($i+1))

			## ajout du sample dans le fichier config de VARANK
			chmod 0777 ${VARANK_FOLDER}/${SAMPLE_GROUP}/${SAMPLE_PROJECT}/configfile
			CONFIG_FILE="${VARANK_FOLDER}/${SAMPLE_GROUP}/${SAMPLE_PROJECT}/configfile"
	
			if [ "$config_cpt" == 0 ] ; 
			then {
				echo "fam1: ${SAMPLE}" >> $CONFIG_FILE
			}
			else {
				tail -n+67 $CONFIG_FILE | awk -F " " '{print $1"\t"$2}' > ${VARANK_FOLDER}/${SAMPLE_GROUP}/${SAMPLE_PROJECT}/config_sample_list.txt
				last_lign_config=$(tail -1 $CONFIG_FILE | awk '{print $1}')
				last_number_config=$(echo $last_lign_config | sed "s/fam//gi" | sed "s/://gi");

				find=0;
				while read famX sample_config 
				do
					if [ ${SAMPLE} == $sample_config ];
					then {
						find=1; 
					}
				fi;
				done < ${VARANK_FOLDER}/${SAMPLE_GROUP}/${SAMPLE_PROJECT}/config_sample_list.txt

				if [ $find == 0 ];
				then {
					let "last_number_config++"
					echo "fam$last_number_config: ${SAMPLE}" >> $CONFIG_FILE;
				}
				fi;
				rm ${VARANK_FOLDER}/${SAMPLE_GROUP}/${SAMPLE_PROJECT}/config_sample_list.txt
			}
		fi;
	fi;
	######################
	done;

	### VARANK ANALYSIS ###
	# need alamut-batch license #
	#######################
	# on garde toujours le dernier résultats généré, au cas où il y a un problème avec le nouveau, d'où la présence des dossiers "old"
	if [ "$VARANK_ANALYSIS" == 1 ] ; then
		for varank_project in `echo ${array[@]} | tr ' ' '\n' | sort -u | tr '\n' ' '`
		do
			if [ -d "${VARANK_FOLDER}/${varank_project}/Alamut" ]
			then
				# if [ -d "${VARANK_FOLDER}/${varank_project}/Alamut_old" ]
				# then
				# 	rm -r ${VARANK_FOLDER}/${varank_project}/Alamut_old
				# fi;
				# mv ${VARANK_FOLDER}/${varank_project}/Alamut ${VARANK_FOLDER}/${varank_project}/Alamut_old
				rm ${VARANK_FOLDER}/${varank_project}/*tsv ## ${VARANK_FOLDER}/${varank_project}/Alamut_old/

			fi;
			unset GZIP; $VARANK/bin/VaRank -vcfdir "${VARANK_FOLDER}/${varank_project}" -alamutHumanDB $ASSEMBLY -SamVa "yes" -AlamutProcesses $THREADS
			sleep 1m
			# copy VaRank results
			if [ ! -d "${FOLDER_REPOSITORY}/${varank_project}/VARANK" ]
			then
				mkdir ${FOLDER_REPOSITORY}/${varank_project}/VARANK
			fi;
			# if [ -d "${FOLDER_REPOSITORY}/${varank_project}/VARANK/Alamut" ]
			# then
			# 	if [ -d "${FOLDER_REPOSITORY}/${varank_project}/VARANK/Alamut_old" ]
			# 	then
			# 		rm -r "${FOLDER_REPOSITORY}/${varank_project}/VARANK/Alamut_old"
			# 	fi;
			# 	#mkdir ${FOLDER_REPOSITORY}/${varank_project}/VARANK/Alamut_old
			# 	mv ${FOLDER_REPOSITORY}/${varank_project}/VARANK/Alamut ${FOLDER_REPOSITORY}/${varank_project}/VARANK/Alamut_old
			# 	mv ${FOLDER_REPOSITORY}/${varank_project}/VARANK/*tsv ${FOLDER_REPOSITORY}/${varank_project}/VARANK/Alamut_old/
			# fi;
			# cp -r ${VARANK_FOLDER}/${varank_project}/Alamut ${VARANK_FOLDER}/${varank_project}/*tsv ${VARANK_FOLDER}/${varank_project}/*vcf* ${FOLDER_REPOSITORY}/${varank_project}/VARANK/
			cp -r ${VARANK_FOLDER}/${varank_project}/*tsv ${VARANK_FOLDER}/${varank_project}/*vcf* ${FOLDER_REPOSITORY}/${varank_project}/VARANK/
		done;
	fi;
	######################


	# Copy files in the run itsef (log, report...)
	if [ "$RESULTS_FOLDER_COPY" != "$RESULTS_FOLDER" ] && [ "$RESULTS_FOLDER_COPY" != "" ] ; then
		for RESULTS_FOLDER_COPY_FOLDER in $RESULTS_FOLDER_COPY;
		do
			echo "["`date '+%Y%m%d-%H%M%S'`"] Copying RUN files '$RUN/*' files from main result folder '$RESULTS_FOLDER/$RUN' to Copy folder '$RESULTS_FOLDER_COPY_FOLDER/$RUN'..."
			chmod $PERMS -R $RESULTS_FOLDER/$RUN 1>/dev/null 2>/dev/null
			mkdir -p $RESULTS_FOLDER_COPY_FOLDER/$RUN
		 	chmod $PERMS -R $RESULTS_FOLDER_COPY_FOLDER/$RUN 1>/dev/null 2>/dev/null
			#find $RESULTS_FOLDER/$RUN -mindepth 1 -maxdepth 1 -type f -exec $COMMAND_COPY -t $RESULTS_FOLDER_COPY_FOLDER/$RUN {} +  1>>$LOGFILE_ANA 2>>$LOGFILE_ANA #&
			while ! nohup find $RESULTS_FOLDER/$RUN -mindepth 1 -maxdepth 1 -type f -exec $COMMAND_COPY -t $RESULTS_FOLDER_COPY_FOLDER/$RUN {} +  1>>$LOGFILE_ANA 2>>$LOGFILE_ANA ; do : ; done;
		 	#$COMMAND_COPY $RESULTS_FOLDER/$RUN/* $RESULTS_FOLDER_COPY_FOLDER/$RUN 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA #&
		done;
	fi;

	# CHANGE PERMS for RUN
	chmod $PERMS -R $RESULTS_FOLDER/$RUN 1>/dev/null 2>/dev/null


done;



###########################
# 6. Database integration #
###########################

# Database integration
if [ $DATABASE_INTEGRATION -gt 0 ] && [ "$DATABASE_INTEGRATION_SCRIPT" != "" ]; then
	echo "["`date '+%Y%m%d-%H%M%S'`"] * Database integration"

#if [ 1 -gt 0 ]; then

#	echo "["`date '+%Y%m%d-%H%M%S'`"] * Integrate into database"

	# 7.1. Send Analysis Log and Report to ANA

	for RUN in $RUNS;
	do

		#echo "["`date '+%Y%m%d-%H%M%S'`"] Integrating RUN '$RUN'... SKIPPED!!!"
		echo "["`date '+%Y%m%d-%H%M%S'`"] Integrating RUN '$RUN'..."
		#echo "$NGS_BIN/scripts/DBintegration.sh '' $RUN 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA"
		#$NGS_BIN/scripts/DBintegration.sh "" $RUN 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA
		$DATABASE_INTEGRATION_SCRIPT $RUN "" "$ENV" 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA

	done;

#fi;
fi;

###############################
# 7. Multi Run PostProcessing #
###############################

if [ $MULTI_RUN_ANALYSIS_PROCESS -gt 0 ]; # No analysis but the rest (report, validation...)
then

	echo "["`date '+%Y%m%d-%H%M%S'`"] * Multi Run PostProcess"

	# 5.1. Multi Report process
	############################

	# Merging Reports
	echo "["`date '+%Y%m%d-%H%M%S'`"] Merging all Reports"
	REPORT_HEADER=$'### Report for RUNs\n'
	for RUN in $RUNS;
	do
		REPORT_HEADER=$REPORT_HEADER"# "$RUN$'\n'
	done;
	#REPORT_HEADER=$REPORT_HEADER$'\n'
	echo "$REPORT_HEADER" >> $FINAL_REPORT
	echo "$REPORT_HEADER" >> $FINAL_REPORT_FULL
	echo "$REPORT_HEADER" >> $FINAL_REPORT_VALIDATION
	for RUN in $RUNS;
	do
		echo "" >> $FINAL_REPORT
		cat $RESULTS_FOLDER/$RUN/run_analysis.V$DATE_MIN.report >> $FINAL_REPORT
		echo "" >> $FINAL_REPORT_FULL
		cat $RESULTS_FOLDER/$RUN/run_analysis.V$DATE_MIN.full.report >> $FINAL_REPORT_FULL
		echo "" >> $FINAL_REPORT_VALIDATION
		cat $RESULTS_FOLDER/$RUN/run_analysis.V$DATE_MIN.validation.report >> $FINAL_REPORT_VALIDATION
	done;


	# PDF Creation
	echo "["`date '+%Y%m%d-%H%M%S'`"] PDF Creation of Final Report '$FINAL_REPORT.pdf'"
	enscript -f Courier5 --header="run_analysis.V$DATE_MIN.report||[`date '+%d/%m/%Y-%H:%M:%S'`]" -p$FINAL_REPORT.ps $FINAL_REPORT 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA
	ps2pdf $FINAL_REPORT.ps $FINAL_REPORT.pdf 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA
	rm $FINAL_REPORT $FINAL_REPORT.ps 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA
	echo "["`date '+%Y%m%d-%H%M%S'`"] PDF Creation of Final Report Full '$FINAL_REPORT_FULL.pdf'"
	enscript -f Courier5 --header="run_analysis.V$DATE_MIN.full.report||[`date '+%d/%m/%Y-%H:%M:%S'`]" -p$FINAL_REPORT_FULL.ps $FINAL_REPORT_FULL 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA
	ps2pdf $FINAL_REPORT_FULL.ps $FINAL_REPORT_FULL.pdf 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA
	rm $FINAL_REPORT_FULL $FINAL_REPORT_FULL.ps 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA
	echo "["`date '+%Y%m%d-%H%M%S'`"] PDF Creation of Final Report Validation '$FINAL_REPORT_VALIDATION.pdf'"
	enscript -f Courier5 --header="run_analysis.V$DATE_MIN.validation.report||[`date '+%d/%m/%Y-%H:%M:%S'`]" -p$FINAL_REPORT_VALIDATION.ps $FINAL_REPORT_VALIDATION 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA
	ps2pdf $FINAL_REPORT_VALIDATION.ps $FINAL_REPORT_VALIDATION.pdf 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA
	rm  $FINAL_REPORT_VALIDATION $FINAL_REPORT_VALIDATION.ps 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA

	# FINAL_REPORT_FULL_VCF
	#FINAL_REPORT_FULL_VCF=$FINAL_REPORT.vcf

	# Multi Reporting from make failed. skipped
	if  [ 0 -eq 1 ];
	then

		echo "["`date '+%Y%m%d-%H%M%S'`"] * Multi Run PostProcess"

		RES_SOURCE=$RESULTS_FOLDER # By default

		# 5.1. Multi Report process (TODO-TEST)
		############################

		echo "["`date '+%Y%m%d-%H%M%S'`"] Report Process (see '$LOGFILE_ANA')..."
		make -j $THREADS -e NGSEnv=$NGS_FOLDER ENV=$ENV PARAM=$MAKEFILE_ANALYSIS VALIDATION=1 OUTDIR=$RES_SOURCE FINAL_REPORT=$FINAL_REPORT FINAL_REPORT_FULL=$FINAL_REPORT_FULL -f $NGS_SCRIPTS/NGSWorkflow.mk 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA

		# 5.2. Multi Validation process (TODO-TEST)
		################################

		echo "["`date '+%Y%m%d-%H%M%S'`"] Validation Process (see '$LOGFILE_VALIDATION')..."
		make -j $THREADS -e NGSEnv=$NGS_FOLDER ENV=$ENV PARAM=$MAKEFILE_ANALYSIS OUTDIR=$RES_SOURCE FINAL_REPORT_VALIDATION=$FINAL_REPORT_VALIDATION -f $NGS_SCRIPTS/NGSWorkflow.validation.mk 1>>$LOGFILE_VALIDATION 2>>$LOGFILE_VALIDATION

	fi;

fi;


##############################
# 8. Copying temporary files # (TODO-TEST)
##############################

if [ $USE_TMP_FOLDER -gt 0 ]; then

	echo "["`date '+%Y%m%d-%H%M%S'`"] * Copying temporary files"

	# 6.3. Send Analysis Log and Report to ANA
	if [ -d $TMP_ANA_DIR ]; then
		echo "["`date '+%Y%m%d-%H%M%S'`"] Copying files from Temporary folder '$TMP_ANA_DIR/*' to Analysis Log and Report folder '$ANALYSIS_FOLDER'..."
		chmod $PERMS -R $TMP_ANA_DIR/*  1>/dev/null 2>/dev/null
	 	#$COMMAND_COPY $TMP_ANA_DIR/* $ANALYSIS_FOLDER 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA #&
	 	while ! nohup $COMMAND_COPY $TMP_ANA_DIR/* $ANALYSIS_FOLDER 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA ; do : ; done;
	fi;

fi;


##########
# 8. End #
##########

echo "["`date '+%Y%m%d-%H%M%S'`"] * Analysis '$DATE_MIN' Complete!"

exit 0;
