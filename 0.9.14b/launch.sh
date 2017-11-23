#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="Launch"
SCRIPT_DESCRIPTION="Launch RUN Analysis"
SCRIPT_RELEASE="0.9.5b"
SCRIPT_DATE="14/04/2016"
SCRIPT_AUTHOR="Antony Le Bechec / Amandine Velt"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

# Realse note
#RELEASE_NOTES="#\n"
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.3b-09/10/2015: Add RUN configuration depending on Group and Project defined in RUN SampleSheet\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.4b-13/10/2015: Add Complete and Running flag files\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.5b-14/04/2016: Add RAW FOLDER and APPLICATION/ENV detection\n";

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
	echo "# -e/--env/--app=<FILE>                       ENV file configuration of the APPLICATION";
	echo "#                                             Must be in the STARK folder if relative path";
	echo "#                                             Default: defined in the RUN SampleSheet, or env.sh if not defined";
	echo "# -r/--runs=<STRING1,STRING2,...>             List of RUNs to analyse, from Illumina sequencers (mandatory).";
	echo "#                                             RUNS will be automatically searched in all ENV RUN folders, if necessary";
	echo "#                                             Format: RUN1,RUN2...";
	echo "# -a/--aligners=<STRING1,STRING2,...>         List of ALIGNERS to use for the analysis";
	echo "#                                             Format: ALIGNER1,ALIGNER2,... (comma separated)";
	echo "#                                             Default: in ENV or 'bwamem' if not defined";
	echo "#                                             e.g. 'bwamem,bwasw'";
	echo "# -c/--callers=<STRING1,STRING2,...>          List of CALLERS to use for the analysis";
	echo "#                                             Format: CALLER1,CALLER2,... (comma separated)";
	echo "#                                             Default: in ENV or 'gatkHC' if not defined";
	echo "#                                             e.g. 'gatkHC,gatkUG,VarScan,samtools'";
	echo "# -k/--annotators=<STRING1,STRING2,...>       List of ANNOTATORS to use for the analysis";
	echo "#                                             Format: ANNOTATOR1,ANNOTATOR2,... (comma separated)";
	echo "#                                             Default: in ENV or 'howard' if not defined";
	echo "#                                             e.g. 'howard,snpeff'";
	echo "# -p/--pipelines=<STRING1,STRING2,...>        List of PIPELINES to use for the analysis";
	echo "#                                             This options has prior over options --aligners, --callers and --annotators";
	echo "#                                             Format: ALIGNER.CALLER.ANNOTATOR,... (comma separated)";
	echo "#                                             Default: in ENV or 'bwamem.gatkHC.howard'";
	echo "#                                             e.g. 'bwamem.gatkHC.howard,bwamem.gatkUG.howard,bwamem.VarScan.howard'";
	echo "# -f/--filter_samples=<STRING1,STRING2,...>   Samples to use for the analysis";
	echo "#                                             Format: Sample1,Sample2...";
	echo "#                                             Default: all samples";
	echo "# -s/--samplesheet=<FILE>                     Illumina SampleSheet.csv file";
	echo "#                                             Default: found in RUN folder";
	echo "# -z/--parallelization                        Parallelization of RUN analyses";
	echo "#                                             Nota Bene: not fully tested.";
	echo "# -b/--by_sample                              Split analysis by SAMPLE, all threads on each sample, one by one";
	#echo "# -v/--remove=<FILE1,FILE2,...>  Files in RES are removed (if any) before the RUN analysis of the run (optional, not fully tested and risked!).";

	#echo "# -w/--applications_infos      Applications informations.";
	#echo "# -y/--pipelines_infos         Pipelines informations.";
	
	echo "# -v/--verbose                                VERBOSE option";
	echo "# -d/--debug                                  DEBUG option";
	echo "# -n/--release                                RELEASE option";
	echo "# -h/--help                                   HELP option";
	echo "#";
}

# header
header;

PARAM=$@

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "e:r:a:c:k:p:f:s:zum:wyvdnh" --long "env:,app:,runs:,aligners:,callers:,annotators:,pipelines:,filter_samples:,samplesheet:,parallelization,by_sample,remove:,application,pipelines_infos,verbose,debug,release,help" -- "$@" 2> /dev/null)
[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && \
	exit 1
eval set -- "$ARGS"
while true
do
	case "$1" in
		-e|--env|--app)
			ENV="$2"
			shift 2 
			;;
		-r|--runs)
			RUNS="$2"
			# transform RUNS list
			RUNS=$(echo $RUNS | tr "," " ")
			shift 2 
			;;
		-a|--aligners)
			ALIGNERS_INPUT="$2"
			ALIGNERS_INPUT=$(echo $ALIGNERS_INPUT | tr "," " ")
			shift 2 
			;;
		-c|--callers)
			CALLERS_INPUT="$2"
			CALLERS_INPUT=$(echo $CALLERS_INPUT | tr "," " ")
			shift 2 
			;;
		-k|--annotators)
			ANNOTATORS_INPUT="$2"
			ANNOTATORS_INPUT=$(echo $ANNOTATORS_INPUT | tr "," " ")
			shift 2 
			;;
		-p|--pipelines)
			PIPELINES_INPUT="$2"
			PIPELINES_INPUT=$(echo $PIPELINES_INPUT | tr "," " ")
			shift 2 
			;;
		-f|--filter_samples)
			FILTER_SAMPLE="$2"
			FILTER_SAMPLE=$(echo $FILTER_SAMPLE | tr "," " ")
			shift 2 
			;;
		-s|--samplesheet)
			SAMPLESHEET="$2"
			shift 2 
			;;
		-z|--parallelization)
			PARALLELIZATION=1
			shift 1
			;;
		-u|--by_sample)
			BY_SAMPLE=1
			shift 1
			;;
		-m|--remove)
			REMOVE="$2"
			REMOVE=$(echo $REMOVE | tr "," " ")
			shift 2 
			;;
		-w|--application_infos)
			APP=1
			shift 1
			;;
		-y|--pipelines_infos)
			PIPELINES_INFOS=1
			shift 1
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
if [ -z $RUNS ]; then # && [ -z $APP ] && [ -z $PIPELINES_INFOS ]; then 
	echo "Option --runs is required. " "Use -h or --help to display the help." && exit 1;
fi;
	

#if [ "$RUNS" == "" ] && [ "${1^^}" == "APP" ] 
#then
#	:
#elif [ "$RUNS" == "" ] && [ "${1^^}" == "APPLICATION" ]
#then
#	:
#elif [ "$RUNS" == "" ]
#then
#	echo "Option --runs or APP is required. " "Use -h or --help to display the help." && exit 1;
#fi
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# RUNS
if [ "$RUNS" == "" ]; then
	echo "#[ERROR] NO RUN defined"
	exit 0;
fi;




# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# RAW LIST
RAW_LIST="";
for ENV_DEF in $SCRIPT_DIR/env.*sh; do
	# SOURCE
	#source $ENV_DEF;
	MISEQ_FOLDER_DEF=$($SCRIPT_DIR/extract_variable_from_env.sh $ENV_DEF "MISEQ_FOLDER");
	# Add RAW
	RAW_LIST=$RAW_LIST" $MISEQ_FOLDER_DEF"
done;
RAW_LIST=$(echo $RAW_LIST | sed "s/ /\n/gi" | sort | uniq | sed "s/\n/ /gi" );

# ENV
ENV_PARAM=0;
if [ -s $ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$ENV;
	ENV_PARAM=1;
elif [ -s $SCRIPT_DIR/$ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$SCRIPT_DIR/$ENV;
	ENV_PARAM=1;
elif [ "$ENV" == "" ] || [ ! -s $ENV ]; then
	if [ -s $SCRIPT_DIR/"env.sh" ]; then
		ENV=$SCRIPT_DIR/"env.sh";
	else
		ENV="";
		echo "#[WARNING] NO ENV defined. Default ENV used."
	fi;
fi;

# ENV
#source $SCRIPT_DIR/env.sh
if [ "$ENV" != "" ]; then
	source $ENV;
fi;

# APP
APP=$(echo $(basename $ENV_DEF) | sed "s/^env.//gi" | sed "s/.sh$//gi" | sed "s/sh$//gi")



# TEST RUNS validity
for RUN in $RUNS;
do
	RUN=$(echo $RUN | tr -d "\"")
	RUN_FOUND=0;
	if [ ! -d "$MISEQ_FOLDER/$RUN" ]; then
		echo "#[WARNING] RUN '$RUN' doesn't exist in the MiSeq Folder '$MISEQ_FOLDER' using ENV '$ENV'"
		for RAW_LIST_ONE in $RAW_LIST; do
			if [ -d $RAW_LIST_ONE/$RUN ]; then
				RUN_FOUND=1;
				echo "#[WARNING] RUN '$RUN' found in '$RAW_LIST_ONE'";
			fi;
		done;
		if [ !$RUN_TEST_FOUND ]; then
			echo "#[WARNING] RUN '$RUN' not found";
			#exit 1;
		fi;
		#exit 0;
	fi;
done;


# MULTI_RUN_ANALYSIS_PROCESS
if [ "$MULTI_RUN_ANALYSIS_PROCESS" == "" ]; then
	MULTI_RUN_ANALYSIS_PROCESS=0 # 0 OR 1 OR 2 for only multireporting
fi;

# FILTER_SAMPLE
if [ "$FILTER_SAMPLE" == "" ]; then
	FILTER_SAMPLE=""
	echo "#[WARNING] NO FILTER_SAMPLE defined. No filter will be applied on SAMPLES."
fi;
# Transform Filter SAMPLE
FILTER_SAMPLE=$(echo $FILTER_SAMPLE | tr -d "\"")

# INTERSEC
if [ "$INTERSEC" == "" ]; then
	#INTERSEC="3"
	INTERSEC=$(echo "scale=0;  ( ( ( "$(echo $ALIGNERS | wc -w )" * "$(echo $CALLERS | wc -w )" * "$(echo $ANNOTATORS | wc -w )" ) + 1 ) / 2 ) + 1" | bc )
	if [ $INTERSEC -lt 2 ]; then INTERSEC=2; fi
	echo "#[WARNING] NO INTERSEC defined. Default INTERSEC '$INTERSEC' will be used."
fi;

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

# PARALLELIZATION
if [ -z $PARALLELIZATION ];then
	PARALLELIZATION=0;
fi;

# BY_SAMPLE
if [ -z $BY_SAMPLE ];then
	BY_SAMPLE=0;
fi;


# OUTPUT
echo "# "
echo "# CONFIGURATION ($APP_NAME)"
echo "################"
echo "# NGS TOOLS Folder:      "$NGS_TOOLS
echo "# MISEQ Folder:          "$MISEQ_FOLDER
echo "# DEMULTIPLEXING Folder: "$DEMULTIPLEXING_FOLDER
echo "# RESULTS Folder:        "$RESULTS_FOLDER
echo "# ANALYSIS Folder:       "$ANALYSIS_FOLDER
echo "# ENV:                   "$(basename $ENV)
echo "# ALIGNERS:              "$ALIGNERS
echo "# CALLERS:               "$CALLERS
echo "# ANNOTATORS:            "$ANNOTATORS
echo "# PIPELINES:             "$PIPELINES
echo "# "

#DEMULTIPLEXING=/media/IRC/RES2/demultiplexing
#DEMULTIPLEXING=$DEMULTIPLEXING_FOLDER
#$RESULTS_FOLDER=/media/miseq/RES
#mkdir -p $DEMULTIPLEXING_FOLDER
#mkdir -p $RESULTS_FOLDER
#mkdir -p $ANALYSIS_FOLDER

# DATABASE_INTEGRATION_SCRIPT
DATABASE_INTEGRATION_SCRIPT=$PEPPER/DBintegration.sh

if (($PARALLELIZATION)); then
	echo "#"
	RUN_LIST="";
	for RUN in $RUNS;
	do
		RUN=$(echo $RUN | tr -d "\"")
		if [ -d "$MISEQ_FOLDER/$RUN" ]; then
			RUN_LIST=$RUN_LIST" "$RUN;
		fi;
	done;
	RELEASE="V"`date '+%Y%m%d-%H%M%S'`
	echo "# RUNs '$RUNS' - Release '$RELEASE'"
	echo "# CONFIGURATION: use default ENV '$ENV' file"
	LOG=$ANALYSIS_FOLDER/run_analysis.$RELEASE.launch.log
	$SCRIPT_DIR/runs_analysis.sh $PARAM -i "$INTERSEC" -u "1" -d "$DATABASE_INTEGRATION_SCRIPT" -m "$MULTI_RUN_ANALYSIS_PROCESS" 1>>$LOG 2>>$LOG
	echo "# DONE @"`date '+%Y%m%d-%H%M%S'`
else
	for RUN in $RUNS;
	do
		echo "#"
		RUN=$(echo $RUN | tr -d "\"")

		# If RUN not in MISEQ_FOLDER define in ENV, found if exists and change MISEQ_FOLDER (last one)
		MISEQ_FOLDER_CHANGED=0;
		if [ ! -d "$MISEQ_FOLDER/$RUN" ]; then
			for RAW_LIST_ONE in $RAW_LIST; do
				if [ -d $RAW_LIST_ONE/$RUN ]; then
					MISEQ_FOLDER_CHANGED=1;
					MISEQ_FOLDER=$RAW_LIST_ONE;
				fi;
			done;
		fi;

		if [ -d "$MISEQ_FOLDER/$RUN" ] ; then

			# GROUP and PROJET and USER
			if [ -z $SAMPLESHEET ]; then
				SAMPLESHEET=$MISEQ_FOLDER/$RUN/SampleSheet.csv
			fi;
			

			# APP_NAME DEFAULT
			APP_NAME=`grep -i 'Investigator Name' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2` 
			# LIST of Sample Project and then env file(s). If Sample_project is empty, we keep the default env, else, we create a list of environment file. 
			# If the $SCRIPT_DIR/env.$RUN_GROUP.$RUN_PROJECT.$RUN_USER.sh file exist, we keep it as environment file
			# else we keep the default environment file
			awk '/Data/{y=1;next}y' $SAMPLESHEET | tr -d '\r' | sed 's/,/\t/g' > ${TMP_FOLDER}/$RUN_SampleSheet.csv_tmp
			# we recover the list of sample project
			SAMPLES_PROJECT_LIST=$( C=1; for i in $(head ${TMP_FOLDER}/$RUN_SampleSheet.csv_tmp -n 1) ; do if [ $i == "Sample_Project" ] ; then break ; else C=$(( $C + 1 )) ; fi ; done ; cut -f $C ${TMP_FOLDER}/$RUN_SampleSheet.csv_tmp | sort -u | sed 's/Sample_Project//' ); rm ${TMP_FOLDER}/$RUN_SampleSheet.csv_tmp
			# if this list is empty, we use the default env file, else we create the list of samples for each environment file
			
			#echo "$SAMPLES_PROJECT_LIST $APP_NAME" | sort -u
			#exit 0;

			#if [[ $SAMPLES_PROJECT_LIST != "" ]]; then
			if [[ $(echo $SAMPLES_PROJECT_LIST | wc -w) -gt 1 ]]; then
				echo "# STARK launches with several applications.";
				echo "# List of applications is : ` echo $SAMPLES_PROJECT_LIST | sed 's/\n/ /g' `";
				for SAMPLES_PROJECT in `echo $SAMPLES_PROJECT_LIST`
				do
					unset LIST_OF_SAMPLES_NEW
					unset LIST_OF_SAMPLES
					
					APP_NAME=$SAMPLES_PROJECT 
					
					RUN_GROUP=`echo $SAMPLES_PROJECT | tr "." "-" | cut -d \- -f 1` #
					RUN_PROJECT=`echo $SAMPLES_PROJECT | tr "." "-" | cut -d \- -f 2`
					RUN_USER=`echo $SAMPLES_PROJECT | tr "." "-" | cut -d \- -f 3`
					
					
					# Define the environment file to use
					RUN_ENV=$ENV
					
					# 3. Find the ENV with the GROUP.PROJECT structure
					if [ -s $SCRIPT_DIR/env.$RUN_GROUP.$RUN_PROJECT.$RUN_USER.sh ]; then
						RUN_ENV=$SCRIPT_DIR/env.$RUN_GROUP.$RUN_PROJECT.$RUN_USER.sh;
					elif [ -s $SCRIPT_DIR/env.$RUN_GROUP.$RUN_PROJECT.sh ]; then
						RUN_ENV=$SCRIPT_DIR/env.$RUN_GROUP.$RUN_PROJECT.sh;
					elif [ -s $SCRIPT_DIR/env.$RUN_GROUP.sh ]; then
						RUN_ENV=$SCRIPT_DIR/env.$RUN_GROUP.sh;
					fi;
					
					# 2. Find the ENV with the ENV file (env*sh))
					if [ -s $SCRIPT_DIR/env.$APP_NAME.sh ]; then
						RUN_ENV=$SCRIPT_DIR/env.$APP_NAME.sh;
					fi;
					
					# 1. Find the ENV with the APP_NAME
					for ENV_DEF in $(ls $SCRIPT_DIR/env.*sh); do
						APP_NAME_DEF=$($SCRIPT_DIR/extract_variable_from_env.sh $ENV_DEF "APP_NAME")
						if [ "$RUN_ENV" == "$APP_NAME_DEF" ] && [ "$RUN_ENV" != "" ] && [ "$APP_NAME_DEF" != "" ]; then
							RUN_ENV=$ENV_DEF
						fi;
					done;
					
					LIST_OF_SAMPLES=($(awk '/Data/{y=1;next}y' $SAMPLESHEET | tr -d '\r' | grep "$SAMPLES_PROJECT" | cut -d"," -f1))
					# echo "for $SAMPLES_PROJECT : list of samples : `for item in ${LIST_OF_SAMPLES[*]}; do printf "%s " $item; done | sed 's/ $//'`"
					# we remove samples which are not in the $FILTER_SAMPLE variable
					if [ "$FILTER_SAMPLE" != "" ]; then
						for item in ${LIST_OF_SAMPLES[*]}
						do
							if grep --quiet $item <<< "$FILTER_SAMPLE"; then
								LIST_OF_SAMPLES_NEW+=($item)
							else
								:
							fi
						done
					else
						for item in ${LIST_OF_SAMPLES[*]}
						do
							LIST_OF_SAMPLES_NEW+=($item)
						done
					fi
					echo "# Running with $SAMPLES_PROJECT ..."
					echo "# List of samples for $SAMPLES_PROJECT: `for item in ${LIST_OF_SAMPLES_NEW[*]}; do printf "%s " $item; done | sed 's/ $//'`"
					echo "#"
					RUN_NGS_TOOLS=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "NGS_TOOLS")
					RUN_MISEQ_FOLDER=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "MISEQ_FOLDER")
					RUN_DEMULTIPLEXING_FOLDER=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "DEMULTIPLEXING_FOLDER")
					RUN_RESULTS_FOLDER=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "RESULTS_FOLDER")
					RUN_ANALYSIS_FOLDER=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "ANALYSIS_FOLDER")
					RUN_ALIGNERS=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "ALIGNERS")
					RUN_CALLERS=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "CALLERS")
					RUN_ANNOTATORS=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "ANNOTATORS")
					RUN_PIPELINES=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "PIPELINES")
					TMP_FOLDER=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "TMP_FOLDER")
					#ANALYSIS_FOLDER=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "ANALYSIS_FOLDER")
					STARK_QUEUED=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "STARK_QUEUED")
					STARK_RUNNING=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "STARK_RUNNING")
					STARK_COMPLETE=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "STARK_COMPLETE")
					STARK_VERSION=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "STARK_VERSION")
					
					RUN_APP=$(echo $(basename $RUN_ENV) | sed "s/^env.//gi" | sed "s/.sh$//gi" | sed "s/sh$//gi")
					if [ "$RUN_APP" == "" ]; then RUN_APP="default"; fi

					# CONFIGURATION RUN
					echo "# CONFIGURATION '$SAMPLES_PROJECT' "
					echo "################################################"
					echo "# RUN:                   $RUN"
					echo "# APP:                   $RUN_APP ($(basename $RUN_ENV))"
					echo "# GROUP:                 $RUN_GROUP"
					echo "# PROJECT:               $RUN_PROJECT"
					echo "# USER:                  $RUN_USER"
					echo "# NGS TOOLS Folder:      "$RUN_NGS_TOOLS
					echo "# MISEQ Folder:          "$RUN_MISEQ_FOLDER
					echo "# DEMULTIPLEXING Folder: "$RUN_DEMULTIPLEXING_FOLDER
					echo "# RESULTS Folder:        "$RUN_RESULTS_FOLDER
					echo "# ANALYSIS Folder:       "$RUN_ANALYSIS_FOLDER
					echo "# DEFAULT ENV:           "$RUN_ENV
					echo "# ALIGNERS:              "$RUN_ALIGNERS
					echo "# CALLERS:               "$RUN_CALLERS
					echo "# ANNOTATORS:            "$RUN_ANNOTATORS
					echo "# PIPELINES:             "$RUN_PIPELINES
					echo "# "
					#echo "# THREADS: $THREADS"
					#echo "# MEMORY:  $MEMORY"
					echo "#"

					# REMOVE?
					if [ "$REMOVE" != "" ]; then
						#for PAT in $REMOVE;
						#do
							echo "# REMOVE FILES '$REMOVE'"
							#echo "ls $RESULTS_FOLDER/$RUN/$REMOVE"
							rm -f $RUN_RESULTS_FOLDER/$RUN/$REMOVE
							rm -f $TMP_FOLDER/RES/ALL/$RUN/$REMOVE
						#done;
					fi;
					RELEASE="V"`date '+%Y%m%d-%H%M%S'`
					echo "# RUN '$RUN' - Release '$RELEASE'"

					#continue;

					LOG=$RUN_ANALYSIS_FOLDER/run_analysis.$RELEASE.launch.log
					if [ "$STARK_QUEUED" == "" ]; then STARK_QUEUED=STARKQueued.txt; fi;
					if [ "$STARK_RUNNING" == "" ]; then  STARK_RUNNING=STARKRunning.txt; fi;
					if [ "$STARK_COMPLETE" == "" ]; then  STARK_COMPLETE=STARKComplete.txt; fi;
					STARK_QUEUED_FILE=$RUN_RESULTS_FOLDER/$RUN/$STARK_QUEUED
					STARK_RUNNING_FILE=$RUN_RESULTS_FOLDER/$RUN/$STARK_RUNNING
					STARK_COMPLETE_FILE=$RUN_RESULTS_FOLDER/$RUN/$STARK_COMPLETE

					# Create directory
					mkdir -p $RUN_RESULTS_FOLDER/$RUN;

#echo "RUN_ENV=$RUN_ENV"; continue;

					# RUNNING
					echo "#["`date '+%Y%m%d-%H%M%S'`"] RUN $RUN running with STARK ($STARK_VERSION)" > $STARK_RUNNING_FILE

					#echo "$SCRIPT_DIR/runs_analysis.sh "$RUN_ENV" "$RUN" "$ALIGNERS_INPUT" "$CALLERS_INPUT" "$ANNOTATORS_INPUT" "$FILTER_SAMPLE" "$INTERSEC" "1" "$DATABASE_INTEGRATION_SCRIPT" "$MULTI_RUN_ANALYSIS_PROCESS" 1>>$LOG 2>>$LOG	"
					#exit 0;
					#echo "for $RUN_ENV, the list of sample have length : ${#LIST_OF_SAMPLES_NEW[@]}"
					if [ ${#LIST_OF_SAMPLES_NEW[@]} -eq 0 ]; then
						echo "#["`date '+%Y%m%d-%H%M%S'`"]  RUN $RUN : No samples for $RUN_ENV - no analysis to do" >> $STARK_COMPLETE_FILE
						:
					elif (($BY_SAMPLE)); then
						for item in ${LIST_OF_SAMPLES_NEW[*]}; do
							echo "# RUN '$RUN' / SAMPLE '$item'"
							$SCRIPT_DIR/runs_analysis.sh $PARAM -e "$RUN_ENV" -f "$item" -i "$INTERSEC" -u "1" -d "$DATABASE_INTEGRATION_SCRIPT" -m "$MULTI_RUN_ANALYSIS_PROCESS" -s "$SAMPLESHEET" 1>>$LOG 2>>$LOG
							(($VERBOSE)) && cat $LOG
						done;
					else
						$SCRIPT_DIR/runs_analysis.sh $PARAM -e "$RUN_ENV" -r "$RUN" -f "`for item in ${LIST_OF_SAMPLES_NEW[*]}; do printf "%s " $item; done | sed 's/ $//'`" -i "$INTERSEC" -u "1" -d "$DATABASE_INTEGRATION_SCRIPT" -m "$MULTI_RUN_ANALYSIS_PROCESS" -s "$SAMPLESHEET" 1>>$LOG 2>>$LOG
						(($VERBOSE)) && cat $LOG
						# COMPLETE
						echo "#["`date '+%Y%m%d-%H%M%S'`"] RUN $RUN analyzed by STARK ($STARK_VERSION)" >> $STARK_COMPLETE_FILE
						grep "ERROR" $LOG >> $STARK_COMPLETE_FILE
						if [ "$(grep "ERROR" $LOG -c)" != "0" ]; then
							echo "# Complete with ERROR "
						fi;
					fi

					# RUNNING stop
					rm -f $STARK_RUNNING_FILE

					# SOURCE ORIGINAL ENV
					#source $ENV;

					echo "# DONE @"`date '+%Y%m%d-%H%M%S'`

				done
			# we launch the run with one environment file
			else
				# APP GROUP and PROJET and USER
				
				APP_NAME=`grep -i 'Investigator Name' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2` 
				RUN_GROUP=`echo $APP_NAME | tr "." "-" | cut -d \- -f 1` 
				RUN_PROJECT=`echo $APP_NAME | tr "." "-"  | cut -d \- -f 2`
				RUN_USER=`echo $APP_NAME | tr "." "-"  | cut -d \- -f 3`

				# RUN ENV
				# ENV-AUTO in ENV allows looking for ENV depending on GROUP/PROJECT/USER
				RUN_ENV=$ENV;
				if ((!$ENV_PARAM)) || (($ENV_AUTO)) || (($MISEQ_FOLDER_CHANGED)); then
					
					# 3. Find the ENV with the GROUP.PROJECT structure
					if [ -s $SCRIPT_DIR/env.$RUN_GROUP.$RUN_PROJECT.$RUN_USER.sh ]; then
						RUN_ENV=$SCRIPT_DIR/env.$RUN_GROUP.$RUN_PROJECT.$RUN_USER.sh;
					elif [ -s $SCRIPT_DIR/env.$RUN_GROUP.$RUN_PROJECT.sh ]; then
						RUN_ENV=$SCRIPT_DIR/env.$RUN_GROUP.$RUN_PROJECT.sh;
					elif [ -s $SCRIPT_DIR/env.$RUN_GROUP.sh ]; then
						RUN_ENV=$SCRIPT_DIR/env.$RUN_GROUP.sh;
					fi;
					
					# 2. Find the ENV with the ENV file (env*sh))
					if [ -s $SCRIPT_DIR/env.$APP_NAME.sh ]; then
						RUN_ENV=$SCRIPT_DIR/env.$APP_NAME.sh;
					fi;
					
					# 1. Find the ENV with the APP_NAME
					for ENV_DEF in $(ls $SCRIPT_DIR/env.*sh); do
						APP_NAME_DEF=$($SCRIPT_DIR/extract_variable_from_env.sh $ENV_DEF "APP_NAME")
						if [ "$RUN_ENV" == "$APP_NAME_DEF" ] && [ "$RUN_ENV" != "" ] && [ "$APP_NAME_DEF" != "" ]; then
							RUN_ENV=$ENV_DEF
						fi;
					done;
					
				fi;
				#echo "$ENV $RUN_ENV"; exit 0;
				# SOUCRE ENV for RUN
				#source $RUN_ENV;
				RUN_APP_NAME=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "APP_NAME")
				RUN_NGS_TOOLS=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "NGS_TOOLS")
				RUN_MISEQ_FOLDER=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "MISEQ_FOLDER")
				RUN_DEMULTIPLEXING_FOLDER=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "DEMULTIPLEXING_FOLDER")
				RUN_RESULTS_FOLDER=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "RESULTS_FOLDER")
				RUN_ANALYSIS_FOLDER=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "ANALYSIS_FOLDER")
				RUN_ALIGNERS=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "ALIGNERS")
				RUN_CALLERS=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "CALLERS")
				RUN_ANNOTATORS=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "ANNOTATORS")
				RUN_PIPELINES=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "PIPELINES")
				RUN_POST_ALIGNMENT=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "POST_ALIGNMENT")
				RUN_POST_ALIGNMENT_STEPS=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "POST_ALIGNMENT_STEPS")
				TMP_FOLDER=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "TMP_FOLDER")
				#ANALYSIS_FOLDER=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "ANALYSIS_FOLDER")
				STARK_QUEUED=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "STARK_QUEUED")
				STARK_RUNNING=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "STARK_RUNNING")
				STARK_COMPLETE=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "STARK_COMPLETE")
				STARK_VERSION=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "STARK_VERSION")
				
				#RUN_APP=$(echo $(basename $RUN_ENV) | sed "s/^env.//gi" | sed "s/.sh$//gi" | sed "s/sh$//gi")
				#if [ "$RUN_APP" == "" ]; then RUN_APP="default"; fi

				# CONFIGURATION RUN
				echo "# CONFIGURATION '$RUN_APP_NAME' : "
				echo "################################################"
				echo "# RUN:                   $RUN"
				echo "# APP:                   $RUN_APP_NAME ($(basename $RUN_ENV))"
				echo "# GROUP:                 $RUN_GROUP"
				echo "# PROJECT:               $RUN_PROJECT"
				echo "# USER:                  $RUN_USER"
				echo "# NGS TOOLS Folder:      "$RUN_NGS_TOOLS
				echo "# MISEQ Folder:          "$RUN_MISEQ_FOLDER
				echo "# DEMULTIPLEXING Folder: "$RUN_DEMULTIPLEXING_FOLDER
				echo "# RESULTS Folder:        "$RUN_RESULTS_FOLDER
				echo "# ANALYSIS Folder:       "$RUN_ANALYSIS_FOLDER
				echo "# ENV:                   "$(basename $RUN_ENV)
				echo "# ALIGNERS:              "$RUN_ALIGNERS
				echo "# CALLERS:               "$RUN_CALLERS
				echo "# ANNOTATORS:            "$RUN_ANNOTATORS
				echo "# PIPELINES:             "$RUN_PIPELINES
				echo "# POST_ALIGNMENT_STEPS:  "$RUN_POST_ALIGNMENT_STEPS
				#echo "# POST_ALIGNMENT:        "$RUN_POST_ALIGNMENT
				echo "# SAMPLES:               "$FILTER_SAMPLE
				#echo "# THREADS: $THREADS"
				#echo "# MEMORY:  $MEMORY"
				echo "#"
				
				# REMOVE?
				if [ "$REMOVE" != "" ]; then
					#for PAT in $REMOVE;
					#do
						echo "# REMOVE FILES '$REMOVE'"
						#echo "ls $RESULTS_FOLDER/$RUN/$REMOVE"
						rm -f $RUN_RESULTS_FOLDER/$RUN/$REMOVE
						rm -f $TMP_FOLDER/RES/ALL/$RUN/$REMOVE
					#done;
				fi;
				RELEASE="V"`date '+%Y%m%d-%H%M%S'`
				echo "# RUN '$RUN' - Release '$RELEASE'"

				#continue;

				LOG=$RUN_ANALYSIS_FOLDER/run_analysis.$RELEASE.launch.log
				if [ "$STARK_QUEUED" == "" ]; then STARK_QUEUED=STARKQueued.txt; fi;
				if [ "$STARK_RUNNING" == "" ]; then  STARK_RUNNING=STARKRunning.txt; fi;
				if [ "$STARK_COMPLETE" == "" ]; then  STARK_COMPLETE=STARKComplete.txt; fi;
				STARK_QUEUED_FILE=$RUN_RESULTS_FOLDER/$RUN/$STARK_QUEUED
				STARK_RUNNING_FILE=$RUN_RESULTS_FOLDER/$RUN/$STARK_RUNNING
				STARK_COMPLETE_FILE=$RUN_RESULTS_FOLDER/$RUN/$STARK_COMPLETE

				# Create directory
				mkdir -p $RUN_RESULTS_FOLDER/$RUN;

				# RUNNING
				echo "#["`date '+%Y%m%d-%H%M%S'`"] RUN $RUN running with STARK ($STARK_VERSION)" > $STARK_RUNNING_FILE

				#echo "$SCRIPT_DIR/runs_analysis.sh "$RUN_ENV" "$RUN" "$ALIGNERS_INPUT" "$CALLERS_INPUT" "$ANNOTATORS_INPUT" "$FILTER_SAMPLE" "$INTERSEC" "1" "$DATABASE_INTEGRATION_SCRIPT" "$MULTI_RUN_ANALYSIS_PROCESS" 1>>$LOG 2>>$LOG	"
				#exit 0;
				if (($BY_SAMPLE)) && [ "$FILTER_SAMPLE" != "" ]; then
					for item in ${FILTER_SAMPLE}; do
						echo "# RUN '$RUN' / SAMPLE '$item'"
						$SCRIPT_DIR/runs_analysis.sh $PARAM -e "$RUN_ENV" -r "$RUN" -f "$item" -i "$INTERSEC" -u "1" -d "$DATABASE_INTEGRATION_SCRIPT" -m "$MULTI_RUN_ANALYSIS_PROCESS" -s "$SAMPLESHEET" 1>>$LOG 2>>$LOG
						(($VERBOSE)) && cat $LOG
					done;
				elif (($BY_SAMPLE)); then
					# list of sample
					NB_SAMPLE=0
					DATA_SECTION_LINE=$(grep "^\[Data\]" $SAMPLESHEET -n | awk -F: '{print $1}')
					SAMPLE_SECTION_FIRST_LINE=$(($DATA_SECTION_LINE+1))
					SAMPLE_LINES=$(tail -n $(($SAMPLE_SECTION_FIRST_LINE-$(grep ^ -c $SAMPLESHEET)))  $SAMPLESHEET)
					for L in $(tail -n $(($SAMPLE_SECTION_FIRST_LINE-$(grep ^ -c $SAMPLESHEET)))  $SAMPLESHEET);
					do
						if [ "$L" == "" ] || [[ $L =~ ^\[.*\] ]]; then
							break;
						else
							((NB_SAMPLE++))
							SAMPLE_ID=$(echo $L | awk -F, '{print $1}')
							echo "# RUN '$RUN' / SAMPLE '$SAMPLE_ID'"
							$SCRIPT_DIR/runs_analysis.sh $PARAM -e "$RUN_ENV" -r "$RUN" -f "$SAMPLE_ID" -i "$INTERSEC" -u "1" -d "$DATABASE_INTEGRATION_SCRIPT" -m "$MULTI_RUN_ANALYSIS_PROCESS" -s "$SAMPLESHEET" 1>>$LOG 2>>$LOG
							(($VERBOSE)) && cat $LOG
						fi;
					done;
				else
				#echo "$SCRIPT_DIR/runs_analysis.sh $PARAM -e "$RUN_ENV" -r "$RUN" -f "$FILTER_SAMPLE" -i "$INTERSEC" -u "1" -d "$DATABASE_INTEGRATION_SCRIPT" -m "$MULTI_RUN_ANALYSIS_PROCESS" -s "$SAMPLESHEET" 1>>$LOG 2>>$LOG"; exit 0;
					$SCRIPT_DIR/runs_analysis.sh $PARAM -e "$RUN_ENV" -r "$RUN" -f "$FILTER_SAMPLE" -i "$INTERSEC" -u "1" -d "$DATABASE_INTEGRATION_SCRIPT" -m "$MULTI_RUN_ANALYSIS_PROCESS" -s "$SAMPLESHEET" 1>>$LOG 2>>$LOG
					(($VERBOSE)) && cat $LOG;
				fi;

				# COMPLETE
				echo "#["`date '+%Y%m%d-%H%M%S'`"] RUN $RUN analyzed by STARK ($STARK_VERSION)" > $STARK_COMPLETE_FILE
				grep "ERROR" $LOG >> $STARK_COMPLETE_FILE

				if [ "$(grep "ERROR" $LOG -c)" != "0" ]; then
					echo "# Complete with ERROR "
				fi;

				# RUNNING stop
				rm -f $STARK_RUNNING_FILE
				rm -f $RUN_RESULTS_FOLDER/$RUN/*/*manifest_by_samples.txt $RUN_RESULTS_FOLDER/$RUN/*/*.genes

				# SOURCE ORIGINAL ENV
				#source $ENV;

				echo "# DONE @"`date '+%Y%m%d-%H%M%S'`
			fi;
		else
			echo "# RUN '$RUN' doesn't exist in the MiSeq Folder '$MISEQ_FOLDER'"
		fi;
	done;
fi;
echo "#"


exit 1;
