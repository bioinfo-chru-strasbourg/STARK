#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARKLaunch"
SCRIPT_DESCRIPTION="STARK launch RUN or FASTQ/BAM analysis"
SCRIPT_RELEASE="0.9b"
SCRIPT_DATE="04/10/2016"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-30/05/2016: Script creation\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

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
	echo "# USAGE: $(basename $0) --runs=<RUN>|--fastq=<FASTQ> [options...]";
	echo "# -e/--env/--app=<FILE>                       ENV file configuration of the APPLICATION.";
	echo "#                                             Must be in the STARK folder if relative path";
	echo "#                                             Default: defined in the RUN SampleSheet, or env.sh if not defined";
	echo "# -f/--fastq=<FILE1,FILE2,...>                List of FASTQ/BAM/SAM/CRAM (mandatory if no --runs option)";
	echo "#                                             Formats: *fastq.gz/*fq.gz/*bam/*ubam/*cram/*ucram/*sam/*usam";
	echo "# -r/--runs=<STRING1,STRING2,...>             List of RUNs to analyse, from Illumina sequencers (mandatory if no --fastq option).";
	echo "#                                             RUNS will be automatically searched in all ENV RUN folders, if necessary";
	echo "#                                             Format: RUN1,RUN2...";
	echo "# -s/--samplesheet=<FILE>                     Illumina SampleSheet.csv file";
	echo "#                                             Default: found in RUN folder";
	

	echo "# -i/--application_infos                      Applications informations.";
	echo "# -y/--pipelines_infos                        Pipelines informations.";
	echo "# -g/--release_infos                          Tools, databases and pipelines release informations.";

	echo "# -v/--verbose                                VERBOSE option";
	echo "# -d/--debug                                  DEBUG option";
	echo "# -n/--release                                RELEASE option";
	echo "# -h/--help                                   HELP option";
	echo "#";
	echo -e "#\n# RUN Analysis\n################";
	$SCRIPT_DIR/launch.sh -h | grep "# [ |-]";
	echo -e "#\n# SAMPLE Analysis\n###################";
	$SCRIPT_DIR/launch.sample.sh -h | grep "# [ |-]";
}

# header
header;

# launch.sh 
# env:,runs:,aligners:,callers:,annotators:,filter_samples:,samplesheet:,parallelization:,remove:,application,help
# e:r:a:c:n:f:s:p:v:wh
# launch.sample.sh
# e:f:q:b:s:r:p:t:a:h
# env:,fastq_R1:,fastq_R2:,bed:,sample:,run:,pipelines:,threads:,adapters:,help

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "e:r:f:w:t:iygvdnh" --long "env:,app:,runs:,fastq:,samplesheet:,application_infos,pipelines_infos,release_infos,verbose,debug,release,help" -- "$@" 2> /dev/null)
# || [ -z $@ ]
if [ $? -ne 0 ]; then
	:
	#echo $?
	#usage;
	#exit;
fi;
PARAM=$@
#PARAM=$(echo $@ | tr "\n" " ")
#echo $PARAM;
#PARAM=$(echo $ARGS | sed s/--//gi);
#exit 0;
 
eval set -- "$ARGS"
while true
do
	#echo "$1=$2"
	#echo "Eval opts";
	case "$1" in
		-e|--env|--app)
			ENV="$2"
			shift 2 
			;;
		-r|--runs)
			RUNS="$2"
			shift 2 
			;;
		-f|--fastq)
			FASTQ=$2 #"$(echo $2 | tr "\n" " ")"
			shift 2 
			;;
		-s|--samplesheet)
			SAMPLESHEET="$2" #"$(echo $2 | tr "\n" " ")"
			shift 2 
			;;
		-i|--application_infos)
			APPLICATIONS_INFOS=1
			shift 1
			;;
		-y|--pipelines_infos)
			PIPELINES_INFOS=1
			shift 1
			;;
		-g|--release_infos)
			RELEASE_INFOS=1
			shift 1
			;;
		-t|--threads)
			THREADS_INPUT="$2"
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
		-w)
			W="$2"
			shift 2 
			;;
		--) shift
			break 
			;;
		*) 	echo "# Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done

#		*) 	echo "# Option $1 is not recognized. " "Use -h or --help to display the help." && \
#			exit 1
#			;;

####################################################################################################################################
# Checking the input parameter
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if [ -z "$FASTQ" ] && [ -z "$RUNS" ] && [ -z "$APPLICATIONS_INFOS" ] && [ -z "$PIPELINES_INFOS" ] && [ -z "$RELEASE_INFOS" ]; then
	echo "Option --runs or --fastq or --applications_infos or --pipelines_infos or --release_infos is required. " "Use -h or --help to display the help." && usage && exit 1;
fi
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#echo "if [ "$FASTQ" == "ollow" ] && [ "$THREADS_INPUT" == "he" ] && [ "$W" == "hite" ] && [ "$RUNS" == "abbit" ];"
if [ "$FASTQ" == "ollow" ] && [ "$THREADS_INPUT" == "he" ] && [ "$W" == "hite" ] && [ "$RUNS" == "abbit" ]; then
	echo "test";
	echo -e "\033[2J\033[?25l"; R=`tput lines` C=`tput cols`;: $[R--] ; while true 
	do ( e=echo\ -e s=sleep j=$[RANDOM%C] d=$[RANDOM%R];for i in `eval $e {1..$R}`;
	do c=`printf '\\\\0%o' $[RANDOM%57+33]` ###  ###
	$e "\033[$[i-1];${j}H\033[32m$c\033[$i;${j}H\033[37m"$c; $s 0.1;if [ $i -ge $d ]
	then $e "\033[$[i-d];${j}H ";fi;done;for i in `eval $e {$[i-d]..$R}`; #[mat!rix]
	do echo -e "\033[$i;${j}f ";$s 0.1;done)& sleep 0.05;done #(c) 2011 -- [ BruXy ]
	exit 0;
fi;
#echo "if [ "$FASTQ" == "ollow" ] && [ "$THREADS_INPUT" == "he" ] && [ "$W" == "hite" ] && [ "$RUNS" == "abbit" ];"
#exit 0;


# ENV
#########

if [ -s $ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$ENV;
elif [ -s $SCRIPT_DIR/$ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$SCRIPT_DIR/$ENV;
elif [ "$ENV" == "" ] || [ ! -s $ENV ]; then
	if [ -s $SCRIPT_DIR/"env.sh" ]; then
		ENV=$SCRIPT_DIR/"env.sh";
	else
		ENV="";
		echo "#[WARNING] NO ENV defined. Default ENV used."
	fi;
fi;
# SOURCE ENV if exists
if [ ! -z $ENV ] && [ -s $ENV ]; then
	source $ENV;
fi;

if [ ! -e $SAMPLESHEET ]; then
	echo "#[ERROR] SampleSheet '$SAMPLESHEET' does NOT exists."
	exit 0;
fi;




# APPLICATIONS_INFOS
#######

if (($APPLICATIONS_INFOS)); then

	SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

	RUNS_TEST=$RUNS

	echo "#"
	echo "######################"
	echo "#### APPLICATIONS ####"
	echo "######################"
	echo "#"
	RAW_LIST="";
	FOLDER_RUN_ORIGINAL=$FOLDER_RUN
	#echo "TRUC "$(echo $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/env.*sh)
	for ENV_DEF in $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/env.*sh; do
		# FIND APP
		APP=$(echo $(basename $ENV_DEF) | sed "s/^env.//gi" | sed "s/.sh$//gi" | sed "s/sh$//gi")
		if [ "$APP" == "" ]; then
			APP="default";
			ENV=$APP;
		fi
		# SOURCE
		FOLDER_RUN=$($SCRIPT_DIR/extract_variable_from_env.sh $ENV_DEF "FOLDER_RUN")
		#RAW_FOLDER=$($SCRIPT_DIR/extract_variable_from_env.sh $ENV_DEF "RAW_FOLDER")
		#MISEQ_FOLDER=$($SCRIPT_DIR/extract_variable_from_env.sh $ENV_DEF "MISEQ_FOLDER");
		FOLDER_RESULTS=$($SCRIPT_DIR/extract_variable_from_env.sh $ENV_DEF "FOLDER_RESULTS")
		PIPELINES=$($SCRIPT_DIR/extract_variable_from_env.sh $ENV_DEF "PIPELINES")
		POST_ALIGNMENT=$($SCRIPT_DIR/extract_variable_from_env.sh $ENV_DEF "POST_ALIGNMENT")
		POST_ALIGNMENT_STEPS=$($SCRIPT_DIR/extract_variable_from_env.sh $ENV_DEF "POST_ALIGNMENT_STEPS")
		# OUTPUT
		#if [ "$RUNS_TEST" == "" ]; then
		echo "### $APP [$(basename $ENV_DEF)]";
		echo "#    RAW DATA $FOLDER_RUN";
		echo "#    RES DATA $MAIN_FOLDER";
		echo "#    PIPELINES "
		echo -e "#       $(echo $PIPELINES | sed "s/ /\n#       /gi")";
		echo "#    POST_ALIGNMENT $POST_ALIGNMENT "
		echo "#    POST_ALIGNMENT_STEPS $POST_ALIGNMENT_STEPS "
		#fi;
		RAW_LIST=$RAW_LIST" $FOLDER_RUN"
	done;
	RAW_LIST=$(echo $RAW_LIST | sed "s/ /\n/gi" | sort | uniq | sed "s/\n/ /gi" );

	if [ "$RUNS_TEST" != "" ]; then

		for RUN_TEST in $(echo $RUNS_TEST | tr "," " "); do 
	
		echo "#";
		echo "#### RUN '$RUN_TEST' ####";
		RUN_TEST_FOUND=0;
		MISEQ_FOLDER_FOUND="";
		if [ -d $FOLDER_RUN_ORIGINAL/$RUN_TEST ]; then
			RUN_TEST_FOUND=1;
			MISEQ_FOLDER_FOUND=$FOLDER_RUN_ORIGINAL;
		fi;

		if ((!$RUN_TEST_FOUND)); then
		for RAW_LIST_ONE in $RAW_LIST; do
			if [ -d $RAW_LIST_ONE/$RUN_TEST ] && ((!$RUN_TEST_FOUND)); then 
				RUN_TEST_FOUND=1;
				MISEQ_FOLDER_FOUND=$RAW_LIST_ONE;
			fi;
		done;
		fi;

		if ((!$RUN_TEST_FOUND)); then
			echo "#[ERROR] RUN not found";
		else
			echo "# RAW FOLDER $MISEQ_FOLDER_FOUND";
		fi;

		# Test SampleSheet to determine RUN_TEST possible Application/ENV
		if [ -z $SAMPLESHEET ]; then
			SAMPLESHEET=$MISEQ_FOLDER_FOUND/$RUN_TEST/SampleSheet.csv
		fi;

		if [ -e $SAMPLESHEET ]; then 
			RUN_GROUP=`grep -i 'Investigator Name' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2 | cut -d \- -f 1`
			RUN_PROJECT=`grep -i 'Investigator Name' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2 | cut -d \- -f 2`
			RUN_USER=`grep -i 'Investigator Name' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2 | cut -d \- -f 3`

			# RUN ENV
			# ENV-AUTO in ENV allows looking for ENV depending on GROUP/PROJECT/USER
			RUN_ENV=$ENV;
				if [ -s $SCRIPT_DIR/env.$RUN_GROUP.$RUN_PROJECT.$RUN_USER.sh ]; then
					RUN_ENV=$SCRIPT_DIR/env.$RUN_GROUP.$RUN_PROJECT.$RUN_USER.sh;
				elif [ -s $SCRIPT_DIR/env.$RUN_GROUP.$RUN_PROJECT.sh ]; then
					RUN_ENV=$SCRIPT_DIR/env.$RUN_GROUP.$RUN_PROJECT.sh;
				elif [ -s $SCRIPT_DIR/env.$RUN_GROUP.sh ]; then
					RUN_ENV=$SCRIPT_DIR/env.$RUN_GROUP.sh;
				fi;
			#fi;
			RUN_APP=$(echo $(basename $RUN_ENV) | sed "s/^env.//gi" | sed "s/.sh$//gi" | sed "s/sh$//gi")
			if [ "$RUN_APP" == "" ]; then RUN_APP="default"; fi
			
			echo "# DEFAULT APPLICATION detection $RUN_APP";
			
			awk '/Data/{y=1;next}y' $SAMPLESHEET | tr -d '\r' | sed 's/,/\t/g' > ${TMP_FOLDER}/$RUN_SampleSheet.csv_tmp
			SAMPLES_PROJECT_LIST=$( C=1; for i in $(head ${TMP_FOLDER}/$RUN_SampleSheet.csv_tmp -n 1) ; do if [ $i == "Sample_Project" ] ; then break ; else C=$(( $C + 1 )) ; fi ; done ; cut -f $C ${TMP_FOLDER}/$RUN_SampleSheet.csv_tmp | sort -u | sed 's/Sample_Project//' );
			rm ${TMP_FOLDER}/$RUN_SampleSheet.csv_tmp
			echo "# List of APPLICATION: "$(echo $SAMPLES_PROJECT_LIST);
		else
			echo "#[WARNNING] SAMPLESHEET '$SAMPLESHEET' not found";
		fi;
		
		done;

	else
		echo "#";
		echo "# USE $(basename $0) --applications_infos --runs=<RUN> to detect RAW FOLDER and APPLICATION of RUN";
	fi;
		echo "#";
	exit 0;
fi;

# PIPELINES
##############

if (($PIPELINES_INFOS)); then
	PIPELINES_INFOS=$TMP_SYS_FOLDER/$RANDOM$RANDOM.pipelines
	MAKEFILE_PARAM=$TMP_SYS_FOLDER/$RANDOM$RANDOM.param
	touch $PIPELINES_INFOS $MAKEFILE_PARAM
	source $ENV;
	make -e -f $SCRIPT_DIR/NGSWorkflow.mk PARAM=$MAKEFILE_PARAM PIPELINES_INFOS=$PIPELINES_INFOS $PIPELINES_INFOS 1>/dev/null 2>/dev/null;
	sort -k1,2 $PIPELINES_INFOS | cut -d: -f1-3 | column -t -s ':'
	#make -e -f $SCRIPT_DIR/NGSWorkflow.mk PARAM=$MAKEFILE_PARAM RELEASE_INFOS=$PIPELINES_INFOS $PIPELINES_INFOS 1>/dev/null 2>/dev/null;
	#cat $PIPELINES_INFOS
	rm $PIPELINES_INFOS $MAKEFILE_PARAM;
	
	exit 0;
fi;


# RELEASE INFOS
#################

if (($RELEASE_INFOS)); then
	RELEASE_INFOS=$TMP_SYS_FOLDER/$RANDOM$RANDOM.release
	MAKEFILE_PARAM=$TMP_SYS_FOLDER/$RANDOM$RANDOM.param
	touch $RELEASE_INFOS $MAKEFILE_PARAM
	source $ENV;
	make -e -f $SCRIPT_DIR/NGSWorkflow.mk PARAM=$MAKEFILE_PARAM RELEASE_INFOS=$RELEASE_INFOS $RELEASE_INFOS 1>/dev/null 2>/dev/null;
	cat $RELEASE_INFOS
	rm $RELEASE_INFOS $MAKEFILE_PARAM;
	
	exit 0;
fi;
#echo "truc"; exit 0;


# FASTQ/BAM/CRAM/SAM SAMPLE
##################################

if [ ! -z "$FASTQ" ]; then # && [ -e "$FASTQ" ]; then
	# TEST FASTQs if exist
	FASTQ_EXISTS=1
	for F in $(echo $FASTQ | tr "," " "); do
		if [ ! -e $F ] || [ -d $F ]; then FASTQ_EXISTS=0; fi;
	done;
	# Launch FASTQs analmysis if all FASTQs exist
	if (($FASTQ_EXISTS)); then
		echo -e "#\n# Launch SAMPLE Analysis '$FASTQ'\n#"; # $PARAM
		$SCRIPT_DIR/launch.sample.sh $PARAM
		exit 0;
	fi;
fi;

# RUN (if no FASTQ files)
##############################


if [ ! -z "$RUNS" ]; then
	echo -e "#\n# Launch RUN Analysis '$RUNS'\n#"; # $PARAM
	$SCRIPT_DIR/launch.sh $PARAM
	exit 0;
fi;

# If nothing...
##################

usage;

exit 0;
