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

# Configuration
ENV_CONFIG=$(find $SCRIPT_DIR/.. -name config.app)
source $ENV_CONFIG


# Header
function header () {
	cat $STARK_FOLDER_DOCS/HEADER
	#echo "#######################################";
	#echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
	#echo "# $SCRIPT_DESCRIPTION ";
	#echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © $SCRIPT_LICENCE";
	#echo "#######################################";
}

# Release
function release () {
	cat $STARK_FOLDER_DOCS/RELEASE_NOTES
	#echo "# RELEASE NOTES:";
	#echo -e $RELEASE_NOTES
}

# Usage
function usage {
	echo "# USAGE: $(basename $0) --runs=<RUN>|--fastq=<FASTQ>|--ana=<FILE> [options...]";
	echo "# -e|--application|--app|--env=<APP|FILE>     APP name or APP file configuration of the APPLICATION.";
	echo "#                                             Must be in the STARK APPS folder if relative path";
	echo "#                                             Default: defined in the RUN SampleSheet, or default.app if not defined";
	echo "# -f|--fastq|--fastq_R1=<FILE1,FILE2,...>     List of FASTQ|BAM|SAM|CRAM (mandatory if no --runs or --config options)";
	echo "#                                             Formats: *fastq.gz|*fq.gz|*bam|*ubam|*cram|*ucram|*sam|*usam";
	echo "# -r|--runs=<STRING1,STRING2,...>             List of RUNs to analyse, from Illumina sequencers (mandatory if no --fastq or --config options).";
	echo "#                                             RUNS will be automatically searched in all ENV RUN folders, if necessary";
	echo "#                                             Format: RUN1,RUN2...";
	echo "# -w|--ana|--analysis|--json=<FILE>           Config file in JSON format listing options (mandatory if no --fastq or --runs options).";
	echo "# -l|--samplesheet=<FILE>                     Illumina SampleSheet.csv file";
	echo "#                                             Default: found in RUN folder";


	echo "# -i|--application_infos                      Applications informations.";
	echo "# -y|--pipelines_infos                        Pipelines informations.";
	echo "# -g|--release_infos                          Pipelines with tools and databases information.";
	echo "# -x|--tools_infos                            Tools release information.";

	echo "# -v|--verbose                                VERBOSE option";
	echo "# -d|--debug                                  DEBUG option";
	echo "# -n|--release                                RELEASE option";
	echo "# -h|--help                                   HELP option";
	echo "#";
	echo -e "#\n# RUN Analysis\n################";
	$STARK_FOLDER_BIN/launch.sh -h | grep "# [ |-]";
	echo -e "#\n# SAMPLE Analysis\n###################";
	$STARK_FOLDER_BIN/launch.sample.sh -h | grep "# [ |-]";
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
ARGS=$(getopt -o "e:r:f:l:w:t:iygxvdnh" --long "env:,app:,application:,runs:,fastq:,fastq_R1:,json:,ana:,analysis:,samplesheet:,application_infos,applications_infos,pipelines_infos,release_infos,tools_infos,verbose,debug,release,help" -- "$@" 2> /dev/null)
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
		-e|--env|--app|--application)
			APP="$2"
			shift 2
			;;
		-r|--runs)
			RUNS="$2"
			shift 2
			;;
		-f|--fastq|--fastq_R1)
			FASTQ=$2 #"$(echo $2 | tr "\n" " ")"
			shift 2
			;;
		-l|--samplesheet)
			SAMPLESHEET="$2" #"$(echo $2 | tr "\n" " ")"
			shift 2
			;;
		-i|--application_infos|--applications_infos)
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
		-x|--tools_infos)
			TOOLS_INFOS=1
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
		-w|--json|--ana|--analysis)
			JSON="$2"
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
if [ -z "$FASTQ" ] && [ -z "$RUNS" ] && [ -z "$JSON" ] && [ -z "$APPLICATIONS_INFOS" ] && [ -z "$PIPELINES_INFOS" ] && [ -z "$RELEASE_INFOS" ] && [ -z "$TOOLS_INFOS" ]; then
	echo "Option --runs or --fastq or --applications_infos or --pipelines_infos or --release_infos or --tools_infos is required. " "Use -h or --help to display the help." && usage && exit 1;
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


# DEBUG
#########
if (($DEBUG)); then
	echo $STARK_FOLDER_BIN
	echo $STARK_FOLDER_APPS
	echo $STARK_FOLDER_DOCS
	echo $STARK_FOLDER_TOOLBOX
	echo $STARK_FOLDER_CONFIG
	echo $STARK_FOLDER
fi;

# ENV
#########

ENV=$(find_app "$APP" "$STARK_FOLDER_APPS")

# SOURCE ENV if exists
if [ ! -z $ENV ] && [ -s $ENV ]; then
	echo "#[INFO] APPLICATION file '"$(echo $ENV | sed s#$STARK_FOLDER_APPS#APPS#)"' found."
	source $ENV;
else
	echo "#[ERROR] NO APPLICATION file '$ENV' found."
	exit 1
fi;

EXIT=0

# CONFIG
##########


if [ -e $JSON ] && [ "$JSON" != "" ]; then

	echo "#[INFO] Config file '$JSON' does exists."
	PARAM_JSON=$(echo "$PARAM " | sed "s/--analysis=[^ |$]*//gi" | sed "s/--json=[^ |$]*//gi" | sed "s/--ana=[^ |$]*//gi" | sed "s/-w [^ |$]*//gi")

	PARAM_STARK=$($STARK_FOLDER_BIN/json_to_options.py $JSON)
	if [ $? -eq 0 ]; then
	    (($DEBUG)) && echo $PARAM_STARK 
	else
	    echo "#[ERROR] Analysis configuration file '$JSON' failed"
	    exit 1;
	fi

	CMD=$(echo $0" "$PARAM_STARK" $PARAM_JSON")
	echo "#[INFO] CMD=$CMD"
	echo ""

	eval $CMD;

	EXIT=1

elif [ "$JSON" != "" ]; then

	echo "[WARNING] Analysis configuration file '$JSON' does NOT exist"

fi;


# APPLICATIONS_INFOS
#######

if (($APPLICATIONS_INFOS)); then

	if [ ! -e $SAMPLESHEET ]; then
		echo "#[ERROR] SampleSheet '$SAMPLESHEET' does NOT exists."
		exit 0;
	fi;


	#SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

	RUNS_TEST=$RUNS

	echo ""
	echo "#[INFO] ######################"
	echo "#[INFO] #### APPLICATIONS ####"
	echo "#[INFO] ######################"
	echo ""
	RAW_LIST="";
	FOLDER_RUN_ORIGINAL=$FOLDER_RUN
	#echo "TRUC "$(echo $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/env.*sh)
	#for ENV_DEF in $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/env.*sh; do

	#for ENV_DEF in $STARK_FOLDER_APPS/*app; do
	for ENV_DEF in $(find $STARK_FOLDER_APPS -name *app); do
		# FIND APP
		#APP=$($STARK_FOLDER_BIN/extract_variable_from_env.sh $ENV_DEF "APP_NAME");
		APP=$(source $ENV_DEF; echo $APP_NAME)
		#APP=$(echo $(basename $ENV_DEF) | sed "s/^env.//gi" | sed "s/.sh$//gi" | sed "s/sh$//gi")
		if [ "$APP" == "" ]; then
			APP="default";
			ENV=$APP;
		fi
		# SOURCE
		FOLDER_RUN=$($STARK_FOLDER_BIN/extract_variable_from_env.sh $ENV_DEF "FOLDER_RUN" 2>/dev/null)
		#RAW_FOLDER=$($STARK_FOLDER_BIN/extract_variable_from_env.sh $ENV_DEF "RAW_FOLDER")
		#MISEQ_FOLDER=$($STARK_FOLDER_BIN/extract_variable_from_env.sh $ENV_DEF "MISEQ_FOLDER");
		FOLDER_RESULTS=$($STARK_FOLDER_BIN/extract_variable_from_env.sh $ENV_DEF "FOLDER_RESULTS" 2>/dev/null)
		PIPELINES=$($STARK_FOLDER_BIN/extract_variable_from_env.sh $ENV_DEF "PIPELINES" 2>/dev/null)
		POST_ALIGNMENT=$($STARK_FOLDER_BIN/extract_variable_from_env.sh $ENV_DEF "POST_ALIGNMENT" 2>/dev/null)
		POST_ALIGNMENT_STEPS=$($STARK_FOLDER_BIN/extract_variable_from_env.sh $ENV_DEF "POST_ALIGNMENT_STEPS" 2>/dev/null)
		# OUTPUT
		#if [ "$RUNS_TEST" == "" ]; then
		echo "#[INFO] ## $APP ["$(echo $ENV_DEF | sed s#$STARK_FOLDER_APPS#APPS#)"]";
		echo "#[INFO]    RAW DATA $FOLDER_RUN";
		echo "#[INFO]    RES DATA $MAIN_FOLDER";
		echo "#[INFO]    PIPELINES "
		echo -e "#[INFO]       $(echo $PIPELINES | sed "s/ /\n#[INFO]       /gi")";
		echo "#[INFO]    POST_ALIGNMENT $POST_ALIGNMENT "
		echo "#[INFO]    POST_ALIGNMENT_STEPS $POST_ALIGNMENT_STEPS "
		#fi;
		RAW_LIST=$RAW_LIST" $FOLDER_RUN"
		
	done;
	RAW_LIST=$(echo $RAW_LIST | sed "s/ /\n/gi" | sort | uniq | sed "s/\n/ /gi" );

	#echo "RAW_LIST=$RAW_LIST"

	if [ "$RUNS_TEST" != "" ]; then

		for RUN_TEST in $(echo $RUNS_TEST | tr "," " "); do

		echo "";
		echo "#[INFO] ### RUN '$RUN_TEST' ####";
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
			echo "#[INFO] RAW FOLDER $MISEQ_FOLDER_FOUND";
		fi;

		# Test SampleSheet to determine RUN_TEST possible Application/ENV
		if [ -z $SAMPLESHEET ]; then
			SAMPLESHEET=$MISEQ_FOLDER_FOUND/$RUN_TEST/SampleSheet.csv
		fi;

		if [ -e $SAMPLESHEET ]; then
			RUN_GROUP=`grep -i 'Investigator Name' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2 | cut -d \- -f 1`
			RUN_PROJECT=`grep -i 'Investigator Name' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2 | cut -d \- -f 2`
			RUN_USER=`grep -i 'Investigator Name' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2 | cut -d \- -f 3`
			RUN_APP_NAME="$RUN_GROUP-$RUN_PROJECT"
			
			# RUN ENV
			# ENV-AUTO in ENV allows looking for ENV depending on GROUP/PROJECT/USER
			RUN_ENV=$ENV;

			# SOURCE ENV if exists
			if [ ! -z $ENV ] && [ -s $ENV ]; then
				echo "#[INFO] APPLICATION file '"$(echo $ENV | sed s#$STARK_FOLDER_APPS#APPS#)"' found."
				source $ENV;
			else
				echo "#[WARNING] NO APPLICATION file '$ENV' found."
			fi;
			
			if ((0)); then
				if [ -s $STARK_FOLDER_BIN/env.$RUN_GROUP.$RUN_PROJECT.$RUN_USER.sh ]; then
					RUN_ENV=$STARK_FOLDER_BIN/env.$RUN_GROUP.$RUN_PROJECT.$RUN_USER.sh;
				elif [ -s $STARK_FOLDER_BIN/env.$RUN_GROUP.$RUN_PROJECT.sh ]; then
					RUN_ENV=$STARK_FOLDER_BIN/env.$RUN_GROUP.$RUN_PROJECT.sh;
				elif [ -s $STARK_FOLDER_BIN/env.$RUN_GROUP.sh ]; then
					RUN_ENV=$STARK_FOLDER_BIN/env.$RUN_GROUP.sh;
				fi;
			fi;

			RUN_ENV=$(find_app "$RUN_APP_NAME" "$STARK_FOLDER_APPS")
			[ "$RUN_ENV" == "" ] && RUN_ENV=$ENV;
			
			RUN_APP=$(source $RUN_ENV; echo $APP_NAME)
			if [ "$RUN_APP" == "" ]; then RUN_APP="default"; fi

			echo "#[INFO] DEFAULT APPLICATION detection $RUN_APP";

			awk '/Data/{y=1;next}y' $SAMPLESHEET | tr -d '\r' | sed 's/,/\t/g' > ${TMP_FOLDER}/$RUN_SampleSheet.csv_tmp
			SAMPLES_PROJECT_LIST=$( C=1; for i in $(head ${TMP_FOLDER}/$RUN_SampleSheet.csv_tmp -n 1) ; do if [ $i == "Sample_Project" ] ; then break ; else C=$(( $C + 1 )) ; fi ; done ; cut -f $C ${TMP_FOLDER}/$RUN_SampleSheet.csv_tmp | sort -u | sed 's/Sample_Project//' );
			rm ${TMP_FOLDER}/$RUN_SampleSheet.csv_tmp
			echo "#[INFO] List of APPLICATION: "$(echo $SAMPLES_PROJECT_LIST);
			for SPL in $SAMPLES_PROJECT_LIST; do
				echo "#[INFO] - '$SPL' ["$(find_app "$SPL" "$STARK_FOLDER_APPS" | sed s#$STARK_FOLDER_APPS#APP#)"]";
			done;
		else
			echo "#[WARNNING] SAMPLESHEET '$SAMPLESHEET' not found";
		fi;

		done;

	else
		echo "";
		echo "#[INFO] USE $(basename $0) --applications_infos --runs=<RUN> to detect RAW FOLDER and APPLICATION of RUN";
	fi;
		echo "";

	EXIT=1
	#exit 0;
fi;

# PIPELINES
##############

if (($PIPELINES_INFOS)); then

	echo ""
	echo "###################"
	echo "#### PIPELINES ####"
	echo "###################"
	echo ""

	PIPELINES_INFOS=$TMP_SYS_FOLDER/$RANDOM$RANDOM.pipelines
	MAKEFILE_PARAM=$TMP_SYS_FOLDER/$RANDOM$RANDOM.param
	touch $PIPELINES_INFOS $MAKEFILE_PARAM
	#source $ENV;
	make -e -f $STARK_FOLDER_BIN/NGSWorkflow.mk PARAM=$MAKEFILE_PARAM PIPELINES_INFOS=$PIPELINES_INFOS $PIPELINES_INFOS 1>/dev/null 2>/dev/null;
	sort -k1,2 $PIPELINES_INFOS | cut -d: -f1-3 | column -t -s ':'
	#make -e -f $STARK_FOLDER_BIN/NGSWorkflow.mk PARAM=$MAKEFILE_PARAM RELEASE_INFOS=$PIPELINES_INFOS $PIPELINES_INFOS 1>/dev/null 2>/dev/null;
	#cat $PIPELINES_INFOS
	rm $PIPELINES_INFOS $MAKEFILE_PARAM;

	EXIT=1
fi;



# TOOLS INFOS
#################

if (($TOOLS_INFOS)); then

	echo ""
	echo "####################"
	echo "#### TOOLS LIST ####"
	echo "####################"
	echo ""
	for T in $TOOLS_LIST; do
		VAR_VERSION=$T"_VERSION"; #echo ${!VAR_VERSION};
		VAR_DESCRIPTION=$T"_DESCRIPTION"; #echo ${!VAR_DESCRIPTION};
		VAR_REF=$T"_REF"; #echo ${!VAR_REF};
		echo "## $T ["${!VAR_VERSION}"]";
		echo "# "${!VAR_DESCRIPTION};
		echo "# "${!VAR_REF};
		echo ""
	done

	EXIT=1
fi;

# RELEASE INFOS
#################

if (($RELEASE_INFOS)); then

	echo ""
	echo "##################"
	echo "#### RELEASES ####"
	echo "##################"
	echo ""

	RELEASE_INFOS=$TMP_SYS_FOLDER/$RANDOM$RANDOM.release
	MAKEFILE_PARAM=$TMP_SYS_FOLDER/$RANDOM$RANDOM.param
	touch $RELEASE_INFOS $MAKEFILE_PARAM
	source $ENV;
	make -e -f $STARK_FOLDER_BIN/NGSWorkflow.mk PARAM=$MAKEFILE_PARAM RELEASE_INFOS=$RELEASE_INFOS $RELEASE_INFOS 1>/dev/null 2>/dev/null;
	cat $RELEASE_INFOS | sed "s/^## /\n## /"

	if ((0)); then
	echo ""
	echo ""
	echo "##############"
	echo "# TOOLS LIST #"
	echo "##############"
	echo ""
	for T in $TOOLS_LIST; do
		VAR_VERSION=$T"_VERSION"; #echo ${!VAR_VERSION};
		VAR_DESCRIPTION=$T"_DESCRIPTION"; #echo ${!VAR_DESCRIPTION};
		VAR_REF=$T"_REF"; #echo ${!VAR_REF};
		echo "## $T ["${!VAR_VERSION}"]";
		echo "# "${!VAR_DESCRIPTION};
		echo "# "${!VAR_REF};
		echo ""
	done
	fi;

	rm $RELEASE_INFOS $MAKEFILE_PARAM;

	EXIT=1
fi;
#echo "truc"; exit 0;



if (($EXIT)); then
	exit 0;
fi;

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
		$STARK_FOLDER_BIN/launch.sample.sh $PARAM
		exit 0;
	fi;
fi;

# RUN (if no FASTQ files)
##############################


if [ ! -z "$RUNS" ]; then
	echo -e "#\n# Launch RUN Analysis '$RUNS'\n#"; # $PARAM
	$STARK_FOLDER_BIN/launch.sh $PARAM
	exit 0;
fi;

# If nothing...
##################

usage;

exit 0;
