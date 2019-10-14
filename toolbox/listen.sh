#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="Listen"
SCRIPT_DESCRIPTION="Listen RUN for Analysis"
SCRIPT_RELEASE="0.9.4.1b"
SCRIPT_DATE="18/04/2016"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"

echo "#######################################";
echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
echo "# $SCRIPT_DESCRIPTION ";
echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © ";
echo "#######################################";

# Realse note
#RELEASE_NOTES="#\n"
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-09/10/2015: Script creation\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.1b-15/10/2015: Add DO option allowing to only show actions to do\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.2b-23/02/2016: Add CONDITIONS option testing if a RUN is ready to be analysed\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.3b-03/03/2016: Add DAYS option limiting analysis to newer RUNs\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.4b-22/03/2016: Bug correction on DAYS option. Changing default DAYS option to 30\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.4.1b-18/04/2016: Change DAYS option to apply to CONDITIONS\n";

# NOTES
if [ "${1^^}" == "RELEASE" ] || [ "${1^^}" == "RELEASE_NOTES" ] || [ "${1^^}" == "NOTES" ]; then
	echo "# RELEASE NOTES:";
	echo -e $RELEASE_NOTES
	exit 0;
fi;


# 1.1. Configuration
#####################

#HELP
if [ "${1^^}" == "HELP" ]; then
	echo "# USAGE:  $0 <ENV> <DO> <CONDITIONS> <DAYS>";
	echo "# ENV          Environnement file (define folders, tools...)";
	echo "# DO           Do actions (launch, remove files...). Default FALSE, i.e. just print info about RUN analysis (DONE, NOT DONE, RUNNING...)";
	echo "# CONDITIONS   Files needed for a RUN to be analysed. Default 'RTAComplete.txt SampleSheet.csv' in RUN folder";
	echo "# DAYS         Consider only RUN (and CONDITIONS) younger than DAYS days. Default 30 (1 month ago)";
	exit 0;
fi;

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# INPUT
ENV=$1
DO=$2
CONDITIONS_INPUT=$3
DAYS=$4

# ENV
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
if [ "$ENV" != "" ]; then
	source $ENV;
fi;

# SHOW
if [ "$DO" == "" ] || [ "$DO" == "0" ]; then
	DO=0;
else
	DO=1;
fi;

# CONDITIONS
# Conditions from input option
if [ "$CONDITIONS_INPUT" != "" ]; then
	CONDITIONS=$CONDITIONS_INPUT;
# Conditions from ENV
elif  [ "$RUN_ANALYSIS_CONDITIONS" ]; then
	CONDITIONS=$RUN_ANALYSIS_CONDITIONS;
# Default Conditions
else
	CONDITIONS="RTAComplete.txt SampleSheet.csv"
fi;

# DAYS
re='^[0-9]+$'
if [[ ! $DAYS =~ $re ]] || [ "$DAYS" == "" ]; then
	DAYS=30
fi;
SECONDS=$((60*60*24*$DAYS))


##########
# LISTEN #
##########

#echo "ls $MISEQ_FOLDER/1510*"
for RUN_FOLDER in $(ls -rd $MISEQ_FOLDER/*); do



	if [ -d $RUN_FOLDER ] && (($(( (`date +%s` - `stat -L --format %Y $RUN_FOLDER `) < $SECONDS )))); then
		RUN=$(basename $RUN_FOLDER)
		echo "# RUN '$RUN'"

		RTA_COMPLETE_FILE=$RUN_FOLDER/RTAComplete.txt
		SAMPLESHEET_FILE=$RUN_FOLDER/SampleSheet.csv
		STARK_QUEUED_FILE=$RESULTS_FOLDER/$RUN/$STARK_QUEUED
		STARK_RUNNING_FILE=$RESULTS_FOLDER/$RUN/$STARK_RUNNING
		STARK_COMPLETE_FILE=$RESULTS_FOLDER/$RUN/$STARK_COMPLETE


		#CONDITIONS="$RTA_COMPLETE_FILE $SAMPLESHEET_FILE"

		# TEST complete RUN CONDITIONS
		CHECK_CONDITION=1
		for CONDITION in $CONDITIONS; do
			echo "# Test condition '$CONDITION'"
			#if [ -e $RUN_FOLDER/$CONDITION ] || [ -e $CONDITION ]; then
			if [ -e $RUN_FOLDER/$CONDITION ] && (($(( (`date +%s` - `stat -L --format %Y $RUN_FOLDER/$CONDITION `) < $SECONDS )))); then
				echo "# $CONDITION OK";
			elif [ -e $CONDITION ] && (($(( (`date +%s` - `stat -L --format %Y $CONDITION `) < $SECONDS )))); then
				echo "# $CONDITION OK";
			else
				echo "# $CONDITION ko due to not existing or to older";
				CHECK_CONDITION=0
			fi;
		done;


		#DONE
		if (($CHECK_CONDITION)) && [ -e $STARK_COMPLETE_FILE ]; then
			echo "#    DONE."
			DONE=1
		fi;
		#TODO
		if (($CHECK_CONDITION)) && [ ! -e $STARK_COMPLETE_FILE ]; then
			echo "#    NOT DONE."
			NOT_DONE=1
		fi;
		#UNQUEUED TEST
		if (($CHECK_CONDITION)) && [ ! -e $STARK_COMPLETE_FILE ]; then
			if [ -e $STARK_RUNNING_FILE ] && [ $($COULSON/running.sh | grep -c "$RUN") -eq 0 ]; then
				echo "#    UNRUNNING!!!"
				UNRUNNING=1
				if (($DO)); then
					rm -f $STARK_RUNNING_FILE;
				fi;
				if [ -e $STARK_QUEUED_FILE ] && [ $($COULSON/queue.sh | grep -c "$(cat $STARK_QUEUED_FILE | tr "[" ".")") -eq 0 ]; then
					# NOT QUEUED in COULCON !!!!
					echo "#    UNQUEUED!!!"
					UNQUEUED=1
					if (($DO)); then
						rm -f $STARK_QUEUED_FILE
					fi;
				fi;
			else
				if  [ ! -e $STARK_RUNNING_FILE ] && [ -e $STARK_QUEUED_FILE ] && [ $($COULSON/queue.sh | grep -c "$(cat $STARK_QUEUED_FILE | tr "[" ".")") -eq 0 ]; then
					# NOT QUEUED in COULCON !!!!
					echo "#    UNQUEUED!!!"
					UNQUEUED=1
					if (($DO)); then
						rm -f $STARK_QUEUED_FILE
					fi;
				fi;
			fi;

		fi;
		#TOLAUNCH
		if (($DO)); then
			if (($UNQUEUED)); then
				echo "#    LAUNCH."
			fi;
		fi
		if (($CHECK_CONDITION)) && [ ! -e $STARK_COMPLETE_FILE ] && [ ! -e $STARK_QUEUED_FILE ]; then
			echo "#    LAUNCH."
			LAUNCH=1
			if (($DO)); then
				mkdir -p $RESULTS_FOLDER/$RUN
				# Create STARKQueued.txt
				STARK_QUEUED_FILE_TEXT="#["`date '+%Y%m%d-%H%M%S'`"] RUN $RUN queued on STARK ($STARK_VERSION) with COULSON ($COULSON_VERSION)"
				echo $STARK_QUEUED_FILE_TEXT > $STARK_QUEUED_FILE
				#sleep 1;
				# Queuing using COULSON with STARK LAUNCH and remove STARKQueued.txt and create STARKComplete.txt
				#$COULSON/show.sh
				CMD_STARK="$SCRIPT_DIR/launch.sh --runs \"$RUN\"; rm $STARK_QUEUED_FILE; #$STARK_QUEUED_FILE_TEXT "
				CMD_COULSON="$COULSON/add.sh '' '"$CMD_STARK"'"
				eval $CMD_COULSON
			fi;
		fi;
		#QUEUED
		if (($CHECK_CONDITION)) && [ ! -e $STARK_COMPLETE_FILE ]; then
			#cat $STARK_RUNNING_FILE;
			if [ -e $STARK_RUNNING_FILE ] && [ $($COULSON/running.sh | grep -c "$RUN") -gt 0 ]; then
				echo "#    RUNNING...";
				RUNNING=1
			elif [ -e $STARK_QUEUED_FILE ] && [ $($COULSON/queue.sh | grep -c "$(cat $STARK_QUEUED_FILE | tr "[" ".")") -gt 0 ]; then # IN QUEUE
				echo "#    QUEUED."
				QUEUED=1
			else
				# NOT IN QUEUE in COULSON !?!
				echo "#    UNQUEUED!!!"
				UNQUEUED=1
				#rm -f $STARK_QUEUED_FILE
			fi;
		fi;
		#KILLED
		#if (($CHECK_CONDITION)) && [ ! -e $STARK_COMPLETE_FILE ] && [ ! -e $STARK_COMPLETE_FILE ]; then
		#	echo "#    NOT DONE."
		#	# Nothing to do
		#fi;
	fi;




done;
