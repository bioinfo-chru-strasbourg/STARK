#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="analysis time"
SCRIPT_DESCRIPTION="Calculates analysis time"
SCRIPT_RELEASE="0.9beta"
SCRIPT_DATE="26/03/2015"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"

echo "#######################################";
echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
echo "# $SCRIPT_DESCRIPTION ";
echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © ";
echo "#######################################";


# 1.1. Configuration
#####################

#HELP
HELP="# USAGE: $0 <LOGS>\n"
HELP=$HELP"# LOGS  List of analysis logs"
if [ "${1^^}" == "HELP" ]; then
	echo -e $HELP;
	exit 0;
fi;

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_DIR/env.sh

#CONFIG=$1
#if [ ! -e $CONFIG ] && [ ! -e $NGS_SCRIPTS/$CONFIG ] || [ "$CONFIG" == "" ]; then
#	CONFIG=$NGS_SCRIPTS/config.ini
#fi;

echo "# "
echo "# NGS Folder:            "$NGS_FOLDER
echo "# NGS BIN Folder:        "$NGS_BIN
echo "# NGS SCRIPTS Folder:    "$NGS_SCRIPTS
echo "# MISEQ Folder:          "$MISEQ_FOLDER
echo "# DEMULTIPLEXING Folder: "$DEMULTIPLEXING_FOLDER
echo "# RESULTS Folder:        "$RESULTS_FOLDER
echo "# ANALYSIS Folder:       "$ANALYSIS_FOLDER
echo "# TMP folder:            "$TMP_FOLDER
echo "# CONFIG ini:            "$CONFIG
echo "# "

# 1.2. INPUT
#############

LOGS=$1
if [ "$LOGS" == "" ]; then
	#echo "# [ERROR] No input log files ";
	#echo -e $HELP
	#exit 1;
	echo "# [WARNING] No input log files ";
	LOGS=$TMP_FOLDER$RESULTS_SUBFOLDER/*/*.log
	echo "# LOGS: $LOGS"
fi;


#echo $LOGS;
for LOG in $LOGS; do
	#echo "$LOG";
	LOG_START="`grep "\[.*\] .* Process for RUN '.*' START$" $LOG`"
	LOG_STOP="`grep "\[.*\] .* Process for RUN '.*' STOP$" $LOG``grep "\[.*\] .* Process for RUN '.*' END$" $LOG`"
	LOG_START_TIME="now"
	LOG_STOP_TIME="now"
	RUNNING="RUNNING"
	if [ "$LOG_START" != "" ]; then
		LOG_START_TIME=`echo $LOG_START | cut  -d' ' -f1-6 | tr -d []`
	fi;
	if [ "$LOG_STOP" != "" ]; then
		LOG_STOP_TIME=`echo $LOG_STOP | cut  -d' ' -f1-6 | tr -d []`
		RUNNING="FINNISHED";
	fi;
	#echo "LOG_START_TIME=$LOG_START_TIME LOG_STOP_TIME=$LOG_STOP_TIME";	
	DSTART=`date -d "$LOG_START_TIME" +"%s"`
	DSTOP=`date -d "$LOG_STOP_TIME" +"%s"`
	#echo "DSTART=$DSTART DSTOP=$DSTOP";
	DIFF=$(( $DSTOP - $DSTART ))
	H=$(($DIFF / 3600))
	M=$((($DIFF - ($H * 3600) )/60))
	S=$((($DIFF - ($H * 3600) - ($M * 60) )))
	TIME=$H"h"$M"m"$S"s"
	#echo "$LOG: Exec Time "$TIME
	printf "# %2sh%2sm%2ss %-10s %50s\n" $H $M $S $RUNNING $LOG
done;



