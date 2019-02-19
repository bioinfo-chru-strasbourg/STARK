#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="MonitorLive"
SCRIPT_DESCRIPTION="Monitor Analysis (TMP folder)"
SCRIPT_RELEASE="0.9.1beta"
SCRIPT_DATE="10/12/2015"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

HEADER="#######################################
# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]
# $SCRIPT_DESCRIPTION 
# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © 
#######################################";
echo "$HEADER";

# Realse note
#RELEASE_NOTES="#\n"
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.3b-09/10/2015: Add RUN configuration depending on Group and Project defined in RUN SampleSheet\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.4b-13/10/2015: Add Complete and Running flag files\n";

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
	echo "# USAGE:  $0 <ENV> <NB_RUN> <NB_LOG> <DETAILS> <TAIL> <FILTER> <SLEEP> <LAST>";
	echo "# ENV             Environnement file (define folders, tools...)";
	echo "# NB_RUN		Number of log to show (default ALL)";
	echo "# NB_LOG		Number of log to show (default 1)";
	echo "# DETAILS		1 for more details. Default 1";
	echo "# TAIL		Number of last log lines to show. default 0";
	echo "# FILTER		Filter RUN on state (RTA_COMPLETE, QUEUED, RUNNING, COMPLETE). default RUNNING";
	echo "# SLEEP		Time to refresh. default 10";
	echo "# LAST		Show LAST COMPLETE run analysis (LAST=NB_RUN, by ENV). default 0";
	exit 0;
fi;

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# INPUT
ENV=$1
NB_RUN=$2
NB_LOG=$3
DETAILS=$4
TAIL=$5
FILTER=$6
SLEEP=$7
LAST=$8

# ENV
: '
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
'

# ENV
#source $SCRIPT_DIR/env.sh
if [ "$ENV" == "" ]; then ENV="ALL"; fi;
if [ "$ENV" != "" ] && [ "$ENV" != "ALL" ] && [ -s $ENV ]; then
	source $ENV;
fi;

# NB_RUN
if [ "$NB_RUN" == "" ]; then
	NB_RUN=100000000000000000
fi;

# NB_LOG
if [ "$NB_LOG" == "" ]; then
	NB_LOG=1
fi;

# DETAILS
if [ "$DETAILS" == "" ]; then
	DETAILS=1
fi;

# TAIL
if [ "$TAIL" == "" ]; then
	TAIL=10
fi;

# FILTER
if [ "$FILTER" == "" ]; then
	FILTER="RUNNING"
fi;

# SLEEP
if [ "$SLEEP" == "" ]; then
	SLEEP=0
fi;

# SLEEP
if [ "$LAST" == "" ]; then
	LAST=3
fi;

# LAST
LAST=0
if [ "$ENV" != "ALL" ]; then LAST=3; else LAST=2; fi;

#$SCRIPT_DIR/monitor.sh "ALL" "1" "1" "-1" "" "COMPLETE" "0"; exit 0

while ((1)); do
	#clear;
	#$SCRIPT_DIR/monitor.sh "$ENV" "$NB_RUN" "$NB_LOG" "$DETAILS" "$FILTER"
	if [ -e $COULSON/queue.sh ]; then 
		COULSON_RUNNING=$($COULSON/running.sh | grep -v "^#" | sed '/^\s*$/d' | grep -v "^$" );
		COULSON_QUEUE=$($COULSON/queue.sh | grep -v "^#" | sed '/^\s*$/d' | grep -v "^$");
	fi;
	if (($LAST)); then TMPERR=$($SCRIPT_DIR/monitor.sh "$ENV" "$LAST" "1" "-1" "" "COMPLETE" "0" 2>/dev/null); fi;
	#TMP=$($SCRIPT_DIR/monitor.sh "ALL" "$NB_RUN" "1" "1" "10" "$FILTER" "0")
	TMP=$($SCRIPT_DIR/monitor.sh "$ENV" "$NB_RUN" "1" "1" "10" "$FILTER" "0")
	clear;
	echo "$HEADER";
	if [ -e $COULSON/queue.sh ]; then echo -e "#\n##### COULSON #####\n#"; echo -ne "$COULSON_RUNNING"; echo -ne "$COULSON_QUEUE"; fi;
	if (($LAST)); then echo -e "#\n##### LAST COMPLETE ANALYSES #####\n#"; echo -e "$TMPERR"; echo "#"; fi;
	echo -e "#\n##### RUNNING ANALYSES #####\n#";
	echo -e "$TMP"
	sleep $SLEEP;
done;

#$NGS_TOOLS/stark/$RELEASE/bin/monitor.live.sh "$ENV" "$NB_RUN" "$NB_LOG" "$DETAILS" "$TAIL" "$FILTER" "0" "$LAST"



