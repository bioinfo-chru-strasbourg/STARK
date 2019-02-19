#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="ExtractVariableFromEnv"
SCRIPT_DESCRIPTION="Extract Variable From Env"
SCRIPT_RELEASE="0.9b"
SCRIPT_DATE="14/04/2016"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

#echo "#######################################";
#echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
#echo "# $SCRIPT_DESCRIPTION ";
#echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © $SCRIPT_LICENCE";
#echo "#######################################";

# Realse note
#RELEASE_NOTES="#\n"
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-14/04/2016: Creation script\n";

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
	echo "# USAGE: $(basename $0) <ENV> <VARIABLE> ";
	echo "# ENV         Environnement file (define folders, tools...)";
	echo "# VARIABLE    Variable name to extract";
	echo "#";
	#echo "# DB_INTEGRATION   integrate data into database";
	#echo "Usage : SELF <CONFIG> <RUNS> <SAMPLES>";

	exit 0;
fi;


# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# INPUT
ENV=$1
RUNS=$2

if [ "$ENV" != "" ] && [ -e $ENV ]; then
	source $1 1>/dev/null 2>/dev/null
	echo ${!2}
fi;


