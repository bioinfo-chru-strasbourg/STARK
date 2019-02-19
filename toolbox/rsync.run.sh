#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="RSYNC"
SCRIPT_DESCRIPTION="RSYNC"
SCRIPT_RELEASE="0.9.5b"
SCRIPT_DATE="27/07/2016"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

echo "#######################################";
echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
echo "# $SCRIPT_DESCRIPTION ";
echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © ";
echo "#######################################";

# Realse note
#RELEASE_NOTES="#\n"
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-21/01/2015: Create script\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.4b-13/10/2015: Add Complete and Running flag files\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.5b-27/07/2016: Add conditions on SOURCE and DEST\n";

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
	echo "# USAGE:  $0 <ENV> <SOURCE> <DEST> <RUN_FILTER> <CONDITIONS> <PERMS>";
	echo "# ENV              Environnement file (define folders, tools...)";
	echo "# SOURCE           Source of runs";
	echo "# DEST             destination of runs";
	echo "# RUN_FILTER       run filter";
	echo "# CONDITIONS       Files to be present. format : 'file1 file2...'";
	echo "# PERMS            Permissions on DEST. format XXX ex 775";
	echo "# EXCLUDE          Exclude patterns";
	exit 0;
fi;

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# INPUT
ENV=$1
SOURCE=$2
DEST=$3
RUN_FILTER=$4
CONDITIONS=$5
PERMS=$6
EXCLUDE=$7

COMMAND_COPY="rsync -auv --progress --acls --no-perms --no-owner --no-group " 

# ENV
if [ -s $ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$ENV;
elif [ -s $APPS/$ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$APPS/$ENV;
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

# ENV
#source $SCRIPT_DIR/env.sh
if [ "$ENV" != "" ]; then
	source $ENV;
fi;


if [ ! -e $SOURCE ] ||  [ "$SOURCE" == "" ]; then
	SOURCE=$MISEQ_FOLDER
	echo "#[WARNING] NO SOURCE folder. Default '$RAW_FOLDER' used ";
fi;
if [ ! -e $SOURCE ] ||  [ "$SOURCE" == "" ]; then
	echo "#[ERROR] NO SOURCE folder";
	exit 1;
fi;

if [ "$DEST" == "" ]; then
	echo "#[ERROR] NO DEST folder";
	exit 1;
fi;

mkdir -p $DEST


echo "# ENV:        $ENV"
echo "# SOURCE:     $SOURCE"
echo "# DEST:       $DEST"
echo "# RUN_FILTER: $RUN_FILTER"
echo "# CONDITIONS: $CONDITIONS"


#echo "ls -d $SOURCE/$RUN_FILTER"
#ls -d $SOURCE/$RUN_FILTER

for RUN_FOLDER in $(ls -d $SOURCE/$RUN_FILTER); do # 1>/dev/null 2>/dev/null 
	echo "#";
	RUN_NAME=$(basename $RUN_FOLDER)
	echo "# RUN '"$(basename $RUN_FOLDER)"'"
	CHECK_CONDITION=1
	for CONDITION in $CONDITIONS; do
		echo "# Test condition '$RUN_FOLDER/$CONDITION'"
		if [ -e $RUN_FOLDER/$CONDITION ]; then 
			echo "# $CONDITION OK";
		else
			echo "# $CONDITION ko";
			CHECK_CONDITION=0
		fi;
	done;
	#echo "	CHECK_CONDITION=$CHECK_CONDITION";
	CHECK_CONDITION_DEST=$(echo $CONDITIONS | wc -w);
	for CONDITION in $CONDITIONS; do
		echo "# Test condition '$DEST/$RUN_NAME/$CONDITION'"
		if [ -e $DEST/$RUN_NAME/$CONDITION ]; then 
			echo "# $CONDITION OK";
		else
			echo "# $CONDITION ko";
			CHECK_CONDITION_DEST=0
		fi;
	done;
	#echo "	CHECK_CONDITION_DEST=$CHECK_CONDITION_DEST";

	if (($CHECK_CONDITION)) && ((! $CHECK_CONDITION_DEST)); then
		echo "# RSYNC $RUN_FOLDER to $DEST ";
		echo "# $COMMAND_COPY $RUN_FOLDER $DEST"
		CONDITIONS_PATTERN=$(for CONDITION in $CONDITIONS; do echo " --exclude=$CONDITION "; done;)
		EXCLUDE_PATTERN=$(for EXC in $EXCLUDE; do echo " --exclude=$EXC "; done;)
		$COMMAND_COPY $RUN_FOLDER $DEST $EXCLUDE_PATTERN $CONDITIONS_PATTERN
		$COMMAND_COPY $RUN_FOLDER $DEST $EXCLUDE_PATTERN
		if [ "$PERMS" != "" ]; then
			echo "# chmod $PERMS $DEST/$RUN_NAME -R"
			chmod $PERMS $DEST/$RUN_NAME -R
		fi;
	fi;

done;

echo "#";

exit 0;


