#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARKRunsCopy"
SCRIPT_DESCRIPTION="STARK Runs Copy"
SCRIPT_RELEASE="0.9b"
SCRIPT_DATE="31/03/2020"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="HUS/CPS"
SCRIPT_LICENCE="GNU GPLA V3"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-31/03/2020: Script creation\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Configuration
ENV_CONFIG=$(find -L $SCRIPT_DIR/.. -name config.app)
source $ENV_CONFIG 1>/dev/null 2>/dev/null


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
	echo "# USAGE: $(basename $0) --sources=<STRING> --dest=<FOLDER> [-h] [options...]";
	echo "#";
	echo "### This script copy runs from multiple source folders to destination folder.";
	echo "#";
	echo "# --sources=<STRING>                         List of sources runs folder (required).";
	echo "#                                            Format: 'folder1,folder2,...'";
	echo "# --dest=<FOLDER>                            Destination runs folder (required)";
	echo "#                                            Format: 'folder'";
	echo "# --days=<INTEGER>                           Run folder modification days";
	echo "#                                            Default: '30'";
	echo "# -v|--verbose                               Verbose mode";
	echo "# -d|--debug                                 Debug mode";
	echo "# -n|--release                               Script Release";
	echo "# -h|--help                                  Help message";
	echo "#";
	#echo -e "#\n# RUN Analysis\n################";
	#$STARK_FOLDER_BIN/launch.sh -h | grep "# [ |-]";
	#echo -e "#\n# SAMPLE Analysis\n###################";
	#$STARK_FOLDER_BIN/launch.sample.sh -h | grep "# [ |-]";
}




####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "vdnh" --long "sources:,dest:,days:,verbose,debug,release,help" -- "$@" 2> /dev/null)
if [ $? -ne 0 ]; then
	:
	echo "#[ERROR] Error in the argument list:";
	echo "#[ERROR] $@"
	echo ""
	usage;
	exit;
fi;




PARAM=$@
DEBUG=0
VERBOSE=0

eval set -- "$ARGS"
while true
do
	#echo "$1=$2"
	#echo "Eval opts";
	case "$1" in
		--sources)
			SOURCES_RUNS=$(echo $2 | tr "," " ")
			shift 2
			;;
		--dest|--destination)
			DEST_RUNS="$2"
			shift 2
			;;
		--days)
			DAYS="$2"
			shift 2
			;;
		-h|--help)
			usage
			exit 0
			;;
		-v|--verbose)
			VERBOSE=1
			shift 1
			;;
		-d|--debug)
			DEBUG=1
			shift 1
			;;
		--) shift
			break
			;;
		*) 	echo "Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done

# header
(($NO_HEADER)) || header;



####################################################################################################################################
# Checking the input parameter
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if [ -z "$SOURCE_RUNS" ] && [ -z "$DEST_RUNS" ] && ((!$DEBUG)); then
	echo "#[ERROR] Required parameter: --sources and --dest. Use --help to display the help." && echo "" && usage && exit 1;
fi
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




# DAYS
########

if [ "$DAYS" == "" ]; then
	DAYS="30"
fi;




# RUNS COPY
#############


for SOURCE_RUNS in $SOURCES_RUNS; do 

	#echo "SOURCE_RUNS=$SOURCE_RUNS"

	if [ -d $SOURCE_RUNS ]; then
		for RUN_SOURCE_FOLDER in $(find -L $SOURCE_RUNS -type d -mindepth 1 -maxdepth 1 -mtime -$DAYS); do
			
			RUN=$(basename $RUN_SOURCE_FOLDER)
			RUN_DEST_FOLDER=$DEST_RUNS/$RUN

			echo $RUN_SOURCE_FOLDER
			echo $RUN_DEST_FOLDER

			EXCLUDE=" --exclude '*Images' "

			(($DEBUG)) && echo "rsync -avz --exclude 'SampleSheet.csv' $EXCLUDE $RUN_SOURCE_FOLDER/ $RUN_DEST_FOLDER/"
			(($DEBUG)) && echo "rsync -avz --include 'SampleSheet.csv' $EXCLUDE $RUN_SOURCE_FOLDER/ $RUN_DEST_FOLDER/"

			if ((1)); then
				#mkdir -p $RUN_DEST_FOLDER
				rsync -avz --exclude 'SampleSheet.csv' $EXCLUDE $RUN_SOURCE_FOLDER/ $RUN_DEST_FOLDER/
				# #ls -l $RUN_SOURCE_FOLDER/RTAComplete.txt
				[ -e $RUN_DEST_FOLDER/RTAComplete.txt ] && rsync -avz --include 'SampleSheet.csv' $EXCLUDE $RUN_SOURCE_FOLDER/ $RUN_DEST_FOLDER/
			fi;

		done;
	fi;

done;


