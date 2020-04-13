#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARKDockerSetup"
SCRIPT_DESCRIPTION="STARK Docker Setup"
SCRIPT_RELEASE="0.9b"
SCRIPT_DATE="12/04/2020"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="HUS/CPS"
SCRIPT_LICENCE="GNU GPLA V3"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-12/04/2020: Script creation\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# GIT_CLONE_DEFAULT
GIT_CLONE_DEFAULT="https://gitlab.bioinfo-diag.fr/Strasbourg/STARK.git"


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
	echo "# USAGE: $(basename $0) [-h] [options...]";
	echo "#";
	echo "### This script setup STARK envronment.";
	echo "### Can be used with curl (e.g. curl https://gitlab.bioinfo-diag.fr/Strasbourg/vision/raw/master/setup.sh | bash)";
	echo "###    (e.g. curl https://gitlab.bioinfo-diag.fr/Strasbourg/STARK/raw/master/setup.sh | bash)";
	echo "#";
	echo "# --git-clone=<STRING>             Download STARK code on GIT ";
	echo "#                                     - 'auto': Will detect '.git' folder to check if GIT clone is needed, if not 'default'";
	echo "#                                     - '0': Will not GIT clone STARK code and consider current directory as STARK code folder ";
	echo "#                                     - 'default': Will DO GIT clone STARK code from default GIT repository (BioInfoDiag GitLab) ";
	echo "#                                     - <URL>: Will DO GIT clone STARK code from GIT URL ";
	echo "#                                  Default: 'auto' (detect '.git' folder to check if GIT clone is needed)";
	echo "#                                  BioInfoDiag GitLab URL: $GIT_CLONE_DEFAULT";
	echo "# -v|--verbose                     Verbose mode";
	echo "# -d|--debug                       Debug mode";
	echo "# -n|--release                     Script Release";
	echo "# -h|--help                        Help message";
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
ARGS=$(getopt -o "vdnh" --long "git-clone:,verbose,debug,release,help" -- "$@" 2> /dev/null)
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
		--git-clone)
			GIT_CLONE="$2"
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
# if [ -z "$SOURCE_RUNS" ] && [ -z "$DEST_RUNS" ] && ((!$DEBUG)); then
# 	echo "#[ERROR] Required parameter: --sources and --dest. Use --help to display the help." && echo "" && usage && exit 1;
# fi
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



## PARAM
##########


## ENV
if [ -z "$GIT_CLONE" ]; then
	GIT_CLONE="auto";
fi;


GIT_CLONE_URL=""
if [ "$GIT_CLONE" == "0" ]; then
	GIT_CLONE_URL=""
elif [ "$GIT_CLONE" == "default" ]; then
	GIT_CLONE_URL=$GIT_CLONE_DEFAULT;
elif [ "$GIT_CLONE" == "auto" ]; then
	if (($(ls -d .git 2>/dev/null | wc -l))); then
		GIT_CLONE_URL=""
	else
		GIT_CLONE_URL=$GIT_CLONE_DEFAULT;
	fi;
else
	GIT_CLONE_URL=$GIT_CLONE;
fi;


(($VERBOSE)) && echo "#[INFO] GIT URL: $GIT_CLONE_URL "


# Git Clone
if [ "$GIT_CLONE_URL" != "" ]; then
	(($VERBOSE)) && echo "#[INFO] GIT Clone "
	git clone $GIT_CLONE_URL
	cd STARK
fi;


# Build
(($VERBOSE)) && echo "#[INFO] Docker Compose Build "
docker-compose build


# Setup
(($VERBOSE)) && echo "#[INFO] Docker Compose Setup "
source .env
mkdir -p $DOCKER_STARK_MAIN_FOLDER
docker-compose up stark-folders
docker-compose up stark-databases


# Start
(($VERBOSE)) && echo "#[INFO] Docker Compose Start "
docker-compose up -d


# Output Message
echo "#[INFO] Open 'http://localhost:$DOCKER_STARK_SERVICE_PORT_PATTERN$DOCKER_STARK_SERVICE_DASHBOARD_PORT' in your browser "
echo "#[INFO] or run 'bin/STARK --help' "


exit 0

