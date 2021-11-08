#!/bin/bash
#################################
##
## STARK Docker Setup 
##
#################################

SCRIPT_NAME="STARKDockerSetup"
SCRIPT_DESCRIPTION="STARK Docker Setup"
SCRIPT_RELEASE="0.9.1.1"
SCRIPT_DATE="28/10/2021"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="HUS/CPS"
SCRIPT_LICENCE="GNU GPLA V3"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-12/04/2020: Script creation\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.1.1-28/10/2021: Script fixes\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# DOCKER_STARK_MAIN_FOLDER_DEFAULT
DOCKER_STARK_MAIN_FOLDER_DEFAULT=${HOME}/STARK

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
	echo "# --env=<FILE>                     Dockerfile environment configuration file ";
	echo "#                                  Default: '.env'";
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
		--env)
			ENV="$2"
			shift 2
			;;
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


(($DEBUG)) && VERBOSE=1


## PARAM
##########



## ENV
if [ -z "$ENV" ]; then
	ENV=".env";
fi;

echo "#[INFO] STARK Docker environment configuration file '$ENV' "

if [ -s $ENV ]; then
	source $ENV;
else
	echo "#[WARNING] STARK Docker environment configuration file '$ENV' does not exist"
	#exit 1;
fi;


# DOCKER_STARK_MAIN_FOLDER
if [ -s $DOCKER_STARK_MAIN_FOLDER ]; then
	DOCKER_STARK_MAIN_FOLDER=$DOCKER_STARK_MAIN_FOLDER_DEFAULT;
	echo "#[INFO] STARK Main folder '$DOCKER_STARK_MAIN_FOLDER' by default"
fi
if mkdir -p $DOCKER_STARK_MAIN_FOLDER; then
	echo "#[INFO] STARK Main folder '$DOCKER_STARK_MAIN_FOLDER' created"
else
	echo "#[ERROR] STARK Main folder '$DOCKER_STARK_MAIN_FOLDER' NOT created "
	exit 1;
fi;

## LOG
DOCKER_STARK_SETUP_LOG="$DOCKER_STARK_MAIN_FOLDER/setup.log"
echo "### STARK Setup" >> $DOCKER_STARK_SETUP_LOG 2>> $DOCKER_STARK_SETUP_LOG
echo -e '# STARK Setup DATE: '`date '+%Y%m%d-%H%M%S'` >> $DOCKER_STARK_SETUP_LOG 2>> $DOCKER_STARK_SETUP_LOG
echo -e '\n' >> $DOCKER_STARK_SETUP_LOG 2>> $DOCKER_STARK_SETUP_LOG
if [ -s $DOCKER_STARK_SETUP_LOG ]; then
	echo "#[INFO] STARK setup log file '$DOCKER_STARK_SETUP_LOG' created "
else
	echo "#[ERROR] STARK setup log file '$DOCKER_STARK_SETUP_LOG' NOT created "
	exit 1;
fi;

(($DEBUG)) && tail $DOCKER_STARK_SETUP_LOG


## GIT CLONE
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


# Git Clone
echo "#[INFO] STARK GIT Clone '$GIT_CLONE_URL'..."
if [ "$GIT_CLONE_URL" != "" ]; then
	#echo "#[INFO] STARK GIT Clone "
	if git clone $GIT_CLONE_URL >> $DOCKER_STARK_SETUP_LOG 2>> $DOCKER_STARK_SETUP_LOG; then
		cd STARK
		echo "#[INFO] STARK GIT Clone '$GIT_CLONE_URL' done."
	else
		echo "#[ERROR] STARK GIT Clone '$GIT_CLONE_URL' failed!"
		exit 1;
	fi;
else
	echo "#[INFO] STARK GIT Clone skipped."
fi;


# Build
echo "#[INFO] STARK Docker Compose Build..."
if docker-compose build >> $DOCKER_STARK_SETUP_LOG 2>> $DOCKER_STARK_SETUP_LOG; then
	echo "#[INFO] STARK Docker Compose Build done."
else
	echo "#[ERROR] STARK Docker Compose Build failed!"
	exit 1;
fi;


# Setup
echo "#[INFO] STARK Docker Compose Setup "
#source .env
#mkdir -p $DOCKER_STARK_MAIN_FOLDER

# Folder creation
echo "#[INFO] STARK Docker Compose Setup - Folders creation..."
if docker-compose --project-name STARK up stark-setup >> $DOCKER_STARK_SETUP_LOG 2>> $DOCKER_STARK_SETUP_LOG; then
	echo "#[INFO] STARK Docker Compose Setup - Folders creation done."
else
	echo "#[ERROR] STARK Docker Compose Setup - Folders creation failed!"
	exit 1;
fi;

# Databases download
echo "#[INFO] STARK Docker Compose Setup - Databases Download..."
if docker-compose --project-name STARK up stark-databases >> $DOCKER_STARK_SETUP_LOG 2>> $DOCKER_STARK_SETUP_LOG; then
	echo "#[INFO] STARK Docker Compose Setup - Databases Download done."
else
	echo "#[ERROR] STARK Docker Compose Setup - Databases Download failed!"
	exit 1;
fi;

# Sources archives
echo "#[INFO] STARK Docker Compose Setup - Sources Archives..."
if docker-compose --project-name STARK up stark-sources-archives >> $DOCKER_STARK_SETUP_LOG 2>> $DOCKER_STARK_SETUP_LOG; then
	echo "#[INFO] STARK Docker Compose Setup - Sources Archives done."
else
	echo "#[ERROR] STARK Docker Compose Setup - Sources Archives failed!"
	exit 1;
fi;

echo "#[INFO] STARK Docker Compose Setup done."


# Start services Build
echo "#[INFO] STARK Docker Compose Build - Services Modules..."
if $SCRIPT_DIR/services/services.sh --modules=* --command=build --verbose >> $DOCKER_STARK_SETUP_LOG 2>> $DOCKER_STARK_SETUP_LOG; then
	echo "#[INFO] STARK Docker Compose Build - Services Modules done."
else
	echo "#[ERROR] STARK Docker Compose Build - Services Modules failed!"
	exit 1;
fi;

# Start services Start
echo "#[INFO] STARK Docker Compose Start - Services Modules..."
if $SCRIPT_DIR/services/services.sh --modules=* --command=up --verbose >> $DOCKER_STARK_SETUP_LOG 2>> $DOCKER_STARK_SETUP_LOG; then
	echo "#[INFO] STARK Docker Compose Start - Services Modules done."
else
	echo "#[ERROR] STARK Docker Compose Start - Services Modules failed!"
	exit 1;
fi;


# Output Message
echo "#[INFO] STARK installed - Read README.md for more information"

exit 0
