
#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARKDockerToSingularity"
SCRIPT_DESCRIPTION="STARK Docker to Singularity"
SCRIPT_RELEASE="0.9b"
SCRIPT_DATE="06/04/2020"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="HUS/CPS"
SCRIPT_LICENCE="GNU GPLA V3"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-06/04/2020: Script creation\n";

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
	echo "# USAGE: $(basename $0) [-h] [options...]";
	echo "#";
	echo "### This script copy runs from multiple source folders to destination folder.";
	echo "#";
	echo "# --docker-env-file=<FILE>                   Docker environment file";
	echo "#                                            Default: '../.env'";
	echo "# --docker-stark-image=<STRING>              Docker STARK image.";
	echo "#                                            Default: $DOCKER_STARK_IMAGE or 'stark:latest'";
	echo "# --output=<FILE>                            Singularity output file.";
	echo "#                                            Default: Docker STARK image 'name_release.simg' in current work directory";
	echo "# -v|--verbose                               Verbose mode";
	echo "# -d|--debug                                 Debug mode";
	echo "# -n|--release                               Script Release";
	echo "# -h|--help                                  Help message";
	echo "#";

}




####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "vdnh" --long "docker-env-file:,docker-stark-image:,output-folder:,output:,verbose,debug,release,help" -- "$@" 2> /dev/null)
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
		--docker-env-file)
			DOCKER_ENV_FILE="$2"
			shift 2
			;;
		--docker-stark-image)
			DOCKER_STARK_IMAGE_INPUT="$2"
			shift 2
			;;
		--output-folder)
			OUTPUT_FOLDER="$2"
			shift 2
			;;
		--output)
			OUTPUT_FILE="$2"
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
if [ -z "$DOCKER_ENV_FILE" ]; then
	if [ -e $SCRIPT_DIR/../.env ]; then
		DOCKER_ENV_FILE=$SCRIPT_DIR/../.env
	fi;
fi;

[ ! -z "$DOCKER_ENV_FILE" ] && [ -e "$DOCKER_ENV_FILE" ] && source $DOCKER_ENV_FILE #&& echo "sourced"


## DOCKER_STARK_IMAGE_INPUT
if [ -z "$DOCKER_STARK_IMAGE_INPUT" ]; then
	DOCKER_STARK_IMAGE_INPUT="stark:latest"
fi;


## OUTPUT_FOLDER
if [ -z "$OUTPUT_FOLDER" ]; then
	OUTPUT_FOLDER=/tmp/$RANDOM
fi;
mkdir -p $OUTPUT_FOLDER


## OUTPUT_FOLDER
if [ -z "$OUTPUT_FILE" ]; then
	OUTPUT_FILE=$(pwd)/$(echo $DOCKER_STARK_IMAGE_INPUT | tr ":" "_").simg
fi;



## PARAM output
################



(($VERBOSE)) && echo "#[INFO] Docker environment file: $DOCKER_ENV_FILE"
(($VERBOSE)) && echo "#[INFO] Docker STARK image: $DOCKER_STARK_IMAGE_INPUT"
(($VERBOSE)) && echo "#[INFO] Singularity output file: $OUTPUT_FILE"



## PARAM test
###############

if [ -z "$DOCKER_STARK_IMAGE_INPUT" ] || [ ! -d "$OUTPUT_FOLDER" ]; then
	echo "#[ERROR] Required parameter: --docker-stark-image and --output. Use --help to display the help." && echo "" && usage && exit 1;
fi;



## APP
########

CMD="
docker run \
-v /var/run/docker.sock:/var/run/docker.sock \
-v $OUTPUT_FOLDER:/output \
--privileged -t --rm \
singularityware/docker2singularity \
$DOCKER_STARK_IMAGE_INPUT
"

(($DEBUG)) && echo $CMD

eval $CMD

mv $(ls $OUTPUT_FOLDER/*simg) $OUTPUT_FILE

rm -r $OUTPUT_FOLDER

exit 0



