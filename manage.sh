#!/bin/bash
#################################
##
## STARK Docker Manage 
##
#################################

SCRIPT_NAME="STARKDockerManage"
SCRIPT_DESCRIPTION="STARK Docker Manage"
SCRIPT_RELEASE="0.9.0.0"
SCRIPT_DATE="19/05/2021"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="HUS/CPS"
SCRIPT_LICENCE="GNU GPLA V3"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.0.0-19/05/2021: Script creation\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


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
	echo "### This script manage STARK envronment.";
	echo "#";
	echo "# --env=<FILE>                     STARK Docker environment configuration file ";
	echo "#                                  Default: '.env'";
	echo "# --env_modules=<FILE>             STARK Docker environment configuration file for modules";
	echo "#                                  Default: --env option";
	echo "# --command=<FILE>                 STARK Docker environment command";
	echo "#                                  Default: 'start'";
	echo "#                                  Options: 'start', 'stop'";
	echo "# --modules=<STRING>               STARK modules as a list";
	echo "# --modules=<STRING>               STARK modules as a list";
	echo "#                                  Format: 'module1,module2,...', use '*' as wildcard";
	echo "#                                  Default: '*' all modules";
	echo "# --submodules=<STRING>            STARK submodules as a list";
	echo "#                                  Format: 'submodule1,submodule2,...'";
	echo "#                                  Default: '' all submodules";
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
ARGS=$(getopt -o "vdnh" --long "env:,env_modules:,command:,modules:,submodules:,verbose,debug,release,help" -- "$@" 2> /dev/null)
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
		--env_modules)
			ENV_MODULES="$2"
			shift 2
			;;
		--command)
			COMMAND="$2"
			shift 2
			;;
		--modules)
			MODULES=$(echo "$2" | tr "," " ")
			shift 2
			;;
    	--submodules)
			SUBMODULES=$(echo "$2" | tr "," " ")
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

DOCKER_COMPOSE_VERBOSITY=""
STARK_SERVICES_VERBOSITY=""
(($VERBOSE)) && DOCKER_COMPOSE_VERBOSITY="$DOCKER_COMPOSE_VERBOSITY --verbose "
(($VERBOSE)) && STARK_SERVICES_VERBOSITY="$STARK_SERVICES_VERBOSITY --verbose "
(($DEBUG)) && STARK_SERVICES_VERBOSITY="$STARK_SERVICES_VERBOSITY --verbose "

# Switch off docker-compose verbosity
DOCKER_COMPOSE_VERBOSITY=""


## PARAM
##########



## ENV
if [ -z "$ENV" ]; then
	if [ -e "$SCRIPT_DIR/.env" ]; then
		ENV="$SCRIPT_DIR/.env";
	elif [ -e ".env" ]; then
		ENV=".env";
	fi;
fi;

(($VERBOSE)) && echo "#[INFO] STARK Docker environment configuration file '$ENV' "

if [ -s "$ENV" ]; then
	source $ENV;
else
	echo "#[ERROR] STARK Docker environment configuration file '$ENV' failed"
	exit 1;
fi;


## ENV
if [ -z "$ENV_MODULES" ]; then
	ENV_MODULES=$ENV;
fi;

(($VERBOSE)) && echo "#[INFO] STARK Docker environment configuration file '$ENV' "


## STARK bin folder (from $ENV)
STARK_BIN_FOLDER=$(dirname $ENV)
if [ -d "$STARK_BIN_FOLDER" ]; then
	(($VERBOSE)) && echo "#[INFO] STARK Docker environment bin folder '$STARK_BIN_FOLDER' "
else
	echo "#[ERROR] STARK Docker environment bin folder '$STARK_BIN_FOLDER' failed"
fi;

## STARK Modules bin folder (from $ENV)
STARK_MODULES_BIN_FOLDER=$(dirname $ENV_MODULES)
if [ -d "$STARK_MODULES_BIN_FOLDER" ]; then
	(($VERBOSE)) && echo "#[INFO] STARK Docker environment bin modules folder '$STARK_MODULES_BIN_FOLDER' "
else
	echo "#[ERROR] STARK Docker environment bin modules folder '$STARK_MODULES_BIN_FOLDER' failed"
fi;

# MODULE
if [ -z "$MODULES" ]  || [ "$MODULES" == "*" ]; then
  MODULES=""
fi;

# SUBMODULES
if [ -z "$SUBMODULES" ]  || [ "$SUBMODULES" == "*" ]; then
  SUBMODULES=""
fi;




## COMMAND
if [ -z "$COMMAND" ]; then
	COMMAND="start";
fi;
(($VERBOSE)) && echo "#[INFO] STARK Docker environment command '$COMMAND' "




if [ "$COMMAND" != "" ]; then

	(($VERBOSE)) && echo "#[INFO] STARK Docker environment - $COMMAND"

	if [ "$COMMAND" == "start" ]; then

		# STARK Core
		(($VERBOSE)) && echo "#[INFO] STARK Docker environment - $COMMAND - STARK Core"
		if (( $(docker-compose --help | grep "\-\-env\-file" -c) )); then
			if docker-compose $DOCKER_COMPOSE_VERBOSITY --file=$STARK_BIN_FOLDER/docker-compose.yml --env-file=$ENV -p STARK up -d; then
				(($VERBOSE)) && echo "#[INFO] STARK Docker environment - $COMMAND - STARK Core - done."
			else
				echo "#[INFO] STARK Docker environment - $COMMAND - STARK Core - failed!!!"
			fi;
		else
			if env $(cat $ENV | grep "#" -v) docker-compose $DOCKER_COMPOSE_VERBOSITY --file=$STARK_BIN_FOLDER/docker-compose.yml -p STARK up -d; then
				(($VERBOSE)) && echo "#[INFO] STARK Docker environment - $COMMAND - STARK Core - done."
			else
				echo "#[INFO] STARK Docker environment - $COMMAND - STARK Core - failed!!!"
			fi;
		fi;

		# STARK Modules
		(($VERBOSE)) && echo "#[INFO] STARK Docker environment - $COMMAND - STARK Modules - STARK"
		if $STARK_MODULES_BIN_FOLDER/services/services.sh --modules="STARK" --command=up $STARK_SERVICES_VERBOSITY; then
			(($VERBOSE)) && echo "#[INFO] STARK Docker environment - $COMMAND - STARK Modules - STARK -  done."
		else
			echo "#[INFO] STARK Docker environment - $COMMAND - STARK Modules - STARK - failed!!!"
		fi;
		(($VERBOSE)) && echo "#[INFO] STARK Docker environment - $COMMAND - STARK Modules"
		if $STARK_MODULES_BIN_FOLDER/services/services.sh --modules="$MODULES" --submodules="$SUBMODULES" --command=up $STARK_SERVICES_VERBOSITY; then
			(($VERBOSE)) && echo "#[INFO] STARK Docker environment - $COMMAND - STARK Modules - done."
		else
			echo "#[INFO] STARK Docker environment - $COMMAND - STARK Modules - failed!!!"
		fi;

	elif [ "$COMMAND" == "stop" ]; then

		# STARK Modules
		(($VERBOSE)) && echo "#[INFO] STARK Docker environment - $COMMAND - STARK Modules"
		if $STARK_MODULES_BIN_FOLDER/services/services.sh --modules="$MODULES" --submodules="$SUBMODULES" --command=down $STARK_SERVICES_VERBOSITY; then
			(($VERBOSE)) && echo "#[INFO] STARK Docker environment - $COMMAND - STARK Modules - done."
		else
			echo "#[INFO] STARK Docker environment - $COMMAND - STARK Modules - failed!!!"
		fi
		if $STARK_MODULES_BIN_FOLDER/services/services.sh --modules="common" --command=down $STARK_SERVICES_VERBOSITY; then
			(($VERBOSE)) && echo "#[INFO] STARK Docker environment - $COMMAND - STARK Modules - 'common' - done."
		else
			echo "#[INFO] STARK Docker environment - $COMMAND - STARK Modules - 'common' - failed!!!"
		fi
		if $STARK_MODULES_BIN_FOLDER/services/services.sh --modules="STARK" --command=down $STARK_SERVICES_VERBOSITY; then
			(($VERBOSE)) && echo "#[INFO] STARK Docker environment - $COMMAND - STARK Modules - STARK - done."
		else
			echo "#[INFO] STARK Docker environment - $COMMAND - STARK Modules - STARK - failed!!!"
		fi

		# STARK Core
		(($VERBOSE)) && echo "#[INFO] STARK Docker environment - $COMMAND - STARK Core"
		if (( $(docker-compose --help | grep "\-\-env\-file" -c) )); then
			if docker-compose $DOCKER_COMPOSE_VERBOSITY --file=$STARK_BIN_FOLDER/docker-compose.yml --env-file=$ENV -p STARK down; then
				(($VERBOSE)) && echo "#[INFO] STARK Docker environment - $COMMAND - STARK Core - done."
			else
				echo "#[INFO] STARK Docker environment - $COMMAND - STARK Core - failed!!!"
			fi;
		else
			if env $(cat $ENV | grep "#" -v) docker-compose $DOCKER_COMPOSE_VERBOSITY --file=$STARK_BIN_FOLDER/docker-compose.yml -p STARK down; then
				(($VERBOSE)) && echo "#[INFO] STARK Docker environment - $COMMAND - STARK Core - done."
			else
				echo "#[INFO] STARK Docker environment - $COMMAND - STARK Core - failed!!!"
			fi;
		fi;

	else

		echo "#[ERROR] STARK Docker environment - $COMMAND NOT found"

	fi;

else

	echo "#[ERROR] STARK Docker environment - NO command found"

fi;


exit 0
