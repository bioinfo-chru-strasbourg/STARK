#!/bin/bash
#################################
##
## STARK Services Modules 
##
#################################

SCRIPT_NAME="STARKServicesModules"
SCRIPT_DESCRIPTION="STARK Services Modules"
SCRIPT_RELEASE="0.9.1b"
SCRIPT_DATE="01/10/2020"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="HUS/CPS"
SCRIPT_LICENCE="GNU GPLA V3"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-8/05/2020: Script creation\n";

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
	echo "# USAGE: $(basename $0) [--help] [options...]";
	echo "#";
	echo "### This script manages STARK modules and services.";
	echo "#";
  	echo "# --env=<FILE>                 STARK Docker Compose environment file ";
	echo "#                              List of files with environment variables";
	echo "#                              Default: '../env,STARK.env'";
  	echo "# --folder_services=<FOLDER>   STARK services modules folder ";
	echo "#                              Default: in env file, or 'services'";
  	echo "# --modules=<STRING>           STARK modules as a list";
	echo "#                              Format: 'module1,module2,...', use '*' as wildcard";
	echo "#                              Default: '*' all modules";
	echo "# --submodules=<STRING>        STARK submodules as a list";
	echo "#                              Format: 'submodule1,submodule2,...'";
	echo "#                              Default: '' all submodules";
	echo "# --common=<STRING>            List of Common modules (subfolders).";
	echo "#                              Format: 'module1,module2,...', use '*' as wildcard";
	echo "#                              These modules will be started before all other modules";
	echo "#                              Use 'none' to not start any common modules";
	echo "#                              Default: 'COMMON' if command is 'up'";
	echo "# --services=<STRING>          STARK modules services as a list from docker-compose";
	echo "#                              Only for command: up, start, stop, restart'";
	echo "#                              Format: 'service1,service2,...'";
	echo "#                              Default: '*' all services";
	echo "# --command=<STRING>           Docker Compose Command ";
	echo "#                              - up: Create and start containers (in daemon mode)";
	echo "#                              - down: Stop and remove containers, networks, images, and volumes";
	echo "#                              - start: Start services";
	echo "#                              - stop: Stop services";
	echo "#                              - restart: Restart services";
	echo "#                              - config: Check config services";
	echo "#                              Default: 'up'";
	echo "# --modules_show               Show modules and services";
	echo "# --command_args=<STRING>      Docker Compose Command args";
	echo "# -v|--verbose                 Verbose mode";
	echo "# -d|--debug                   Debug mode";
	echo "# -n|--release                 Script Release";
	echo "# -h|--help                    Help message";
	echo "#";

}




####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "vdnh" --long "env:,folder_services:,modules:,submodules:,modules_show,common:,services:,command:,command_args:,verbose,debug,release,help" -- "$@" 2> /dev/null)
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
			ENV=$(echo "$2" | tr "," " ")
			shift 2
			;;
    	--folder_services)
			FOLDER_SERVICES="$2"
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
    	--modules_show)
			MODULES_SHOW=1
			shift 1
			;;
    	--common)
			MODULE_COMMON=$(echo "$2" | tr "," " ")
			shift 2
			;;
  		--services)
			SERVICES=$(echo "$2" | tr "," " ")
			shift 2
			;;
    	--command)
			COMMAND="$2"
			shift 2
			;;
  		--command_args)
			COMMAND_ARGS="$2"
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

# Copy Command
#COMMAND_COPY="rsync -auczqAXhi --no-links --no-perms --no-owner --no-group"
COMMAND_COPY="rsync -az --update"


# TMP
TMP_FOLDER=/tmp/STARK/services.$RANDOM
mkdir -p $TMP_FOLDER


# DOCKER
DOCKER_VERSION=$(docker --version)
(($VERBOSE)) && echo "#[INFO] STARK Module Docker version '$DOCKER_VERSION'"


# ENV
if [ -z "$ENV" ]; then
	ENV="$SCRIPT_DIR/../.env $SCRIPT_DIR/STARK.env"
fi;

ENV_LIST="";
for ENV_FILE in $ENV; do
	if [ -e "$ENV_FILE" ]; then
		ENV_LIST=$ENV_LIST" "$ENV_FILE;
	fi;
done;


if [ ! -z "$ENV_LIST" ]; then
	for ENV_FILE in $ENV_LIST; do
		(($VERBOSE)) && echo "#[INFO] STARK Module Docker Compose environment file '$ENV_FILE'"
		source $ENV_FILE
	done
else
	echo "#[ERROR] STARK Module Docker Compose environment file '$ENV' NOT exist"
	usage;
	exit 1;
fi;



# DOCKER_STARK_MAIN_FOLDER
if [ -z "$DOCKER_STARK_MAIN_FOLDER" ]; then
  if [ "$DOCKER_STARK_MAIN_FOLDER" != "" ] && [ -d $DOCKER_STARK_MAIN_FOLDER ]; then
    DOCKER_STARK_MAIN_FOLDER=$DOCKER_STARK_MAIN_FOLDER
  elif [ -d "~/STARK" ]; then
    DOCKER_STARK_MAIN_FOLDER="~/STARK"
  else
    echo "#[ERROR] STARK Module Main folder"
    exit 1;
  fi;
fi;

(($VERBOSE)) && echo "#[INFO] STARK Module main folder '$DOCKER_STARK_MAIN_FOLDER'"


# SERVICES
if [ -z "$FOLDER_SERVICES" ]; then
  if [ "$DOCKER_STARK_FOLDER_SERVICES" != "" ]; then
    FOLDER_SERVICES=$DOCKER_STARK_MAIN_FOLDER/$DOCKER_STARK_FOLDER_SERVICES
  else
  	FOLDER_SERVICES=$DOCKER_STARK_MAIN_FOLDER/services
  fi;
fi;

mkdir -p $FOLDER_SERVICES
chmod a+x $FOLDER_SERVICES  2>/dev/null

(($VERBOSE)) && echo "#[INFO] STARK Module services modules folder '$FOLDER_SERVICES'"


# MODULE
if [ -z "$MODULES" ]  || [ "$MODULES" == "*" ]; then
  MODULES="*"
fi;


# SERVICES
if [ -z "$SERVICES" ]; then
  SERVICES=''
fi;


# COMMAND
if [ -z "$COMMAND" ]; then
  COMMAND='up'
fi;

[ "$COMMAND" == "up" ] && COMMAND="up -d"
[ "$COMMAND" == "config" ] && SERVICES=""

(($VERBOSE)) && echo "#[INFO] STARK Module Docker Compose Command '$COMMAND'"


# MODULE COMMON
if [ -z "$MODULE_COMMON" ] && [ "$COMMAND" == "up -d" ]; then
	MODULE_COMMON="common"
fi;

(($VERBOSE)) && echo "#[INFO] STARK Module Common services '$MODULE_COMMON'"


# Main STARK prefix

MAIN_MODULE_PREFIX="STARK"



### FUNCTIONS


# list_include_item "10 11 12" "2"
function list_include_item {
  local list="$1"
  local item="$2"
  if [[ $list =~ (^|[[:space:]])"$item"($|[[:space:]]) ]] ; then
    # yes, list include item
    result=0
  else
    result=1
  fi
  return $result
}



### Find services
MODULE_checked=""

for service_module in \
	$( if [ ! -z "$MODULE_COMMON" ]; then echo "$MODULE_COMMON" | sed 's#[^ ]* *#'$SCRIPT_DIR'/&#g' | tr " " "\n"; fi ) \
	$( if [ ! -z "$MODULES" ]; then echo "$MODULES" | sed 's#[^ ]* *#'$SCRIPT_DIR'/&#g' | tr " " "\n"; fi );
	do

	if [ -d $service_module ]; then

		# Module name
		module_name=$(basename $service_module)
		
		# Module check
		if $(list_include_item "$MODULE_checked" "$module_name"); then
			(($DEBUG)) && echo "#[INFO] STARK Module '$module_name' already checked"
			continue
		fi;
		MODULE_checked=$MODULE_checked" $module_name"

		# Module common test
		#(($VERBOSE)) && echo "" && echo "#[INFO] STARK Module '$module_name'"
		(($VERBOSE)) && echo ""
		echo "#[INFO] STARK Module '$module_name'"
		if $(list_include_item "$MODULE_COMMON" "$module_name"); then
			(($VERBOSE)) && echo "#[INFO] STARK Module '$module_name' is a 'common' module"
		fi;
	
		# Services list
		list_services="$MAIN_MODULE_PREFIX"
		# for services_full_path in $(ls $service_module/$MAIN_MODULE_PREFIX.*.docker-compose.yml $service_module/*/$MAIN_MODULE_PREFIX.*.docker-compose.yml 2>/dev/null); do
		# 	#list_services=$list_services" "$(echo $services_full_path | xargs basename | sed 's/.docker-compose.yml$//' )
		# 	list_services=$list_services" "$(echo $services_full_path | sed 's#'$service_module'/##' | sed 's/.docker-compose.yml$//' )
		# done;
		for services_full_path in $(ls $service_module/*/$MAIN_MODULE_PREFIX*.docker-compose.yml 2>/dev/null); do
			#list_services=$list_services" "$(echo $services_full_path | xargs basename | sed 's/.docker-compose.yml$//' )
			list_services=$list_services" "$(echo $services_full_path | sed 's#'$service_module'/##' | sed 's/.docker-compose.yml$//' )
		done;

	
		#for service_prefix in "STARK" $list_services; do
		for service_prefix in $list_services; do


			# Service name
			service_name=$(echo $service_prefix | xargs dirname | sed 's/^'$MAIN_MODULE_PREFIX'//' | sed 's/^\.//')
			service_folder=$service_name
			[[ $service_name == "" ]] && service_name="main" && service_folder=""
			#(($VERBOSE)) && echo "#[INFO] Module '$module_name' Service '$service_prefix/$service_name'"
			#(($VERBOSE)) && echo "#[INFO] STARK Module '$module_name' Service '$service_name'"

			if $(list_include_item "$SUBMODULES" "$service_name") || [ "$SUBMODULES" == "" ]; then

			# Find YML file
			service_module_yml="";
			if [ -e $service_module"/"$service_prefix".docker-compose.yml" ]; then
				service_module_yml=$service_module"/"$service_prefix".docker-compose.yml";
			# elif [ -e $service_module"/"$service_prefix".yml" ]; then
			# 	service_module_yml=$service_module"/"$service_prefix".yml";
			# elif [ -e $service_module"/docker-compose.yml" ]; then
			# 	service_module_yml=$service_module"/docker-compose.yml";
			else
				(($DEBUG)) && echo "#[WARNING] STARK Module '$module_name' folder does not have docker configuration file"
			fi;


			# Find ENV file
			service_module_env="";
			if [ -e $service_module"/"$service_prefix".env" ]; then
				service_module_env=$service_module"/"$service_prefix".env";
			# elif [ -e $service_module"/.env" ]; then
			# 	service_module_env=$service_module"/.env";
			else
				(($DEBUG)) && echo "#[WARNING] STARK Module '$module_name' folder does not have environment file"
			fi;


			# Find MODULE file
			service_module_module="";
			if [ -e $service_module"/"$service_prefix".module" ]; then
				service_module_module=$service_module"/"$service_prefix".module";
			# elif [ -e $service_module"/.module" ]; then
			# 	service_module_module=$service_module"/.module";
			else
				(($DEBUG)) && echo "#[WARNING] STARK Module '$module_name' folder does not have module configuration file"
			fi;


			# Find README file
			service_module_readme="";
			if [ -e $service_module"/"$service_prefix".README.md" ]; then
				service_module_readme=$service_module"/"$service_prefix".README.md";
			# elif [ -e $service_module"/README.md" ]; then
			# 	service_module_readme=$service_module"/README.md";
			else
				(($DEBUG)) && echo "#[WARNING] STARK Module '$module_name' folder does not have README file"
			fi;


			### LAUNCH modules/services
			if [ "$service_module_yml" != "" ] && [ "$service_module_env" != "" ] && [ "$service_module_module" != "" ]; then

				### SHOW modules/services
				if (($MODULES_SHOW)); then

					echo "#[INFO] STARK Module '$module_name' SubModule '$service_name'"

				else

					#(($VERBOSE)) && echo "#[INFO] STARK Module '$module_name'"
					(($VERBOSE)) && echo "#[INFO] STARK Module '$module_name' SubModule '$service_name'"
					(($DEBUG)) && echo "#[INFO] STARK Module '$module_name' SubModule '$service_name' docker configuration file: '$service_module_yml'"
					(($DEBUG)) && echo "#[INFO] STARK Module '$module_name' SubModule '$service_name' environment file: '$service_module_env'"
					(($DEBUG)) && echo "#[INFO] STARK Module '$module_name' SubModule '$service_name' module configuration file: '$service_module_module'"
					(($DEBUG)) && echo "#[INFO] STARK Module '$module_name' SubModule '$service_name' README file: '$service_module_readme'"
					(($DEBUG)) && echo "#[INFO] STARK Module '$module_name' SubModule '$service_name' command '$COMMAND'"

					# Alternative ENV
					(($DEBUG)) && echo "#[INFO] STARK Module '$module_name' - Alternate environment file creation"
					cat $ENV_LIST 1> $TMP_FOLDER/.env 2>$TMP_FOLDER/.env.err
					cat $service_module"/"$MAIN_MODULE_PREFIX".env" 1>> $TMP_FOLDER/.env 2>>$TMP_FOLDER/.env.err
					cat $service_module_env 1>> $TMP_FOLDER/.env 2>>$TMP_FOLDER/.env.err

					# Create folder
					mkdir -p $FOLDER_SERVICES/$module_name/$service_folder
					chmod o+rwx $FOLDER_SERVICES/$module_name 2>/dev/null
					chmod o+rwx $FOLDER_SERVICES/$module_name/$service_folder 2>/dev/null

					# Module configuration file
					(($DEBUG)) && echo "#[INFO] STARK Module '$module_name' - Module configuration file copy"
					$COMMAND_COPY $service_module_module $FOLDER_SERVICES/$module_name/$service_folder 2>/dev/null;

					# Readme file
					(($DEBUG)) && echo "#[INFO] STARK Module '$module_name' - Readme file copy"
					if [ "$service_module_readme" != "" ]; then
						$COMMAND_COPY $service_module_readme $FOLDER_SERVICES/$module_name/$service_folder 2>/dev/null;
					fi;
					
					# Files permissions
					chmod o+rxw -R $FOLDER_SERVICES/$module_name/* 2>/dev/null

					# DEBUG
					(($DEBUG)) && cat $TMP_FOLDER/.env $TMP_FOLDER/.env.err
					(($DEBUG)) && echo "    docker-compose --file $service_module_yml --env-file $TMP_FOLDER/.env -p $module_name $COMMAND"


					(($DEBUG)) && docker-compose --file $service_module_yml --env-file $TMP_FOLDER/.env -p $module_name config

					

					# patch for --env-file unavailable
					if (( $(docker-compose --help | grep "\-\-env\-file" -c) )); then
						
						# Command
						if docker-compose --file $service_module_yml --env-file $TMP_FOLDER/.env -p $module_name config 1>$TMP_FOLDER/docker-compose.log 2>$TMP_FOLDER/docker-compose.err; then
							> $TMP_FOLDER/docker-compose.err

							if [ "$(cat $TMP_FOLDER/docker-compose.log)" == "services: {}" ]; then

								echo "#[WARNING] docker-compose 'empty yml'"
							
							else

								if docker-compose --file $service_module_yml --env-file $TMP_FOLDER/.env -p $module_name $COMMAND $COMMAND_ARGS $SERVICES 2>>$TMP_FOLDER/docker-compose.err; then
									(($VERBOSE)) && cat $TMP_FOLDER/docker-compose.err | grep "Found orphan containers" -v
								else
									if (($(cat $TMP_FOLDER/docker-compose.err | grep "no service selected" -c))); then
										echo "#[WARNING] docker-compose 'no service selected'";
									fi;
									if (($(cat $TMP_FOLDER/docker-compose.err | grep " level=warning " -c))); then
										echo "#[WARNING] docker-compose 'level=warning'";
									fi;
									if (($(cat $TMP_FOLDER/docker-compose.err | grep -e "no service selected" -e " level=warning " -vc))); then
										echo "#[ERROR] docker-compose error1";
										cat $TMP_FOLDER/docker-compose.err
										exit 1;
									fi;
								fi;

							fi;

						else
							echo "#[ERROR] docker-compose error2";
							cat $TMP_FOLDER/docker-compose.log $TMP_FOLDER/docker-compose.err;
							exit 1;
						fi;
						
					else

						(($VERBOSE)) &&  echo "#[WARNING] STARK Module Docker version does not provide '--env-file' parameter. STARK module is patched to work, but please update your docker version."

						# Command
						if env $(cat $TMP_FOLDER/.env | grep "#" -v) docker-compose --file $service_module_yml -p $module_name config 1>$TMP_FOLDER/docker-compose.log 2>$TMP_FOLDER/docker-compose.err; then
							> $TMP_FOLDER/docker-compose.err
							if env $(cat $TMP_FOLDER/.env | grep "#" -v) docker-compose --file $service_module_yml -p $module_name $COMMAND $COMMAND_ARGS $SERVICES 2>>$TMP_FOLDER/docker-compose.out ; then
								(($VERBOSE)) && cat $TMP_FOLDER/docker-compose.err | grep "Found orphan containers" -v
							else
								# echo "#[ERROR] docker-compose error";
								# cat $TMP_FOLDER/docker-compose.err
								# exit 1;
								if (($(cat $TMP_FOLDER/docker-compose.err | grep "no service selected" -c))); then
									echo "#[WARNING] docker-compose 'no service selected'";
								fi;
								if (($(cat $TMP_FOLDER/docker-compose.err | grep " level=warning " -c))); then
									echo "#[WARNING] docker-compose 'level=warning'";
								fi;
								if (($(cat $TMP_FOLDER/docker-compose.err | grep -e "no service selected" -e " level=warning " -vc))); then
									echo "#[ERROR] docker-compose error1";
									cat $TMP_FOLDER/docker-compose.err
									exit 1;
								fi;
							fi;
						else
							echo "#[ERROR] docker-compose error";
							cat $TMP_FOLDER/docker-compose.log $TMP_FOLDER/docker-compose.err;
							exit 1;
						fi;

					fi;

				fi;

			else
				(($VERBOSE)) && echo "#[WARNING] STARK Module '$module_name' is not a docker service module (no STARK*.docker-compose.yml, STARK*.env and STARK*.module)"
			fi;

		else

			(($VERBOSE)) && echo "#[INFO] STARK Module '$module_name' SubModule '$service_name' not selected"

		fi;


		done;

	fi;

done;

# TMP
rm -rf $TMP_FOLDER


exit 0
