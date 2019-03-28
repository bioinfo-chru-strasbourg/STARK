#!/bin/bash
#################################
## STARK environment
#################################


source_app () {
# source an env file from an application definition
# $1: APP definition
# $2: APPS folder
# return APP env file or default APP env file or null

	local APP_LIST=$(echo $1 | tr "," " ") FOLDER_APPS=$2

	if [ "$FOLDER_APPS" == "" ]; then FOLDER_APPS="$STARK_FOLDER_APPS"; fi
	if [ "$FOLDER_APPS" == "" ]; then FOLDER_APPS=".."; fi

	#echo "opt $APP_LIST F $FOLDER_APPS"

	[ "$CONFIG_HEADER" != "" ] && [ -e "$CONFIG_HEADER" ] && source $CONFIG_HEADER #1>/dev/null 2>/dev/null;

	#source $(find_app "$APP_LIST" "$FOLDER_APPS")
	local ENV=$(find_app "$APP_LIST" "$FOLDER_APPS")
	#echo "ENV=$ENV"
	for ENV_ONE in $ENV; do
		#echo "ENV_ONE=$ENV_ONE"
		local APP_FOLDER=$(dirname $ENV_ONE)
		#echo "APP_FOLDER=$APP_FOLDER"
		#echo "STARK_FOLDER_APPS $STARK_FOLDER_APPS"
		source $ENV_ONE #1>/dev/null 2>/dev/null;
		#echo "E=$E OK"
	done
	#source $ENV
	#echo "ENV=$ENV loaded"

	#for APP in $APP_LIST; do
	#	local ENV=$(find_app "$APP" "$FOLDER_APPS")
	#	[ "$ENV" != "" ] && [ -e "$ENV" ] && source $ENV;
	#	#source $ENV
	#
	#done;

	[ "$CONFIG_FOOTER" != "" ] && [ -e "$CONFIG_FOOTER" ] && source $CONFIG_FOOTER #1>/dev/null 2>/dev/null;

} # source_app


name_app () {
# source an env file from an application definition
# $1: APP definition
# $2: APPS folder
# return APP env file or default APP env file or null

	local APP_LIST=$(echo $1 | tr "," " ") FOLDER_APPS=$2

	#source $(find_app "$APP_LIST" "$FOLDER_APPS")
	local ENV=$(find_app "$APP_LIST" "$FOLDER_APPS")
	APP_LIST_NAME=""
	for E in $ENV; do
		APP_NAME_DEF=$(source_app "$E"; echo $APP_NAME)
		APP_RELEASE_DEF=$(source_app "$E"; echo $APP_RELEASE)
		APP_LIST_NAME=$APP_LIST_NAME" "$APP_NAME_DEF":"$APP_RELEASE_DEF
	done

	echo $APP_LIST_NAME
	return 0;

} # source_app



find_app () {
# find a env file from an application definition
# $1: APP definition
# $2: APPS folder
# return APP env file or default APP env file or null
# Usage find_app "APP" "FOLDER_APPS"

	local APP_LIST=$(echo $1 | tr "," " ") FOLDER_APPS=$2

	# default APP

	if [ "$APP_LIST" == "" ]; then

		return 0

	elif [ $(echo "$APP_LIST" | wc -w) -gt 1 ]; then


		APP_TO_RETURN=""
		for APP_ONE in $APP_LIST; do
			#echo "APP_ONE=$APP_ONE"
			APP_TO_RETURN=$APP_TO_RETURN" "$(find_app "$APP_ONE" "$FOLDER_APPS")
		done;
		echo $APP_TO_RETURN #| sed "s/^,//";
		return 0;

	else
		local APP=$APP_LIST

		# APP is a app file
		[ -f $APP ] && echo $APP && return 0

		local APP_DIRNAME=$(dirname $APP 2>/dev/null | tr -d ".")
		local APP_BASENAME=$(basename $APP 2>/dev/null)

		if [ "$FOLDER_APPS" == "" ] && [ "$APP_DIRNAME" == "" ]; then FOLDER_APPS=".."; fi

		local APP_tr=$(echo $APP | tr "-" ".")
		local APP_BASENAME_tr=$(echo $APP_BASENAME | tr "-" ".")

		local FIND_APP_FILE=$(find $FOLDER_APPS/$APP_DIRNAME -name $APP_BASENAME -type f)
		[ ! -z $FIND_APP_FILE ] && echo $FIND_APP_FILE && return 0

		local FIND_APP_SUB_FILE=$(find $FOLDER_APPS/$APP_DIRNAME -name $APP_BASENAME.app -or -name $APP_BASENAME_tr.app -name $APP_BASENAME.plugapp -or -name $APP_BASENAME_tr.plugapp | sort -t'/' -k2.2r -k2.1 | head -n1 )
		[ ! -z "$FIND_APP_SUB_FILE" ] && echo $FIND_APP_SUB_FILE && return 0

		for APP_TEST in $(find $FOLDER_APPS -name *.app -or -name *.plugapp); do

			local APP_NAME=$(source $APP_TEST 2>/dev/null; echo $APP_NAME);
			[[ $APP = $APP_NAME ]] && echo $APP_TEST && return 0
			[[ $APP_tr = $APP_NAME ]] && echo $APP_TEST && return 0
			[[ $APP_BASENAME = $APP_NAME ]] && echo $APP_TEST && return 0
			[[ $APP_BASENAME_tr = $APP_NAME ]] && echo $APP_TEST && return 0

		done;

	fi;

	#echo "not found"
	#echo $APP_DEFAULT && return 0
	return 1

} #function find_app


convertsecs() {
 ((h=${1}/3600))
 ((m=(${1}%3600)/60))
 ((s=${1}%60))
 printf "%02d"h"%02d"m"%02d"s"\n" $h $m $s
}