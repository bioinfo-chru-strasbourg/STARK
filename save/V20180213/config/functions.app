#!/bin/bash
#################################
## STARK environment
#################################


find_app () {
# find a env file from an application definition
# $1: APP definition
# $2: APPS folder
# return APP env file or default APP env file or null

	local APP=$1
	local APP_tr=$(echo $APP | tr "-" ".")
	local FOLDER_APPS=$2
	if [ ! -d $FOLDER_APPS ] || [ "$FOLDER_APPS" == "" ]; then FOLDER_APPS=".."; fi
	
	# default APP
	#local APP_DEFAULT=$(find $FOLDER_APPS -name default.app | head -n1)

	if [ "$APP" == "" ]; then
	
		echo $CONFIG_DEFAULT_APP && return 0
		#[ -e "$APP_DEFAULT" ] && echo $APP_DEFAULT && return 0
	
	else

		local FIND_APP_FILE=$(find $FOLDER_APPS -name $APP)
		[ ! -z $FIND_APP_FILE ] && echo $FIND_APP_FILE && return 0
	
		local FIND_APP_SUB_FILE=$(find $FOLDER_APPS -name $APP.app -or -name $APP_tr.app | head -n1)
		[ ! -z $FIND_APP_SUB_FILE ] && echo $FIND_APP_SUB_FILE && return 0

		for APP_TEST in $(find $FOLDER_APPS -name *app); do
			#local APP_NAME=$($STARK_FOLDER_BIN/extract_variable_from_env.sh $APP_TEST "APP_NAME" 2>/dev/null);
			local APP_NAME=$(source $APP_TEST 2>/dev/null; echo $APP_NAME);
			[[ $APP = $APP_NAME ]] && echo $APP_TEST && return 0
			[[ $APP_tr = $APP_NAME ]] && echo $APP_TEST && return 0
		done;
	
	fi;
	
	#echo $APP_DEFAULT && return 0
	return 1
	
} #function find_app


