#!/bin/bash
#################################
## STARK environment
#################################


# Timestamp
function timestamp() {
    date +"%Y-%m-%d %H:%M:%S"
}

# Log prefix function
# $1: LEVEL (INFO, WARNING, ERROR, DEBUG)
# $2: message
# $3: TYPE (FULL, NOCOLOR, or anything else for colored level only; default NOCOLOR)
log_prefix () {

	local LEVEL=$1 LOG=$2 TYPE=$3
	local LABEL COLOR="" NC=$'\033[0m' TS PAD

	# Short label and color per level
	case "$LEVEL" in
		INFO)    LABEL="INFO";  COLOR=$'\033[0;32m' ;;
		WARNING) LABEL="WARN";  COLOR=$'\033[0;33m' ;;
		ERROR)   LABEL="ERROR"; COLOR=$'\033[0;31m' ;;
		DEBUG)   LABEL="DEBUG"; COLOR=$'\033[0;35m' ;;
		*)       LABEL="$LEVEL"; NC="" ;;
	esac

	TS=$(timestamp)
	printf -v PAD '%5s' "$LABEL"   # right-aligned on 5 chars

	case "$TYPE" in
		FULL)    printf '#[%s] %s|%s %s%s\n' "$TS" "$PAD" "$COLOR" "$LOG" "$NC" ;;
		NOCOLOR) printf '#[%s] %s| %s\n'     "$TS" "$PAD" "$LOG" ;;
		*)       printf '#[%s] %s%s%s| %s\n' "$TS" "$COLOR" "$PAD" "$NC" "$LOG" ;;
	esac

}

# function in_array
# input: $element $array
in_array () 
{ 
	local param=$1;
	shift;
	for elem in "$@";
	do
		[[ "$param" = "$elem" ]] && return 0;
	done;
	return 1
}


get_cores() {
# Get the number of cores available for the current process, taking into account cgroup limits if applicable.

	compute_cores() {
        local cores=$(((quota + period - 1) / period))

        # Clamp minimum to 1
        [ "$cores" -lt 1 ] && cores=1

        echo "$cores"
    }

    # --- Cgroups v2 ---
    if [ -f /sys/fs/cgroup/cpu.max ]; then

        read quota period < /sys/fs/cgroup/cpu.max

        if [ "$quota" != "max" ]; then
            compute_cores
            return
        fi

        nproc
        return

    fi

    # --- Cgroups v1 ---
    if [ -f /sys/fs/cgroup/cpu/cpu.cfs_quota_us ]; then

        local quota=$(cat /sys/fs/cgroup/cpu/cpu.cfs_quota_us)
        local period=$(cat /sys/fs/cgroup/cpu/cpu.cfs_period_us)

        if [ "${quota:-0}" -gt 0 ] 2>/dev/null; then
			compute_cores
			return
        fi

        nproc
        return
    fi

    # --- No file available ---
    nproc
}

extract_tag () {
# Extract APP variable from STARK TAG List
# $1 : STARK TAG List
# format: "TYPE#TAG[#TAG]{ !}[TYPE#TAG[#TAG]]".
# TYPE can be null. TAG null will not be considered. TAG can be separated by " " (:space:) or "!". Words without "#" will not be considered as TAG. TAG can be cumulative for a TYPE (TYPE#TAG1#TAG2)
# e.g. "TYPE1#TAG1 TAG_TO_FIND#TAG_FOUND TAG_TO_FIND#TAG_FOUND2 TAG_TO_FIND2#TAG_FOUND3#TAG_FOUND4 #TAG2 TAG_NOT_CONSIDERED# Other word not considered as TAG")
# $2 : TAGs to find. default "" ALL TAGs.
# format: "TYPE{ ,}[TYPE]"
# e.g. "TYPE1 TYPE2", "TYPE1,TYPE2"
# $3 : Extract mode type. default "" or "TAG" only TAG.
# TAG: Extract only TAG (without TYPE, alla TAG separated)
# #TAG: Extract only TAG (without TYPE, alla TAG separated)
# TYPE#TAG: Extract TYPE#TAG TYPE#TAG...
# TYPE#TAG#TAG: Extract TYPE#TAG#TAG...

	local TAG_STARK=$(echo $1 | tr "!" " ")
	local TAG_FILTER=$(echo $2 | tr "," " ")
	local MODE=$3
	local TAG_LIST=""
	local TAG_TYPE=""
	local TAG_VALUE=""

	[ "$MODE" == "" ] && MODE="TAG"

	for T in $TAG_STARK; do

		if (( $(echo $T | grep "#" -c) )); then
			TYPE=$(echo "$T" | awk -F"#" '{print $1}' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
			TAG=$(echo "$T" | cut -d"#" -f2- | tr "#" " " | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')

			# MODE
			TAG_TO_EXTRACT=""
			[ "$MODE" == "TAG" ] && TAG_TO_EXTRACT="$TAG"
			[ "$MODE" == "#TAG" ] && TAG_TO_EXTRACT="#$TAG"
			[ "$MODE" == "TYPE#TAG" ] && for T in $TAG; do TAG_TO_EXTRACT="$TAG_TO_EXTRACT$TYPE#$T "; done;
			[ "$MODE" == "TYPE#TAG#TAG" ] && TAG_TO_EXTRACT=$(echo "$TAG_TO_EXTRACT$TYPE#$TAG" | tr " " "#")

			if [ "$TAG" != "" ]; then

				if [ "$TAG_FILTER" == "" ]; then
					TAG_LIST="$TAG_LIST $TAG_TO_EXTRACT"
				else
					for V in $TAG_FILTER; do
						if [ "$V" == "$TYPE" ]; then
							TAG_LIST="$TAG_LIST $TAG_TO_EXTRACT"
						fi;
					done;
				fi;
			fi;
		fi;
	done;
	echo $TAG_LIST | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//'
	return 0;

}


source_app () {
# source an env file from an application definition
# $1: APP definition
# $2: APPS folder
# return APP env file or default APP env file or null

	[ -z $TMP_SYS_FOLDER ] && TMP_SYS_FOLDER="/tmp"

	local APP_LIST=$(echo $1 | tr "," " " | tr "+" " ") FOLDER_APPS=$2 VERBOSE=$3

	local TMP_VERBOSE=$TMP_SYS_FOLDER/$RANDOM

	if [ "$FOLDER_APPS" == "" ]; then FOLDER_APPS="$STARK_FOLDER_APPS"; fi
	if [ "$FOLDER_APPS" == "" ]; then FOLDER_APPS=".."; fi

	if [ "$CONFIG_HEADER" != "" ]; then

		# Source config header
		[[ -f "$CONFIG_HEADER" ]] || { echo "File not found: $CONFIG_HEADER" >&2; exit 1; }	
		source "$CONFIG_HEADER"

	fi

	local ENV=$(find_app "$APP_LIST" "$FOLDER_APPS")

	for ENV_ONE in $ENV; do

		# Add app folder to allow specific rules defined in apps
		local APP_FOLDER=$(dirname "$ENV_ONE")

		# Source app env file
		[[ -f "$ENV_ONE" ]] || { echo "File not found: $ENV_ONE" >&2; exit 1; }
 		source "$ENV_ONE"

	done

	# FOOTER
	if [ "$CONFIG_FOOTER" != "" ]; then

		# Source config footer
		[[ -f "$CONFIG_FOOTER" ]] || { echo "File not found: $CONFIG_FOOTER" >&2; exit 1; }
 		source "$CONFIG_FOOTER"

	fi

} # source_app

name_app () {
# source an env file from an application definition
# $1: APP definition
# $2: APPS folder
# return APP env file or default APP env file or null

	local APP_LIST=$(echo $1 | tr "," " ") FOLDER_APPS=$2

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


description_app () {
# source an env file from an application definition
# $1: APP definition
# $2: APPS folder
# return APP env file or default APP env file or null

	local APP_LIST=$(echo $1 | tr "," " ") FOLDER_APPS=$2

	local ENV=$(find_app "$APP_LIST" "$FOLDER_APPS")
	APP_LIST_DESCRIPTION=""
	for E in $ENV; do
		APP_DESCRIPTION=$(source_app "$E"; echo $APP_DESCRIPTION)
		APP_LIST_DESCRIPTION=$APP_LIST_DESCRIPTION" "$APP_DESCRIPTION". "
	done

	echo $APP_LIST_DESCRIPTION
	return 0;

} # source_app

check_app () {
    APP_TEST="$1"

    # Environment isolation (no source in the parent shell)
    APP_NAME=$(
        cd "$(dirname "$APP_TEST")" &&
        source_app "$APP_TEST" >/dev/null 2>&1
        echo "$APP_NAME"
    )

    [[ "$APP" = "$APP_NAME" ]] && echo "$APP_TEST"
    [[ "$APP_tr" = "$APP_NAME" ]] && echo "$APP_TEST"
    [[ "$APP_BASENAME" = "$APP_NAME" ]] && echo "$APP_TEST"
    [[ "$APP_BASENAME_tr" = "$APP_NAME" ]] && echo "$APP_TEST"
}
export -f check_app

find_app () {
# find -L a env file from an application definition
# $1: APP definition
# $2: APPS folder
# return APP env file or default APP env file or null
# Usage find_app "APP" "FOLDER_APPS"

	local APP_LIST=$(echo $1 | tr "," " "  | tr "+" " ") FOLDER_APPS=$2

	[ -z $TMP_SYS_FOLDER ] && TMP_SYS_FOLDER="/tmp"

	# default APP

	if [ "$APP_LIST" == "" ]; then

		return 0

	elif [ $(echo "$APP_LIST" | wc -w) -gt 1 ]; then


		APP_TO_RETURN=""
		for APP_ONE in $APP_LIST; do
			APP_TO_RETURN=$APP_TO_RETURN" "$(find_app "$APP_ONE" "$FOLDER_APPS")
		done;
		echo $APP_TO_RETURN
		return 0;

	else
		local APP=$APP_LIST

		# APP is a app file
		[ -f $APP ] && echo $APP && return 0

		local APP_DIRNAME=$(dirname $APP  | tr -d ".")
		local APP_BASENAME=$(basename $APP )

		if [ "$FOLDER_APPS" == "" ] && [ "$APP_DIRNAME" == "" ]; then FOLDER_APPS=".."; fi

		local APP_tr=$(echo $APP | tr "-" ".")
		local APP_BASENAME_tr=$(echo $APP_BASENAME | tr "-" ".")

		local FIND_APP_FILE=$(find -L $FOLDER_APPS/$APP_DIRNAME -name $APP_BASENAME -type f)
		[ ! -z $FIND_APP_FILE ] && echo $FIND_APP_FILE && return 0

		local FIND_APP_SUB_FILE=$(find -L $FOLDER_APPS/$APP_DIRNAME -name $APP_BASENAME.app -or -name $APP_BASENAME_tr.app -name $APP_BASENAME.plugapp -or -name $APP_BASENAME_tr.plugapp | sort -t'/' -k2.2r -k2.1 | head -n1 )
		[ ! -z "$FIND_APP_SUB_FILE" ] && echo $FIND_APP_SUB_FILE && return 0

		[[ "${THREADS:-}" =~ ^[1-9][0-9]*$ ]] || THREADS=1

		# Find all .app and .plugapp files in the FOLDER_APPS directory and check if they match the APP using check_app function in parallel.
		find -L "$FOLDER_APPS" \
			-name '*.app' -or -name '*.plugapp' |
		xargs -P$THREADS -I{} bash -c 'check_app "$@"' _ {} |
		sed -n '1p'

	fi;

	return 1

} #function find_app


convertsecs() {
 ((h=${1}/3600))
 ((m=(${1}%3600)/60))
 ((s=${1}%60))
 printf "%02d"h"%02d"m"%02d"s"\n" $h $m $s
}

list_include_item() {
# list_include_item "$LIST" "$ITEM"
	local list=$(echo "$1" | tr "," " ")
	local item="$2"
	if [[ $list =~ (^|[[:space:]])"$item"($|[[:space:]]) ]] ; then
		# yes, list include item
		result=0
	else
		result=1
	fi
	echo $return;
	return $result;
}
