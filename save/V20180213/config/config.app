#!/bin/bash
#################################
## STARK environment
#################################

STARK_FOLDER_CONFIG="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export CONFIG_HEADER=$STARK_FOLDER_CONFIG/header.app
export CONFIG_FOOTER=$STARK_FOLDER_CONFIG/footer.app
export CONFIG_FUNCTIONS=$STARK_FOLDER_CONFIG/functions.app
export CONFIG_DATABASES=$STARK_FOLDER_CONFIG/databases.app
export CONFIG_TOOLS=$STARK_FOLDER_CONFIG/tools.app

source $CONFIG_HEADER 1>/dev/null 2>/dev/null
source $CONFIG_FOOTER 1>/dev/null 2>/dev/null
source $CONFIG_FUNCTIONS 1>/dev/null 2>/dev/null

[ -e $STARK_FOLDER_CONFIG/default.app ] && export CONFIG_DEFAULT_APP=$STARK_FOLDER_CONFIG/default.app
[ -e $STARK_FOLDER_APPS/default.app ] && export CONFIG_DEFAULT_APP=$STARK_FOLDER_APPS/default.app


#echo $RULES;
#exit 0

