#!/bin/bash
#################################
## STARK environment
#################################

export ENV_NAME="STARK"
export ENV_DESCRIPTION="Stellar Tools for variants Analysis and RanKing"
export ENV_RELEASE="0.9.18.1"
export ENV_DATE="02/06/2020"
export ENV_AUTHOR="Antony Le Bechec/Amandine Velt/Sinthuja Pachchek/Vincent Zilliox/Samuel Nicaise/Jean Muller"
export ENV_COPYRIGHT="HUS"
export ENV_LICENCE="GNU GPLA V3"

export STARK_FOLDER_CONFIG="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export CONFIG_HEADER=$STARK_FOLDER_CONFIG/header.app
export CONFIG_FOOTER=$STARK_FOLDER_CONFIG/footer.app
export CONFIG_FUNCTIONS=$STARK_FOLDER_CONFIG/functions.app
export CONFIG_DATABASES=$STARK_FOLDER_CONFIG/databases.app
export CONFIG_TOOLS=$STARK_FOLDER_CONFIG/tools.app
export CONFIG_CONFIG=$STARK_FOLDER_CONFIG/config.app

source $CONFIG_HEADER 1>/dev/null 2>/dev/null
source $CONFIG_FOOTER 1>/dev/null 2>/dev/null
source $CONFIG_FUNCTIONS 1>/dev/null 2>/dev/null
source $CONFIG_TOOLS 1>/dev/null 2>/dev/null
source $CONFIG_DATABASES 1>/dev/null 2>/dev/null

[ -e $STARK_FOLDER_CONFIG/default.app ] && export CONFIG_DEFAULT_APP=$STARK_FOLDER_CONFIG/default.app
[ -e $STARK_FOLDER_APPS/default.app ] && export CONFIG_DEFAULT_APP=$STARK_FOLDER_APPS/default.app
