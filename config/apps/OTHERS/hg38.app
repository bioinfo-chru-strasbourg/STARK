#!/bin/bash
## STARK env for HUSHEMATO RUNs

# DEFAULT ENV
######################
#source $STARK_FOLDER_APPS/default.app
#source $CONFIG_DEFAULT_APP

# APPLICATION INFOS
#####################

# Name of the APP
# Used to autodetect APP in RUNS (SampleSheet/investigator Name field)
# Usefull to include specific rules if exists (in $STARK/$APP_NAME.rules.mk/*rules.mk)
# AUTO detect with ENV file name:
APP_NAME="hg38"

# Release of the APP
APP_RELEASE="1.0"

# GROUP and PROJECT Associated with the APP
# Use to structure data in the repository folder
# AUTO detect with ENV file name:
APP_GROUP=""
APP_PROJECT=""


# ASSEMBLY
#############

ASSEMBLY=hg38


# RULES
##########
#echo "truc"
#SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )"
#echo $SCRIPT_DIR
#APP_RULES="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )"/hg38.rules.mk/*rules.mk

# FOOTER
###########
#source $CONFIG_FOOTER



