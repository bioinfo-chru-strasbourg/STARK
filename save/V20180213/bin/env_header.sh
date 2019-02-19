#!/bin/bash
#################################
## STARK environment
#################################

export ENV_NAME="STARK"
export ENV_DESCRIPTION="Stellar Tools for variants Analysis and RanKing"
export ENV_RELEASE="0.9.18d"
export ENV_DATE="14/01/2019"
export ENV_AUTHOR="Antony Le Bechec/Amandine Velt/Sinthuja PACHCHEK/Vincent ZILLIOX"
export ENV_COPYRIGHT="IRC"
export ENV_LICENCE="GNU-GPL-CeCILL"

# SCRIPT DIR
##############

# Folder of the ENV
export STARK_FOLDER_CONFIG="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )"


# VIARIABLE INITIALISATION
############################

# RULES for APP
RULES_APP=""

# FILE to provide in folder repository
RESULTS_SUBFOLDER_ROOT_FILE_PATTERNS=""

# PIPELINES
PIPELINES=""
ALIGNERS=""
CALLERS=""
ANNOTATORS=""

# PANEL/MANIFEST/BED
MANIFEST=""
BED=""


# HOWARD FILTER DEFAULT
HOWARD_PRIORITIZATION_DEFAULT="default"
ANNOVAR_DATABASES=""
SNPEFF_DATABASES=""

