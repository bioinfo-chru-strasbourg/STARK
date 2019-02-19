#!/bin/bash
## STARK env for CPSGEN.MASTR RUNs

# DEFAULT ENV
######################
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env.ONCO.sh

# APPLICATION INFOS
#####################

APP_NAME="HUSDIAGGEN"
APP_RELEASE="1.0"

GROUP="HUSDIAGGEN"
PROJECT="UNKNOWN"


# ANALYSIS PARAMETERS
#######################

# ALIGNERS
ALIGNERS="bwamem"

# CALLERS
CALLERS="gatkHC_HUSDIAGGEN gatkUG_HUSDIAGGEN canoes"

# ANNOTATORS
ANNOTATORS="howard"

# PIPELINES
PIPELINES=""


# POST_ALIGNMENT (Capture)
POST_ALIGNMENT_STEPS="sorting markduplicates realignment compress"


# HOWARD FILTER
################
#HOWARD_FILTER="DIAGGEN,GERMLINE"
HOWARD_PRIORITIZATION="DIAGGEN,GERMLINE"
HOWARD_PRIORITIZATION_MINIMAL="DIAGGEN,GERMLINE"
HOWARD_PRIORITIZATION_REPORT="DIAGGEN,GERMLINE"


# FOOTER
##########
#source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env_footer.sh
source $STARK_FOLDER/env_footer.sh


