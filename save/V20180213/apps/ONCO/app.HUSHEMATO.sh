#!/bin/bash
## STARK env for HUSHEMATO RUNs

# DEFAULT ENV
######################
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env.ONCO.sh

# APPLICATION INFOS
#####################

APP_NAME="HUSHEMATO"
APP_RELEASE="1.0"

GROUP="HUSHEMATO"
PROJECT="UNKNOWN"


# ANALYSIS PARAMETERS
#######################

# ALIGNERS
ALIGNERS="bwamem"

# CALLERS
CALLERS="gatkUG_HUSHEMATO VarScan_HUSHEMATO samtools_HUSHEMATO itdseek canoes"

# ANNOTATORS
ANNOTATORS="howard"

# PIPELINES
PIPELINES=""


# HOWARD FILTER
################
#HOWARD_FILTER="HEMATO,HEMATOLOGY"
HOWARD_PRIORITIZATION="HEMATO,HEMATOLOGY"
HOWARD_PRIORITIZATION_MINIMAL="HEMATO,HEMATOLOGY"
HOWARD_PRIORITIZATION_REPORT="HEMATO,HEMATOLOGY"


# FOOTER
###########
#source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env_footer.sh
source $STARK_FOLDER/env_footer.sh



