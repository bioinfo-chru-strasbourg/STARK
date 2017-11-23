#!/bin/bash
## STARK env for HEMATOLOGY

# DEFAULT ENV
######################
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env.sh

# APPLICATION INFOS
#####################
APP_NAME=HEMATOLOGY
GROUP=HEMATOLOGY
PROJECT=HEMATOLOGY

# ANALYSIS PARAMETERS
#######################

# ALIGNERS CALLERS ANNOTATORS
ALIGNERS="bwamem"
CALLERS="gatkUG_HEMATOLOGY VarScan_HEMATOLOGY"
ANNOTATORS="howard"
# PIPELINES
PIPELINES="bwamem.gatkUG_HEMATOLOGY.howard bwamem.VarScan_HEMATOLOGY.howard"

# NB_BASES_AROUND
NB_BASES_AROUND=0

# HOWARD FILTER
HOWARD_FILTER=HEMATOLOGY

# FOOTER
###########

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env_footer.sh


