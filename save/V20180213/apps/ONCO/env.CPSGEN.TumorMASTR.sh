#!/bin/bash
## STARK env for HUSDIAGGEN.MASTR RUNs

# DEFAULT ENV
######################
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env.HUSDIAGGEN.sh

# APPLICATION INFOS
#####################

# Name of the APP. usefull to include specific rules if exists (in $STARK/$APP_NAME.rules.mk/*rules.mk)
# AUTO detect: $(basename ${BASH_SOURCE[0]} | sed 's/^env//' | sed 's/\.sh$//gi' | cut -d. -f2-)
APP_NAME="HUSDIAGGEN.TumorMASTR"
APP_RELEASE="1.0"

GROUP="HUSDIAGGEN"
PROJECT="TumorMASTR"


# ANALYSIS PARAMETERS
#######################

# CALLERS
CALLERS="gatkHC_CPSGEN_MASTR gatkUG_CPSGEN_MASTR"

POST_ALIGNMENT_STEPS="sorting realignment recalibration clipping compress"

# HOWARD FILTER
################
HOWARD_PRIORITIZATION="CPSGEN.MASTR.HARD,CPSGEN.MASTR,CPSGEN.MASTR.SOFT,GERMLINE"
HOWARD_PRIORITIZATION_MINIMAL="CPSGEN.MASTR.HARD,CPSGEN.MASTR,CPSGEN.MASTR.SOFT,GERMLINE"
HOWARD_PRIORITIZATION_REPORT="CPSGEN.MASTR.HARD,CPSGEN.MASTR,CPSGEN.MASTR.SOFT,GERMLINE"


# PIPELINES
PIPELINES="bwamem.gatkHC_CPSGEN_MASTR.howard bwamem.gatkUG_CPSGEN_MASTR.howard canoe"


# INTERVAL_PADDING (default 0)
# Add some “padding” to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=25

# FOOTER
###########
#source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env_footer.sh
source $STARK_FOLDER/env_footer.sh

