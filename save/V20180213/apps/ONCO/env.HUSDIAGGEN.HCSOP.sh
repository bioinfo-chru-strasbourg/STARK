#!/bin/bash
## STARK env for HUSDIAGGEN.HCSOP RUNs

# DEFAULT ENV
######################
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env.HUSDIAGGEN.sh

# APPLICATION INFOS
#####################

# Name of the APP. usefull to include specific rules if exists (in $STARK/$APP_NAME.rules.mk/*rules.mk)
# AUTO detect: $(basename ${BASH_SOURCE[0]} | sed 's/^env//' | sed 's/\.sh$//gi' | cut -d. -f2-)
APP_NAME="HUSDIAGGEN.HCSOP"
APP_RELEASE="1.0"

GROUP="HUSDIAGGEN"
PROJECT="HCSOP"


# ANALYSIS PARAMETERS
#######################

# CALLERS
CALLERS="gatkHC_CPSGEN_MASTR gatkUG_CPSGEN_MASTR"

# HOWARD FILTER
################
#HOWARD_FILTER="CPSGEN.MASTR.HARD,CPSGEN.MASTR,CPSGEN.MASTR.SOFT,GERMLINE"
HOWARD_PRIORITIZATION="HUSDIAGGEN.MASTR.HARD,HUSDIAGGEN.MASTR,HUSDIAGGEN.MASTR.SOFT,GERMLINE"
HOWARD_PRIORITIZATION_MINIMAL="HUSDIAGGEN.MASTR.HARD,HUSDIAGGEN.MASTR,HUSDIAGGEN.MASTR.SOFT,GERMLINE"
HOWARD_PRIORITIZATION_REPORT="HUSDIAGGEN.MASTR.HARD,HUSDIAGGEN.MASTR,HUSDIAGGEN.MASTR.SOFT,GERMLINE"

# COVERAGE CRITERIA (default "1,30")
COVERAGE_CRITERIA="1,100,300"

# PIPELINES
PIPELINES="bwamem.gatkHC_CPSGEN_MASTR.howard bwamem.gatkUG_CPSGEN_MASTR.howard"


# POST_ALIGNMENT (Capture)
POST_ALIGNMENT_STEPS="sorting markduplicates realignment compress"

# INTERVAL_PADDING (default 0)
# Add some “padding” to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=25

# FOOTER
###########
#source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env_footer.sh
source $STARK_FOLDER/env_footer.sh






