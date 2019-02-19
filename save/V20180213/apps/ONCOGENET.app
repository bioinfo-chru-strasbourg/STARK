#!/bin/bash

# DEFAULT ENV
######################
source $CONFIG_DEFAULT_APP

# APPLICATION INFOS
#####################

# Name of the APP. usefull to include specific rules if exists (in $STARK/$APP_NAME.rules.mk/*rules.mk)
# AUTO detect: $(basename ${BASH_SOURCE[0]} | sed 's/^env//' | sed 's/\.sh$//gi' | cut -d. -f2-)
APP_NAME="ONCOGENET"

GROUP="UNKNOWN"
PROJECT="UNKNOWN"

# ANALYSIS PARAMETERS
#######################

# ALIGNERS
ALIGNERS="bwamem"

# CALLERS
CALLERS="gatkHC_ONCOGENET gatkUG_ONCOGENET canoes"

# ANNOTATORS
ANNOTATORS="howard"

# PIPELINES
PIPELINES=""

# POST_ALIGNMENT (Capture)
POST_ALIGNMENT_STEPS="sorting markduplicates realignment compress"


# FOOTER
###########
source $CONFIG_FOOTER






