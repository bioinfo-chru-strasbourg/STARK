#!/bin/bash
## STARK env for HUSHEMATO RUNs

# DEFAULT ENV
######################
source $CONFIG_DEFAULT_APP

# APPLICATION INFOS
#####################

APP_NAME="ITD"

GROUP="UNKNOWN"
PROJECT="UNKNOWN"


# ANALYSIS PARAMETERS
#######################

# ALIGNERS
ALIGNERS="bwamem"

# CALLERS
CALLERS="itdseek"

# ANNOTATORS
ANNOTATORS="howard"

# PIPELINES
PIPELINES=""

# POST ALIGNEMENT STEPS (default "sorting realignment clipping compress")
POST_ALIGNMENT_STEPS="sorting realignment recalibration clipping compress"


# FOOTER
###########
source $CONFIG_FOOTER



