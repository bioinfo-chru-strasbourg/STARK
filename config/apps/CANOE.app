#!/bin/bash
## STARK env for HUSHEMATO RUNs

# DEFAULT ENV
######################
source $CONFIG_DEFAULT_APP

# FOLDER
#####################

FOLDER_REPOSITORY=/home1/DIAG/DATA/NGS/REP

# APPLICATION INFOS
#####################

APP_NAME="CANOES"

GROUP="UNKNOWN"
PROJECT="UNKNOWN"


# ANALYSIS PARAMETERS
#######################

# ALIGNERS
ALIGNERS="bwamem"

# CALLERS
CALLERS="canoes"

# ANNOTATORS
ANNOTATORS="howard"

# PIPELINES
PIPELINES=""

# POST ALIGNEMENT STEPS (default "sorting realignment clipping compress")
POST_ALIGNMENT_STEPS="sorting markduplicates realignment compress"


# FOOTER
###########
source $CONFIG_FOOTER



