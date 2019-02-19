#!/bin/bash
## STARK env for GERMLINE Analysis

# DEFAULT ENV
######################
source $CONFIG_DEFAULT_APP

# APPLICATION INFOS
#####################
APP_NAME="GERMLINE"
GROUP=GERMLINE
PROJECT=GERMLINE

# ANALYSIS PARAMETERS
#######################

# INTERVAL_PADDING / add some padding to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=100

# ALIGNERS CALLERS ANNOTATORS
ALIGNERS="bwamem"
CALLERS="gatkHC_GERMLINE gatkUG_GERMLINE canoes"
ANNOTATORS="howard"

# PIPELINES
PIPELINES="bwamem.gatkHC_GERMLINE.howard bwamem.gatkUG_GERMLINE.howard"

# POST_ALIGNEMNT
POST_ALIGNMENT_STEPS="sorting markduplicates realignment compress"

# FOOTER
###########
source $CONFIG_FOOTER

