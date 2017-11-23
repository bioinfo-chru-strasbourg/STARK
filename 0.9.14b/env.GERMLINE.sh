#!/bin/bash
## STARK env for GERMLINE Analysis

# DEFAULT ENV
######################
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env.sh

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
CALLERS="gatkHC_GERMLINE gatkUG_GERMLINE"
ANNOTATORS="howard"

# PIPELINES
PIPELINES="bwamem.gatkHC_GERMLINE.howard bwamem.gatkUG_GERMLINE.howard"

# POST_ALIGNEMNT
POST_ALIGNMENT_STEPS="sorting markduplicates realignment compress"

# HOWARD FILTER
HOWARD_FILTER=GERMLINE

# FOOTER
###########
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env_footer.sh

