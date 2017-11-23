#!/bin/bash
## STARK env for GENOME Analysis

# DEFAULT ENV
######################
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env.sh

# APPLICATION INFOS
#####################
APP_NAME=GENOME
GROUP=GENOME
PROJECT=GENOME

# ANALYSIS PARAMETERS
#######################

# INTERVAL_PADDING / add some padding to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=100

# ALIGNERS CALLERS ANNOTATORS
ALIGNERS="bwamem"
CALLERS="gatkHC_GENOME gatkUG_GENOME"
ANNOTATORS="howard"
# PIPELINES
PIPELINES="bwamem.gatkHC_GENOME.howard bwamem.gatkUG_GENOME.howard"

# POST_ALIGNMENT
POST_ALIGNMENT_STEPS="sorting markduplicates realignment compress"

# HOWARD FILTER
HOWARD_FILTER=GENOME

# FOOTER
###########
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env_footer.sh
