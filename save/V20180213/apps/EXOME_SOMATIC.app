#!/bin/bash
## STARK env for GERMLINE Analysis

# DEFAULT ENV
######################
source $STARK_FOLDER_APPS/EXOME.app

# APPLICATION INFOS
#####################
APP_NAME="EXOME_SOMATIC"
GROUP=EXOME_SOMATIC
PROJECT=EXOME_SOMATIC

# ANALYSIS PARAMETERS
#######################

# INTERVAL_PADDING / add some padding to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=100

# ALIGNERS CALLERS ANNOTATORS
ALIGNERS="bwamem"
CALLERS="gatkHC_EXOME_SOMATIC gatkUG_EXOME_SOMATIC VarScan_EXOME_SOMATIC"
ANNOTATORS="howard"

# PIPELINES
PIPELINES="bwamem.gatkHC_EXOME_SOMATIC.howard bwamem.gatkUG_EXOME_SOMATIC.howard bwamem.VarScan_EXOME_SOMATIC.howard"

# POST_ALIGNEMNT
POST_ALIGNMENT_STEPS="sorting markduplicates realignment compress"


# FOOTER
###########
source $CONFIG_FOOTER

