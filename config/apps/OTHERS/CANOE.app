#!/bin/bash
## STARK application CANOES

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP

# APPLICATION INFOS
#####################
APP_NAME="CANOES"
APP_RELEASE="1.0"
APP_DESCRIPTION="CNV detection (CANOES tool), for capture sequencing"
APP_GROUP=""
APP_PROJECT=""


# ANALYSIS PARAMETERS
#######################

# PIPELINES
PIPELINES="bwamem.canoes.howard"

# POST ALIGNEMENT STEPS (default "sorting realignment clipping compress")
POST_ALIGNMENT_STEPS="sorting markduplicates realignment compress"
