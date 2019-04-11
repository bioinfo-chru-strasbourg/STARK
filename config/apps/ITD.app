#!/bin/bash
## STARK application ITD

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP,AMPLICON

# APPLICATION INFOS
#####################
APP_NAME="ITD"
APP_RELEASE="1.0"
APP_DESCRIPTION="ITD-FMT3 mutation detection (ITDSeek tool), for amplicon sequencing"
APP_GROUP=""
APP_PROJECT=""

# ANALYSIS PARAMETERS
#######################

# PIPELINES
PIPELINES="bwamem.itdseek.howard"

# POST ALIGNEMENT STEPS (default "sorting realignment clipping compress")
#POST_ALIGNMENT_STEPS="sorting realignment recalibration clipping compress"
