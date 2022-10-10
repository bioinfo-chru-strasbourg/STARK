#!/bin/bash
## STARK application LONG_INDELS

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP

# APPLICATION INFOS
#####################
APP_NAME="LONG_INDELS"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect long indels in gene panel sequencing data"
APP_GROUP="GENETIC"
APP_PROJECT="LONG_INDELS"

# ANALYSIS PARAMETERS
#######################

# PIPELINES
PIPELINES="bwamem.gatkHC_LONG_INDELS.howard bwamem.gatkUG_LONG_INDELS.howard"

# INTERVAL_PADDING / add some padding to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=100
