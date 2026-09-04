#!/bin/bash
## STARK application EXOME

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP

# APPLICATION INFOS
#####################
APP_NAME="EXOME"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect germline mutations in exome sequencing data"
APP_GROUP="GENETIC"
APP_PROJECT="EXOME"

# ANALYSIS PARAMETERS
#######################

# PIPELINES
PIPELINES="bwamem.gatkHC_EXOME.howard bwamem.gatkUG_EXOME.howard"

# INTERVAL_PADDING / add some padding to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=100
