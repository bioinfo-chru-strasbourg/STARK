#!/bin/bash
## STARK application GENOME

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP

# APPLICATION INFOS
#####################
APP_NAME="GENOME"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect germline mutations in genome sequencing data"
APP_GROUP="GENETIC"
APP_PROJECT="GENOME"

# ANALYSIS PARAMETERS
#######################

# PIPELINES
PIPELINES="bwamem.gatkHC_GENOME.howard bwamem.gatkUG_GENOME.howard"

# INTERVAL_PADDING / add some padding to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=100

