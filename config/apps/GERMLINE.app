#!/bin/bash
## STARK application GERMLINE

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP

# APPLICATION INFOS
#####################
APP_NAME="GERMLINE"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect germline mutations in gene panel sequencing data"
APP_GROUP="GENETIC"
APP_PROJECT="GERMLINE"

# ANALYSIS PARAMETERS
#######################

# PIPELINES
PIPELINES="bwamem.gatkHC_GERMLINE.howard bwamem.gatkUG_GERMLINE.howard"

# INTERVAL_PADDING / add some padding to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=100
