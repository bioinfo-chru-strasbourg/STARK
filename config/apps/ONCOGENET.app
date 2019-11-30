#!/bin/bash
## STARK application ONCOGENET

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP

# APPLICATION INFOS
#####################
APP_NAME="ONCOGENET"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect germline mutations in oncology gene panel sequencing data"
APP_GROUP="GENETIC"
APP_PROJECT="ONCOGENETIC"

# ANALYSIS PARAMETERS
#######################

# PIPELINES
PIPELINES="bwamem.gatkHC_ONCOGENET.howard bwamem.gatkUG_ONCOGENET.howard"
