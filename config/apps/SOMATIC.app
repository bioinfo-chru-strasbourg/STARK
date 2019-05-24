#!/bin/bash
## STARK application SOLIDTUMOR

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP,SOMATIC_PARAMETERS

# APPLICATION INFOS
#####################
APP_NAME="SOMATIC"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect somatic mutations in solid tumor gene panel sequencing data"
APP_GROUP=""
APP_PROJECT=""

# ANALYSIS PARAMETERS
#######################

# PIPELINES
PIPELINES="bwamem.gatkUG_SOMATIC.howard bwamem.gatkHC_SOMATIC.howard bwamem.VarScan_SOMATIC.howard"