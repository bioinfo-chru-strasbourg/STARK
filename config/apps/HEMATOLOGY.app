#!/bin/bash
## STARK application HEMATOLOGY

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP,SOMATIC_PARAMETERS

# APPLICATION INFOS
#####################
APP_NAME="HEMATOLOGY"
APP_RELEASE="1.1"
APP_DESCRIPTION="Application to detect somatic mutations in hematology gene panel sequencing data"
APP_GROUP=""
APP_PROJECT=""

# ANALYSIS PARAMETERS
#######################

# PIPELINES
PIPELINES="bwamem.gatkUG_HEMATOLOGY.howard bwamem.VarScan_HEMATOLOGY.howard bwamem.itdseek.howard" #
