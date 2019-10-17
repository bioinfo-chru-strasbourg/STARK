#!/bin/bash
## STARK application MUTECT

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP,SOMATIC_PARAMETERS

# APPLICATION INFOS
#####################
APP_NAME="MUTECT"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect somatic mutations in solid tumor gene panel sequencing data with MuTect"
APP_GROUP=""
APP_PROJECT=""

# ANALYSIS PARAMETERS
#######################

export MUTECT_INTERVAL_PADDING=50

# PIPELINES
PIPELINES="bwamem.MuTect2.howard bwamem.GATK3_MuTect2.howard bwamem.MuTect.howard"
#PIPELINES="bwamem.mutect2.howard"
