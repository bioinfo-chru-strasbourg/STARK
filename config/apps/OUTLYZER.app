#!/bin/bash
## STARK application OUTLYZER

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP,SOMATIC_PARAMETERS

# APPLICATION INFOS
#####################
APP_NAME="OUTLYZER"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect somatic mutations in solid tumor gene panel sequencing data with OutLyzer"
APP_GROUP=""
APP_PROJECT=""

# ANALYSIS PARAMETERS
#######################


# PIPELINES
PIPELINES="bwamem.outLyzer.howard"
