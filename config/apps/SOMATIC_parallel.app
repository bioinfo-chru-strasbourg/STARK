#!/bin/bash
## STARK application SOMATIC parallel

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP,SOMATIC_PARAMETERS

# APPLICATION INFOS
#####################
APP_NAME="SOMATIC_parallel"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect somatic mutations in solid tumor gene panel sequencing data"
APP_GROUP="SOMATIC"
APP_PROJECT="UNKNOWN"

# ANALYSIS PARAMETERS
#######################

# PIPELINES
PIPELINES="bwamem.gatkUG_SOMATIC.howard bwamem.gatkHC_SOMATIC.howard bwamem.MuTect2_parallel_stringent.howard bwamem.VarScan_parallel_SOMATIC.howard bwamem.outLyzer_parallel_filtered.howard "

