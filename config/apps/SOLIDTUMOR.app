#!/bin/bash
## STARK application SOLIDTUMOR

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP,SOMATIC_PARAMETERS

# APPLICATION INFOS
#####################
APP_NAME="SOLIDTUMOR"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect somatic mutations in solid tumor gene panel sequencing data"
APP_GROUP=""
APP_PROJECT=""

# ANALYSIS PARAMETERS
#######################

FOLDER_DEMULTIPLEXING=/STARK/output/demuSOLIDTUM

# PIPELINES
PIPELINES="bwamem.gatkUG_SOLIDTUMOR.howard bwamem.gatkHC_SOLIDTUMOR.howard bwamem.VarScan_SOLIDTUMOR.howard"
