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
APP_GROUP="SOMATIC"
APP_PROJECT="UNKNOWN"

# ANALYSIS PARAMETERS
#######################

# PIPELINES
PIPELINES="bwamem.gatkUG_SOMATIC.howard bwamem.gatkHC_SOMATIC.howard bwamem.MuTect2_stringent.howard bwamem.VarScan_SOMATIC.howard bwamem.outLyzer_filtered.howard "
#PIPELINES="bwamem.gatkUG_SOMATIC bwamem.gatkHC_SOMATIC bwamem.MuTect2_stringent bwamem.VarScan_SOMATIC bwamem.outLyzer_filtered"
#POST_ALIGNMENT_STEPS="sorting markduplicates realignment recalibration compress"