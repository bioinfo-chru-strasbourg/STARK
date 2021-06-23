#!/bin/bash
## STARK application SOLIDTUMOR

# DEFAULT ENV
######################
#source_app $CONFIG_DEFAULT_APP,SOMATIC_PARAMETERS
source_app SOMATIC

# APPLICATION INFOS
#####################
APP_NAME="SOLIDTUMOR"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect somatic mutations in solid tumor gene panel sequencing data"
APP_GROUP="SOMATIC"
APP_PROJECT="SOLIDTUMOR"

# ANALYSIS PARAMETERS
#######################

# PIPELINES
#PIPELINES="bwamem.gatkUG_SOMATIC.howard bwamem.gatkHC_SOMATIC.howard bwamem.MuTect2.howard bwamem.VarScan.howard bwamem.outLyzer.howard "
#PIPELINES="bwamem.gatkUG_SOMATIC.howard bwamem.gatkHC_SOMATIC.howard bwamem.MuTect2_stringent.howard bwamem.MuTect2.howard bwamem.VarScan.howard bwamem.outLyzer.howard "
#PIPELINES="bwamem.MuTect2.howard bwamem.MuTect2_filtered.howard bwamem.MuTect2_stringent.howard bwamem.outLyzer_filtered.howard"
#PIPELINES="bwamem.gatkUG_SOMATIC.howard bwamem.gatkHC_SOMATIC.howard bwamem.MuTect2_stringent.howard bwamem.VarScan_SOMATIC.howard bwamem.outLyzer_filtered.howard "

