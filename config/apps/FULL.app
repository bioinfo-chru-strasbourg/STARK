#!/bin/bash
## STARK application HEMATOLOGY

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP

# APPLICATION INFOS
#####################
APP_NAME="FULL"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application with all available pipelines"
APP_GROUP="UNLNOWN"
APP_PROJECT="UNKNOWN"

# ANALYSIS PARAMETERS
#######################

# PIPELINES
#PIPELINES="bowtie.gatkHC.howard bwamem.gatkHC.howard" #
PIPELINES="bowtie.gatk4HC.howard" #
#PIPELINES="bwamem.gatkHC.howard" #
#PIPELINES="bwamem.gatkHC.howard bwamem.gatkUG.howard bwamem.outLyzer.howard bwamem.MuTect2.howard bwamem.itdseek.howard bwamem.VarScan.howard " 
#PIPELINES="bowtie2.gatkHC.howard bowtie2.gatk4HC.howard bowtie2.gatkUG.howard bowtie2.outLyzer.howard bowtie2.MuTect2.howard bowtie2.itdseek.howard bowtie2.VarScan.howard bwamem.gatkHC.howard bwamem.gatk4HC.howard bwamem.gatkUG.howard bwamem.outLyzer.howard bwamem.MuTect2.howard bwamem.itdseek.howard bwamem.VarScan.howard " 
