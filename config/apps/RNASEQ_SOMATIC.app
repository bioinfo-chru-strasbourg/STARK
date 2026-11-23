#!/bin/bash
## STARK application RNASEQ

# DEFAULT ENV
######################
# Depends on the RNASEQ application
source_app RNASEQ


# APPLICATION INFOS
#####################
APP_NAME="RNASEQ"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect mutations in RNA-Seq data"
APP_GROUP="UNKNOWN"
APP_PROJECT="UNKNOWN"


# ANALYSIS PARAMETERS
#######################

# PIPELINES
# Use STAR to align RNA-Seq data and then different callers to detect mutations (Fusions and SNVs/Indels)
# Add callers for specific somatic mutation detection
PIPELINES=$PIPELINES" star.gatkUG.howard star.outLyzer.howard star.MuTect2.howard star.VarScan.howard "

