#!/bin/bash
## STARK application EXOME_SOMATIC

# DEFAULT ENV
######################
source_app SOMATIC,CAPTURE

# APPLICATION INFOS
#####################
APP_NAME="EXOME_SOMATIC"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect somatic mutations in exome sequencing data"
APP_GROUP="SOMATIC"
APP_PROJECT="EXOME"

# ANALYSIS PARAMETERS
#######################

# PIPELINES
PIPELINES="bwamem.gatkHC_EXOME_SOMATIC.howard bwamem.gatkUG_EXOME_SOMATIC.howard bwamem.outLyzer_filtered.howard bwamem.VarScan_EXOME_SOMATIC.howard bwamem.MuTect2_stringent.howard"

# INTERVAL_PADDING / add some padding to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=100
