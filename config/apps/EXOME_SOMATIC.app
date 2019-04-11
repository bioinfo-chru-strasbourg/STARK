#!/bin/bash
## STARK application EXOME_SOMATIC

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP,CAPTURE,SOMATIC_PARAMETERS

# APPLICATION INFOS
#####################
APP_NAME="EXOME_SOMATIC"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect somatic mutations in exome sequencing data"
APP_GROUP=""
APP_PROJECT=""

# ANALYSIS PARAMETERS
#######################

# PIPELINES
PIPELINES="bwamem.gatkHC_EXOME_SOMATIC.howard bwamem.gatkUG_EXOME_SOMATIC.howard bwamem.VarScan_EXOME_SOMATIC.howard"

# INTERVAL_PADDING / add some padding to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=100

# HOWARD ANNOTATION/PRIOTITIZATION/TRANSLATION CONFIGURATION

# ANNOTATION
# Default annotation with HOWARD for intermediate VCF (for each caller) used by default with annotation rule "howard"
#ANNOTATION_TYPE="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs" "core,symbol,location,outcome,hgvs,snpeff,snpeff_hgvs,snpeff_split"
HOWARD_ANNOTATION=""
# Default annotation with HOWARD for minimal VCF annotation (rule howard_minimal)
HOWARD_ANNOTATION_MINIMAL=""
# Default annotation with HOWARD for report
HOWARD_ANNOTATION_REPORT="core,frequency,score,annotation,prediction,snpeff_split"
