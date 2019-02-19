#!/bin/bash
## STARK env for GENOME Analysis

# DEFAULT ENV
######################
source $CONFIG_DEFAULT_APP

# APPLICATION INFOS
#####################
APP_NAME=GENOME
GROUP=GENOME
PROJECT=GENOME

# ANALYSIS PARAMETERS
#######################

# INTERVAL_PADDING / add some padding to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=100

# ALIGNERS CALLERS ANNOTATORS
ALIGNERS="bwamem"
CALLERS="gatkHC_GENOME gatkUG_GENOME"
ANNOTATORS="howard"
# PIPELINES
PIPELINES="bwamem.gatkHC_GENOME.howard bwamem.gatkUG_GENOME.howard"

# POST_ALIGNMENT (Capture)
POST_ALIGNMENT_STEPS="sorting markduplicates realignment compress"


# HOWARD ANNOTATION/PRIOTITIZATION/TRANSLATION CONFIGURATION

# ANNOTATION
# Default annotation with HOWARD for intermediate VCF (for each caller) used by default with annotation rule "howard"
#ANNOTATION_TYPE="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs" "core,symbol,location,outcome,hgvs,snpeff,snpeff_hgvs,snpeff_split"
HOWARD_ANNOTATION=""
# Default annotation with HOWARD for minimal VCF annotation (rule howard_minimal)
HOWARD_ANNOTATION_MINIMAL=""
# Default annotation with HOWARD for report
HOWARD_ANNOTATION_REPORT="core,frequency,score,annotation,prediction,snpeff_split"


# CALCULATION
# Default calculation with HOWARD for all VCF/pipelines
HOWARD_CALCULATION="VARTYPE,NOMEN"
# Default minimal calculation with HOWARD for final VCF report
HOWARD_CALCULATION_MINIMAL="VARTYPE,NOMEN"
# Default calculation with HOWARD for final VCF report
HOWARD_CALCULATION_REPORT="FindByPipelines,GenotypeConcordance,VAF_STATS,CALLING_QUALITY,CALLING_QUALITY_EXPLODE,VARTYPE,NOMEN"


# PRIORITIZATION
# Default filter to prioritize/rank variant.
# This option create ranking scores in VCF and comment in TXT (after translation).
# Scores can be used to sort variant in the TXT
# HOWARD_FILTER_DEFAULT="dfault" # in env_header.sh 
# Default calculation with HOWARD 
HOWARD_PRIORITIZATION=$HOWARD_PRIORITIZATION_DEFAULT # "default"
# Minimal calculation with HOWARD 
HOWARD_PRIORITIZATION_MINIMAL=$HOWARD_PRIORITIZATION_DEFAULT # "default"
# Default calculation with HOWARD for Report (full/final VCF)
HOWARD_PRIORITIZATION=$HOWARD_PRIORITIZATION_DEFAULT # "default"



# FOOTER
###########
source $CONFIG_FOOTER
