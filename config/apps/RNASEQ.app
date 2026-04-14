#!/bin/bash
## STARK application RNASEQ

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP


# APPLICATION INFOS
#####################
APP_NAME="RNASEQ"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect mutations in RNA-Seq data"
APP_GROUP="UNKNOWN"
APP_PROJECT="UNKNOWN"


# ANALYSIS PARAMETERS
#######################

# Option to keep raw BAM file from STAR (with splitNcigar for SNV calling). This option is useful to keep BAM without post-alignment steps (e.g. realignment, recalibration...).
# Default: False (0)
#STAR_KEEP_RAW_BAM=0

# POST ALIGNMENT STEPS 
# No need to realign because of STAR alignemnt
# POST_ALIGNMENT_STEPS="sorting markduplicates realignment recalibration compress" # DEFAULT
POST_ALIGNMENT_STEPS="sorting splitncigar markduplicates recalibration compress"

# PIPELINES
# Use STAR to align RNA-Seq data and then different callers to detect mutations (Fusions and SNVs/Indels)
PIPELINES="star.Arriba.howard star.STARFusion.howard star.gatkHC.howard"

# POST_CALLING_MERGING_STEPS
# No variant filtration, because of dealing only with SNV and InDels (not Fusions)
POST_CALLING_MERGING_STEPS="sorting normalization variantrecalibration"
