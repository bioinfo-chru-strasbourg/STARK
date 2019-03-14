#!/bin/bash
## STARK application ONCO

# DEFAULT ENV
######################
source_app "ONCO"

# APPLICATION INFOS
#####################
APP_NAME="HUSDIAGGEN"
APP_RELEASE="1.0"
APP_DESCRIPTION="Applications des Hopitaux universitaire de Strasbourg - Département d'onco-génétique - Laboratoire de diagnostic onco-génétique"
APP_GROUP="HUSDIAGGEN"
APP_PROJECT=""


# ANALYSIS PARAMETERS
#######################


# Rules
# Add specific rules to load.
# These files will be added to the list of rules files from APPS folder
# example: "MYGROUP/*.rules.mk"
# for inheritance, use RULES_APP="$RULES_APP MYGROUP/*.rules.mk"
#RULES_APP="$RULES_APP HUS/ONCO/HUSDIAGGEN.rules.mk/*.rules.mk"
RULES_APP="$RULES_APP $APP_FOLDER/HUSDIAGGEN.rules.mk/*.rules.mk"


# PIPELINES
PIPELINES="bwamem.gatkHC_HUSDIAGGEN.howard bwamem.gatkUG_HUSDIAGGEN.howard bwamem.canoes.howard"


# POST_ALIGNMENT (Capture)
POST_ALIGNMENT_STEPS="sorting markduplicates realignment compress"


# HOWARD FILTER
################
HOWARD_PRIORITIZATION="DIAGGEN,GERMLINE"
HOWARD_PRIORITIZATION_MINIMAL="DIAGGEN,GERMLINE"
HOWARD_PRIORITIZATION_REPORT="DIAGGEN,GERMLINE"
