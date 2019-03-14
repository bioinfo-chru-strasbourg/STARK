#!/bin/bash
## STARK application ONCO

# DEFAULT ENV
######################
source_app "ONCO"

# APPLICATION INFOS
#####################
APP_NAME="HUSHEMATO"
APP_RELEASE="1.0"
APP_DESCRIPTION="Applications des Hopitaux universitaire de Strasbourg - Département d'onco-génétique - Laboratoire d'hématologie"
APP_GROUP="HUSHEMATO"
APP_PROJECT=""


# ANALYSIS PARAMETERS
#######################


# Rules
# Add specific rules to load.
# These files will be added to the list of rules files from APPS folder
# example: "MYGROUP/*.rules.mk"
# for inheritance, use RULES_APP="$RULES_APP MYGROUP/*.rules.mk"
#RULES_APP="$RULES_APP HUS/ONCO/HUSHEMATO.rules.mk/*.rules.mk"
RULES_APP="$RULES_APP $APP_FOLDER/HUSHEMATO.rules.mk/*.rules.mk"


# PIPELINES
PIPELINES="bwamem.gatkUG_HUSHEMATO.howard bwamem.VarScan_HUSHEMATO.howard bwamem.samtools_HUSHEMATO.howard bwamem.itdseek.howard bwamem.canoes.howard"


# HOWARD FILTER
################
HOWARD_PRIORITIZATION="HEMATO,HEMATOLOGY"
HOWARD_PRIORITIZATION_MINIMAL="HEMATO,HEMATOLOGY"
HOWARD_PRIORITIZATION_REPORT="HEMATO,HEMATOLOGY"

# COVERAGE CRITERIA (default "1,30")
COVERAGE_CRITERIA="1,300,1000"
