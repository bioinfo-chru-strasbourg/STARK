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

# COVERAGE
##############

# COVERAGE CRITERIA (default "1,30")
# For gene coverage metrics
# the criteria to calculate the percent of bases over $COVERAGE_CRITERIA X (eg 30 for 30X)
#COVERAGE_CRITERIA="1,30"
COVERAGE_CRITERIA="1,5,10,20,30,50,100,200,300,500,1000"

# COVERAGE DP THRESHOLD (default "30" "100" "1")
# For gene coverage metrics
# the criteria test if genes failed (or just warning) the coverage threshold
SEQUENCING_DEPTH="1" # Sequencing depth threshold
SEQUENCING_COVERAGE_THRESHOLD="1" # Sequencing coverage threshold
MINIMUM_DEPTH="100" # fail DP threshold (default 30X)
EXPECTED_DEPTH="300" # warn DP threshold (default 100X)
DEPTH_COVERAGE_THRESHOLD="1" # threshold percentage of bases over the DP threshold
