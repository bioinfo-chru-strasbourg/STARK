#!/bin/bash
## STARK application ONCO

# DEFAULT ENV
######################
source_app "ONCO"

# APPLICATION INFOS
#####################
APP_NAME="HUSTUMSOL"
APP_RELEASE="1.0"
APP_DESCRIPTION="Applications des Hopitaux universitaire de Strasbourg - Département d'onco-génétique - Laboratoire biologie moléculaire"
APP_GROUP="HUSTUMSOL"
APP_PROJECT=""


# ANALYSIS PARAMETERS
#######################


# Rules
# Add specific rules to load.
# These files will be added to the list of rules files from APPS folder
# example: "MYGROUP/*.rules.mk"
# for inheritance, use RULES_APP="$RULES_APP MYGROUP/*.rules.mk"
#RULES_APP="$RULES_APP HUS/ONCO/HUSTUMSOL.rules.mk/*.rules.mk"
RULES_APP="$RULES_APP $APP_FOLDER/HUSTUMSOL.rules.mk/*.rules.mk"



# PIPELINES
PIPELINES="bwamem.gatkHC_HUSTUMSOL.howard bwamem.gatkUG_HUSTUMSOL.howard bwamem.VarScan_HUSTUMSOL.howard bwamem.canoes.howard"


# HOWARD FILTER
################
HOWARD_PRIORITIZATION="TUMSOL,SOLIDTUMOR"
HOWARD_PRIORITIZATION_MINIMAL="TUMSOL,SOLIDTUMOR"
HOWARD_PRIORITIZATION_REPORT="TUMSOL,SOLIDTUMOR"

# COVERAGE CRITERIA (default "1,30")
COVERAGE_CRITERIA="1,300,1000"

# COVERAGE DP THRESHOLD (default "30" "100" "1")
# For gene coverage metrics
# the criteria test if genes failed (or just warning) the coverage threshold
DP_FAIL="300"	# fail DP threshold (default 30X)
DP_WARN="1000"	# warn DP threshold (default 100X)
DP_THRESHOLD="1" # threshold percentage of bases over the DP threshold
