#!/bin/bash
## STARK application ONCO

# DEFAULT ENV
######################
source_app "HUSDIAGGEN"

# APPLICATION INFOS
#####################
APP_NAME="HUSDIAGGEN.HCSOP"
APP_RELEASE="1.0"
APP_GROUP="HUSDIAGGEN"
APP_PROJECT="HCSOP"


# ANALYSIS PARAMETERS
#######################

# Rules
RULES_APP="$RULES_APP HUS/ONCO/HUSDIAGGEN.rules.mk/*.rules.mk"

# PIPELINES
PIPELINES="bwamem.gatkHC_CPSGEN_MASTR.howard bwamem.gatkUG_CPSGEN_MASTR.howard bwamem.canoes.howard"

# POST_ALIGNMENT (Capture)
POST_ALIGNMENT_STEPS="sorting markduplicates realignment compress"

# HOWARD FILTER
################
HOWARD_PRIORITIZATION="HUSDIAGGEN.MASTR.HARD,HUSDIAGGEN.MASTR,HUSDIAGGEN.MASTR.SOFT,GERMLINE"
HOWARD_PRIORITIZATION_MINIMAL="HUSDIAGGEN.MASTR.HARD,HUSDIAGGEN.MASTR,HUSDIAGGEN.MASTR.SOFT,GERMLINE"
HOWARD_PRIORITIZATION_REPORT="HUSDIAGGEN.MASTR.HARD,HUSDIAGGEN.MASTR,HUSDIAGGEN.MASTR.SOFT,GERMLINE"

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


# INTERVAL_PADDING (default 0)
# Add some “padding” to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=25
