#!/bin/bash

# DEFAULT ENV
######################
source_app SOMATIC


# APPLICATION INFOS
#####################
APP_NAME="HEMATOLOGY"
APP_RELEASE="1.1"
APP_DESCRIPTION="Application to detect somatic mutations in hematology gene panel sequencing data"
APP_GROUP="SOMATIC"
APP_PROJECT="HEMATOLOGY"


# ANALYSIS PARAMETERS
#######################

# PIPELINES
# Add the ITDSeek pipeline to the existing pipelines
PIPELINES=$PIPELINES" bwamem.itdseek.howard"
