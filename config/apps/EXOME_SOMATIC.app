#!/bin/bash

# DEFAULT ENV
######################
# Depend on SOMATIC App and CAPTURE PlugApp
source_app SOMATIC,CAPTURE


# APPLICATION INFOS
#####################
APP_NAME="EXOME_SOMATIC"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect somatic mutations in exome sequencing data"
APP_GROUP="SOMATIC"
APP_PROJECT="EXOME"


# METRICS
##########

# INTERVAL_PADDING / add some padding to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=100
