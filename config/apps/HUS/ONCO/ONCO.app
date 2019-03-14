#!/bin/bash
## STARK application ONCO

# DEFAULT ENV
######################
source_app "HUS"

# APPLICATION INFOS
#####################
APP_NAME="ONCO"
APP_RELEASE="1.0"
APP_DESCRIPTION="Applications des Hopitaux universitaire de Strasbourg - Département d'onco-génétique"
APP_GROUP="ONCO"
APP_PROJECT=""


# FOLDERS
###########

# TOOLS FOLDER
# tools: All tools needed for STARK, and more, including STARK
FOLDER_TOOLS=/home1/TOOLS/tools

# DATABASES FOLDER
# genomes: All references genomes. format: $FOLDER_GENOMES/$ASSEMBLY/$ASSEMBLY.fa (with ASSEMBLY=hg19, hg38, mmu19...)
# db: Folder with mandatory databases: dbSNP database for variant calling, such as "dbsnp_138.hg19.vcf.gz" (mandatory, depending on ASSEMBLY), VCF databases for recalibration such as "dbsnp_137.hg19.vcf" (mandatory, depending on ASSEMBLY)
FOLDER_DATABASES=/home1/TOOLS/databases

# RUN FOLDER
# Illumina Sequencer repository Folder. Subfolder as runs
FOLDER_RUN=/home1/IRC/DATA/RAW/MSR

# MANIFEST FOLDER
# Illumina Manifests repository.
# Files to provide in the SampleSheet of each run
FOLDER_MANIFEST=/home1/IRC/DATA/RAW/Manifests


# RESULTS FOLDER
# All results will be generated in this folder :
# RES: RUN files such as BAM, VCF, metrics
# DEM: Demultiplexing folder
# LOG: log files
# TMP: temporary files
FOLDER_RESULTS=/home1/IRC/DATA

# REPOSITORY FOLDER
# Results data can be copy in a repository folder. leave it blank for no copy
FOLDER_REPOSITORY=/home1/L


# PARAMETERS
######################

# Rules
# Add specific rules to load.
# These files will be added to the list of rules files from APPS folder
# example: "MYGROUP/*.rules.mk" "$APP_FOLDER/*.rules.mk"
# for inheritance, use RULES_APP="$RULES_APP MYGROUP/*.rules.mk"
#RULES_APP="$RULES_APP HUS/ONCO/*.rules.mk"
RULES_APP="$RULES_APP $APP_FOLDER/*.rules.mk"


# POST ALIGNEMENT STEPS (Amplicon)
POST_ALIGNMENT_STEPS="sorting realignment clipping compress"


# HOWARD Configuration files for Annotation and Prioritization
# Use APP_FOLDER if necessary
# Example: HOWARD_CONFIG_ANNOTATION=$APP_FOLDER/config.annotation.stark.ini
# Example: HOWARD_CONFIG_ANNOTATION=$STARK_FOLDER_APPS/MY_APP_GROUP/config.annotation.stark.ini
#HOWARD_CONFIG_ANNOTATION=$APP_FOLDER/config.annotation.stark.ini
#HOWARD_CONFIG_PRIORITIZATION=$APP_FOLDER/config.prioritization.stark.ini
