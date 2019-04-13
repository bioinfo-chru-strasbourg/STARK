#!/bin/bash
#################################
## STARK environment
#################################

# APPLICATION INFOS
#####################
APP_NAME=""
APP_RELEASE=""
APP_DESCRIPTION=""
APP_GROUP=""
APP_PROJECT=""

# FOLDERS
###########

# TOOLS FOLDER
# tools: All tools needed for STARK, and more, including STARK
FOLDER_TOOLS=

# DATABASES FOLDER
# genomes: All references genomes. format: $FOLDER_GENOMES/$ASSEMBLY/$ASSEMBLY.fa (with ASSEMBLY=hg19, hg38, mmu19...)
# db: Folder with mandatory databases: dbSNP database for variant calling, such as "dbsnp_138.hg19.vcf.gz" (mandatory, depending on ASSEMBLY), VCF databases for recalibration such as "dbsnp_137.hg19.vcf" (mandatory, depending on ASSEMBLY)
FOLDER_DATABASES=/STARK/databases

# RUN FOLDER
# Illumina Sequencer repository Folder. Subfolder as runs
FOLDER_RUN=

# MANIFEST FOLDER
# Illumina Manifests repository.
# Files to provide in the SampleSheet of each run
FOLDER_MANIFEST=/STARK/manifests

# RESULTS FOLDER
# All results will be generated in this folder :
# RES: RUN files such as BAM, VCF, metrics
# DEM: Demultiplexing folder
# LOG: log files
# TMP: temporary files
FOLDER_RESULTS=/STARK/results

# REPOSITORY FOLDER
# Results data can be copy in a repository folder. leave it blank for no copy
FOLDER_REPOSITORY=/STARK/repository

# ARCHIVE FOLDER
# Results data can be copy in a archive folder. leave it blank for no copy
FOLDER_ARCHIVE=/STARK/archive



# VARIABLE INITIALISATION
############################

APP_NAME=""
APP_GROUP=""
APP_PROJECT=""

# RULES for APP
RULES_APP=""
RULES=""

# FILE to provide in folder repository
RESULTS_SUBFOLDER_ROOT_FILE_PATTERNS=""

# PIPELINES
PIPELINES=""
ALIGNERS=""
CALLERS=""
ANNOTATORS=""

# PANEL/MANIFEST/BED
MANIFEST=""
BED=""


# HOWARD FILTER DEFAULT
HOWARD_PRIORITIZATION_DEFAULT="default"
HOWARD_PRIORITIZATION_REPORT="default"
HOWARD_PRIORITIZATION_MINIMUM="default"
#ANNOVAR_DATABASES=""
#SNPEFF_DATABASES=""

# SYSTEM Config
export LANG=C
