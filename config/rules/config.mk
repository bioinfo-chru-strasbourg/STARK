##############################
# CONFIGURATION Rules
# Release: 0.9.2
# Date: 02/02/2015
# Author: Antony Le Bechec
##############################


#MK_PATH=$(abspath $(lastword $(MAKEFILE_LIST)))
#MK_DIR_PATH=$(shell dirname $(MK_PATH))

#NGS_SCRIPTS?=$(MK_DIR_PATH)
NGSscripts?=$(NGS_SCRIPTS)
NGSEnv?=/tool/
STARK_FOLDER_ROOT?=/tool
STARK_FOLDER_CONFIG?=$(STARK_FOLDER_ROOT)/config
STARK_FOLDER_BIN?=$(STARK_FOLDER_ROOT)/bin
STARK_FOLDER_APPS?=$(STARK_FOLDER_CONFIG)/apps
STARK_FOLDER_RULES?=$(STARK_FOLDER_CONFIG)/rules
ENV?=$(STARK_FOLDER_CONFIG)/default.app
#NGSscripts?=$(NGSEnv)/scripts
#PWD=.


# RELEASE FILE
RELEASE_FILE?=$(shell source $(STARK_FOLDER_CONFIG)/config.app; source_app $(ENV); echo $$RELEASE_FILE)
THREADS?=$(shell source $(STARK_FOLDER_CONFIG)/config.app; source_app $(ENV); echo $$THREADS)

# THREAD Calculation
#THREADS?=$(shell ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w) # DYNAMIC
#THREADS?=4
#THREADS_BY_SAMPLE?=$(THREADS)


# FOLDERS
MISEQDIR?=$(MISEQ_FOLDER)
INPUTDIR?=$(DEMULTIPLEXING_FOLDER)
OUTDIR?=$(RESULTS_FOLDER)
TMP_FOLDER_TMP?=/tmp				# NGS Temporary folder
TMP_SYS_FOLDER?=/tmp				# System Temporary folder
ASSEMBLY?=hg19					# Default assembly
GENOMES?=genomes				# Genomes folder
#REF?=$(GENOMES)/$(ASSEMBLY)/$(ASSEMBLY).fa	# Default Reference genome FASTA file
REF?=$(GENOMES)/current/$(ASSEMBLY).fa	# Default Reference genome FASTA file

# OPTIONS
JAVA_MEMORY?=4
JAVA_MEMORY_BY_SAMPLE?=4
JAVA_FLAGS_TMP_FOLDER?=-Dorg.xerial.snappy.tempdir=$(TMP_FOLDER_TMP) -Djava.io.tmpdir=$(TMP_FOLDER_TMP)
JAVA_FLAGS_OTHER_PARAM?=-Dsnappy.disable=true
JAVA_FLAGS?= -Xmx$(JAVA_MEMORY)g $(JAVA_FLAGS_TMP_FOLDER) $(JAVA_FLAGS_OTHER_PARAM)
JAVA_FLAGS_BY_SAMPLE?= -Xmx$(JAVA_MEMORY_BY_SAMPLE)g $(JAVA_FLAGS_TMP_FOLDER) $(JAVA_FLAGS_OTHER_PARAM)

# PICARD
PICARD_FLAGS?=-SORT_ORDER coordinate -RGLB 001 -RGPL ILLUMINA -RGPU PU -VALIDATION_STRINGENCY SILENT
PICARD_UNALIGNED_FLAGS?=-COMPRESSION_LEVEL 1 -MAX_RECORDS_IN_RAM 500000
PICARD_UNALIGNED_NAME_FLAGS?=-LIBRARY_NAME 001 -PLATFORM ILLUMINA -PLATFORM_UNIT PU -READ_GROUP_NAME A

# DATABASES
DBFOLDER?=/media/IRCV2/NGSEnv/annovar_sources
VCFDBSNP?=$(DBFOLDER)/snp138.vcf.gz
VCFDBSNP137VCF?=$(DBFOLDER)/dbsnp_137.hg19.vcf
VCF1000G?=$(DBFOLDER)/1000G_phase1.indels.hg19.vcf
VCFMILLS1000G?=$(DBFOLDER)/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
COSMIC?=$(DBFOLDER)/COSMIC.CodingMuts.vcf
KNOWN_ALLELES?=$(VCFMILLS1000G)
HAPMAP?=$(DBFOLDER)/hapmap_3.3.hg19.sites.vcf
OMNI?=$(DBFOLDER)/1000G_omni2.5.hg19.vcf
PHASE1_1000G?=$(DBFOLDER)/1000G_phase1.snps.high_confidence.hg19.sites.vcf


# BAM Check option
BAM_CHECK_STEPS?=1

# BAM METRICS
BAM_METRICS?=1

# BAM compression
BAM_COMPRESSION?=5

# POST_SEQUENCING
POST_SEQUENCING?=

# POST_ALIGNMENT
POST_ALIGNMENT?=
#.sorting

# POST_CALLING
POST_CALLING?=
#.normalization.sorting
#.filtration.recalibration

# POST_ANNOTATION
POST_ANNOTATION?=
#.normalization.sorting
#.filtration.recalibration


# HEADER
#RELEASE_CMD := $(shell echo "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#" >> $(RELEASE_INFOS) )
#RELEASE_CMD := $(shell echo "\# $(ENV_NAME) - $(ENV_DESCRIPTION)" >> $(RELEASE_INFOS) )
#RELEASE_CMD := $(shell echo "\# RELEASE $(ENV_RELEASE) - $(ENV_DATE) - COPYRIGHT © $(ENV_COPYRIGHT) - $(ENV_AUTHOR) ($(ENV_LICENCE) licence)" >> $(RELEASE_INFOS) )
#RELEASE_CMD := $(shell echo "\# CONFIG  Scripts '$(NGSscripts)' - ENV '$(ENV)' - TOOLS '$(NGSEnv)'" >> $(RELEASE_INFOS) )
#RELEASE_CMD := $(shell echo "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#" >> $(RELEASE_INFOS) )

RELEASE_CMD := $(shell echo "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\# $(ENV_NAME) - $(ENV_DESCRIPTION)" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\# RELEASE $(ENV_RELEASE) [$(ENV_DATE)]" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\# AUTHORS $(ENV_AUTHOR)" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\# COPYRIGHT $(ENV_COPYRIGHT)" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\# LICENCE $(ENV_LICENCE)" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#" >> $(RELEASE_INFOS) )


# RELEASE
RELEASE_CMD := $(shell echo "\#\# STARK_RELEASE: $(ENV_RELEASE)" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\# APP_NAME: $(APP_NAME)" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\# APP_RELEASE: $(APP_RELEASE)" >> $(RELEASE_INFOS) )

# GROUP/PROJECT
RELEASE_CMD := $(shell echo "\#\# GROUP: $(GROUP)" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\# PROJECT: $(PROJECT)" >> $(RELEASE_INFOS) )

# INFRASTRUCTURE
RELEASE_CMD := $(shell echo "\#\# INFRASTRUCTURE: THREADS=$(THREADS), THREADS_BY_SAMPLE=$(THREADS_BY_SAMPLE), MEMORY=$(MEMORY)" >> $(RELEASE_INFOS) )

RELEASE_CMD := $(shell echo "" >> $(RELEASE_INFOS) )

# TOOLS
#RELEASE_CMD := $(shell echo "\#\# TOOLS INFOS" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell $(STARK_FOLDER_BIN)/STARK --tools_infos --app=$(ENV) >> $(RELEASE_INFOS) )

# DATABASES
#RELEASE_CMD := $(shell echo "\#\# DATABASES INFOS" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell $(STARK_FOLDER_BIN)/STARK --databases_infos --app=$(ENV) >> $(RELEASE_INFOS) )

# OTHER INFORMATIONS
RELEASE_CMD := $(shell echo "" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#\#\# OTHER INFORMATIONS \#\#\#\#" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#" >> $(RELEASE_INFOS) )


# DATABASES
#RELEASE_COMMENT := "\#\# DATABASES: VCFDBSNP='$(VCFDBSNP)', VCFDBSNP137VCF='$(VCFDBSNP137VCF)', VCF1000G='$(VCF1000G)', VCFMILLS1000G='$(VCFMILLS1000G)', KNOWN_ALLELES='$(KNOWN_ALLELES)', COSMIC='$(COSMIC)'"
#RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

# TOOLS
#RELEASE_COMMENT := "\#\# TOOLS: JAVA=$(JAVA_VERSION)\&$(JAVA6_VERSION), BWA=$(BWA_VERSION), GATK=$(GATK_VERSION), SAMTOOLS=$(SAMTOOLS_VERSION), IGV=$(IGV_VERSION), IGVTOOLS=$(IGVTOOLS_VERSION), PICARD=$(PICARD_VERSION), TABIX=$(TABIX_VERSION), BCFTOOLS=$(BCFTOOLS_VERSION), FASTQC=$(FASTQC_VERSION), CAP=$(CAP_VERSION), CAP_ManifestToBED=$(CAP_ManifestToBED_VERSION), CAP_BEDToOPTION=$(CAP_BEDToOPTION_VERSION)"
#RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CONFIG FILES: CONFIG=$(CONFIG), FUNCTIONS=$(FUNCTIONS), PARAM=$(PARAM)"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

# OTHER INFORMATIONS
RELEASE_CMD := $(shell echo "" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#\#\# INFORMATIONS ON RULES \#\#\#\#" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#" >> $(RELEASE_INFOS) )
