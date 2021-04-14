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




# RELEASE FILE
RELEASE_FILE?=$(shell source $(STARK_FOLDER_CONFIG)/config.app; source_app $(ENV); echo $$RELEASE_FILE)
THREADS?=$(shell source $(STARK_FOLDER_CONFIG)/config.app; source_app $(ENV); echo $$THREADS)

# THREAD Calculation
THREADS?=$(shell ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w) # DYNAMIC
THREADS_BY_SAMPLE?=$(THREADS)
THREADS_BY_ALIGNER?=$(THREADS)
THREADS_BY_CALLER?=$(THREADS)
THREADS_BWA?=$(THREADS_BY_SAMPLE)
THREADS_SAMTOOLS?=$(THREADS_BY_SAMPLE)



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



### TOOLS
BWA?=bwa
SAMTOOLS?=samtools
GZ?=gzip
GUNZIP?=gzip -d
GATK?=GenomeAnalysisTK.jar
IGVTOOLS?=igvtools.jar
JAVA?=java
TABIX?=tabix
FASTQC?=fastqc
BGZIP?=bgzip
GZ?=gzip
HOWARD?=HOWARD
SNPEFF?=snpEff.jar
VARSCAN?=$(NGSbin)/varscan.jar
CAP?=CAP
CAP_SOFTCLIPTOQ0?=CAP.SoftClipToQ0.pl



# JAVA OPTIONS
JAVA_MEMORY?=4
JAVA_MEMORY_BY_SAMPLE?=4
JAVA_FLAGS_TMP_FOLDER?=-Dorg.xerial.snappy.tempdir=$(TMP_FOLDER_TMP) -Djava.io.tmpdir=$(TMP_FOLDER_TMP)
JAVA_FLAGS_OTHER_PARAM?=-Dsnappy.disable=true
JAVA_FLAGS?= -Xmx$(JAVA_MEMORY)g $(JAVA_FLAGS_TMP_FOLDER) $(JAVA_FLAGS_OTHER_PARAM)
JAVA_FLAGS_BY_SAMPLE?= -Xmx$(JAVA_MEMORY_BY_SAMPLE)g $(JAVA_FLAGS_TMP_FOLDER) $(JAVA_FLAGS_OTHER_PARAM)


# PICARD
PICARD?=picard.jar
PICARD_FLAGS?=-SORT_ORDER coordinate -RGLB 001 -RGPL ILLUMINA -RGPU PU -VALIDATION_STRINGENCY SILENT
PICARD_UNALIGNED_FLAGS?=-COMPRESSION_LEVEL 1 -MAX_RECORDS_IN_RAM 500000
PICARD_UNALIGNED_NAME_FLAGS?=-LIBRARY_NAME 001 -PLATFORM ILLUMINA -PLATFORM_UNIT PU -READ_GROUP_NAME A


# HOWARD OPTIONS
HOWARD_ANNOTATION?=$(HOWARD_FOLDER)/VCFannotation.pl
HOWARD_PRIORITIZATION?=$(HOWARD_FOLDER)/VCFprioritization.pl
HOWARD_TRANSLATION?=$(HOWARD_FOLDER)/VCFtranslation.pl



# DATABASES
DBFOLDER?=/STARK/databases
VCFDBSNP?=$(DBFOLDER)/snp138.vcf.gz
VCFDBSNP137VCF?=$(DBFOLDER)/dbsnp_137.hg19.vcf
VCF1000G?=$(DBFOLDER)/1000G_phase1.indels.hg19.vcf
VCFMILLS1000G?=$(DBFOLDER)/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
COSMIC?=$(DBFOLDER)/COSMIC.CodingMuts.vcf
KNOWN_ALLELES?=$(VCFMILLS1000G)
HAPMAP?=$(DBFOLDER)/hapmap_3.3.hg19.sites.vcf
OMNI?=$(DBFOLDER)/1000G_omni2.5.hg19.vcf
PHASE1_1000G?=$(DBFOLDER)/1000G_phase1.snps.high_confidence.hg19.sites.vcf

# HOWARD
HOWARD_CONFIG_OPTIONS?=--config=$(HOWARD_CONFIG) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION) --config_annotation=$(HOWARD_CONFIG_ANNOTATION) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --multithreading --threads=$(THREADS) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS)

# BAM Check option
BAM_CHECK_STEPS?=1

# BAM METRICS
BAM_METRICS?=1

# BAM compression
BAM_COMPRESSION?=5

# CRAM OPTIONS
CRAM_OPTIONS?=version=3.0,level=9

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



# COLORS RGB
PASS_COLOR_RGB?=0,255,0
WARN_COLOR_RGB?=252,161,6
FAIL_COLOR_RGB?=255,0,0
MISS_COLOR_RGB?=139,0,0
FILTERED_COLOR_RGB?=252,161,6



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
#RELEASE_CMD := $(shell echo "\#\# APP_RELEASE: $(APP_RELEASE)" >> $(RELEASE_INFOS) )

# GROUP/PROJECT
#RELEASE_CMD := $(shell echo "\#\# GROUP: $(GROUP)" >> $(RELEASE_INFOS) )
#RELEASE_CMD := $(shell echo "\#\# PROJECT: $(PROJECT)" >> $(RELEASE_INFOS) )

# INFRASTRUCTURE
RELEASE_CMD := $(shell echo "\#\# INFRASTRUCTURE:" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#    THREADS=$(THREADS)" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#    THREADS_BY_SAMPLE=$(THREADS_BY_SAMPLE)" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#    THREADS_BY_PIPELINE=$(THREADS_BY_PIPELINE)" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#    THREADS_BY_ALIGNER=$(THREADS_BY_ALIGNER)" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#    THREADS_BY_CALLER=$(THREADS_BY_CALLER)" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#    MEMORY=$(MEMORY)" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\# CONFIG FILES:" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#    CONFIG=$(CONFIG)" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#    FUNCTIONS=$(FUNCTIONS)" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#    PARAM=$(PARAM)" >> $(RELEASE_INFOS) )


# RELEASE_CMD := $(shell echo "" >> $(RELEASE_INFOS) )

# # APPLICATION
# RELEASE_CMD := $(shell $(STARK_FOLDER_BIN)/STARK --applications_infos --app="$(ENV)" >> $(RELEASE_INFOS) )

# # TOOLS
# #RELEASE_CMD := $(shell echo "\#\# TOOLS INFOS" >> $(RELEASE_INFOS) )
# RELEASE_CMD := $(shell $(STARK_FOLDER_BIN)/STARK --tools_infos --app="$(ENV)" >> $(RELEASE_INFOS) )

# # DATABASES
# #RELEASE_CMD := $(shell echo "\#\# DATABASES INFOS" >> $(RELEASE_INFOS) )
# RELEASE_CMD := $(shell $(STARK_FOLDER_BIN)/STARK --databases_infos --app="$(ENV)" >> $(RELEASE_INFOS) )
# RELEASE_CMD := $(shell echo "" >> $(RELEASE_INFOS) )

# # APPLICATION ALL
# RELEASE_CMD := $(shell $(STARK_FOLDER_BIN)/STARK --applications_infos_all --app="$(ENV)" >> $(RELEASE_INFOS) )


# OTHER INFORMATIONS
RELEASE_CMD := $(shell echo "" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#\#\# INFORMATIONS ON RULES \#\#\#\#" >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#" >> $(RELEASE_INFOS) )
