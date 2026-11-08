############################
# Clipping Rules
# Release: 1.0.0
# Date: 05/11/2026
# Author: Antony Le Bechec
############################
MK_RELEASE="1.0.0"
MK_DATE="05/11/2026"

## Release note
# 25/07/2014: DEV into PROD. Creation V0.9. Splitting clipping. Option file from BED file. Checking Option file empty. Merge resulting SAM
# 01/06/2015: New release of FATBAM
# 16/06/2015: Add option CLIPPING to skip clipping step
# 29/09/2016: New clipping script including multithreading
# 27/08/2019: Changing FATBAM to CAP tool
# 03/02/2023: Extract clipping rule
# 05/11/2026: Change compression to 1, change threads to THREADS_BY_ALIGNER


# TOOLS & OPTION
CAP_THREADS?=$(THREADS_BY_ALIGNER)
CAP_OPTIONS?= 	#--verbose # Default options + verbose
CAP_TMP_FOLDER?=$(TMP_FOLDER_TMP) # TMP_SYS_FOLDER ???


# Soft Clipping primer bam file from a manifest
%.bam: %.clipping.bam %.clipping.bam.bai %.manifest
	# Clipping
	+$(CAP) --function=clipping --env=$(CONFIG_TOOLS) --ref=$(GENOME) --bam=$< --manifest=$*.manifest --output=$@ --compress=1 --threads=$(CAP_THREADS) --bedtools=$(BEDTOOLS) --samtools=$(SAMTOOLS) --picard=$(PICARD) --tmp=$(CAP_TMP_FOLDER) $(CAP_OPTIONS)
	-rm $*.clipping.*



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CLIPPING: CAP tool removes primers described in BED file from a BAM file. It parallelizes over chromosomes. Options: CAP='$(CAP)', CAP_THREADS=$(CAP_THREADS), CAP_OPTIONS='$(CAP_OPTIONS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "POST_ALIGNMENT:clipping:Primers Soft-Clipping for Amplicon sequencing"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
