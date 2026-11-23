############################
# Clipping Rules
# Author: Antony Le Bechec
############################
MK_RELEASE=1.0.0
MK_DATE=05/11/2026

## Release note
# 25/07/2014: DEV into PROD. Creation V0.9. Splitting clipping. Option file from BED file. Checking Option file empty. Merge resulting SAM
# 01/06/2015: New release of FATBAM
# 16/06/2015: Add option CLIPPING to skip clipping step
# 29/09/2016: New clipping script including multithreading
# 27/08/2019: Changing FATBAM to CAP tool
# 03/02/2023: Extract clipping rule
# 05/11/2026: Change compression to 1, change threads to THREADS_BY_ALIGNER


############
# Clipping #
############

# Threads and temporary folder for CAP tool
CAP_THREADS?=$(THREADS_BY_ALIGNER)
CAP_TMP_FOLDER?=$(TMP_FOLDER_TMP)

# CAP Option
# This variable defines additional options to pass to the CAP tool during clipping.
# It can be set directly in the Makefile or inherited from the environment.
CAP_OPTIONS?=


# Rule: Soft Clipping primer BAM file from a manifest
%.bam: %.clipping.bam %.clipping.bam.bai %.manifest
	# Clipping
	+$(CAP) --function=clipping --env=$(CONFIG_TOOLS) --ref=$(GENOME) --bam=$< --manifest=$*.manifest --output=$@ --compress=1 --threads=$(CAP_THREADS) --bedtools=$(BEDTOOLS) --samtools=$(SAMTOOLS) --picard=$(PICARD) --tmp=$(CAP_TMP_FOLDER) $(CAP_OPTIONS)
	# Clean
	-rm $*.clipping.*


# Documentation
$(file >>$(RELEASE_INFOS),## POST_ALIGNMENT clipping '$(MK_RELEASE)': Clipping tool for Amplicon sequencing, to remove primers described in a manifest file from BAM files by soft-clipping. Parameters CAP_OPTIONS='$(CAP_OPTIONS)')
$(file >>$(PIPELINES_INFOS),POST_ALIGNMENT:clipping:Primers Soft-Clipping for Amplicon sequencing:CAP_OPTIONS='$(CAP_OPTIONS)')
