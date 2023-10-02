############################
# Clipping Rules
# Release: 0.9.4
# Date: 27/09/2019
# Author: Antony Le Bechec
############################
MK_RELEASE="0.9.4"
MK_DATE="27/09/2019"

## Release note
# 25/07/2014: DEV into PROD. Creation V0.9. Splitting clipping. Option file from BED file. Checking Option file empty. Merge resulting SAM
# 01/06/2015: New release of FATBAM
# 16/06/2015: Add option CLIPPING to skip clipping step
# 29/09/2016: New clipping script including multithreading
# 27/08/2019: Changing FATBAM to CAP tool
# 03/02/2023: Extract softclippedtoq0 rule

# TOOLS & OPTION
CAP_SOFTCLIPTOQ0?=$(CAP_FOLDER)/CAP.SoftClipToQ0.pl


# Change SoftClip to Quality 0
%.softclippedtoq0.sam: %.bam %.bam.bai %.manifest
	# Clipping
	$(SAMTOOLS) view $< -h | perl $(CAP_SOFTCLIPTOQ0) -v1 > $@



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# SOFTCLIPPEDTOQ0: CAP tool change softclipped base quality to 0."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )
