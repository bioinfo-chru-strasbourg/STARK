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

# TOOLS & OPTION
CAP_SOFTCLIPTOQ0?=$(CAP_FOLDER)/CAP.SoftClipToQ0.pl
CAP_THREADS?=$(THREADS)
CAP_OPTIONS?=--verbose # Default options + verbose
CAP_TMP_FOLDER?=$(TMP_FOLDER_TMP) # TMP_SYS_FOLDER ???


# Soft Clipping primer bam file from a manifest
%.bam: %.clipping.bam %.clipping.bam.bai %.manifest %.genome
	# Clipping
	+$(CAP) --function=clipping --env=$(CONFIG_TOOLS) --ref=`cat $*.genome` --bam=$< --manifest=$*.manifest --output=$@ --threads=$(CAP_THREADS) --bedtools=$(BEDTOOLS) --tmp=$(CAP_TMP_FOLDER) $(CAP_OPTIONS)
	#-rm $*.unclipped.*
	-rm $*.clipping.*



# Change SoftClip to Quality 0
%.softclippedtoq0.sam: %.bam %.bam.bai %.manifest %.genome
	# Clipping
	$(SAMTOOLS) view $< -h | perl $(CAP_SOFTCLIPTOQ0) -v1 > $@



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CLIPPING: CAP tool removes primers described in BED file from a BAM file. It parallelizes over chromosomes. Options: CAP='$(CAP)', CAP_THREADS=$(CAP_THREADS), CAP_OPTIONS='$(CAP_OPTIONS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# SOFTCLIPPEDTOQ0: CAP tool change softclipped base quality to 0."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "POST_ALIGNMENT:clipping:Primers Soft-Clipping for Amplicon sequencing"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
