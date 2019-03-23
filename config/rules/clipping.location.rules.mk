############################
# FATBAM Clipping Rules
# Release: 0.9.3
# Date: 01/06/2016
# Author: Antony Le Bechec
############################
MK_RELEASE="0.9.3beta"
MK_DATE="29/09/2016"

## Release note
# 25/07/2014: DEV into PROD. Creation V0.9. Splitting clipping. Option file from BED file. Checking Option file empty. Merge resulting SAM
# 01/06/2015: New release of FATBAM
# 16/06/2015: Add option CLIPPING to skip clipping step
# 29/09/2016: New clipping script including multithreading

# TOOLS
#FATBAM?=$(NGSscripts)
FATBAM_CLIPPING?=$(FATBAM)/FATBAM.clipping.sh
FATBAM_SOFTCLIPTOQ0?=$(FATBAM)/FATBAM.SoftClipToQ0.pl
CLIPPING?=1

# OPTIONS
THREADS_FATBAM?=$(THREADS_BY_SAMPLE)
FATBAM_OPTIONS?=--verbose # Default options + verbose
FATBAM_TMP_FOLDER?=$(TMP_FOLDER_TMP) # TMP_SYS_FOLDER ???

#%.bam: %.unclipped.bam %.unclipped.bam.bai %.manifest
%.bam: %.clipping.bam %.clipping.bam.bai %.manifest %.genome
	# Clipping
	#+$(FATBAM_CLIPPING) --env=$(ENV) --bam=$< --manifest=$*.manifest --output=$@ --verbose --multithreading --threads=$(THREADS) --tmp=$(FATBAM_TMP_FOLDER)
	+$(FATBAM_CLIPPING) --env=$(CONFIG_TOOLS) --ref=`cat $*.genome` --bam=$< --manifest=$*.manifest --output=$@ --verbose --multithreading --threads=$(THREADS) --tmp=$(FATBAM_TMP_FOLDER)
	#-rm $*.unclipped.*
	-rm $*.clipping.*



# Change SoftClip to Quality 0
%.softclippedtoq0.sam: %.bam %.bam.bai %.manifest %.genome
	# Clipping
	#+$(FATBAM_CLIPPING) --env=$(ENV) --bam=$< --manifest=$*.manifest --output=$@ --verbose --multithreading --threads=$(THREADS) --tmp=$(FATBAM_TMP_FOLDER)
	$(SAMTOOLS) view $< -h | perl $(FATBAM_SOFTCLIPTOQ0) -v1 > $@



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CLIPPING: FATBAM tool removes primers described in BED file from a BAM file. It parallelizes over chromosomes. Options: FATBAM='$(FATBAM)', FATBAM_CLIPPING='$(FATBAM_CLIPPING)', FATBAM_COVERAGE='$(FATBAM_COVERAGE)', FATBAM_OPTIONS='$(FATBAM_OPTIONS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# SOFTCLIPPEDTOQ0: FATBAM tool change softclipped base quality to 0."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "POST_ALIGNMENT:clipping:Primers Soft-Clipping for Amplicon sequencing"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
