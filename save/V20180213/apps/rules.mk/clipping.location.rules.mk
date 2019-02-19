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
CLIPPING?=1

# OPTIONS
THREADS_FATBAM?=$(THREADS_BY_SAMPLE)
FATBAM_OPTIONS?=--verbose # Default options + verbose
FATBAM_TMP_FOLDER?=$(TMP_FOLDER_TMP) # TMP_SYS_FOLDER ???

#%.bam: %.unclipped.bam %.unclipped.bam.bai %.manifest
%.bam: %.clipping.bam %.clipping.bam.bai %.manifest
	# Clipping
	+$(FATBAM_CLIPPING) --env=$(ENV) --bam=$< --manifest=$*.manifest --output=$@ --verbose --multithreading --threads=$(THREADS) --tmp=$(FATBAM_TMP_FOLDER)
	# CHECK NUMBER of READS in BAM
	if (($(BAM_CHECK_STEPS))); then \
		#echo $$($(SAMTOOLS) view -c -@ $(THREADS_SAMTOOLS) -F 0x0100 $<)" READS for $<"; \
		#echo $$($(SAMTOOLS) view -c -@ $(THREADS_SAMTOOLS) -F 0x0100 $@)" READS for $@"; \
		#$(SAMTOOLS) idxstats $< | awk '{SUM+=$$3+$$4} END {print SUM}' \
		#if [ ! -e $@.idx ]; then $(SAMTOOLS) index $@; fi; \
		#if [ "$(SAMTOOLS) idxstats $< | awk '{SUM+=$$3+$$4} END {print SUM}'" != "$(SAMTOOLS) idxstats $@ | awk '{SUM+=$$3+$$4} END {print SUM}'" ]; then \
		if [ "$$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $<)" != "$$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $@)" ]; then \
			echo "# ERROR in Number of reads between $< and $@ !!!"; \
			echo "# BCFTOOLS STATS for $<";  \
			$(SAMTOOLS) index $<; \
			$(SAMTOOLS) stats $< | grep ^SN; \
			$(SAMTOOLS) idxstats $<; \
			echo "# BCFTOOLS STATS for $@";  \
			$(SAMTOOLS) index $@; \
			$(SAMTOOLS) stats $@ | grep SN; \
			$(SAMTOOLS) idxstats $@; \
			exit 1; \
		else \
			echo "# Number of reads OK between $< and $@"; \
		fi; \
	fi;
	#-rm $*.unclipped.*
	-rm $*.clipping.*
	


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CLIPPING: FATBAM tool removes primers described in BED file from a BAM file. It parallelizes over chromosomes. Options: FATBAM='$(FATBAM)', FATBAM_CLIPPING='$(FATBAM_CLIPPING)', FATBAM_COVERAGE='$(FATBAM_COVERAGE)', FATBAM_OPTIONS='$(FATBAM_OPTIONS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "POST_ALIGNMENT:clipping:Prime Clipping for Amplicon sequencing"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


