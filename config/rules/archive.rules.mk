############################
# CRAM Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.1b"
MK_DATE="04/04/2017"

# Release note
# 12/10/2016: Creation

# TOOLS
SAMTOOLS?=$(NGSbin)/samtools

# OPTIONS
THREADS_BY_SAMPLE?=1
MEMORY?=1


CRAM_OPTIONS?=


## FASTQ from ILLUMINA ##


%.archive.cram: %.bams.list %.genome %.R1.fastq.gz %.R2.fastq.gz $(REF_CACHE_FOLDER)
	# Archive aligned BAM only if all original reads present. otherwise, FASTQ compressed 
	mkdir -p $@.metrics;
	$(STARK_FOLDER_BIN)/STARK.archive --fastq="$*.R1.fastq.gz $*.R2.fastq.gz" --bam="$$(cat $< | grep "$*")" --cram=$@ --cram_read_integrity_checksum --cram_options=$(CRAM_OPTIONS) --threads=$(THREADS_BY_ALIGNER) --stats=$@.metrics/stats.tsv --verbose

	


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CRAM '$(MK_RELEASE)': CRAM tool generate *.archive.cram files."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )
