############################
# CRAM Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.2.1"
MK_DATE="20/09/2021"

# Release note
# 0.9.1b-12/10/2016: Creation
# 0.9.2.0-12/04/2021: Archive from script STARK.archive
# 0.9.2.1-20/09/2021: Add option for STARK.archive script, no --cram_read_integrity_checksum if KEEP_ALIGNMENT


STARK_ARCHIVE_OPTIONS?=$(shell if ! (( $(KEEP_ALIGNMENT) )); then echo "--cram_read_integrity_checksum"; else echo ""; fi;)

## FASTQ from ILLUMINA ##

%.archive.cram: %.bams.list %.genome %.R1.fastq.gz %.R2.fastq.gz $(REF_CACHE_FOLDER)
	# Archive aligned BAM only if all original reads present. otherwise, FASTQ compressed 
	# KEEP_ALIGNMENT: $(KEEP_ALIGNMENT)
	mkdir -p $@.metrics;
	#$(STARK_FOLDER_BIN)/STARK.archive --fastq="$*.R1.fastq.gz $*.R2.fastq.gz" --bam="$$(cat $< | grep "$*")" --cram=$@ --cram_read_integrity_checksum --cram_options=$(CRAM_OPTIONS) --threads=$(THREADS_BY_ALIGNER) --sort --stats=$@.metrics/stats.tsv --reheader --verbose
	#$(STARK_FOLDER_BIN)/STARK.archive --fastq="$*.R1.fastq.gz $*.R2.fastq.gz" --bam="$$(cat $< | grep "$*")" --cram=$@ --cram_read_integrity_checksum --cram_options=$(CRAM_OPTIONS) --threads=$(THREADS_BY_ALIGNER) --stats=$@.metrics/stats.tsv --reheader --verbose
	$(STARK_FOLDER_BIN)/STARK.archive --fastq="$*.R1.fastq.gz $*.R2.fastq.gz" --bam="$$(cat $< | grep "$*")" --cram=$@ $(STARK_ARCHIVE_OPTIONS) --cram_options=$(CRAM_OPTIONS) --threads=$(THREADS_BY_ALIGNER) --stats=$@.metrics/stats.tsv --reheader --remove_tags=$(CRAM_REMOVE_TAGS) --verbose

	


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CRAM '$(MK_RELEASE)': CRAM tool generate *.archive.cram files."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )
