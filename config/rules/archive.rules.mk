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
# 0.9.3.0-18/09/2026: Add BAM list from sample to ensure only BAM files with corresponding BAI indexes are included, select and list BAM in a comma separated format


STARK_ARCHIVE_OPTIONS?=$(shell if ! (( $(KEEP_ALIGNMENT) )); then echo "--cram_read_integrity_checksum"; else echo ""; fi;)


%.archive.cram: %.bams_sample.list %.R1.fastq.gz %.R2.fastq.gz
	# Archive aligned BAM only if all original reads present. otherwise, FASTQ compressed 
	mkdir -p $@.metrics;
	$(STARK_FOLDER_BIN)/STARK.archive --fastq="$*.R1.fastq.gz $*.R2.fastq.gz" --bam="$$(cat $*.bams_sample.list | grep "$*" | tr "\n" "," | tr " " ",")" --cram=$@ $(STARK_ARCHIVE_OPTIONS) --cram_options=$(CRAM_OPTIONS) --threads=$(THREADS_BY_ALIGNER) --stats=$@.metrics/stats.tsv --reheader --remove_tags=$(CRAM_REMOVE_TAGS) --verbose

	


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CRAM '$(MK_RELEASE)': CRAM tool generate *.archive.cram files."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )
