############################
# BWA Aligner Rules
# Release: 0.9.4.0
# Date: 27/10/2025
# Author: Antony Le Bechec
############################

# Release note
# 25/07/20140.2b: add clipping step (".unclipped" on targets)
# 10/03/20150.9.1b: change genome reference location, in the file %.genome
# 29/09/2016-0.9.2b: Cleaning, PICARD new release picard.jar
# 13/04/2021-0.9.3.0: Cleaning, removing old BWA alignment release
# 23/05/2021-0.9.3.1: Remove samtools view step
# 27/10/2025-0.9.4.0: Add function to manage memory for BWA MEM



###################################
# BWA-MEM By Default (FROM FASTQ) #
###################################

## BWA MEM (Last powerful algorithm, including SW, HMM...)
## A lockfile is used to check BWA ressources. This is due to occasional BWA failures on some systems (e.g. memory issue).

# Options

# BWA MEM flags to set specific parameters (see BWA help)
BWAMEM_FLAGS?=-C -M #-t $(THREADS_BY_ALIGNER)

# Add flag to post alignment samtools view to remove specific reads (e.g. secondary, supplementary)
BWAMEM_SAMTOOLS_FILTER_FLAG?=

# Memory management
# Use a minimum of 6 Go per BWA MEM process
# The calculation of MAX_CONCURRENT_ALIGNMENTS_BWAMEM is done to avoid overloading the system memory
BWAMEM_MIN_MEM?=6
MAX_CONCURRENT_ALIGNMENTS_BWAMEM?=$(shell if [ $(shell echo "$(MEMTOTAL_IN_GO)/$(BWAMEM_MIN_MEM)/$(NB_ALIGNERS)" | bc) -lt 1 ]; then echo 1; else echo "$(MEMTOTAL_IN_GO)/$(BWAMEM_MIN_MEM)/$(NB_ALIGNERS)" | bc; fi)

# Threads per BWA MEM process
THREADS_BWAMEM?=$(shell echo " if ($(MAX_CONCURRENT_ALIGNMENTS_BWAMEM)<$(NB_SAMPLE)) ($(THREADS)/$(MAX_CONCURRENT_ALIGNMENTS_BWAMEM)) else ($(THREADS)/$(NB_SAMPLE))" | bc)


%.bwamem$(POST_ALIGNMENT).bam: %.R1$(POST_SEQUENCING).fastq.gz %.R2$(POST_SEQUENCING).fastq.gz
	# List of FASTQs
	if (($$(zcat $*.R2.fastq.gz | head -n 1 | wc -l))); then \
		echo "$*.R1$(POST_SEQUENCING).fastq.gz $*.R2$(POST_SEQUENCING).fastq.gz" > $@.fastq_list; \
	else \
		echo "$*.R1$(POST_SEQUENCING).fastq.gz" > $@.fastq_list; \
	fi;
	# Alignment
	$(PYTHON3) $(STARK_FOLDER_BIN)/functions.py launch \
		--cmd "$(BWA) mem $(BWAMEM_FLAGS) -t $(THREADS_BWAMEM) -R '@RG\tID:1\tPL:ILLUMINA\tPU:PU\tLB:001\tSM:$(*F)' $(GENOME) $$(cat $@.fastq_list) -o $@.sam" \
		--lockfile_prefix $$(echo $@ | xargs -0 dirname | xargs -0 dirname)/lockfile.bwamem. \
		--target $@ \
		--max_jobs $(MAX_CONCURRENT_ALIGNMENTS_BWAMEM);
	# Sorting
	echo "#[INFO] Sorting BAM file for $*:"
	$(SAMTOOLS) view -h $(BWAMEM_SAMTOOLS_FILTER_FLAG) $@.sam -@ $(THREADS_SAMTOOLS) | $(SAMTOOLS) sort -l 1 -O BAM -o $@.tmp -T $@.SAMTOOLS_PREFIX -@ $(THREADS_SAMTOOLS)
	rm $@.sam
	# AddOrReplaceReadGroups
	if (($$($(SAMTOOLS) view $@.tmp -H | grep "^@RG" -c))); then \
		echo "#[INFO] BAM $@.tmp with read group"; \
		mv $@.tmp $@; \
	else \
		echo "#[INFO] BAM $@.tmp without read group"; \
		$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) AddOrReplaceReadGroups $(PICARD_FLAGS) -I $@.tmp O=$@ -COMPRESSION_LEVEL 1 -RGSM $(*F); \
	fi;
	-rm $@.tmp $@.RG $@.fastq_list


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# BWA ALIGNMENT '$(MK_RELEASE)': BWA generates an aligned BAM file from FASTQ file, and ask for post alignment processes 'sorting', 'realignment', 'clipping' \(if needed\) and 'recalibration'. PICARD TOOL is used to Add Or Replace Read Groups and modified the BAM header. Options: BWA='$(BWA)', BWAMEM_FLAGS='$(BWAMEM_FLAGS)', PICARD='$(PICARD)', PICARD_FLAGS='$(PICARD_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


# PIPELINES INFOS
PIPELINES_COMMENT := "ALIGNER:bwamem:BWA MEM - Last powerful algorithm. From FASTQ files."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

