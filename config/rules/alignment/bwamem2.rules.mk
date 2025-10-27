############################
# BWA Aligner Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.4.1"
MK_DATE="27/10/2025"

# Release note
# 25/07/20140.2b: add clipping step (".unclipped" on targets)
# 10/03/20150.9.1b: change genome reference location, in the file %.genome
# 29/09/2016-0.9.2b: Cleaning, PICARD new release picard.jar
# 13/04/2021-0.9.3.0: Cleaning, removing old BWA alignment release
# 23/05/2021-0.9.3.1: Remove samtools view step
# 12/05/2023-0.9.4.0: BWA MEM2
# 27/10/2025-0.9.4.1: Add checks for BWA2 MEM success



####################################
# BWA-MEM2 By Default (FROM FASTQ) #
####################################

## BWA MEM2 (Last powerful algorithm, including SW, HMM...)
## A lockfile is used to check BWA ressources. This is due to occasional BWA failures on some systems (e.g. memory issue).

# Options
BWAMEM2_FLAGS?= mem -C -M -t $(THREADS_BWA)
BWAMEM2_MIN_MEM?=12
MAX_CONCURRENT_ALIGNMENTS_BWAMEM2?=$(shell if [ $(shell echo "$(MEMTOTAL_IN_GO)/$(BWAMEM2_MIN_MEM)" | bc) -lt 1 ]; then echo 1; else echo "$(MEMTOTAL_IN_GO)/$(BWAMEM2_MIN_MEM)" | bc; fi)


%.bwamem2$(POST_ALIGNMENT).bam: %.R1$(POST_SEQUENCING).fastq.gz %.R2$(POST_SEQUENCING).fastq.gz
	# List of FASTQs
	if (($$(zcat $*.R2.fastq.gz | head -n 1 | wc -l))); then \
		echo "$*.R1$(POST_SEQUENCING).fastq.gz $*.R2$(POST_SEQUENCING).fastq.gz" > $@.fastq_list; \
	else \
		echo "$*.R1$(POST_SEQUENCING).fastq.gz" > $@.fastq_list; \
	fi;
	# Alignment
	$(PYTHON3) $(STARK_FOLDER_BIN)/functions.py launch \
		--cmd "$(BWA2) $(BWAMEM2_FLAGS) -R '@RG\tID:1\tPL:ILLUMINA\tPU:PU\tLB:001\tSM:$(*F)' $(GENOME) $$(cat $@.fastq_list) -o $@.sam 2> $@.bwa.log" \
		--lockfile_prefix $$(echo $@ | xargs -0 dirname | xargs -0 dirname)/lockfile.bwamem2. \
		--target $@ \
		--max_jobs $(MAX_CONCURRENT_ALIGNMENTS_BWAMEM2);
	# Sorting
	echo "#[INFO] Sorting BAM file for $*:"
	head -n50 $@.sam
	$(SAMTOOLS) sort $@.sam -l 1 -O BAM -o $@.tmp -T $@.SAMTOOLS_PREFIX -@ $(THREADS_SAMTOOLS)
	rm $@.sam
	# AddOrReplaceReadGroups
	if (($$($(SAMTOOLS) view $@.tmp -H | grep "^@RG" -c))); then \
		echo "#[INFO] BAM $@.tmp with read group"; \
		mv $@.tmp $@; \
	else \
		echo "#[INFO] BAM $@.tmp without read group"; \
		$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) AddOrReplaceReadGroups $(PICARD_FLAGS) -I $@.tmp O=$@ -COMPRESSION_LEVEL 1 -RGSM $(*F); \
	fi;
	-rm $@.tmp $@.RG $@.fastq_list $@.bwa.log



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# BWA2 ALIGNMENT '$(MK_RELEASE)': BWA2 generates an aligned BAM file from FASTQ file, and ask for post alignment processes 'sorting', 'realignment', 'clipping' \(if needed\) and 'recalibration'. PICARD TOOL is used to Add Or Replace Read Groups and modified the BAM header. Options: BWA2='$(BWA2)', BWAMEM2_FLAGS='$(BWAMEM2_FLAGS)', PICARD='$(PICARD)', PICARD_FLAGS='$(PICARD_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


# PIPELINES INFOS
PIPELINES_COMMENT := "ALIGNER:bwamem:BWA2 MEM - Last powerful algorithm. From FASTQ files."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

