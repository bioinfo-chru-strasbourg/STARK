############################
# BWA Aligner Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.3.1"
MK_DATE="23/05/2021"

# Release note
# 25/07/20140.2b: add clipping step (".unclipped" on targets)
# 10/03/20150.9.1b: change genome reference location, in the file %.genome
# 29/09/2016-0.9.2b: Cleaning, PICARD new release picard.jar
# 13/04/2021-0.9.3.0: Cleaning, removing old BWA alignment release
# 23/05/2021-0.9.3.1: Remove samtools view step



###################################
# BWA-MEM By Default (FROM FASTQ) #
###################################

## BWA MEM (Last powerful algorithm, including SW, HMM...)

# Options
BWAMEM_FLAGS?= mem -C -M -t $(THREADS_BWA)

%.bwamem$(POST_ALIGNMENT).bam: %.R1$(POST_SEQUENCING).fastq.gz %.R2$(POST_SEQUENCING).fastq.gz %.genome #check in the code
	# Read group
	echo "@RG\tID:1\tPL:ILLUMINA\tPU:PU\tLB:001\tSM:$(*F)" > $@.RG
	if [ "`cat $@.RG`" != "" ]; then echo " -R "`cat $@.RG` > $@.RG; fi;
	# Alignment
	if (($$(zcat $*.R2.fastq.gz | head -n 1 | wc -l))); then \
		echo "BWA MEM Paired-End"; \
		#$(BWA) $(BWAMEM_FLAGS) $$(cat $@.RG) $$(cat $*.genome) $*.R1$(POST_SEQUENCING).fastq.gz $*.R2$(POST_SEQUENCING).fastq.gz | sed 's/[1-2]\:N\:0\:\([A-Z]*\)/BC\:Z\:\1/' | $(SAMTOOLS) view -b -1 -S -T $$(cat $*.genome) - -@ $(THREADS_SAMTOOLS) | $(SAMTOOLS) sort - -l 1 -O BAM -o $@.tmp -T $@.SAMTOOLS_PREFIX -@ $(THREADS_SAMTOOLS); \
		#$(BWA) $(BWAMEM_FLAGS) $$(cat $@.RG) $$(cat $*.genome) $*.R1$(POST_SEQUENCING).fastq.gz $*.R2$(POST_SEQUENCING).fastq.gz | sed 's/[1-9]*\:N\:0\:\([A-Z]*\)/BC\:Z\:\1/' | $(SAMTOOLS) sort - -l 1 -O BAM -o $@.tmp -T $@.SAMTOOLS_PREFIX -@ $(THREADS_SAMTOOLS); \
		$(BWA) $(BWAMEM_FLAGS) $$(cat $@.RG) $$(cat $*.genome) $*.R1$(POST_SEQUENCING).fastq.gz $*.R2$(POST_SEQUENCING).fastq.gz | $(SAMTOOLS) sort - -l 1 -O BAM -o $@.tmp -T $@.SAMTOOLS_PREFIX -@ $(THREADS_SAMTOOLS); \
	else \
		echo "BWA MEM Single-End"; \
		#$(BWA) $(BWAMEM_FLAGS) $$(cat $@.RG) $$(cat $*.genome) $*.R1$(POST_SEQUENCING).fastq.gz | sed 's/[1-2]\:N\:0\:\([A-Z]*\)/BC\:Z\:\1/' | $(SAMTOOLS) view -b -1 -S -T $$(cat $*.genome) - -@ $(THREADS_SAMTOOLS) | $(SAMTOOLS) sort - -l 1 -O BAM -o $@.tmp -T $@.SAMTOOLS_PREFIX -@ $(THREADS_SAMTOOLS); \
		#$(BWA) $(BWAMEM_FLAGS) $$(cat $@.RG) $$(cat $*.genome) $*.R1$(POST_SEQUENCING).fastq.gz | sed 's/[1-9]*\:N\:0\:\([A-Z]*\)/BC\:Z\:\1/' | $(SAMTOOLS) sort - -l 1 -O BAM -o $@.tmp -T $@.SAMTOOLS_PREFIX -@ $(THREADS_SAMTOOLS); \
		$(BWA) $(BWAMEM_FLAGS) $$(cat $@.RG) $$(cat $*.genome) $*.R1$(POST_SEQUENCING).fastq.gz | $(SAMTOOLS) sort - -l 1 -O BAM -o $@.tmp -T $@.SAMTOOLS_PREFIX -@ $(THREADS_SAMTOOLS); \
	fi;
	# AddOrReplaceReadGroups
	if (($$($(SAMTOOLS) view $@.tmp -H | grep "^@RG" -c))); then \
		echo "# BAM $@.tmp with read group"; \
		mv $@.tmp $@; \
	else \
		echo "# BAM $@.tmp without read group"; \
		$(JAVA11) $(JAVA_FLAGS) -jar $(PICARD) AddOrReplaceReadGroups $(PICARD_FLAGS) -I $@.tmp O=$@ -COMPRESSION_LEVEL 1 -RGSM $(*F); \
	fi;
	-rm $@.tmp $@.RG



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# BWA ALIGNMENT '$(MK_RELEASE)': BWA generates an aligned BAM file from FASTQ file, and ask for post alignment processes 'sorting', 'realignment', 'clipping' \(if needed\) and 'recalibration'. PICARD TOOL is used to Add Or Replace Read Groups and modified the BAM header. Options: BWA='$(BWA)', BWAMEM_FLAGS='$(BWAMEM_FLAGS)', PICARD='$(PICARD)', PICARD_FLAGS='$(PICARD_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


# PIPELINES INFOS
PIPELINES_COMMENT := "ALIGNER:bwamem:BWA MEM - Last powerful algorithm. From FASTQ files."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

