############################
# BOWTIE Aligner Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.2b"
MK_DATE="29/09/2016"

# Release note
# 25/07/20140.2b: add clipping step (".unclipped" on targets)
# 10/03/20150.9.1b: change genome reference location, in the file %.genome
# 29/09/2016-0.9.2b: Cleaning, PICARD new release picard.jar

# TOOLS
JAVA?=java
BWA?=$(NGSbin)/bwa
PICARDLIB?=$(NGSbin)/picard-tools
# OPTIONS
JAVA_FLAGS?= -Xmx16g
PICARD_FLAGS=SORT_ORDER=coordinate RGLB=001 RGPL=Illumina RGPU=A3 VALIDATION_STRINGENCY=SILENT
THREADS_BWA?=$(THREADS_BY_SAMPLE)
POST_ALIGNMENT?=.unrecalibrated.unclipped.unrealigned.unsorted
PICARD_UNALIGNED_FLAGS?=COMPRESSION_LEVEL=1 MAX_RECORDS_IN_RAM=500000

%.bowtie$(POST_ALIGNMENT).sam: %.unaligned.bam %.genome 
	# SAM TO FASTQ
	$(JAVA) -jar $(PICARD) SamToFastq I=$< FASTQ=$@.R1.fastq SECOND_END_FASTQ=$@.R2.fastq UNPAIRED_FASTQ=$@.RU.fastq
	# ALIGNMENT
	#$(BOWTIE) -x $$(cat $*.genome | sed -e 's/\.fa$$//gi') -1 $@.R1.fastq -2 $@.R2.fastq -S $@.aligned.sam -p $(THREADS_BY_SAMPLE)
	#$(BOWTIE) -x $$(cat $*.genome | sed -e 's/\.fa$$//gi') -12 $@.R1.fastq,$@.R2.fastq,$@.RU.fastq -S $@.aligned.sam -p $(THREADS_BY_SAMPLE)
	$(BOWTIE) -x $$(cat $*.genome | sed -e 's/\.fa$$//gi') -1 $@.R1.fastq -2 $@.R2.fastq -U $@.RU.fastq -S $@.aligned.sam -p $(THREADS_BY_SAMPLE)
	-rm -f  $@.R1.fastq $@.R2.fastq $@.RU.fastq
	# UNPAIRED READS
	#$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) FastqToSam $(PICARD_UNALIGNED_FLAGS) FASTQ=$@.RU.fastq OUTPUT=$@.RU.sam SAMPLE_NAME=$(*F) PLATFORM=PL
	#-rm -f $@.RU.fastq
	# ALIGNED SAM
	#$(SAMTOOLS) view $@.RU.sam >> $@.aligned.sam
	#-rm -f $@.RU.sam
	# AddOrReplaceReadGroups
	$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) AddOrReplaceReadGroups $(PICARD_FLAGS) I=$@.aligned.sam O=$@ RGSM=$(*F)
	-rm -f $@.aligned.sam
	
	


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# ALIGNMENT '$(MK_RELEASE)': BOWTIE generates an aligned BAM file from a unaligned BAM file, and ask for post alignment processes 'sorting', 'realignment', 'clipping' \(if needed\) and 'recalibration'."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

# PIPELINES INFOS
PIPELINES_COMMENT := "ALIGNER:bowtie:BOWTIE - Bowtie alignment algorithm"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

#PIPELINES_COMMENT := "\#\#STEP_TYPE	STEP_NAME	STEP_DESCRIPTION2"



