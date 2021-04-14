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



%.bowtie$(POST_ALIGNMENT).sam: %.R1$(POST_SEQUENCING).fastq.gz %.R2$(POST_SEQUENCING).fastq.gz %.genome 
	# SAM TO FASTQ
	zcat $*.R1$(POST_SEQUENCING).fastq.gz > $*.R1$(POST_SEQUENCING).for_bowtie.fastq
	zcat $*.R2$(POST_SEQUENCING).fastq.gz > $*.R2$(POST_SEQUENCING).for_bowtie.fastq
	$(BOWTIE) -x $$(cat $*.genome | sed -e 's/\.fa$$//gi') -1 $*.R1$(POST_SEQUENCING).for_bowtie.fastq -2 $*.R2$(POST_SEQUENCING).for_bowtie.fastq -S $@.aligned.sam -p $(THREADS_BY_SAMPLE)
	-rm -f $*.R1$(POST_SEQUENCING).for_bowtie.fastq $*.R2$(POST_SEQUENCING).for_bowtie.fastq
	# AddOrReplaceReadGroups
	if (($$($(SAMTOOLS) view $@.aligned.sam -H | grep "^@RG" -c))); then \
		echo "# BAM $@.aligned.sam with read group"; \
		mv $@.aligned.sam $@; \
	else \
		echo "# BAM $@.aligned.sam without read group"; \
		$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) AddOrReplaceReadGroups $(PICARD_FLAGS) -I $@.aligned.sam -O $@ -COMPRESSION_LEVEL 1 -RGSM $(*F); \
	fi;
	-rm -f $@.aligned.sam



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# ALIGNMENT '$(MK_RELEASE)': BOWTIE generates an aligned BAM file from FASTQ file, and ask for post alignment processes 'sorting', 'realignment', 'clipping' \(if needed\) and 'recalibration'."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

# PIPELINES INFOS
PIPELINES_COMMENT := "ALIGNER:bowtie:BOWTIE - Bowtie alignment algorithm"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
