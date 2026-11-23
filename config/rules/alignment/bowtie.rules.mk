############################
# BOWTIE Aligner Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE=0.9.3
MK_DATE=11/11/2026

# Release note
# 25/07/20140.2b: add clipping step (".unclipped" on targets)
# 10/03/20150.9.1b: change genome reference location, in the file %.genome
# 29/09/2016-0.9.2b: Cleaning, PICARD new release picard.jar
# 11/11/2026-0.9.3: Change number of threads


##########
# BOWTIE #
##########

# Flags
BOWTIE_FLAGS?=


%.bowtie$(POST_ALIGNMENT).sam: %.R1$(POST_SEQUENCING).fastq.gz %.R2$(POST_SEQUENCING).fastq.gz
	# Uncompress FASTQ
	$(UNGZ) -p $(THREADS_BY_SAMPLE) -c $*.R1$(POST_SEQUENCING).fastq.gz > $*.R1$(POST_SEQUENCING).for_bowtie.fastq
	$(UNGZ) -p $(THREADS_BY_SAMPLE) -c $*.R2$(POST_SEQUENCING).fastq.gz > $*.R2$(POST_SEQUENCING).for_bowtie.fastq
	$(BOWTIE) -x $$($(GENOME) | sed -e 's/\.fa$$//gi') -1 $*.R1$(POST_SEQUENCING).for_bowtie.fastq -2 $*.R2$(POST_SEQUENCING).for_bowtie.fastq -S $@.aligned.sam -p $(THREADS_BY_ALIGNER) --sam-append-comment $(BOWTIE_FLAGS)
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


# Documentation
$(file >>$(RELEASE_INFOS),## ALIGNER bowtie '$(MK_RELEASE)': BOWTIE generates an aligned BAM file from FASTQ file. Parameters BOWTIE_FLAGS='$(BOWTIE_FLAGS)')
$(file >>$(PIPELINES_INFOS),ALIGNER:bowtie:BOWTIE - Bowtie alignment algorithm:BOWTIE_FLAGS='$(BOWTIE_FLAGS)')
