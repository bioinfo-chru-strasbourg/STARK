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

## FASTQ from ILLUMINA ##



#%.archive.cram: %.bams.list %.unaligned.bam %.genome
#%.archive.bam: %.bams.list %.unaligned.bam %.genome
#	#echo "BAMS_LIST: $< $$(cat $<)"
#	if (($$(head -n 1 $< | wc -l))); then \
#		echo "# BAM file archived : $$(basename $$(head -n 1 $< | wc -l)))" ; \
#		$(SAMTOOLS) sort $$(head -n 1 $<) -@ $(THREADS_BY_SAMPLE) | $(SAMTOOLS) view -o $@ -O BAM -S -T `cat $*.genome` - -@ $(THREADS_BY_SAMPLE); \
#	else \
#		echo "# BAM file archived : $$(basename $*.unaligned.bam)" ; \
#		$(SAMTOOLS) sort $*.unaligned.bam -@ $(THREADS_BY_SAMPLE) | $(SAMTOOLS) view -o $@ -O BAM -S -T `cat $*.genome` - -@ $(THREADS_BY_SAMPLE); \
#	fi;

%.archive.cram: %.bams.list %.genome %.R1.fastq.gz %.R2.fastq.gz
	# Archive aligned BAM only if all original reads present. otherwise, FASTQ compressed into uBAM is archived
	#echo "BAMS_LIST: $< $$(cat $<)"
	if (($$(cat $< | wc -l))); then \
		for b in $$(cat $< | grep "$*"); do \
			echo "# BAM file to archived ? $$b" ; \
			if [ ! -s $@ ]; then \
				if [ "$$(zcat $*.R1.fastq.gz $*.R2.fastq.gz | echo $$((`wc -l`/4)) )" == "$$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $$b)" ]; then \
					echo "# BAM file archived : $$b" ; \
					#$(SAMTOOLS) sort -@ $(THREADS_BY_SAMPLE) -l $(BAM_COMPRESSION) -m $(MEMORY)G -T $@.SAMTOOLS_PREFIX -o $@ -O CRAM $$b ; \
					$(SAMTOOLS) sort $$b -@ $(THREADS_BY_SAMPLE) -T $@.SAMTOOLS_PREFIX | $(SAMTOOLS) view -o $@ -O CRAM -S -T `cat $*.genome` - -@ $(THREADS_BY_SAMPLE) -h; \
				fi; \
			fi; \
		done; \
	fi;
	if [ ! -s $@ ] && (($$(zcat $*.R1.fastq.gz $*.R2.fastq.gz | head -n 1 | wc -l))); then \
	#if ((1)); then \
		# FASTQ to BAM \
		#if [ -s $*.R2.fastq.gz ]; then \
		if (($$(zcat $*.R2.fastq.gz | head -n 1 | wc -l))); then \
			echo "PAIRED-END" ; \
			$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) FastqToSam $(PICARD_UNALIGNED_FLAGS) $(PICARD_UNALIGNED_NAME_FLAGS) FASTQ=$*.R1.fastq.gz FASTQ2=$*.R2.fastq.gz OUTPUT=$@.tmp SAMPLE_NAME=$(*F); \
		else \
			echo "SIGLE-END (NO reads in $*.R2.fastq.gz)" ; \
			$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) FastqToSam $(PICARD_UNALIGNED_FLAGS) $(PICARD_UNALIGNED_NAME_FLAGS) FASTQ=$*.R1.fastq.gz OUTPUT=$@.tmp  SAMPLE_NAME=$(*F); \
		fi; \
		# Fix Mate Information BAM $@.tmp \
		$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) FixMateInformation  $(PICARD_UNALIGNED_FLAGS) INPUT=$@.tmp ASSUME_SORTED=true VALIDATION_STRINGENCY=STRICT ; \
		# BAM Sorting and Compression \
		$(SAMTOOLS) sort $@.tmp -@ $(THREADS_BY_SAMPLE) -T $@.SAMTOOLS_PREFIX | $(SAMTOOLS) view -o $@ -O CRAM -S -T `cat $*.genome` - -@ $(THREADS_BY_SAMPLE) -h; \
		# Validation BAM $@.tmp \
		$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) ValidateSamFile $(PICARD_UNALIGNED_FLAGS)  VALIDATE_INDEX=true IGNORE_WARNINGS=true INDEX_VALIDATION_STRINGENCY=EXHAUSTIVE I=$@ > $@.validation; \
		if [ $$(grep "^ERROR" $@.validation -c) -gt 0 ]; then \
			echo "[ERROR] Input file error. Generated Archive file '$@' malformed!"; \
			exit 0; \
		fi;  \
		rm $@.tmp* $@.validation; \
	fi;
	if  [ ! -s $@ ]; then \
		echo "[ERROR] Error in generation of Archive file '$@'!"; \
		exit 1; \
	fi;
	

%.archiveOLD2.cram: %.bams.list %.unaligned.bam %.genome
	# Archive aligned BAM only if all original reads present. otherwise, unaligned BAM is archived
	#echo "BAMS_LIST: $< $$(cat $<)"
	if (($$(cat $< | wc -l))); then \
		for b in $$(cat $<); do \
			echo "# BAM file to archived ? $$b" ; \
			if [ ! -s $@ ]; then \
				if [ "$$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $*.unaligned.bam)" == "$$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $$b)" ]; then \
					echo "# BAM file archived : $$b" ; \
					#$(SAMTOOLS) sort -@ $(THREADS_BY_SAMPLE) -l $(BAM_COMPRESSION) -m $(MEMORY)G -T $@.SAMTOOLS_PREFIX -o $@ -O CRAM $$b ; \
					$(SAMTOOLS) sort $$b -@ $(THREADS_BY_SAMPLE) -T $@.SAMTOOLS_PREFIX | $(SAMTOOLS) view -o $@ -O CRAM -S -T `cat $*.genome` - -@ $(THREADS_BY_SAMPLE); \
				fi; \
			fi; \
		done; \
	else \
		echo "# BAM file archived : $*.unaligned.bam" ; \
		$(SAMTOOLS) sort $*.unaligned.bam -@ $(THREADS_BY_SAMPLE) -T $@.SAMTOOLS_PREFIX | $(SAMTOOLS) view -o $@ -O CRAM -S -T `cat $*.genome` - -@ $(THREADS_BY_SAMPLE); \
		#$(SAMTOOLS) sort $*.unaligned.bam -@ $(THREADS_BY_SAMPLE) | $(SAMTOOLS) view -o $@ -O CRAM -S -T `cat $*.genome` - -@ $(THREADS_BY_SAMPLE); \
	fi;
# -l $(BAM_COMPRESSION)

	
%.archiveOLD.cram: %.bams.list %.unaligned.bam %.genome
	# Archive aligned BAM only if all original reads present. otherwise, unaligned BAM is archived
	#echo "BAMS_LIST: $< $$(cat $<)"
	if (($$(head -n 1 $< | wc -l))); then \
		echo "# BAM file archived : $$(basename $$(head -n 1 $< | wc -l)))" ; \
		$(SAMTOOLS) sort $$(head -n 1 $<) -@ $(THREADS_BY_SAMPLE) | $(SAMTOOLS) view -o $@ -O CRAM -S -T `cat $*.genome` - -@ $(THREADS_BY_SAMPLE); \
	else \
		echo "# BAM file archived : $$(basename $*.unaligned.bam)" ; \
		$(SAMTOOLS) sort $*.unaligned.bam -@ $(THREADS_BY_SAMPLE) | $(SAMTOOLS) view -o $@ -O CRAM -S -T `cat $*.genome` - -@ $(THREADS_BY_SAMPLE); \
	fi;


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CRAM '$(MK_RELEASE)': CRAM tool generate *.archive.cram files."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )
