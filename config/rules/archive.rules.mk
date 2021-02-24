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


%.archive.cram: %.bams.list %.genome %.R1.fastq.gz %.R2.fastq.gz $(REF_CACHE_FOLDER) #%.R1.fastq.gz.format
	# Archive aligned BAM only if all original reads present. otherwise, FASTQ compressed into uBAM is archived
	# echo "BAMS_LIST: $< $$(cat $<)"
	# SEPARATOR (to find aligner name), FILES_TRIED and VALID_CRAM_MSG are used for metrics file
	# The first valid cram is kept.
	FILES_TRIED=""; \
	VALID_CRAM_MSG="NO CRAM could pass the validation"; \
	if (($$(cat $< | wc -l))); then \
		for b in $$(cat $< | grep "$*"); do \
			if [ ! -s $@ ]; then \
				echo "# BAM file to archived ? $$b" ; \
				FILES_TRIED="$$FILES_TRIED $$b"; \
				####Cram creation and test, from each bam in list. Cram is removed if test is failed \
				#1) CRAM creation \
				$(SAMTOOLS) view -h -@ $(THREADS_BY_SAMPLE) -O CRAM -T `cat $*.genome` $$b | $(SAMTOOLS) sort -l 9 -O CRAM -@ $(THREADS_BY_SAMPLE) -T $@.SAMTOOLS_PREFIX -o $*.archive.cram; \
				echo "# BAM file archived : $$b" ; \
				#2) Generate FASTQs from the new CRAM \
				$(SAMTOOLS) bam2fq -@ $(THREADS_BY_SAMPLE) -O --reference `cat $*.genome` -1 $*.R1.fromCram.fastq -2 $*.R2.fromCram.fastq $*.archive.cram > $*.R0.fromCram.fastq; \
				# If R0 not empty \
				[ -s $*.R0.fromCram.fastq ] && echo "[WARNING] $*.R0.fromCram.fastq is not empty"; \
				# IF Single End \
				! (( $$($(UNGZ) -c $*.R2.fastq.gz | head -n 1 | wc -l) )) && echo "[WARNING] Single End detected ($*.R2.fastq.gz empty) " && cat $*.R2.fromCram.fastq $*.R0.fromCram.fastq >> $*.R1.fromCram.fastq && > $*.R2.fromCram.fastq; \
				#3) sort original and fromCram FASTQs, md5 on each pair to check if identical content \
				if [ "`$(UNGZ) -c $*.R1.fastq.gz | cut -f1 -d " " | paste - - - - | LC_COLLATE=C sort -S 100% -T $@.SORT_PREFIX -n -k1,1 -t " " | tr "\t" "\n" | md5sum`" == "`cat $*.R1.fromCram.fastq | cut -f1 -d " " | paste - - - - | LC_COLLATE=C sort -S 100% -T $@.SORT_PREFIX -n -k1,1 -t " " | tr "\t" "\n" | md5sum`" ] && [ "`$(UNGZ) -c $*.R2.fastq.gz | cut -f1 -d " " | paste - - - - | LC_COLLATE=C sort -S 100% -T $@.SORT_PREFIX -n -k1,1 -t " " | tr "\t" "\n" | md5sum`" == "`cat $*.R2.fromCram.fastq | cut -f1 -d " " | paste - - - - | LC_COLLATE=C sort -S 100% -T $@.SORT_PREFIX -n -k1,1 -t " " | tr "\t" "\n" | md5sum`" ]; then \
					echo "$*.archive.cram contains as expected the same reads as $*.R1.fastq.gz and $*.R2.fastq.gz"; \
					VALID_CRAM_MSG="CRAM created from $$b SUCCESSFULLY passed validation: contained the same reads as $*.R1.fastq.gz, $*.R2.fastq.gz"; \
				else \
					echo "$*.archive.cram does not contain the exact same reads as $*.R1.fastq.gz and $*.R2.fastq.gz and will be removed"; \
					rm -f  $*.archive.cram; \
				fi; \
				# Cleaning fromCram \
				rm $*.R1.fromCram.fastq $*.R2.fromCram.fastq $*.R0.fromCram.fastq; \
				#### \
			fi; \
		done; \
	fi; \
	if [ ! -s $@ ] && (($$($(UNGZ) -c $*.R1.fastq.gz $*.R2.fastq.gz | head -n 1 | wc -l))); then \
		FILES_TRIED="$$FILES_TRIED $*.R1.fastq.gz $*.R2.fastq.gz"; \
	#if ((1)); then \
		# FASTQ to BAM \
		#if [ -s $*.R2.fastq.gz ]; then \
		if (($$($(UNGZ) -c $*.R2.fastq.gz | head -n 1 | wc -l))); then \
			echo "PAIRED-END" ; \
			$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) FastqToSam $(PICARD_UNALIGNED_FLAGS) $(PICARD_UNALIGNED_NAME_FLAGS) -FASTQ $*.R1.fastq.gz -FASTQ2 $*.R2.fastq.gz -OUTPUT $@.tmp -SAMPLE_NAME $(*F); \
		else \
			echo "SINGLE-END (NO reads in $*.R2.fastq.gz)" ; \
			$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) FastqToSam $(PICARD_UNALIGNED_FLAGS) $(PICARD_UNALIGNED_NAME_FLAGS) -FASTQ $*.R1.fastq.gz -OUTPUT $@.tmp -SAMPLE_NAME $(*F); \
		fi; \
		# Fix Mate Information BAM $@.tmp \
		$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) FixMateInformation  $(PICARD_UNALIGNED_FLAGS) -INPUT $@.tmp -ASSUME_SORTED true -VALIDATION_STRINGENCY STRICT ; \
		# BAM Sorting and Compression \
		$(SAMTOOLS) sort $@.tmp -@ $(THREADS_BY_SAMPLE) -T $@.SAMTOOLS_PREFIX | $(SAMTOOLS) view -o $@ -O CRAM -S -T `cat $*.genome` - -@ $(THREADS_BY_SAMPLE) -h; \
		# Validation BAM $@.tmp \
		$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) ValidateSamFile $(PICARD_UNALIGNED_FLAGS) -VALIDATE_INDEX true -IGNORE_WARNINGS true -INDEX_VALIDATION_STRINGENCY EXHAUSTIVE -I $@ > $@.validation; \
		if [ $$(grep "^ERROR" $@.validation -c) -gt 0 ]; then \
			echo "[ERROR] Input file error. Generated Archive file '$@' malformed!"; \
			exit 0; \
		fi;  \
		rm $@.tmp* $@.validation; \
		####Same cram test as before. Cram is removed if test is failed \
		$(SAMTOOLS) bam2fq -@ $(THREADS_BY_SAMPLE) --reference `cat $*.genome` -1 $*.R1.fromCram.fastq -2 $*.R2.fromCram.fastq $*.archive.cram > $*.R0.fromCram.fastq ; \
		# If R0 not empty \
		[ -s $*.R0.fromCram.fastq ] && echo "[WARNING] $*.R0.fromCram.fastq is not empty"; \
		# IF Single End \
		! (( $$($(UNGZ) -c $*.R2.fastq.gz | head -n 1 | wc -l) )) && echo "[WARNING] Single End detected ($*.R2.fastq.gz empty) " && cat $*.R2.fromCram.fastq $*.R0.fromCram.fastq >> $*.R1.fromCram.fastq && > $*.R2.fromCram.fastq; \
		# sort original and fromCram FASTQs, md5 on each pair to check if identical content \
		if [ "`$(UNGZ) -c $*.R1.fastq.gz | cut -f1 -d " " | paste - - - - | LC_COLLATE=C sort -S 100% -T $@.SORT_PREFIX -n -k1,1 -t " " | tr "\t" "\n" | md5sum`" == "`cat $*.R1.fromCram.fastq | cut -f1 -d " " | paste - - - - | LC_COLLATE=C sort -S 100% -T $@.SORT_PREFIX -n -k1,1 -t " " | tr "\t" "\n" | md5sum`" ] && [ "`$(UNGZ) -c $*.R2.fastq.gz | cut -f1 -d " " | paste - - - - | LC_COLLATE=C sort -S 100% -T $@.SORT_PREFIX -n -k1,1 -t " " | tr "\t" "\n" | md5sum`" == "`cat $*.R2.fromCram.fastq | cut -f1 -d " " | paste - - - - | LC_COLLATE=C sort -S 100% -T $@.SORT_PREFIX -n -k1,1 -t " " | tr "\t" "\n" | md5sum`" ]; then \
			echo "$*.archive.cram contains as expected the same reads as $*.R1.fastq.gz and $*.R2.fastq.gz"; \
			VALID_CRAM_MSG="CRAM created from $*.R1.fastq.gz, $*.R2.fastq.gz SUCCESSFULLY passed validation: contained the same reads as $*.R1.fastq.gz, $*.R2.fastq.gz"; \
		else \
			echo "$*.archive.cram does not contain the exact same reads as $*.R1.fastq.gz and $*.R2.fastq.gz and will be removed"; \
			rm -f  $*.archive.cram; \
		fi; \
		rm $*.R1.fromCram.fastq $*.R2.fromCram.fastq $*.R0.fromCram.fastq; \
		#### \
	fi; \
	if  [ ! -s $@ ]; then \
		echo "[ERROR] Error in generation of Archive file '$@'!"; \
		exit 1; \
	fi; \
	#### METRICS \
	mkdir -p "$*.archive.cram.metrics"; \
	echo "Tried to generate CRAM from: $$FILES_TRIED" > $*.archive.cram.metrics/metrics; \
	echo "$$VALID_CRAM_MSG" >> $*.archive.cram.metrics/metrics;

	#if (($$($(UNGZ) -c $*.R2.fastq.gz | head -n 1 | wc -l))); then


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CRAM '$(MK_RELEASE)': CRAM tool generate *.archive.cram files."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )
