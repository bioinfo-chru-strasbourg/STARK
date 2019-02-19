############################
# Unaligned BAM Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.7b"
MK_DATE="23/09/2016"

# Release note
# 10/03/2015: change genome reference location, in the file %.genome
# 22/05/2015: Deal with empty FASTQ files
# 04/12/2015: Change rules for bcl2fastq demultiplexing parameters
# 25/08/2016: Change fastq file generation to include R1 and R2. Simplification of code. Add $PICARD_FLAGS with compression level 0. Test number of reads using samtools. Add compression 9 with samtools
# 23/09/2016: Change PICARD LIB to BIN due to new release and PICARD use (unique JAR file)

# TOOLS
SICKLE?=$(NGSbin)/sickle
SAMTOOLS?=$(NGSbin)/samtools

# OPTIONS
FASTQ_FILTER=0 # 0 OR 1
FASTQ_R1?=
FASTQ_R2?=

GUNZIP?=gunzip

PICARD_FLAGS?=SORT_ORDER=coordinate RGLB=001 RGPL=ILLUMINA RGPU=PU VALIDATION_STRINGENCY=SILENT
PICARD_UNALIGNED_FLAGS?=COMPRESSION_LEVEL=1 MAX_RECORDS_IN_RAM=500000
PICARD_UNALIGNED_NAME_FLAGS?=LIBRARY_NAME=001 PLATFORM=ILLUMINA PLATFORM_UNIT=PU READ_GROUP_NAME=A

THREADS_BY_SAMPLE?=1

BAM_COMPRESSION?=5

GZ?=gzip

## FASTQ from ILLUMINA ##

%.fastq.gz: %.R1.fastq.gz %.R2.fastq.gz
	# Create directory
	-mkdir -p $(@D)
	# Contatenate all fastq.gz files
	-zcat $^ | $(GZ) - --fast -f -q > $@;
	# Test FASTQ content and archive file if necessary
	# Other command with BEDTOOLS: bamToFastq [OPTIONS] -i <BAM> -fq <FASTQ> -fq2 <FASTQR2> -tags
	if [ $$(zcat $@ | head -n 1 | wc -l) -lt 1 ] && [ -s $*.archive.cram ]; then \
		$(SAMTOOLS) bam2fq $*.archive.cram > $*.fastq; \
		$(GZ) $*.fastq --fast -f -q; \
	fi;



%.R1.fastq.gz: $(NEEDED)
	# Create directory
	-mkdir -p $(@D)
	# Contatenate all fastq.gz files
	-cat $(INPUTDIR)/$$(echo $$(basename $$(dirname $(@D))))/$(*F)_S*_R1_*.fastq.gz $(INPUTDIR)/$$(echo $$(basename $$(dirname $(@D))))/*/$(*F)_S*_R1_*.fastq.gz > $@;
	# if no reads
	if [ $$(zcat $@ | head -n 1 | wc -l) -lt 1 ] && [ -s $*.unaligned.bam ]; then \
		$(SAMTOOLS) bam2fq $*.unaligned.bam -1 $*.R1.fastq -2 $*.R2_1.fastq >> $*.R0_1.fastq; \
		if [ -s $*.R0_1.fastq ]; then \
			cat $*.R2_1.fastq $*.R0_1.fastq >> $*.R1.fastq; \
			rm -f $*.R2_1.fastq $*.R0_1.fastq; \
		fi; \
		$(GZ) $*.R1.fastq --fast -f -q; \
	fi;
	# if no reads
	if [ $$(zcat $@ | head -n 1 | wc -l) -lt 1 ] && [ -s $*.archive.cram ]; then \
		$(SAMTOOLS) bam2fq $*.archive.cram -1 $*.R1.fastq -2 $*.R2_1.fastq >> $*.R0_1.fastq; \
		if [ -s $*.R0_1.fastq ]; then \
			cat $*.R2_1.fastq $*.R0_1.fastq >> $*.R1.fastq; \
			rm -f $*.R2_1.fastq $*.R0_1.fastq; \
		fi; \
		$(GZ) $*.R1.fastq --fast -f -q; \
	fi;

%.R2.fastq.gz: $(NEEDED)
	# Create directory
	-mkdir -p $(@D)
	# Contatenate all fastq.gz files
	-cat $(INPUTDIR)/$$(echo $$(basename $$(dirname $(@D))))/$(*F)_S*_R2_*.fastq.gz $(INPUTDIR)/$$(echo $$(basename $$(dirname $(@D))))/*/$(*F)_S*_R2_*.fastq.gz > $@;
	# if no reads
	if [ $$(zcat $@ | head -n 1 | wc -l) -lt 1 ] && [ -s $*.unaligned.bam ]; then \
		$(SAMTOOLS) bam2fq $*.unaligned.bam -1 $*.R1_2.fastq -2 $*.R2.fastq >> $*.R0_2.fastq; \
		if [ -s $*.R0_1.fastq ]; then \
			> $*.R2.fastq; \
			rm -f $*.R1_2.fastq $*.R0_2.fastq; \
		fi; \
		$(GZ) $*.R2.fastq --fast -f -q; \
	fi;
	# if no reads
	if [ $$(zcat $@ | head -n 1 | wc -l) -lt 1 ] && [ -s $*.archive.cram ]; then \
		$(SAMTOOLS) bam2fq $*.archive.cram -1 $*.R1_2.fastq -2 $*.R2.fastq >> $*.R0_2.fastq; \
		if [ -s $*.R0_1.fastq ]; then \
			> $*.R2.fastq; \
			rm -f $*.R1_2.fastq $*.R0_2.fastq; \
		fi; \
		$(GZ) $*.R2.fastq --fast -f -q; \
	fi;

## UNALIGNED BAM ###

%.unaligned.bam: %.R1.fastq.gz %.R2.fastq.gz 
	# Creation of the output folder $(@D)
	@mkdir -p $(@D)
	
	# Test FASTQ content and archive file if necessary
	# Other command with BEDTOOLS: bamToFastq [OPTIONS] -i <BAM> -fq <FASTQ> -fq2 <FASTQR2> -tags
	if [ $$(zcat $^ | head -n 1 | wc -l) -lt 1 ] && [ -s $*.archive.cram ]; then \
		$(SAMTOOLS) bam2fq $*.archive.cram -1 $*.R1.fastq -2 $*.R2.fastq >$*.R0.fastq; \
		if [ -s $*.R0.fastq ]; then \
			cat $*.R2.fastq $*.R0.fastq >> $*.R1.fastq; \
			> $*.R2.fastq; \
			rm -f $*.R0.fastq; \
		fi; \
		$(GZ) $*.R1.fastq --fast -f -q; \
		$(GZ) $*.R2.fastq --fast -f -q; \
	fi;
	
	# If fastq not empty
	if (($$(zcat $^ | head -n 1 | wc -l))); then \
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
		$(SAMTOOLS) sort -o $@ -l $(BAM_COMPRESSION) -@ $(THREADS_BY_SAMPLE) -T $@.SAMTOOLS.SORT $@.tmp; \
		# Validation BAM $@.tmp \
		$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) ValidateSamFile $(PICARD_UNALIGNED_FLAGS)  VALIDATE_INDEX=true IGNORE_WARNINGS=true INDEX_VALIDATION_STRINGENCY=EXHAUSTIVE I=$@ > $@.validation; \
		if [ $$(grep "^ERROR" $@.validation -c) -gt 0 ]; then \
			echo "[ERROR] Input file error. Generated uBAM file '$@' malformed!"; \
			exit 0; \
		fi;  \
		#rm $@.tmp* $@.validation; \
	else \
		echo "[ERROR] Input file error. No Fastq files and no archive files. Error in generation of uBAM file '$@'!"; \
		exit 1; \
	fi;
	
	# check if there are same number of reads in fastq and unaligned.bam files
	#if [ "$$($(SAMTOOLS) view $@ -c)" == "$$(zcat $*.fastq.R1.gz $*.fastq.R2.gz | awk 'END{print NR/4}')" ]; then \
	#	echo "# Conversion of  '$<' to '$@' OK"; \
	#else \
	#	echo "# Error : not the same number of reads in '$<' and '$@' - EXIT STARK "; \
	#	exit 1; \
	#fi;
	#if [ ! -s $@ ]; then echo "# ERROR in $@ generation"; exit 1; fi;
	
	# cleaning
	-rm $@.tmp* $@.validation



%.unaligned.OK.bam: %.R1.fastq.gz %.R2.fastq.gz #%.fastq.gz
	# Creation of the output folder $(@D)
	@mkdir -p $(@D)
	#ls -l $(@D)
	# Create FASTQ
	#$(GUNZIP) $(INPUTDIR)/`echo $$(basename $$(dirname $(@D)))`/$(*F)_S*_R*_*.fastq.gz -c > $@.fastq;
	# FASTQ to BAM
	if [ -s $*.R2.fastq.gz ]; then \
		echo "PAIRED-END" ; \
		$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) FastqToSam $(PICARD_UNALIGNED_FLAGS) $(PICARD_UNALIGNED_NAME_FLAGS) FASTQ=$*.R1.fastq.gz FASTQ2=$*.R2.fastq.gz OUTPUT=$@.tmp SAMPLE_NAME=$(*F); \
	else \
		echo "SIGLE-END (NO reads in $*.R2.fastq.gz)" ; \
		$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) FastqToSam $(PICARD_UNALIGNED_FLAGS) $(PICARD_UNALIGNED_NAME_FLAGS) FASTQ=$*.R1.fastq.gz OUTPUT=$@.tmp  SAMPLE_NAME=$(*F); \
	fi;
	# AddOrReplaceGROUP
	#$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) AddOrReplaceReadGroups $(PICARD_FLAGS) I=$@.tmp O=$@.tmp2 RGSM=$(*F) 
	# Fix Mate Information BAM $@.tmp
	#$(JAVA) $(JAVA_FLAGS) -jar $(PICARDLIB)/FixMateInformation.jar  $(PICARD_UNALIGNED_FLAGS) INPUT=$@.tmp;
	$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) FixMateInformation  $(PICARD_UNALIGNED_FLAGS) INPUT=$@.tmp ASSUME_SORTED=true VALIDATION_STRINGENCY=STRICT ;
	
	# VALIDATION_STRINGENCY=STRICT ???
	# BAM Sorting and Compression
	$(SAMTOOLS) sort -o $@ -l $(BAM_COMPRESSION) -@ $(THREADS_BY_SAMPLE) -T $@.SAMTOOLS.SORT $@.tmp;
	
	# Validation BAM $@.tmp
	#-$(JAVA) $(JAVA_FLAGS) -jar $(PICARDLIB)/ValidateSamFile.jar  $(PICARD_UNALIGNED_FLAGS) VALIDATE_INDEX=true I=$@.tmp ;
	$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) ValidateSamFile $(PICARD_UNALIGNED_FLAGS)  VALIDATE_INDEX=true IGNORE_WARNINGS=true INDEX_VALIDATION_STRINGENCY=EXHAUSTIVE I=$@ > $@.validation;
	if [ $$(grep "^ERROR" $@.validation -c) -gt 0 ]; then \
		echo "[ERROR] Input file error. Generated uBAM file '$@' malformed!"; \
		exit 0; \
	fi;
	
	# check if there are same number of reads in fastq and unaligned.bam files
	#if [ "$$($(SAMTOOLS) view $@ -c)" == "$$(zcat $*.fastq.R1.gz $*.fastq.R2.gz | awk 'END{print NR/4}')" ]; then \
	#	echo "# Conversion of  '$<' to '$@' OK"; \
	#else \
	#	echo "# Error : not the same number of reads in '$<' and '$@' - EXIT STARK "; \
	#	exit 1; \
	#fi;
	#if [ ! -s $@ ]; then echo "# ERROR in $@ generation"; exit 1; fi;
	rm $@.tmp* $@.validation


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# UNALIGNED BAM '$(MK_RELEASE)': PICARD tool generate a *.unaligned.bam file from the two FASTQ Read1 and Read2 files."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


