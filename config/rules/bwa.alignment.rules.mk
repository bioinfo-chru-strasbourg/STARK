############################
# BWA Aligner Rules
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
PICARD_FLAGS?=SORT_ORDER=coordinate RGLB=001 RGPL=ILLUMINA RGPU=PU VALIDATION_STRINGENCY=SILENT
THREADS_BWA?=$(THREADS_BY_SAMPLE)
POST_ALIGNMENT?=.unrecalibrated.unclipped.unrealigned.unsorted



###################################
# BWA-MEM By Default (FROM FASTQ) #
###################################

## BWA MEM (Last powerful algorithm, including SW, HMM...)

# Options
BWAMEM_FLAGS?= mem -C -M -t $(THREADS_BWA)
#BWAMEM_FLAGS?= mem -C -a -M -t $(THREADS_BWA)
#BWAMEM_FLAGS?= mem -M -t $(THREADS_BWA)
#-a

%.bwamem$(POST_ALIGNMENT).bam: %.R1.fastq.gz %.R2.fastq.gz %.genome #check in the code
	# Read group
	echo "@RG\tID:1\tPL:ILLUMINA\tPU:PU\tLB:001\tSM:$(*F)" > $@.RG
	if [ "`cat $@.RG`" != "" ]; then echo " -R "`cat $@.RG` > $@.RG; fi;
	# Alignment
	if (($$(zcat $*.R2.fastq.gz | head -n 1 | wc -l))); then \
		echo "BWA MEM Paired-End"; \
		echo "$(BWA) $(BWAMEM_FLAGS) `cat $@.RG` `cat $*.genome` $*.R1.fastq.gz $*.R2.fastq.gz | sed 's/[1-2]\:N\:0\:\([A-Z]*\)/BC\:Z\:\1/' | $(SAMTOOLS) view -o $@.tmp -b -1 -S -T `cat $*.genome` - -@ $(THREADS_SAMTOOLS);"; \
		$(BWA) $(BWAMEM_FLAGS) `cat $@.RG` `cat $*.genome` $*.R1.fastq.gz $*.R2.fastq.gz | sed 's/[1-2]\:N\:0\:\([A-Z]*\)/BC\:Z\:\1/' | $(SAMTOOLS) view -b -1 -S -T `cat $*.genome` - -@ $(THREADS_SAMTOOLS) | $(SAMTOOLS) sort - -l 1 -O BAM -o $@.tmp -T $@.SAMTOOLS_PREFIX -@ $(THREADS_SAMTOOLS); \
	else \
		echo "BWA MEM Single-End"; \
		$(BWA) $(BWAMEM_FLAGS) `cat $@.RG` `cat $*.genome` $*.R1.fastq.gz | sed 's/[1-2]\:N\:0\:\([A-Z]*\)/BC\:Z\:\1/' | $(SAMTOOLS) view -b -1 -S -T `cat $*.genome` - -@ $(THREADS_SAMTOOLS) | $(SAMTOOLS) sort - -l 1 -O BAM -o $@.tmp -T $@.SAMTOOLS_PREFIX -@ $(THREADS_SAMTOOLS); \
	fi;
	# AddOrReplaceReadGroups
	if (($$($(SAMTOOLS) view $@.tmp -H | grep "^@RG" -c))); then \
		echo "# BAM $@.tmp with read group"; \
		mv $@.tmp $@; \
	else \
		echo "# BAM $@.tmp without read group"; \
		$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) AddOrReplaceReadGroups $(PICARD_FLAGS) I=$@.tmp O=$@  COMPRESSION_LEVEL=1 RGSM=$(*F); \
	fi;
	-rm $@.tmp $@.RG




###########
# BWA-MEM #
###########

## BWA MEM (Last powerful algorithm, including SW, HMM...)

# Options
BWAMEM_FromUBAM_FLAGS= mem -C -a -Mp -t $(THREADS_BWA)
#-a

%.bwamem_FromUBAM$(POST_ALIGNMENT).sam: %.unaligned.bam %.genome #check in the code
	# Read group
	$(SAMTOOLS) view $< -H | grep '@RG' | head -n 1 | sed 's/\t/\\t/gi' > $@.RG
	if [ "`cat $@.RG`" != "" ]; then echo " -R "`cat $@.RG` > $@.RG; fi;
	# Alignment
	-$(SAMTOOLS) bam2fq $< | $(BWA) $(BWAMEM_FromUBAM_FLAGS) `cat $*.genome` `cat $@.RG` - > $@.tmp
	# AddOrReplaceReadGroups
	#echo "# AddOrReplaceReadGroups %.bwamem.unrecalibrated.unclipped.unrealigned.unsorted.sam"
	#$(JAVA) $(JAVA_FLAGS) -jar $(PICARDLIB)/AddOrReplaceReadGroups.jar $(PICARD_FLAGS) I=$@.tmp O=$@ RGSM=$(*F)
	if (($$($(SAMTOOLS) view $@.tmp -H | grep "^@RG" -c))); then \
		echo "# BAM $@.tmp with read group"; \
		mv $@.tmp $@; \
	else \
		echo "# BAM $@.tmp without read group"; \
		$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) AddOrReplaceReadGroups $(PICARD_FLAGS) I=$@.tmp O=$@  COMPRESSION_LEVEL=1 RGSM=$(*F); \
	fi;
	# clean
	-rm $@.tmp $@.RG


######################
# BWA-MEM FROM FASTQ #
######################

## BWA MEM (Last powerful algorithm, including SW, HMM...)

# Options
BWAMEM_FromFASTQ_FLAGS= mem -C -a -M -t $(THREADS_BWA)
#-a

%.bwamem_FromFASTQ$(POST_ALIGNMENT).bam: %.R1.fastq.gz %.R2.fastq.gz %.genome #check in the code
	# Read group
	#$(SAMTOOLS) view $< -H | grep '@RG' | head -n 1 | sed 's/\t/\\t/gi' > $@.RG
	echo "@RG\tID:1\tPL:ILLUMINA\tPU:PU\tLB:001\tSM:$(*F)" > $@.RG
	if [ "`cat $@.RG`" != "" ]; then echo " -R "`cat $@.RG` > $@.RG; fi;
	# Alignment
	if (($$(zcat $*.R2.fastq.gz | head -n 1 | wc -l))); then \
		echo "BWA MEM Paired-End"; \
		#echo "(BWA) $(BWAMEM_FromFASTQ_FLAGS) -M `cat $@.RG` `cat $*.genome` $*.R1.fastq.gz $*.R2.fastq.gz | $(SAMTOOLS) view -o $@.tmp -b -1 -S -T `cat $*.genome` - -@ $(THREADS_SAMTOOLS)"; exit 0; \
		$(BWA) $(BWAMEM_FromFASTQ_FLAGS) `cat $@.RG` `cat $*.genome` $*.R1.fastq.gz $*.R2.fastq.gz | sed 's/[1-2]\:N\:0\:\([A-Z]*\)/BC\:Z\:\1/' | $(SAMTOOLS) view -o $@.tmp -b -1 -S -T `cat $*.genome` - -@ $(THREADS_SAMTOOLS); \
	else \
		echo "BWA MEM Single-End"; \
		$(BWA) $(BWAMEM_FromFASTQ_FLAGS) `cat $@.RG` `cat $*.genome` $*.R1.fastq.gz | sed 's/[1-2]\:N\:0\:\([A-Z]*\)/BC\:Z\:\1/' | $(SAMTOOLS) view -o $@.tmp -b -1 -S -T `cat $*.genome` - -@ $(THREADS_SAMTOOLS); \
	fi;
	# AddOrReplaceReadGroups
	#echo "# AddOrReplaceReadGroups %.bwamem.unrecalibrated.unclipped.unrealigned.unsorted.sam"
	#$(JAVA) $(JAVA_FLAGS) -jar $(PICARDLIB)/AddOrReplaceReadGroups.jar $(PICARD_FLAGS) I=$@.tmp O=$@ RGSM=$(*F)
	if (($$($(SAMTOOLS) view $@.tmp -H | grep "^@RG" -c))); then \
		echo "# BAM $@.tmp with read group"; \
		mv $@.tmp $@; \
	else \
		echo "# BAM $@.tmp without read group"; \
		$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) AddOrReplaceReadGroups $(PICARD_FLAGS) I=$@.tmp O=$@  COMPRESSION_LEVEL=1 RGSM=$(*F); \
	fi;
	#mv $@.tmp $@
	# CHECK NUMBER of READS in BAM
	if ((0)); then \
		#echo $$($(SAMTOOLS) view -c -@ $(THREADS_SAMTOOLS) -F 0x0100 $<)" READS for $<"; \
		#echo $$($(SAMTOOLS) view -c -@ $(THREADS_SAMTOOLS) -F 0x0100 $@)" READS for $@"; \
		if [ "$$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $<)" != "$$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $@)" ]; then \
			echo "# ERROR in Number of reads between $< and $@ !!!"; \
			echo "# BCFTOOLS STATS for $<";  \
			$(SAMTOOLS) index $<; \
			$(SAMTOOLS) stats $< | grep ^SN; \
			$(SAMTOOLS) idxstats $<; \
			echo "# BCFTOOLS STATS for $@";  \
			$(SAMTOOLS) index $@; \
			$(SAMTOOLS) stats $@ | grep SN; \
			$(SAMTOOLS) idxstats $@; \
			exit 1; \
		else \
			echo "# Number of reads OK between $< and $@"; \
		fi; \
	fi;
	-rm $@.tmp $@.RG




###########
# BWA-ALN #
###########

## BWA ALN (First BWA algorithm)

# OPTIONS
BWASAMPE_FLAGS=sampe -r '@RG\tID:$(*F)\tSM:$(*F)\tPL:Illumina' -a 600
BWAALN_FLAGS=aln -t $(THREADS_BWA)

%.bwaaln$(POST_ALIGNMENT).sam: %.unaligned.bam %.genome
	# Alignment
	$(BWA) $(BWAALN_FLAGS) `cat $*.genome` -b1 $*.unaligned.bam > $*.unaligned.bam.R1.sai
	$(BWA) $(BWAALN_FLAGS) `cat $*.genome` -b2 $*.unaligned.bam > $*.unaligned.bam.R2.sai
	$(BWA) $(BWASAMPE_FLAGS) `cat $*.genome` $*.unaligned.bam.R1.sai $*.unaligned.bam.R2.sai $*.unaligned.bam $*.unaligned.bam > $@
	# Clean
	-rm -Rf $*.unaligned.bam.R1.sai $*.unaligned.bam.R2.sai


##########
# BWA-SW #
##########

## BWA SW (Smith Watermann algorithm)

# Options
BWASW_FLAGS=bwasw -t $(THREADS_BWA)

%.bwasw$(POST_ALIGNMENT).sam: %.unaligned.bam %.genome #%.unaligned.bam #check in the code
	# Read group
	$(SAMTOOLS) view $< -H | grep '@RG' | head -n 1 | sed 's/\t/\\t/gi' > $@.RG
	if [ "`cat $@.RG`" != "" ]; then echo " -R "`cat $@.RG` > $@.RG; fi;
	# Alignment
	-$(SAMTOOLS) bam2fq $*.unaligned.bam | $(BWA) $(BWASW_FLAGS) `cat $*.genome` `cat $@.RG` - > $@.tmp
	#$(JAVA) $(JAVA_FLAGS) -jar $(PICARDLIB)/AddOrReplaceReadGroups.jar $(PICARD_FLAGS) I=$@.tmp O=$@ RGSM=$(*F)
	$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) AddOrReplaceReadGroups $(PICARD_FLAGS) I=$@.tmp O=$@  COMPRESSION_LEVEL=1 RGSM=$(*F)
	# Clean
	-rm $@.tmp $@.RG




## UNCLIPPED
# Skip the clipping step. Need a sorting
%.bwaalnUnclipped.unrecalibrated.unsorted.bam: %.bwaaln.unrecalibrated.clipped.bam
	cp $< $@
%.bwaswUnclipped.unrecalibrated.unsorted.bam: %.bwasw.unrecalibrated.clipped.bam
	cp $< $@
%.bwamemUnclipped.unrecalibrated.unsorted.bam: %.bwamem.unrecalibrated.clipped.bam
	cp $< $@

## CLEAN
#CLEAN += *.unrecalibrated.unclipped.bam

# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# ALIGNMENT '$(MK_RELEASE)': BWA generates an aligned BAM file from a unaligned BAM file, and ask for post alignment processes 'sorting', 'realignment', 'clipping' \(if needed\) and 'recalibration'. Three options a available: 'ALN', 'SW' and 'MEM', generating files '*.bwaaln.bam', '*.bwasw.bam' or '*.bwamem.bam', respectively. A special Unclipped BAM file is also available for each options, generating files such as '*.bwaalnUnclipped.bam'. PICARD TOOL is used to Add Or Replace Read Groups and modified the BAM header. Options: BWA='$(BWA)', BWA_VERSION='$(BWA_VERSION)', BWASAMPE_FLAGS='$(BWASAMPE_FLAGS)', BWAALN_FLAGS='$(BWAALN_FLAGS)', BWASW_FLAGS='$(BWASW_FLAGS)', BWAMEM_FLAGS='$(BWAMEM_FLAGS)', BWAMEM_FLAGS='$(BWAMEM_FromFASTQ_FLAGS)', PICARD='$(PICARD)', PICARD_VERSION='$(PICARD_VERSION)', PICARD_FLAGS='$(PICARD_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


# PIPELINES INFOS
PIPELINES_COMMENT := "ALIGNER:bwamem:BWA MEM - Last powerful algorithm. From FASTQ files."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
PIPELINES_COMMENT := "ALIGNER:bwamem_FromUBAM:BWA MEM - Last powerful algorithm from unaligned BAM files"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
PIPELINES_COMMENT := "ALIGNER:bwamem_FromFASTQ:BWA MEM - Last powerful algorithm from FASTQ files."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
PIPELINES_COMMENT := "ALIGNER:bwaaln:BWA ALN - First BWA algorithm"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
PIPELINES_COMMENT := "ALIGNER:bwasw:BWA SW - Smith Watermann algorithm"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

#PIPELINES_COMMENT := "\#\#STEP_TYPE	STEP_NAME	STEP_DESCRIPTION2"
