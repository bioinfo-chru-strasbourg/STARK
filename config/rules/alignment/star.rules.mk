############################
# STAR Aligner Rules
# Release: 0.9.4.8
# Date: 24/08/2022
# Author: Samuel Nicaise, Thomas Lavaux
############################
# Release
MK_RELEASE="1.0.0"
MK_DATE="19/03/2026"

# Release note
# 0.9.3.1-24/08/2022: First rules for STAR aligner, with memory management for STAR processes. The memory management is based on a lockfile system to avoid multiple STAR processes to run at the same time if not enough memory is available. The lockfile is created in the parent folder of the target (e.g. postalignment) to avoid multiple STAR processes to run at the same time if not enough memory is available.
# 1.0.0-19/03/2026: Updated STAR aligner rules with improved memory management and lockfile system, and raw bam file handling.


###################################
# STAR By Default (FROM FASTQ) #
###################################

# Memory management
# Use a minimum of 36 Go per STAR process
# The calculation of MAX_CONCURRENT_ALIGNMENTS_STAR is done to avoid overloading the system memory
STAR_MIN_MEM?=36
MAX_CONCURRENT_ALIGNMENTS_STAR?=$(shell if [ $(shell echo "$(MEMTOTAL_IN_GO)/$(STAR_MIN_MEM)/$(NB_ALIGNERS)" | bc) -lt 1 ]; then echo 1; else echo "$(MEMTOTAL_IN_GO)/$(STAR_MIN_MEM)/$(NB_ALIGNERS)" | bc; fi)

# Threads per STAR process
THREADS_STAR?=$(shell echo " if ($(MAX_CONCURRENT_ALIGNMENTS_STAR)<$(NB_SAMPLE)) ($(THREADS)/$(MAX_CONCURRENT_ALIGNMENTS_STAR)) else ($(THREADS)/$(NB_SAMPLE))" | bc)

# Option to keep raw BAM file from STAR (with splitNcigar for SNV calling). This option is useful to keep BAM without post-alignment steps (e.g. realignment, recalibration...).
STAR_KEEP_RAW_BAM?=0

# STAR Flags to set specific parameters (see STAR help)
STAR_FLAGS?=--outSAMtype BAM SortedByCoordinate --chimOutJunctionFormat 1 --outSAMunmapped Within --outBAMcompression 0 --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 10 --chimOutType Junctions WithinBAM --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --quantTranscriptomeSAMoutput BanSingleEnd

# Raw alignemnt with STAR
# File star_raw.bam is an intermediate file
%.star.star_raw.bam: %.R1$(POST_SEQUENCING).fastq.gz %.R2$(POST_SEQUENCING).fastq.gz
	# DEVEL
	echo "DEVEL FOR STAR ALIGNMENT"
	echo "STAR_MIN_MEM=$(STAR_MIN_MEM)"
	echo "MAX_CONCURRENT_ALIGNMENTS_STAR=$(MAX_CONCURRENT_ALIGNMENTS_STAR)"
	echo "THREADS_STAR=$(THREADS_STAR)"
	echo "MEMTOTAL_IN_GO=$(MEMTOTAL_IN_GO) NB_ALIGNERS=$(NB_ALIGNERS) STAR_MIN_MEM=$(STAR_MIN_MEM) THREADS=$(THREADS) NB_SAMPLE=$(NB_SAMPLE)"
	# Create metrics folder
	mkdir -p $*.star.bam.metrics/$(*F).star_raw.bam.metrics;
	# Launch STAR alignment with lockfile to manage memory and avoid multiple STAR processes to run at the same time if not enough memory is available. The lockfile is created in the parent folder of the target (e.g. postalignment) to avoid multiple STAR processes to run at the same time if not enough memory is available.
	$(PYTHON3) $(STARK_FOLDER_BIN)/functions.py launch \
		--cmd "$(STAR) --genomeDir $(GENOME_RNA).star.idx \
			--runThreadN $(THREADS_STAR) \
			--readFilesIn $*.R1$(POST_SEQUENCING).fastq.gz $*.R2$(POST_SEQUENCING).fastq.gz \
			--readFilesCommand zcat \
			--outFileNamePrefix $*.star.bam.metrics/$(*F).star_raw.bam.metrics/$(*F).star_raw. \
			--outSAMattrRGline ID:1 PL:ILLUMINA PU:PU LB:001 \"SM:$(*F)\" $(STAR_FLAGS)" \
		--lockfile_prefix $$(echo $@ | xargs -0 dirname | xargs -0 dirname)/lockfile. \
		--target $@ \
		--max_jobs $(MAX_CONCURRENT_ALIGNMENTS_STAR);
	# Catch output bam file
	if (( $(STAR_KEEP_RAW_BAM) )); then \
		cp $*.star.bam.metrics/$(*F).star_raw.bam.metrics/$(*F).star_raw.Aligned.sortedByCoord.out.bam $@; \
		cp $*.star.bam.metrics/$(*F).star_raw.bam.metrics/$(*F).star_raw.Chimeric.out.junction $*.star.junction; \
	else \
		mv $*.star.bam.metrics/$(*F).star_raw.bam.metrics/$(*F).star_raw.Aligned.sortedByCoord.out.bam $@; \
		cp $*.star.bam.metrics/$(*F).star_raw.bam.metrics/$(*F).star_raw.Chimeric.out.junction $*.star.junction; \
		rm -rf $*.star.bam.metrics/$(*F).star_raw.bam.metrics/*.bam; \
	fi;
	# Copy junction file from STAR alignment (Fail if not exists)
	cp $*.star.bam.metrics/$(*F).star_raw.bam.metrics/$(*F).star_raw.Chimeric.out.junction $*.star.junction;
	# Clean
	-rm -rf $*.star.bam.metrics/$(*F).star_raw.bam.metrics/$(*F).star_raw._STARgenome $*.star.bam.metrics/$(*F).star_raw.bam.metrics/$(*F).star_raw._STARpass1;


# Alignement with STAR from raw alignement (with post alignment for SNV calling)
%.star$(POST_ALIGNMENT).bam: %.star.star_raw.bam
	ln -s $< $@

# POST ALIGNMENT STEPS

# DEVEL: Integrated into post alignment rules (folder postalignment)
# # Post alignment spécific for STAR: we need to use the bam with splitNcigar for SNV calling. However, splitNcigar is forbiden for fusion detection tools (Arriba and STARFusion) to work properly (see STARFusion.rules.mk and Arriba.rules.mk for explanation).
# %.bam: %.splitncigar.bam %.splitncigar.bam.bai
# 	$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) SplitNCigarReads -R $(GENOME) -I $< -O $@


# DEVEL: integrated into star alignment rules (see above)
# # Copy junction file from STAR alignments if exist, otherwise create empty file to avoid error in STARFusion rules
# %.star.junction: %.star.star_raw.bam
# 	if [ -e $@ ]; then touch $@; else cp $<.Chimeric.out.junction $@; fi;



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# STAR ALIGNMENT '$(MK_RELEASE)': STAR generates an aligned BAM file from FASTQ files."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


# PIPELINES INFOS
PIPELINES_COMMENT := "ALIGNER:STAR:STAR - Excellent for RNA-Seq data. From FASTQ files."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )