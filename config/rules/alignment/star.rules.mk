############################
# STAR Aligner Rules
# Release: 0.9.4.8
# Date: 24/08/2022
# Author: Samuel Nicaise, Thomas Lavaux
############################
# Release
MK_RELEASE="0.9.3.1"
MK_DATE="24/08/2022"

# Release note

###################################
# STAR By Default (FROM FASTQ) #
###################################

MAX_CONCURRENT_ALIGNMENTS_STAR?=1
STAR_KEEP_RAW_BAM?=0
#STAR_FLAGS?=--outSAMtype BAM SortedByCoordinate --chimOutJunctionFormat 1 --outSAMunmapped Within --outBAMcompression 0 --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 10 --chimOutType Junctions WithinBAM --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --quantTranscriptomeBan Singleend
STAR_FLAGS?=--outSAMtype BAM SortedByCoordinate --chimOutJunctionFormat 1 --outSAMunmapped Within --outBAMcompression 0 --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 10 --chimOutType Junctions WithinBAM --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --quantTranscriptomeSAMoutput BanSingleEnd

# Raw alignemnt with STAR
%.star.star_raw.bam: %.R1$(POST_SEQUENCING).fastq.gz %.R2$(POST_SEQUENCING).fastq.gz
	# Create metrics folder
	mkdir -p $@.metrics;
	# Create read group file for STAR
	#echo "ID:1\tPL:ILLUMINA\tPU:PU\tLB:001\tSM:$(*F)" > $@.RG_STAR;
	$(PYTHON3) $(STARK_FOLDER_BIN)/functions.py launch \
		--cmd "$(STAR) --genomeDir $(GENOME_RNA).star.idx \
			--runThreadN $(THREADS_BY_SAMPLE) \
			--readFilesIn $*.R1$(POST_SEQUENCING).fastq.gz $*.R2$(POST_SEQUENCING).fastq.gz \
			--readFilesCommand zcat \
			--outFileNamePrefix $@.metrics/star. \
			--outSAMattrRGline ID:1 PL:ILLUMINA PU:PU LB:001 \"SM:$(*F)\" $(STAR_FLAGS)" \
		--lockfile_prefix $$(echo $@ | xargs -0 dirname | xargs -0 dirname)/lockfile. \
		--target $@ \
		--max_jobs $(MAX_CONCURRENT_ALIGNMENTS_STAR);
	# Rename output bam file
	if (( $(STAR_KEEP_RAW_BAM) )); then \
		cp $@.metrics/star.Aligned.sortedByCoord.out.bam $@; \
		cp $@.metrics/star.Chimeric.out.junction $*.star.junction; \
	else \
		mv $@.metrics/star.Aligned.sortedByCoord.out.bam $@; \
		cp $@.metrics/star.Chimeric.out.junction $*.star.junction; \
	fi;
	# Copy junction file from STAR alignments if exist, otherwise create empty file to avoid error in STARFusion rules
	cp $@.metrics/star.Chimeric.out.junction $*.star.junction;
	# Clean
	-rm -rf $@.Aligned.sortedByCoord.out.bam $@.RG_STAR $@._STARgenome $@._STARpass1 $@.rg_args;

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