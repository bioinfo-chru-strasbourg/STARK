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

STAR_FLAGS?=--outSAMtype BAM SortedByCoordinate --chimOutJunctionFormat 1 --outSAMunmapped Within --outBAMcompression 0 --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 10 --chimOutType Junctions WithinBAM --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --quantTranscriptomeBan Singleend

%.star_raw.bam: %.R1$(POST_SEQUENCING).fastq.gz %.R2$(POST_SEQUENCING).fastq.gz %.genome
	echo "ID:1\tPL:ILLUMINA\tPU:PU\tLB:001\tSM:$(*F)" > $@.RG_STAR
	STAR --genomeDir $$(cat $*.genome | xargs -0 dirname)/ref_genome.fa.star.idx/ --runThreadN 4 --readFilesIn $*.R1$(POST_SEQUENCING).fastq.gz $*.R2$(POST_SEQUENCING).fastq.gz --readFilesCommand zcat --outFileNamePrefix $@. --outSAMattrRGline $$(cat $@.RG_STAR) $(STAR_FLAGS)
	# fix issue with base recalibration and rename output
	$(JAVA11) $(JAVA_FLAGS) -jar $(PICARD) AddOrReplaceReadGroups $(PICARD_FLAGS) -I $@.Aligned.sortedByCoord.out.bam -O $@ -COMPRESSION_LEVEL 1 -RGSM $(*F);
	-rm -rf $@.Aligned.sortedByCoord.out.bam $@.RG_STAR $@._STARgenome $@._STARpass1

%.bam: %.splitncigar.bam %.splitncigar.bam.bai %.genome
	$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) SplitNCigarReads -R $$(cat $*.genome) -I $< -O $@

%.star.junction: %.star_raw.bam
	mv $<.Chimeric.out.junction $@

%.star$(POST_ALIGNMENT).bam: %.star_raw.bam
	ln -s $< $@

# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# STAR ALIGNMENT '$(MK_RELEASE)': STAR generates an aligned BAM file from FASTQ file
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


# PIPELINES INFOS
PIPELINES_COMMENT := "ALIGNER:STAR:STAR - Excellent for RNA-Seq data. From FASTQ files."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

