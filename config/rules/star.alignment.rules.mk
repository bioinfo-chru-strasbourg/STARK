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

%.star$(POST_ALIGNMENT).bam: %.R1$(POST_SEQUENCING).fastq.gz %.R2$(POST_SEQUENCING).fastq.gz %.genome
	STAR --genomeDir $$(cat $*.genome | dirname)/ref_genome.fa.star.idx/ --runThreadN 4 --readFilesIn $*.R1$(POST_SEQUENCING).fastq.gz $*.R2$(POST_SEQUENCING).fastq.gz --readFilesCommand zcat --outFileNamePrefix $@. --outSAMattrRGline "ID:1 PL:ILLUMINA PU:PU LB:001 SM:"$$(basename $@ | cut -f1 -d ".") --outSAMtype BAM SortedByCoordinate --chimOutJunctionFormat 1 --outSAMunmapped Within --outBAMcompression 0 --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --quantTranscriptomeBan Singleend
	mv $@.Aligned.sortedByCoord.out.bam $@

%.bam: %.splitncigar.bam %.splitncigar.bam.bai %.genome
	$(JAVA11) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) SplitNCigarReads -R $$(cat $*.genome) -I $< -O $@


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# STAR ALIGNMENT '$(MK_RELEASE)': STAR generates an aligned BAM file from FASTQ file
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


# PIPELINES INFOS
PIPELINES_COMMENT := "ALIGNER:STAR:STAR - Excellent for RNA-Seq data. From FASTQ files."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

