############################
# GATK Recalibration Rules
# Release: 0.9.1
# Date: 19/01/2015
# Author: Antony Le Bechec
############################

# TOOLS
JAVA?=java

# OPTIONS
GATK?=$(NGSbin)/GenomeAnalysisTK.jar
#GATK_BQSR= -dcov 10000 -nt 2


%.bam.grp: %.bam %.bam.bai %.from_manifest.intervals %.genome
	# Generate BaseRecalibrator grp file for recalibration
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T BaseRecalibrator -R `cat $*.genome` -knownSites $(VCFDBSNP) -I $< -o $@ -L $*.from_manifest.intervals -nct $(THREADS_BY_SAMPLE) -U -compress 0

%.bam: %.recalibration.bam %.genome %.recalibration.from_manifest.intervals %.recalibration.bam.bai %.recalibration.bam.grp #%.recalibration.bam.bai %.from_manifest.intervals
	# Recalibrate BAM with BaseRecalibrator grp file
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T PrintReads -R `cat $*.genome` -I $< -BQSR $*.recalibration.bam.grp -o $@ -L $*.recalibration.from_manifest.intervals -nct $(THREADS_BY_SAMPLE) -U
	# clean
	-rm -f $<;

RELEASE_COMMENT := "\#\# RECALIBRATION: GATK BaseRecalibrator and PrintReads are used to recalibrate BAM files."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_ALIGNMENT:recalibration:BaseRecalibrator of reads in BAM. Warning: step BAM destructive, i.e. remove reads"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
