############################
# GATK Recalibration Rules
# Release: 0.9.2
# Date: 29/07/2022
# Author: Antony Le Bechec
############################

# Release note
# 10/03/2015-0.9.0: Creation, BAM recalibration and variant recalibration
# 29/07/2022-0.9.2: Remove variant recalibration


# BAM RECALIBRATION

%.bam.grp: %.bam %.bam.bai %.from_manifest.interval_list
	# Generate BaseRecalibrator grp file for recalibration
	$(JAVA8) $(JAVA_FLAGS) -jar $(GATK3) -T BaseRecalibrator -R $(GENOME) -knownSites $(VCFDBSNP) -I $< -o $@ -L $*.from_manifest.interval_list -nct $(THREADS_BY_SAMPLE) -U -compress 0

%.bam: %.recalibration.bam %.recalibration.bam.bai %.recalibration.bam.grp 
	# Recalibrate BAM with BaseRecalibrator grp file
	$(JAVA8) $(JAVA_FLAGS) -jar $(GATK3) -T PrintReads -R $(GENOME) -I $< -BQSR $*.recalibration.bam.grp -o $@ -nct $(THREADS_BY_SAMPLE) -U -EOQ
	-rm -f $*.recalibration.*;



RELEASE_COMMENT := "\#\# BAM RECALIBRATION: GATK BaseRecalibrator and PrintReads are used to recalibrate BAM files."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_ALIGNMENT:recalibration:BaseRecalibrator of reads in BAM. Warning: step BAM destructive, i.e. remove reads"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

