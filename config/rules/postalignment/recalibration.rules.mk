############################
# GATK Recalibration Rules
# Release: 0.9.3
# Date: 31/10/2025
# Author: Antony Le Bechec
############################

# Release note
# 10/03/2015-0.9.0: Creation, BAM recalibration and variant recalibration
# 29/07/2022-0.9.2: Remove variant recalibration
# 31/10/2025-0.9.3: Change GATK3 to GATK4 for recalibration


# BAM RECALIBRATION

%.bam.grp: %.bam %.bam.bai %.from_manifest.interval_list
	# Generate BaseRecalibrator grp file for recalibration
	$(JAVA) $(JAVA_FLAGS_GATK4) -XX:ParallelGCThreads=2 -jar $(GATK4) \
		BaseRecalibrator \
		-I $< \
		-R $(GENOME) \
		$(GATK_RECALIBRATION_KNOWN_OPTIONS) \
		--use-original-qualities \
		-L $*.from_manifest.interval_list \
		-O $@


%.bam: %.recalibration.bam %.recalibration.bam.bai %.recalibration.bam.grp 
	# Recalibrate BAM with BaseRecalibrator grp file
	$(JAVA) $(JAVA_FLAGS_GATK4) -XX:ParallelGCThreads=$(THREADS_BY_SAMPLE) -jar $(GATK4) \
		ApplyBQSR \
		-R $(GENOME) \
		-I $< \
		--bqsr-recal-file $*.recalibration.bam.grp \
		--use-original-qualities \
		-O $@
	-rm -f $*.recalibration.*;



RELEASE_COMMENT := "\#\# BAM RECALIBRATION: GATK BaseRecalibrator and PrintReads are used to recalibrate BAM files."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_ALIGNMENT:recalibration:BaseRecalibrator of reads in BAM. Warning: step BAM destructive, i.e. remove reads"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

