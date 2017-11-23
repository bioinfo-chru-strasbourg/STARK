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
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T BaseRecalibrator -R `cat $*.genome` -knownSites $(VCFDBSNP) -I $< -o $@ -L $*.from_manifest.intervals -nct $(THREADS_BY_SAMPLE) -U -compress 0

%.bam: %.recalibration.bam %.genome %.recalibration.from_manifest.intervals %.recalibration.bam.bai %.recalibration.bam.grp #%.recalibration.bam.bai %.from_manifest.intervals 
	#$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T BaseRecalibrator -R $(REF) -knownSites $(VCFDBSNP) -I $< -o $<.grp
	#$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T BaseRecalibrator -R `cat $*.genome` -knownSites $(VCFDBSNP) -I $< -o $*.recalibration.bam.grp -L $*.recalibration.from_manifest.intervals -nct $(THREADS_BY_SAMPLE) -U -compress 0
	#$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T PrintReads -R $(REF) -I $< -BQSR $<.grp -o $*.bam
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T PrintReads -R `cat $*.genome` -I $< -BQSR $*.recalibration.bam.grp -o $@ -L $*.recalibration.from_manifest.intervals -nct $(THREADS_BY_SAMPLE) -U
	# clean
	#-rm $<.grp
	#echo "# No recalibration!!! Targeted genes project may not have enough data..."
	#mv $< $@
	# CHECK NUMBER of READS in BAM
	#if (($(BAM_CHECK_STEPS))); then
	if (($(BAM_CHECK_STEPS))); then \
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
	# CHECK NUMBER of READS in BAM
	#if [ "$$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $<)" != "$$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $@)" ]; then echo "# ERROR in Number of reads between $< and $@ !!!"; exit 1; else echo "# Number of reads OK between $< and $@"; fi
	#-mv $<.bai $@.bai

RELEASE_COMMENT := "\#\# RECALIBRATION: GATK BaseRecalibrator and PrintReads are used to recalibrate BAM files."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_ALIGNMENT:recalibration:BaseRecalibrator of reads in BAM. Warning: step BAM destructive, i.e. remove reads"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


## Recalibration
#%.recalibrated.bam: %.bam
#	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T BaseRecalibrator $(GATK_BQSR) -R $(REF) -knownSites $(VCFDBSNP) -I $< -o $<.grp
#	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T PrintReads -R $(REF) -I $< -BQSR $<.grp -o $@

# CONFIG/RELEASE
#RELEASE_COMMENT := "\#\# RECALIBRATION: GATK BaseRecalibrator and PrintReads are used to recalibrate BAM files. Currently switched off"
#RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

#PIPELINES_COMMENT := "POST_ALIGNMENT:unrecalibrated:BaseRecalibrator of reads in BAM"
#PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )



