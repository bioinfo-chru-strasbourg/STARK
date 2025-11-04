############################
# GATK Realignment Rules
# Release: 0.9.1
# Date: 31/10/2025
# Author: Antony Le Bechec
############################

# Release note
# 10/03/2015-0.9.0: change genome reference location, in the file %.genome
# 31/10/2025-0.9.1: Reduce IndelRealigner by filtering intervals on chromosomes


## INTERVALS
##############

# FLAGS and Options
THREADS_RTC?=$(THREADS_BY_SAMPLE)
GATKRealignerTargetCreatorFLAGS= -nt $(THREADS_RTC) #-nt 8
GATKRealignerTargetCreatorOptions= -known $(VCFDBSNP) -allowPotentiallyMisencodedQuals
GATKIndelRealignerFLAGS=
GATKIndelRealignerOptions= -known $(VCFDBSNP) --LODThresholdForCleaning 2.0 -compress 1 --maxReadsForRealignment 50000 --maxReadsForConsensuses 120 --maxReadsInMemory 2000000 --maxConsensuses 30 -model USE_READS -allowPotentiallyMisencodedQuals -dfrac 1

%.bam: %.realignment.bam %.realignment.bam.bai %.realignment.design.bed
	# RealignerTargetCreator 
	$(JAVA8) $(JAVA_FLAGS) -jar $(GATK3) $(GATKRealignerTargetCreatorFLAGS) $(GATKRealignerTargetCreatorOptions) \
			-T RealignerTargetCreator \
			-R $(GENOME) \
			-I $< \
			-o $*.for_realignment.RealignerTargetCreator.intervals \
			$$(if (($$(grep ^ -c $*.realignment.design.bed))); then echo "-L $*.realignment.design.bed"; fi;);
	# IF READS
	rm -f $*.realignment*.mk;
	+if (($$($(SAMTOOLS) idxstats $< | awk '{SUM+=$$3+$$4} END {print SUM}'))); then \
		echo "$*.for_realignment.unmapped.bam: $*.realignment.bam" >> $*.realignment1.mk; \
		echo "	$(SAMTOOLS) view -b -f 12 $*.realignment.bam > $*.for_realignment.unmapped.bam;" >> $*.realignment1.mk; \
		echo -n " $*.for_realignment.unmapped.bam " > $*.realignment2.mk; \
		for chr in $$($(SAMTOOLS) idxstats $< | grep -v "\*" | awk '{ if ($$3+$$4>0) print $$1 }'); do \
			grep "^$$chr:" $*.for_realignment.RealignerTargetCreator.intervals > $*.for_realignment.RealignerTargetCreator.$$chr.intervals; \
			if [ -s $*.for_realignment.RealignerTargetCreator.$$chr.intervals ]; then \
				echo "$*.for_realignment.$$chr.bam: $*.realignment.bam" >> $*.realignment1.mk; \
				echo "	$(JAVA8) $(JAVA_FLAGS) -jar $(GATK3) $(GATKIndelRealignerFLAGS) --analysis_type IndelRealigner --reference_sequence $(GENOME) --input_file $*.realignment.bam --out $*.for_realignment.$$chr.bam --interval_padding $(INTERVAL_PADDING) --targetIntervals $*.for_realignment.RealignerTargetCreator.$$chr.intervals $(GATKIndelRealignerOptions) --intervals $$chr" >> $*.realignment1.mk; \
				echo -n " $*.for_realignment.$$chr.bam " >> $*.realignment2.mk; \
			else \
				echo "#[INFO] No intervals to realign on chromosome $$chr for $*:"; \
				continue; \
			fi; \
		done; \
		echo -n "$@: " | cat - $*.realignment2.mk > $*.realignment3.mk; \
		echo ""  >> $*.realignment3.mk; \
		echo "	$(SAMTOOLS) merge -f $@ $$(cat $*.realignment2.mk) -@ $(THREADS_BY_SAMPLE)" >> $*.realignment3.mk; \
		cat $*.realignment1.mk $*.realignment3.mk >> $*.realignment.mk; \
		cat $*.realignment.mk; \
		make -f $*.realignment.mk $@; \
	else \
		cp $< $@; \
	fi;
	# clean
	-rm -f $*.realignment.bam $*.realignment.bam.bai $*.realignment*.mk $*.for_realignment.*


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# LOCAL REALIGNEMNT '$(MK_RELEASE)': Using GATK IndelRealigner and SAMTOOLS, BAM file is locally realigned. Options: GATKIndelRealignerFLAGS='$(GATKIndelRealignerFLAGS)', GATKIndelRealignerOptions='$(GATKIndelRealignerOptions)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "POST_ALIGNMENT:realignment:Local Realignment of reads in BAM"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
