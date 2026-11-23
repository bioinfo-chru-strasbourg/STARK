############################
# GATK Realignment Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE=1.0.0
MK_DATE=12/11/2026

# Release note
# 10/03/2015-0.9.0: change genome reference location, in the file %.genome
# 31/10/2025-0.9.1: Reduce IndelRealigner by filtering intervals on chromosomes
# 18/09/2026-0.9.2: Fix chromosomes not in design/intervals
# 12/11/2026-1.0.0: Add GATK_INDELREALIGNER_FLAGS and GATK_REALIGNER_TARGET_CREATOR_FLAGS as configurable options


####################
# GATK Realignment #
####################

# JAVA MEMORY REALIGNMENT
JAVA_MEMORY_REALIGNMENT_MAX?=4
JAVA_MEMORY_REALIGNMENT?=$(shell echo " if ($(JAVA_MEMORY_REALIGNMENT_MAX)>$(JAVA_MEMORY)) ($(JAVA_MEMORY_REALIGNMENT_MAX)) else ($(JAVA_MEMORY))" | bc)
JAVA_FLAGS_REALIGNMENT=-Xmx$(JAVA_MEMORY_REALIGNMENT)g $(JAVA_FLAGS_OTHER_PARAM) $(JAVA_FLAGS_TMP_FOLDER)

# Flags RealignerTargetCreator
GATK_REALIGNERTARGETCREATOR_FLAGS?= \
	$(GATK_REALIGNMENT_KNOWN_OPTIONS) \
	-allowPotentiallyMisencodedQuals

# Flags IndelRealigner
GATK_INDELREALIGNER_FLAGS?= \
	$(GATK_REALIGNMENT_KNOWN_OPTIONS) \
	--interval_padding $(INTERVAL_PADDING) \
	--LODThresholdForCleaning 2.0 \
	-compress 1 \
	--maxReadsForRealignment 50000 \
	--maxReadsForConsensuses 120 \
	--maxReadsInMemory 2000000 \
	--maxConsensuses 30 \
	-model USE_READS \
	-allowPotentiallyMisencodedQuals \
	-dfrac 1 \
	--filter_reads_with_N_cigar


# Rule
%.bam: %.realignment.bam %.realignment.bam.bai %.realignment.design.bed
	# RealignerTargetCreator 
	$(JAVA8) $(JAVA_FLAGS_REALIGNMENT) -jar $(GATK3) \
			-T RealignerTargetCreator \
			-R $(GENOME) \
			-I $< \
			-o $*.for_realignment.RealignerTargetCreator.intervals \
			$(GATK_REALIGNERTARGETCREATOR_FLAGS) \
			-nt $(THREADS_RTC) \
			$$(if (($$(grep ^ -c $*.realignment.design.bed))); then echo "-L $*.realignment.design.bed"; fi;);
	# IF READS
	rm -f $*.realignment*.mk;
	+if (($$($(SAMTOOLS) idxstats $< | awk '{SUM+=$$3+$$4} END {print SUM}'))); then \
		echo "$*.for_realignment.unmapped.bam: $*.realignment.bam" >> $*.realignment1.mk; \
		echo "	$(SAMTOOLS) view --output-fmt-option level=1 -b $*.realignment.bam '*' > $*.for_realignment.unmapped.bam;" >> $*.realignment1.mk; \
		echo -n " $*.for_realignment.unmapped.bam " > $*.realignment2.mk; \
		for chr in $$($(SAMTOOLS) idxstats $< | grep -v "\*" | awk '$$3+$$4>0 {print $$1, $$3+$$4}' | sort -k2,2nr | awk '{print $$1}'); do \
			grep "^$$chr:" $*.for_realignment.RealignerTargetCreator.intervals > $*.for_realignment.RealignerTargetCreator.$$chr.intervals; \
			if [ -s $*.for_realignment.RealignerTargetCreator.$$chr.intervals ]; then \
				echo "$*.for_realignment.$$chr.bam: $*.realignment.bam" >> $*.realignment1.mk; \
				echo "	$(JAVA8) $(JAVA_FLAGS_REALIGNMENT) -Dsamjdk.compression_level=1 -jar $(GATK3) $(GATK_INDELREALIGNER_FLAGS) --analysis_type IndelRealigner --reference_sequence $(GENOME) --input_file $*.realignment.bam --out $*.for_realignment.$$chr.bam --targetIntervals $*.for_realignment.RealignerTargetCreator.$$chr.intervals --intervals $$chr" >> $*.realignment1.mk; \
				echo -n " $*.for_realignment.$$chr.bam " >> $*.realignment2.mk; \
			else \
				echo "#[INFO] No intervals to realign on chromosome $$chr for $*: reads kept as is"; \
				echo "$*.for_realignment.$$chr.bam: $*.realignment.bam" >> $*.realignment1.mk; \
				echo "	$(SAMTOOLS) view --output-fmt-option level=1 -b $*.realignment.bam '$$chr' > $*.for_realignment.$$chr.bam" >> $*.realignment1.mk; \
				echo -n " $*.for_realignment.$$chr.bam " >> $*.realignment2.mk; \
			fi; \
		done; \
		echo -n "$@: " | cat - $*.realignment2.mk > $*.realignment3.mk; \
		echo ""  >> $*.realignment3.mk; \
		echo "	$(SAMTOOLS) merge --output-fmt-option level=1 -f $@ $$(cat $*.realignment2.mk) -@ $(THREADS_BY_ALIGNER)" >> $*.realignment3.mk; \
		cat $*.realignment1.mk $*.realignment3.mk >> $*.realignment.mk; \
		make -f $*.realignment.mk $@; \
	else \
		cp $< $@; \
	fi;
	# clean
	-rm -f $*.realignment.bam $*.realignment.bam.bai $*.realignment*.mk $*.for_realignment.*


# Documentation
$(file >>$(RELEASE_INFOS),## POST_ALIGNMENT realignment '$(MK_RELEASE)': Local realignment of reads in BAM with GATK IndelRealigner. Parameters GATK_REALIGNERTARGETCREATOR_FLAGS='$(GATK_REALIGNERTARGETCREATOR_FLAGS)', GATK_INDELREALIGNER_FLAGS='$(GATK_INDELREALIGNER_FLAGS)')
$(file >>$(PIPELINES_INFOS),POST_ALIGNMENT:realignment:Local realignment of reads with GATK IndelRealigner:GATK_REALIGNERTARGETCREATOR_FLAGS='$(GATK_REALIGNERTARGETCREATOR_FLAGS)', GATK_INDELREALIGNER_FLAGS='$(GATK_INDELREALIGNER_FLAGS)')
