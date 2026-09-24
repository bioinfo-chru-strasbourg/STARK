############################
# GATK Recalibration Rules
# Release: 0.9.4
# Date: 18/09/2026
# Author: Antony Le Bechec
############################

# Release note
# 10/03/2015-0.9.0: Creation, BAM recalibration and variant recalibration
# 29/07/2022-0.9.2: Remove variant recalibration
# 31/10/2025-0.9.3: Change GATK3 to GATK4 for recalibration
# 18/09/2026-0.9.4: Fix option --emit-original-quals for GATK4


# BAM RECALIBRATION

# BQSR calculation with gatk in one line
%.bam.grp: %.bam %.bam.bai %.from_manifest.interval_list
	# Generate BaseRecalibrator grp file for recalibration
	$(JAVA) $(JAVA_FLAGS_GATK4) -XX:ParallelGCThreads=$(THREADS_BY_ALIGNER) -jar $(GATK4) \
		BaseRecalibrator \
		-I $< \
		-R $(GENOME) \
		$(GATK_RECALIBRATION_KNOWN_OPTIONS) \
		--use-original-qualities \
		-L $*.from_manifest.interval_list \
		-O $@


# Recalibration with gatk using ApplyBQSR with BQSR, through parallelization
%.bam: %.recalibration.bam %.recalibration.bam.bai %.recalibration.bam.grp
	# If both the recalibration group file and the BAM file have valid content, proceed with recalibration
	rm -f $*.recalibration*.mk;
	#+if (($$($(SAMTOOLS) idxstats $< | awk '{SUM+=$$3+$$4} END {print SUM}'))); then
	#+if ((grp_valid && reads_valid)); then
	+if awk -f $(STARK_FOLDER_BIN)/grp_valid.awk "$*.recalibration.bam.grp" && (( $$( $(SAMTOOLS) idxstats "$<" | awk '{SUM += $$3 + $$4} END {print SUM + 0}') > 0 )); then \
		echo "$*.for_recalibration.unmapped.bam: $*.recalibration.bam" >> $*.recalibration1.mk; \
		echo "	$(SAMTOOLS) view --output-fmt-option level=1 -b $*.recalibration.bam '*' > $*.for_recalibration.unmapped.bam;" >> $*.recalibration1.mk; \
		echo -n " $*.for_recalibration.unmapped.bam " > $*.recalibration2.mk; \
		cat $*.recalibration.bam.grp > /tmp/recalibration.bam.grp; \
		for chr in $$($(SAMTOOLS) idxstats $< | grep -v "\*" | awk '$$3+$$4>0 {print $$1, $$3+$$4}' | sort -k2,2nr | awk '{print $$1}'); do \
			echo "$*.for_recalibration.$$chr.bam: $*.recalibration.bam" >> $*.recalibration1.mk; \
			echo "	$(JAVA) $(JAVA_FLAGS_GATK4) -Dsamjdk.compression_level=1 -jar $(GATK4) ApplyBQSR -R $(GENOME) -I $*.recalibration.bam --bqsr-recal-file $*.recalibration.bam.grp --use-original-qualities --emit-original-quals -L $$chr -O $*.for_recalibration.$$chr.bam" >> $*.recalibration1.mk; \
			echo -n " $*.for_recalibration.$$chr.bam " >> $*.recalibration2.mk; \
		done; \
		echo -n "$@: " | cat - $*.recalibration2.mk > $*.recalibration3.mk; \
		echo ""  >> $*.recalibration3.mk; \
		echo "	$(SAMTOOLS) merge --output-fmt-option level=1 -f $@ $$(cat $*.recalibration2.mk) -@ $(THREADS_BY_ALIGNER)" >> $*.recalibration3.mk; \
		cat $*.recalibration1.mk $*.recalibration3.mk >> $*.recalibration.mk; \
		make -f $*.recalibration.mk $@; \
	else \
		cp $< $@; \
	fi;
	# clean
	-rm -f $*.recalibration.bam $*.recalibration.bam.bai $*.recalibration*.mk $*.for_recalibration.*


# BQSR calculation with gatk and splitting interval lists for parallelization
# As multiple JVM are used, performances are not optimal, and can cause issues (e.g. with gathering)
# %.bam.grp: %.bam %.bam.bai %.from_manifest.interval_list
# 	# Scatter the interval list
# 	rm -rf $*.recalibration.scatter $*.recalibration.mk
# 	mkdir -p $*.recalibration.scatter
# 	$(JAVA) $(JAVA_FLAGS_GATK4) -XX:ParallelGCThreads=$(THREADS_BY_ALIGNER) -jar $(GATK4) \
# 		SplitIntervals \
# 		-R $(GENOME) \
# 		-L $*.from_manifest.interval_list \
# 		--scatter-count $(THREADS) \
# 		--subdivision-mode INTERVAL_SUBDIVISION \
# 		-O $*.recalibration.scatter
# 	# Generate the sub-makefile
# 	@rm -f $*.recalibration.mk; \
# 	reports=""; \
# 	report_args=""; \
# 	for interval in $*.recalibration.scatter/*-scattered.interval_list; do \
# 		shard=$$(basename "$$interval" -scattered.interval_list); \
# 		report="$*.recalibration.$$shard.grp"; \
# 		echo "$$report: $< $<.bai $$interval" >> $*.recalibration.mk; \
# 		echo "	$(JAVA) $(JAVA_FLAGS_GATK4) -XX:ParallelGCThreads=1 -jar $(GATK4) \\" >> $*.recalibration.mk; \
# 		echo "		BaseRecalibrator \\" >> $*.recalibration.mk; \
# 		echo "			-I $< \\" >> $*.recalibration.mk; \
# 		echo "			-R $(GENOME) \\" >> $*.recalibration.mk; \
# 		echo "			$(GATK_RECALIBRATION_KNOWN_OPTIONS) \\" >> $*.recalibration.mk; \
# 		echo "			--use-original-qualities \\" >> $*.recalibration.mk; \
# 		echo "			-L $$interval \\" >> $*.recalibration.mk; \
# 		echo "			-O $$report" >> $*.recalibration.mk; \
# 		reports="$$reports $$report"; \
# 		report_args="$$report_args -I $$report"; \
# 	done; \
# 	echo "" >> $*.recalibration.mk; \
# 	echo "$@:$$reports" >> $*.recalibration.mk; \
# 	echo "	$(JAVA) $(JAVA_FLAGS_GATK4) -XX:ParallelGCThreads=$(THREADS_BY_ALIGNER) -jar $(GATK4) \\" >> $*.recalibration.mk; \
# 	echo "		GatherBQSRReports \\" >> $*.recalibration.mk; \
# 	echo "			$$report_args \\" >> $*.recalibration.mk; \
# 	echo "			-O $@" >> $*.recalibration.mk
# 	# Debug
# 	@cat $*.recalibration.mk
# 	# Run the scattered BaseRecalibrator jobs in parallel, then GatherBQSRReports
# 	+make -f $*.recalibration.mk -j $(THREADS_BY_ALIGNER) $@
# 	-rm -f $*.recalibration.*.grp
# 	-rm -rf $*.recalibration.scatter
# 	-rm -f $*.recalibration.mk


# Recalibration with gatk using ApplyBQSR with BQSR, in one line
# %.bam: %.recalibration.bam %.recalibration.bam.bai %.recalibration.bam.grp 
# 	# Recalibrate BAM with BaseRecalibrator grp file
# 	if ! $(JAVA) $(JAVA_FLAGS_GATK4) -XX:ParallelGCThreads=$(THREADS_BY_SAMPLE) -jar $(GATK4) \
# 		ApplyBQSR \
# 		-R $(GENOME) \
# 		-I $< \
# 		--bqsr-recal-file $*.recalibration.bam.grp \
# 		--use-original-qualities \
# 		--emit-original-quals \
# 		-O $@; \
# 	then \
# 		echo "Error: GATK ApplyBQSR failed for $<"; \
# 		mv $< $@; \
# 	fi
# 	-rm -f $*.recalibration.*;


RELEASE_COMMENT := "\#\# BAM RECALIBRATION: GATK BaseRecalibrator and PrintReads are used to recalibrate BAM files."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_ALIGNMENT:recalibration:BaseRecalibrator of reads in BAM. Warning: step BAM destructive, i.e. remove reads"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

