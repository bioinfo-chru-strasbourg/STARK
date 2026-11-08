############################
# GATK Recalibration Rules
# Release: 1.0.0
# Date: 23/10/2026
# Author: Antony Le Bechec
############################

# Release note
# 10/03/2015-0.9.0: Creation, BAM recalibration and variant recalibration
# 29/07/2022-0.9.2: Remove variant recalibration
# 31/10/2025-0.9.3: Change GATK3 to GATK4 for recalibration
# 18/09/2026-0.9.4: Fix option --emit-original-quals for GATK4
# 23/10/2026-1.0.0: Add filter on chromosomes using manifest/design 


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
%.bam: %.recalibration.bam %.recalibration.bam.bai %.recalibration.bam.grp %.recalibration.from_manifest.interval_list
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
			if (( $$(awk -F'\t' '$$1=="'"$$chr"'"{print $$0}' "$*.recalibration.from_manifest.interval_list" | wc -l) )); then \
				echo "$*.for_recalibration.$$chr.bam: $*.recalibration.bam" >> $*.recalibration1.mk; \
				echo "	$(JAVA) $(JAVA_FLAGS_GATK4) -Dsamjdk.compression_level=1 -jar $(GATK4) ApplyBQSR -R $(GENOME) -I $*.recalibration.bam --bqsr-recal-file $*.recalibration.bam.grp --use-original-qualities --emit-original-quals -L $$chr -O $*.for_recalibration.$$chr.bam" >> $*.recalibration1.mk; \
				echo -n " $*.for_recalibration.$$chr.bam " >> $*.recalibration2.mk; \
			else \
				echo "$*.for_recalibration.$$chr.bam: $*.recalibration.bam" >> $*.recalibration1.mk; \
				echo "	$(SAMTOOLS) view --output-fmt-option level=1 -b $*.recalibration.bam '$$chr' > $*.for_recalibration.$$chr.bam;" >> $*.recalibration1.mk; \
				echo -n " $*.for_recalibration.$$chr.bam " >> $*.recalibration2.mk; \
			fi; \
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


RELEASE_COMMENT := "\#\# BAM RECALIBRATION: GATK BaseRecalibrator and PrintReads are used to recalibrate BAM files."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_ALIGNMENT:recalibration:BaseRecalibrator of reads in BAM. Warning: step BAM destructive, i.e. remove reads"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

