############################
# Metrics Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.5b"
MK_DATE="17/03/2019"

# Release note
# 18/12/2015 - 0.9.2b : Force gzip metrics files
# 23/09/2016 - 0.9.3b : Add BAM Check metrics.bam_check. Change PICARD version
# 29/09/2016 - 0.9.4b : Add Amplicon coverage metrics.amplicon_coverage
# 29/09/2016 - 0.9.5b : Chenge metrics.genes rule


# TOOLS
SAMTOOLS?=$(NGSbin)/samtools
JAVA?=java
PICARDLIB?=$(NGSbin)/picard-tools
FASTQC?=$(NGSbin)/fastqc
NGSscripts?=$(NGS_SCRIPTS)
GZ?=gzip

# OPTIONS
BED?=
PRIMER_BED?=
BAM_METRICS?=1
FULL_COVERAGE?=0

FATBAM_TMP_FOLDER?=$(TMP_FOLDER_TMP)

METRICS_SNPEFF?=0


# Coverage
MINIMUM_DEPTH?=30
EXPECTED_DEPTH?=100
DEPTH_COVERAGE_THRESHOLD?=1
COVERAGE_CRITERIA?=1,5,10,20,30,50,100,200,300



# GLOBAL METRICS VARIABLES
METRICS_MINIMUM_MAPPING_QUALITY?=20
METRICS_MINIMUM_BASE_QUALITY?=20
CLIP_OVERLAPPING_READS?=1
METRICS_FLAGS=UNMAP,SECONDARY,QCFAIL,DUP

# GATK Guidelines
# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.1/picard_analysis_directed_CollectHsMetrics.php


# SAMTOOLS
SAMTOOLS_METRICS_DEPTH_base_q?=$(METRICS_MINIMUM_BASE_QUALITY)
SAMTOOLS_METRICS_DEPTH_map_q?=$(METRICS_MINIMUM_MAPPING_QUALITY)
SAMTOOLS_METRICS_DEPTH_PARAM?= -d $(SAMTOOLS_METRICS_DEPTH_base_q) -q $(SAMTOOLS_METRICS_DEPTH_map_q)
SAMTOOLS_METRICS_FLAG_PARAM?= -F 0x4 -F 0x100 -F 0x200 -F 0x400

SAMTOOLS_METRICS_VIEW_PARAM?= $(SAMTOOLS_METRICS_FLAG_PARAM) -q $(SAMTOOLS_METRICS_DEPTH_map_q)

SAMTOOLS_METRICS_MPILEUP_DEPTH_base_q?=$(METRICS_MINIMUM_BASE_QUALITY)
SAMTOOLS_METRICS_MPILEUP_DEPTH_map_q?=$(METRICS_MINIMUM_MAPPING_QUALITY)
SAMTOOLS_METRICS_MPILEUP_DEPTH_max_depth?=100000000
SAMTOOLS_METRICS_MPILEUP_DEPTH_adjust_MQ?=0
SAMTOOLS_METRICS_MPILEUP_FLAGS?=$(METRICS_FLAGS)
SAMTOOLS_METRICS_MPILEUP_DEPTH_excl_flags?=$(shell if [ "$(SAMTOOLS_METRICS_MPILEUP_FLAGS)" == "" ]; then echo ""; else echo "--excl-flags $(SAMTOOLS_METRICS_MPILEUP_FLAGS)"; fi;)
SAMTOOLS_METRICS_MPILEUP_DEPTH_count_orphans?=--count-orphans
SAMTOOLS_METRICS_MPILEUP_PARAM?= $(SAMTOOLS_METRICS_MPILEUP_DEPTH_excl_flags) --min-BQ $(SAMTOOLS_METRICS_MPILEUP_DEPTH_base_q) --min-MQ $(SAMTOOLS_METRICS_MPILEUP_DEPTH_map_q) --max-depth $(SAMTOOLS_METRICS_MPILEUP_DEPTH_max_depth) --adjust-MQ $(SAMTOOLS_METRICS_MPILEUP_DEPTH_adjust_MQ) $(SAMTOOLS_METRICS_MPILEUP_DEPTH_count_orphans) $(shell if ! (( $(CLIP_OVERLAPPING_READS) )); then echo " -x "; fi )

# PICARD
PICARD_CollectHsMetrics_MINIMUM_MAPPING_QUALITY?=$(METRICS_MINIMUM_MAPPING_QUALITY)
PICARD_CollectHsMetrics_MINIMUM_BASE_QUALITY?=$(METRICS_MINIMUM_BASE_QUALITY)
PICARD_CollectHsMetrics_PARAM?=MINIMUM_MAPPING_QUALITY=$(PICARD_CollectHsMetrics_MINIMUM_MAPPING_QUALITY) MINIMUM_BASE_QUALITY=$(PICARD_CollectHsMetrics_MINIMUM_BASE_QUALITY) CLIP_OVERLAPPING_READS=$(CLIP_OVERLAPPING_READS)



################################
## BED / INTERVALS for Metrics #
################################


# BED from BAM
################

%.bam.bed: %.bam %.bam.bai
	#BAM.BED from BAM
	# samtools view P1335.bwamem.bam -b | /STARK/tools/bedtools/current/bin/bedtools genomecov -ibam stdin -bg | /STARK/tools/bedtools/current/bin/bedtools merge -i stdin
	-+if ((1)); then \
	if (($$($(SAMTOOLS) idxstats $< | awk '{SUM+=$$3+$$4} END {print SUM}'))); then \
		rm -f $<.genomeCoverageBed.mk $<.genomeCoverageBed1.mk $<.genomeCoverageBed2.mk $<.genomeCoverageBed3.mk; \
		for chr in $$($(SAMTOOLS) idxstats $< | grep -v "\*" | awk '{ if ($$3+$$4>0) print $$1 }'); do \
			echo "$<.genomeCoverageBed.$$chr.bed: $<" >> $<.genomeCoverageBed1.mk; \
			echo "	$(SAMTOOLS) view $< -b $$chr | $(BEDTOOLS)/genomeCoverageBed -ibam stdin -bg | $(BEDTOOLS)/mergeBed -i - > $<.genomeCoverageBed.$$chr.bed " >> $<.genomeCoverageBed1.mk; \
			echo -n " $<.genomeCoverageBed.$$chr.bed" >> $<.genomeCoverageBed2.mk; \
		done; \
		echo -n "$@: " | cat - $<.genomeCoverageBed2.mk > $<.genomeCoverageBed3.mk; \
		echo ""  >> $<.genomeCoverageBed3.mk; \
		echo "	cat $$^ > $@ " >> $<.genomeCoverageBed3.mk; \
		echo "	-rm -f $$^ " >> $<.genomeCoverageBed3.mk; \
		#echo "	-rm -f \$^ " >> $<.genomeCoverageBed3.mk; \
		cat $<.genomeCoverageBed1.mk $<.genomeCoverageBed3.mk >> $<.genomeCoverageBed.mk; \
		make -j $(THREADS) -i -f $<.genomeCoverageBed.mk $@ 1>/dev/null 2>/dev/null; \
		rm $<.genomeCoverageBed*.mk; \
	fi; \
	fi;



# DESIGN BED
##############
# BED for metrics (either a provided BED or generated from BAM)
# without header

%.design.bed: %.withoutheader.for_metrics_bed
	cp $< $@;



# BAM for metrics
###################

%.metrics.bam: %.bam
	-ln -P $< $@;
	if [ ! -e $@ ]; then cp $< $@; fi; # if ln des not work



# BED for metrics
###################
# BED for metrics (either a provided BED or generated from BAM)
# Check if .metrics.region_clipped.bed and metrics.bed exists
# Else generate BED from BAM

%.for_metrics_bed: %.bam %.bam.bai %.metrics.from_manifest.intervals %.metrics.bed %.metrics.region_clipped.bed #%.bam.bed
	if [ "`grep ^ -c $*.metrics.bed`" == "0" ] && [ "`grep ^ -c $*.metrics.region_clipped.bed`" == "0" ]; then \
		echo "# No BED file provided... BED from the BAM file"; \
		if (($$($(SAMTOOLS) idxstats $< | awk '{SUM+=$$3+$$4} END {print SUM}'))); then \
			rm -f $<.genomeCoverageBed.mk $<.genomeCoverageBed1.mk $<.genomeCoverageBed2.mk $<.genomeCoverageBed3.mk; \
			for chr in $$($(SAMTOOLS) idxstats $< | grep -v "\*" | awk '{ if ($$3+$$4>0) print $$1 }'); do \
				echo "$<.genomeCoverageBed.$$chr.bed: $<" >> $<.genomeCoverageBed1.mk; \
				echo "	$(SAMTOOLS) view $< -b $$chr | $(BEDTOOLS)/genomeCoverageBed -ibam stdin -bg | $(BEDTOOLS)/mergeBed -i - > $<.genomeCoverageBed.$$chr.bed " >> $<.genomeCoverageBed1.mk; \
				echo -n " $<.genomeCoverageBed.$$chr.bed" >> $<.genomeCoverageBed2.mk; \
			done; \
			echo -n "$@: " | cat - $<.genomeCoverageBed2.mk > $<.genomeCoverageBed3.mk; \
			echo ""  >> $<.genomeCoverageBed3.mk; \
			echo "	cat $$^ > $@ " >> $<.genomeCoverageBed3.mk; \
			echo "	-rm -f $$^ " >> $<.genomeCoverageBed3.mk; \
			cat $<.genomeCoverageBed1.mk $<.genomeCoverageBed3.mk >> $<.genomeCoverageBed.mk; \
			make -j $(THREADS) -i -f $<.genomeCoverageBed.mk $@ 1>/dev/null 2>/dev/null; \
			rm $<.genomeCoverageBed*.mk; \
		fi; \
	else \
		echo "# BED file provided..."; \
		echo "# Extract BAM Header"; \
		$(SAMTOOLS) view -H $< > $@; \
		if [ -e $*.metrics.region_clipped.bed ] && [ "`grep ^ -c $*.metrics.region_clipped.bed`" != "0" ]; then \
			echo "# BED file from clipped BED"; \
			grep -v ^@ $*.metrics.region_clipped.bed >> $@; \
		else \
			if [ -e $*.metrics.bed ]; then \
				echo "# BED FILE from BED"; \
				grep -v ^@ $*.metrics.bed >> $@; \
			fi; \
		fi; \
	fi; \
	if [ ! -e $@ ]; then touch $@; fi



# BED without header for metrics
##################################

%.withoutheader.for_metrics_bed: %.for_metrics_bed
	if [ -s $< ]; then \
		echo "[INFO] File '$@' generated from '$<'"; \
		grep -v ^@ $< > $@; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi



# BED with only 3 fields for Picards
######################################

%.3fields.for_metrics_bed: %.for_metrics_bed
	if [ -s $< ]; then \
		echo "[INFO] File '$@' generated from '$<'"; \
		grep ^@ $< > $@; \
		grep -v ^@ $< > $*.withoutheader.for_metrics_bed.3fields.bed; \
		cat $*.withoutheader.for_metrics_bed.3fields.bed | awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t+\t"$$1":"$$2"-"$$3}' >> $@; \
		rm $*.withoutheader.for_metrics_bed.3fields.bed; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi



#############
## METRICS ##
#############

# ALL METRICS
###############

%.bam.metrics/metrics: %.bam.metrics/metrics.design %.bam.metrics/metrics.gatk %.bam.metrics/metrics.picard %.bam.metrics/metrics.samtools %.bam.metrics/metrics.regions_coverage %.bam.metrics/metrics.amplicon_coverage #%.bam.metrics/metrics.bam_check
	#Create directory
	mkdir -p $(@D)
	cat $^ > $@
	echo "#[$$(date)] BAM Metrics done" >> $@
	-rm -f $*.for_metrics_bed



# Design
##########

%.bam.metrics/metrics.design: %.list.genes %.design.bed
	# Create directory
	mkdir -p $(@D)
	touch $@
	# Foreach design and genes files
	for one_bed in $$(cat $*.list.genes) $*.design.bed; do \
		if [ -s $$one_bed ]; then \
			cp $$one_bed $(@D)/$$(basename $$one_bed); \
			if [ -s $(@D)/$$(basename $$one_bed) ]; then \
				echo "#[INFO] COPY of design/genes '"$$(basename $$one_bed)"' done. See '$(@D)/$$(basename $$one_bed)'. " >> $@; \
			else \
				echo "#[INFO] COPY of design/genes '"$$(basename $$one_bed)"' failed. " >> $@; \
			fi; \
		fi; \
	done;
	[ ! -z $@ ] && echo "#[INFO] COPY of design/genes not done because not bed/genes files. " >> $@;



# MarkDuplicates metrics
##########################

%.bam.metrics/metrics.markDuplicates: %.bam %.bam.bai
	cp -f $**.markduplicates.bam.metrics/*.markDuplicates.metrics $(@D)/;
	rm -rf $**.markduplicates.bam.metrics;



# Amplicon metrics
####################
# From FATBAM

%.bam.metrics/metrics.amplicon_coverage: %.bam %.bam.bai %.manifest %.genome
	mkdir -p $(@D) ;
	-+$(FATBAM_COVERAGE) --env=$(CONFIG_TOOLS) --ref=$$(cat $*.genome)  --bam=$< --output=$(@D)/$(*F).amplicon_coverage --manifest=$*.manifest --multithreading --threads=$(THREADS) -v --tmp=$(FATBAM_TMP_FOLDER) --verbose 1>$(@D)/$(*F).amplicon_coverage.log 2>$(@D)/$(*F).amplicon_coverage.err;
	cat $(@D)/$(*F).amplicon_coverage.log $(@D)/$(*F).amplicon_coverage.err;
	echo "#[$$(date)] BAM Amplicon Coverage Metrics done" > $@;



# GATK metrics
################

GATKDOC_FLAGS= -rf BadCigar -allowPotentiallyMisencodedQuals
%.bam.metrics/metrics.gatk: %.bam %.bam.bai %.genome %.for_metrics_bed %.3fields.for_metrics_bed
	# TODO: speed up ! Too loog for exome/genome...
	# use split algorithm with makefile and $(SAMTOOLS) view $< -b $$chr | $(BEDTOOLS)/genomeCoverageBed -ibam stdin -bg
	# see rule %.bam.bed
	# Create directory
	mkdir -p $(@D);
	if (($(BAM_METRICS))) && ((1)); then \
		grep -v ^@ $*.for_metrics_bed > $*.withoutheader.for_metrics_bed.gatk.bed ; \
		# GATK DepthOfCoverage needs BED without HEADER!!! ; \
		if [ ! -e $(@D)/$(*F) ]; then \
			$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKDOC_FLAGS) \
				-T DepthOfCoverage \
				-R $$(cat $*.genome) \
				-o $(@D)/$(*F) \
				-I $< \
				-L $*.withoutheader.for_metrics_bed.gatk.bed; \
		fi; \
		rm $*.withoutheader.for_metrics_bed.gatk.bed; \
		echo "#[$$(date)] BAM GATK Metrics done" > $@; \
	else \
		echo "#[$$(date)] BAM GATK not done because BAM_METRICS=0" > $@; \
	fi;





# BAM Validation
##################

%.validation.bam: %.bam %.bam.bai #%.list.genes %.design.bed
	# Create directory ;
	mkdir -p $(@D);
	# BAM Validation
	$(SAMTOOLS) view $(SAMTOOLS_METRICS_VIEW_PARAM) -h $< -1 -@ $(THREADS) > $@ ;



# Empty HsMetrics file
########################

%.empty.HsMetrics:
	echo "## METRICS CLASS	picard.analysis.directed.HsMetrics" > $@
	echo "BAIT_SET	GENOME_SIZE	BAIT_TERRITORY	TARGET_TERRITORY	BAIT_DESIGN_EFFICIENCY	TOTAL_READS	PF_READS	PF_UNIQUE_READS	PCT_PF_READS	PCT_PF_UQ_READS	PF_UQ_READS_ALIGNED	PCT_PF_UQ_READS_ALIGNED	PF_UQ_BASES_ALIGNED	ON_BAIT_BASES	NEAR_BAIT_BASES	OFF_BAIT_BASES	ON_TARGET_BASES	PCT_SELECTED_BASES	PCT_OFF_BAIT	ON_BAIT_VS_SELECTED	MEAN_BAIT_COVERAGE	MEAN_TARGET_COVERAGE	PCT_USABLE_BASES_ON_BAIT	PCT_USABLE_BASES_ON_TARGET	FOLD_ENRICHMENT	ZERO_CVG_TARGETS_PCT	FOLD_80_BASE_PENALTY	PCT_TARGET_BASES_2X	PCT_TARGET_BASES_10X	PCT_TARGET_BASES_20X	PCT_TARGET_BASES_30X	PCT_TARGET_BASES_40X	PCT_TARGET_BASES_50X	PCT_TARGET_BASES_100X	HS_LIBRARY_SIZE	HS_PENALTY_10X	HS_PENALTY_20X	HS_PENALTY_30X	HS_PENALTY_40X	HS_PENALTY_50X	HS_PENALTY_100X	AT_DROPOUT	GC_DROPOUT	SAMPLE	LIBRARY	READ_GROUP" >> $@
	echo $$(echo $$(basename $@) | awk -F"." '{print $$1}')"	?	?	?	?	0	0	0	?	?	0	?	0	0	0	0	0	?	?	?	0	?	?	?	?	1	?	0	0	0	0	0	0	0	?	0	0	0	0	0	0	0	0" >> $@



# PICARD Metrics
##################

%.bam.metrics/metrics.picard: %.validation.bam %.validation.bam.bai %.empty.HsMetrics %.genome %.list.genes %.design.bed %.dict
	# Create directory
	mkdir -p $(@D)
	touch $@
	# Picard Metrics needs BED with Header and 3+3 fields!!!
	for one_bed in $$(cat $*.list.genes) $*.design.bed; do \
		if [ -s $$one_bed ]; then \
			$(JAVA) $(JAVA_FLAGS_BY_SAMPLE) -jar $(PICARD) BedToIntervalList I=$$one_bed O=$(@D)/$(*F).$$(basename $$one_bed).interval SD=$$(cat $*.dict) ; \
			$(JAVA) $(JAVA_FLAGS_BY_SAMPLE) -jar $(PICARD) CollectHsMetrics INPUT=$*.validation.bam OUTPUT=$(@D)/$(*F).$$(basename $$one_bed).HsMetrics R=$$(cat $*.genome) BAIT_INTERVALS=$(@D)/$(*F).$$(basename $$one_bed).interval TARGET_INTERVALS=$(@D)/$(*F).$$(basename $$one_bed).interval PER_TARGET_COVERAGE=$(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_target_coverage $(PICARD_CollectHsMetrics_PARAM) VALIDATION_STRINGENCY=SILENT 2>$(@D)/$(*F).$$(basename $$one_bed).HsMetrics.err ; \
			rm $(@D)/$(*F).$$(basename $$one_bed).interval; \
		fi; \
		if [ ! -s $(@D)/$(*F).$$(basename $$one_bed).HsMetrics ]; then \
			cp $*.empty.HsMetrics $(@D)/$(*F).$$(basename $$one_bed).HsMetrics; \
			echo "#[ERROR] BAM PICARD Metrics failed. Empty HsMetrics file generated. See '$(@D)/$(*F).$$(basename $$one_bed).HsMetrics.err'" >> $@; \
		else \
			echo "#[INFO] BAM PICARD Metrics done. HsMetrics file generated. See '$(@D)/$(*F).$$(basename $$one_bed).HsMetrics'" >> $@; \
		fi \
	done;
	echo "#[INFO] BAM PICARD Metrics done" >> $@




# SAMTOOLS metrics
####################


%.bam.metrics/metrics.samtools: %.bam.metrics/metrics.samtools.flagstat %.bam.metrics/metrics.samtools.idxstats %.bam.metrics/metrics.samtools.genomeCoverage %.bam.metrics/metrics.samtools.depth \
	%.bam.metrics/metrics.samtools.depthbed.coverage %.bam.metrics/metrics.samtools.on.target %.bam.metrics/metrics.samtools.off.target %.bam.metrics/metrics.samtools.genomeCoverageBed
	# Create directory ;
	mkdir -p $(@D);
	cat $^ > $@
	echo "#[$$(date)] SAMTOOLS Metrics done" >> $@


%.bam.metrics/metrics.samtools.flagstat: %.bam %.bam.bai
	# Create directory ;
	mkdir -p $(@D);
	> $@;
	$(SAMTOOLS) flagstat $< > $(@D)/$(*F).flagstat ;
	if [ -s $(@D)/$(*F).flagstat ]; then \
		echo "#[INFO] SAMTOOLS flagstat done. See '$(@D)/$(*F).flagstat'. " >> $@; \
	else \
		echo "#[INFO] SAMTOOLS flagstat failed. " >> $@; \
	fi;


%.bam.metrics/metrics.samtools.idxstats: %.bam %.bam.bai
	# Create directory ;
	mkdir -p $(@D);
	> $@;
	$(SAMTOOLS) idxstats $< > $(@D)/$(*F).idxstats ;
	if [ -z $(@D)/$(*F).idxstats ]; then \
		echo "#[INFO] SAMTOOLS idxstats done. See '$(@D)/$(*F).idxstats'. " >> $@; \
	else \
		echo "#[INFO] SAMTOOLS idxstats failed. " >> $@; \
	fi;


%.bam.metrics/metrics.samtools.genomeCoverage: %.bam %.bam.bai
	# Create directory ;
	mkdir -p $(@D);
	> $@;
	if (($(BAM_METRICS))); then \
		if (($(FULL_COVERAGE))); then \
			$(BEDTOOLS)/genomeCoverageBed -ibam $< -dz > $(@D)/$(*F).genomeCoverage; \
			$(GZ) --best -f $(@D)/$(*F).genomeCoverage; \
			if [ -s $(@D)/$(*F).genomeCoverage.gz ]; then \
				echo "#[INFO] BEDTOOLS genomeCoverage done. See '$(@D)/$(*F).genomeCoverage.gz'. " >> $@; \
			else \
				echo "#[INFO] SAMTOOLS genomeCoverage failed. " >> $@; \
			fi; \
		else \
			echo "#[INFO] SAMTOOLS genomeCoverage not done because FULL_COVERAGE=0. " >> $@; \
		fi; \
	else \
		echo "#[INFO] SAMTOOLS genomeCoverage not done because BAM_METRICS=0. " >> $@; \
	fi;


%.bam.metrics/metrics.samtools.depth: %.bam %.bam.bai
	# Create directory ;
	mkdir -p $(@D);
	> $@;
	if (($(BAM_METRICS))); then \
		if (($(FULL_COVERAGE))); then \
			$(SAMTOOLS) mpileup $(SAMTOOLS_METRICS_MPILEUP_PARAM) $< | cut -f1,2,4 > $(@D)/$(*F).depth; \
			$(GZ) --best -f $(@D)/$(*F).depth; \
			if [ -s $(@D)/$(*F).idxstats.gz ]; then \
				echo "#[INFO] SAMTOOLS depth done. See '$(@D)/$(*F).depth.gz'. " >> $@; \
			else \
				echo "#[INFO] SAMTOOLS depth failed. " >> $@; \
			fi; \
		else \
			echo "#[INFO] SAMTOOLS depth not done because FULL_COVERAGE=0. " >> $@; \
		fi; \
	else \
		echo "#[INFO] SAMTOOLS depth not done because BAM_METRICS=0. " >> $@; \
	fi;


%.bam.metrics/metrics.samtools.depthbed.coverage: %.validation.bam %.validation.bam.bai %.list.genes %.design.bed
	# Create directory ;
	mkdir -p $(@D);
	> $@;
	# For each BED GENES
	for one_bed in $$(cat $*.list.genes) $*.design.bed; do \
		# Test if BED exists \
		if [ -e $$one_bed ]; then \
			# depthbed generation \
			cat $$one_bed | awk -F"\t" '{print $$1"\t"$$2-1"\t"$$3"\t"$$4"\t"$$5}' > $$one_bed.0_based; \
			$(SAMTOOLS) mpileup $(SAMTOOLS_METRICS_MPILEUP_PARAM) -aa -l $$one_bed.0_based $< | cut -f1,2,4 > $(@D)/$(*F).$$(basename $$one_bed).depthbed; \
			rm -f $$one_bed.0_based; \
			if [ -s $(@D)/$(*F).$$(basename $$one_bed).depthbed ]; then \
				echo "#[INFO] SAMTOOLS depthbed with '"$$(basename $$one_bed)"' done. See '$(@D)/$(*F).$$(basename $$one_bed).depthbed'. " >> $@; \
			else \
				echo "#[INFO] SAMTOOLS depthbed with '"$$(basename $$one_bed)"' failed. " >> $@; \
			fi; \
			# coverage generation \
			cat $(@D)/$(*F).$$(basename $$one_bed).depthbed | \
			awk -v MDP=$$(echo $(COVERAGE_CRITERIA) | tr "," "\n" | sort -n | tail -n 1)  '{SUM++} { if ($$3>MDP) {DP[MDP]++} else {DP[$$3]++} } END { for (i=MDP; i>=0; i-=1) {print i" "DP[i]" SUM"SUM}}' | \
			sort -g -r | awk -v COVERAGE_CRITERIA=$(COVERAGE_CRITERIA) '{SUM+=$$2} {CUM[$$1]=SUM} {split(COVERAGE_CRITERIA,C,",")} END { for (j in C) {print C[j]" "CUM[C[j]]" SUM "(CUM[C[j]]/SUM)} }' | \
			sort -g > $(@D)/$(*F).$$(basename $$one_bed).coverage; \
			echo "#[INFO] SAMTOOLS depthbed and coverage with '"$$(basename $$one_bed)"' done. See '$(@D)/$(*F).$$(basename $$one_bed).depthbed' and  '$(@D)/$(*F).$$(basename $$one_bed).coverage'. " >> $@; \
			if [ -s $(@D)/$(*F).$$(basename $$one_bed).coverage ]; then \
				echo "#[INFO] SAMTOOLS coverage with '"$$(basename $$one_bed)"' done. See '$(@D)/$(*F).$$(basename $$one_bed).coverage'. " >> $@; \
			else \
				echo "#[INFO] SAMTOOLS coverage with '"$$(basename $$one_bed)"' failed. " >> $@; \
			fi; \
		else \
			echo "#[INFO] SAMTOOLS depthbed and coverage with '"$$(basename $$one_bed)"' failed, because no '"$$(basename $$one_bed)"'." >> $@; \
		fi; \
	done;
	[ ! -z $@ ] && echo "#[INFO] SAMTOOLS depthbed and coverage not done because not bed/genes files. " >> $@;


%.bam.metrics/metrics.samtools.on.target: %.validation.bam %.validation.bam.bai %.list.genes %.design.bed
	# Create directory ;
	mkdir -p $(@D);
	> $@;
	# For each BED GENES
	for one_bed in $$(cat $*.list.genes) $*.design.bed; do \
		# Test if BED exists \
		if [ -e $$one_bed ]; then \
			$(SAMTOOLS) view -c $*.validation.bam -L $$one_bed > $(@D)/$(*F).$$(basename $$one_bed).on.target; \
			if [ -s $(@D)/$(*F).$$(basename $$one_bed).coverage ]; then \
				echo "#[INFO] ON TARGET with '"$$(basename $$one_bed)"' done. Number of reads aligned within the regions defined in the BED. See '$(@D)/$(*F).$$(basename $$one_bed).on.target'. " >> $@; \
			else \
				echo "#[INFO] ON TARGET with '"$$(basename $$one_bed)"' failed. " >> $@; \
			fi; \
		else \
			echo "#[INFO] ON TARGET with '"$$(basename $$one_bed)"' failed, because no '"$$(basename $$one_bed)"'." >> $@; \
		fi; \
	done;
	[ ! -z $@ ] && echo "#[INFO] ON TARGET not done because not bed/genes files. " >> $@;


%.bam.metrics/metrics.samtools.off.target: %.validation.bam %.validation.bam.bai %.list.genes %.design.bed
	# Create directory ;
	mkdir -p $(@D);
	> $@;
	# For each BED GENES
	for one_bed in $$(cat $*.list.genes) $*.design.bed; do \
		# Test if BED exists \
		if [ -e $$one_bed ]; then \
			$(BEDTOOLS)/intersectBed -abam $*.validation.bam -b $$one_bed -v | $(SAMTOOLS) view -c - > $(@D)/$(*F).$$(basename $$one_bed).off.nbreads; \
			if [ -s $(@D)/$(*F).$$(basename $$one_bed).coverage ]; then \
				echo "#[INFO] OFF TARGET with '"$$(basename $$one_bed)"' done. Number of reads aligned not in the regions defined in the BED. See '$(@D)/$(*F).$$(basename $$one_bed).off.target'. " >> $@; \
			else \
				echo "#[INFO] OFF TARGET with '"$$(basename $$one_bed)"' failed." >> $@; \
			fi; \
		else \
			echo "#[INFO] OFF TARGET with '"$$(basename $$one_bed)"' failed, because no '"$$(basename $$one_bed)"'." >> $@; \
		fi; \
	done;
	[ ! -z $@ ] && echo "#[INFO] OFF TARGET not done because not bed/genes files. " >> $@;



%.bam.metrics/metrics.samtools.genomeCoverageBed: %.validation.bam %.validation.bam.bai %.list.genes %.design.bed
	# Create directory ;
	mkdir -p $(@D);
	> $@;
	if (($(BAM_METRICS))); then \
		# For each BED GENES \
		for one_bed in $$(cat $*.list.genes) $*.design.bed; do \
			# Test if BED exists \
			if [ -e $$one_bed ]; then \
				$(SAMTOOLS) view -b $*.validation.bam -L $$one_bed | $(SAMTOOLS) sort | $(BEDTOOLS)/genomeCoverageBed -ibam stdin -dz > $(@D)/$(*F).$$(basename $$one_bed).genomeCoverageBedbed; \
				if [ -s $(@D)/$(*F).$$(basename $$one_bed).coverage ]; then \
					echo "#[INFO] BEDTOOLS genomeCoverage with '"$$(basename $$one_bed)"' done. See '$(@D)/$(*F).$$(basename $$one_bed).genomeCoverageBedbed'. " >> $@; \
				else \
					echo "#[INFO] BEDTOOLS genomeCoverage with '"$$(basename $$one_bed)"' failed" >> $@; \
				fi; \
			else \
				echo "#[INFO] BEDTOOLS genomeCoverage with '"$$(basename $$one_bed)"' failed, because no '"$$(basename $$one_bed)"'. " >> $@; \
			fi; \
		done; \
		[ ! -z $@ ] && echo "#[INFO] BEDTOOLS genomeCoverage not done because not bed/genes files. " >> $@; \
	else \
		echo "#[INFO] BEDTOOLS genomeCoverage not done because BAM_METRICS=0. " >> $@; \
	fi;





# FASTQ metrics
#################
# FatsQC metrics and counts metrics

%.sequencing/metrics: %.sequencing/metrics.fastqc %.sequencing/metrics.counts %.sequencing/metrics.Q30  %.sequencing/metrics.infos
	# create directory
	mkdir -p $(@D)
	cat $^ > $@
	-rm -f $^



# FastQC metrics
##################

#%.fastqc/metrics.fastqc: %.fastq.gz
%.sequencing/metrics.fastqc: %.fastq.gz
	# create directory
	mkdir -p $(@D)
	# create link
	ln -s $< $(@D)/$(*F).fastq.gz
	# touch target
	touch $@;
	# FASTQC
	-if (($$(zcat $< | head -n 1 | wc -l))); then \
		#$(FASTQC) $< --outdir=$(@D) --casava --extract; \
		$(FASTQC) $(@D)/$(*F).fastq.gz --outdir=$(@D) --casava --extract --threads $(THREADS_BY_SAMPLE) ; \
		cp $(@D)/$(*F)_fastqc/fastqc_data.txt $(@D)/metrics.fastqc.txt; \
		echo "#[INFO] FASTQC done. See 'metrics.fastqc.txt' file." >> $@; \
	else \
		echo "#[ERROR] FASTQC can't be launched. No Reads in FASTQ file '$<'" >> $@; \
	fi;
	# create link
	-rm -f $(@D)/$(*F).fastq.gz



# Q30
########

%.sequencing/metrics.Q30: %.sequencing/metrics.fastqc
	# create directory
	mkdir -p $(@D)
	cat $(@D)/metrics.fastqc.txt | awk '/>>Per sequence quality scores/,/>>END_MODULE/' | head -n -1 | tail -n+3 | awk '{s+=$$2}END{print s}' > $@.Q30_ALL;
	cat $(@D)/metrics.fastqc.txt | awk '/>>Per sequence quality scores/,/>>END_MODULE/' | awk '/>>Per sequence quality scores/,/30\t/' | head -n -1 | tail -n+3 | awk '{s+=$$2}END{print s}' > $@.Q30_UNTIL;
	echo -e "Q30\t"$$(echo "scale = 4; (($$(cat $@.Q30_ALL) - $$(cat $@.Q30_UNTIL)) / $$(cat $@.Q30_ALL) * 100)" | bc | sed s/[0]*$$//) > $@.txt;
	echo "#[INFO] Q30 calculation done. See 'metrics.Q30.txt' file." > $@;
	-rm -f $@.Q30_ALL $@.Q30_UNTIL



# Reads Count
###############

#%.fastqc/metrics.counts: %.fastq.gz #%.fastqc/metrics.fastqc
%.sequencing/metrics.counts: %.fastq.gz #%.fastqc/metrics.fastqc
	# create directory
	mkdir -p $(@D)
	-if (($$(zcat $< | head -n 1 | wc -l))); then \
		# header \
		echo "total unique %unique maxRead count_maxRead %count_maxRead" > $@.tmp; \
		# Counts \
		if [ $$(zcat $< | head -n 1 | wc -l) ]; then \
			zcat $< | awk '{unique=0} ((NR-2)%4==0){read=$$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}' >> $@.tmp; \
		else \
			echo "0 0 - - - -" >> $@.tmp; \
		fi; \
		# transposition \
		awk '{ for (i=1; i<=NF; i++)  {a[NR,i] = $$i} } NF>p { p = NF } END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++){ str=str" "a[i,j]; } print str } }' $@.tmp | tr " " "\t" > $@.txt; \
		# nb bases \
		echo -e "count_bases\t"$$(zcat $< | paste - - - - | cut -f4 | wc -c) >> $@.txt; \
		# script \
		# echo "" > $*.fastqc/metrics.counts.txt; \
		echo "#[INFO] FASTQ Counts done. See 'metrics.counts.txt' file." > $@; \
	else \
		echo "#[ERROR] FASTQ Counts ERROR. No reads in '$<'..." > $@; \
	fi;
	-rm -f $@.tmp



# techno
########

%.sequencing/metrics.infos: %.manifest %.R2.fastq.gz
	# create directory
	mkdir -p $(@D)
	touch $@.txt
	echo -e $$((($$(zcat $*.R2.fastq.gz | head -n1 | wc -l))) && echo "mode\tPaired-End" || echo "mode\tSingle-End") >> $@.txt
	echo -e $$((($$(grep "Amplicon Start" $*.manifest | grep -c "Upstream Probe Length"))) && echo "technology\tAmplicon" || echo "technology\tCapture") >> $@.txt
	echo "#[INFO] SEQUENCING INFOS done. See 'metrics.infos.txt' file." > $@;
	#-rm -f $@.*




# UnAligned BAM metrics (deprecated)
######################################

%.from_unaligned_bam.fastqc/metrics: %.unaligned.bam
	mkdir -p $(@D)
	touch $@;
	-if [ "$$($(SAMTOOLS) view $< | grep ^ -c)" == "0" ]; then \
		echo "#FASTQGC can't be launched. No Reads in the unaligned BAM file '$<'" >> $@; \
	else \
		#$(FASTQC) $< --outdir=$(@D) --casava --extract; \
		$(FASTQC) $< --outdir=$(@D) --casava --extract --threads $(THREADS) ; \
		echo "#FASTQGC done" >> $@; \
	fi;



# VCF METRICS
###############
# SNPEFF and BCFTOOLS metrics

%.vcf.metrics/metrics: %.vcf.metrics/metrics.snpeff %.vcf.metrics/metrics.bcftools %.vcf.metrics/metrics.info_field %.vcf.metrics/metrics.genes
	cat $^ > $@



# SNPEFF Metrics
##################
# SNPEFF metrics through HOWARD

%.vcf.metrics/metrics.snpeff: %.vcf
	mkdir -p $(@D);
	touch $@;
	if (($(METRICS_SNPEFF))); then \
		+$(HOWARD) --input=$< --output=$@.vcf --snpeff_stats=$@.html --annotation=null --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS)  --force; \
		echo "#[INFO] snpEff metrics done. See '$@.html' file." > $@; \
	else \
		echo "#[INFO] snpEff metrics NOT done." > $@; \
	fi;



# BCFTOOLS metrics
####################
# Stats from BCFTOOLS

%.vcf.metrics/metrics.bcftools: %.vcf
	mkdir -p $(@D);
	touch $@;
	-$(BCFTOOLS) stats $< > $@.stats
	echo "#[INFO] BCFTOOLS metrics done. See '$@.stats' file." > $@;



# INFO stats
####################
# Stats of INFO field

%.vcf.metrics/metrics.info_field: %.vcf
	mkdir -p $(@D);
	touch $@;
	grep -v "^#" $< | cut -f8 | tr ";" "\n" | sort | uniq -c | sed "s/^      / /gi" | tr "=" " " | awk '{print $$2"\t"$$3"\t"$$1} {a[$$2]+=$$1} END { for (key in a) { print "#\t" key "\t" a[key] } }' | sort > $@.stats
	echo "#[INFO] INFO field stats done. See '$@.stats' file." > $@;



# INFO stats
####################
# Stats of INFO field

%.vcf.metrics/metrics.genes: %.vcf.gz %.vcf.gz.tbi %.list.genes
	mkdir -p $(@D);
	> $@;
	echo "LIST_GENES:"
	cat $*.list.genes;
	+for one_bed in $$(cat $*.list.genes); do \
		cat "ONE_BED: "$$one_bed; \
		$(BCFTOOLS) view $*.vcf.gz -R $$one_bed > $@.$$(basename $$one_bed).vcf; \
		$(BCFTOOLS) stats $@.$$(basename $$one_bed).vcf > $@.$$(basename $$one_bed).vcf.bcftools.stats; \
		grep -v "^#" $@.$$(basename $$one_bed).vcf | cut -f8 | tr ";" "\n" | sort | uniq -c | sed "s/^      / /gi" | tr "=" " " | awk '{print $$2"\t"$$3"\t"$$1} {a[$$2]+=$$1} END { for (key in a) { print "#\t" key "\t" a[key] } }' | sort > $@.$$(basename $$one_bed).vcf.info_field.stats; \
		(($(METRICS_SNPEFF))) && $(HOWARD) --input=$@.$$(basename $$one_bed).vcf --output=$@.$$(basename $$one_bed).vcf.snpeff.vcf --snpeff_stats=$@.$$(basename $$one_bed).vcf.snpeff.html --annotation=null --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS)  --force ; \
		echo "#[INFO] VCF filtered by '$$one_bed' done. See '$@.$$(basename $$one_bed).vcf' file, and stats '$@.$$(basename $$one_bed).vcf.*'" >> $@; \
	done;




# List of Genes BED
#####################


%.list.genes: %.bed %.manifest
	mkdir -p $(@D);
	#touch $@;
	+if [ ! -s $@ ]; then \
		if [ -s $$(sample=$$(basename $* | cut -d. -f1 ); echo "$(*D)/$$sample.list.genes") ] ; then \
			echo "BEDFILE_GENES for $* from SAMPLE.list.genes in same directory"; \
			for g in $$(sample=$$(basename $* | cut -d. -f1 ); cat "$(*D)/$$sample.list.genes"); do \
				echo "BEDFILE_GENES for $* from SAMPLE.list.genes: $(*D)/$$g "; \
				bedfile_genes_list="$$bedfile_genes_list $(*D)/$$g"; \
			done; \
		elif [ -s $$(sample=$$(basename $* | cut -d. -f1 ); echo "$$(dirname $(*D))/$$sample.list.genes") ] ; then \
			echo "BEDFILE_GENES for $* from SAMPLE.list.genes from previous directory"; \
			for g in $$(sample=$$(basename $* | cut -d. -f1 ); cat "$$(dirname $(*D))/$$sample.list.genes"); do \
				echo "BEDFILE_GENES for $* from SAMPLE.list.genes: $$(dirname $(*D))/$$g "; \
				bedfile_genes_list="$$bedfile_genes_list $$(dirname $(*D))/$$g"; \
			done; \
		elif [ -s `file=$$( echo $* | cut -d. -f1 ); echo "$$file.genes"` ] ; then \
			echo "BEDFILE_GENES for $* from SAMPLE.genes"; \
			bedfile_genes_list=`file=$$( echo $* | cut -d. -f1 ); echo "$$file.genes"`; \
		#elif [ "$(BEDFILE_GENES)" != "" ]; then \
		elif (( $$(echo $(BEDFILE_GENES) | wc -w) )); then \
			echo "BEDFILE_GENES for $* from BEDFILE_GENES variable. Use all files to generate metrics."; \
			bedfile_genes_list=""; \
			for bedfile_genes in $$(echo $(BEDFILE_GENES) | tr "+" " " | tr " " "\n" | sort -u); \
			do \
				bedfile_name=$$( basename $$bedfile_genes | sed "s/\.genes$$//" ); \
				cp $$bedfile_genes $(@D)/$(*F)_$${bedfile_name}.bed; \
				bedfile_genes_list="$$bedfile_genes_list $(@D)/$(*F)_$${bedfile_name}.bed"; \
			done; \
		else \
			echo "BEDFILE_GENES for $* generated from SAMPLE.manifest"; \
			bedfile_genes_list=`file=$$( echo $* | cut -d. -f1 ); echo "$$file.from_design.genes"`; \
			# ici on va faire une intersection avec le manifest et ce sera notre bed par defaut sil nexiste pas; \
			#bed=$$(sample=$$(basename $* | cut -d. -f1 ); echo "$(*D)/$$sample.bed"); \
			#if [ ! -s $$(echo "$$bed") ] ; then \
			#	bed=$$(sample=$$(basename $* | cut -d. -f1 ); echo "$$(dirname $(*D))/$$sample.bed"); \
			#fi; \
			#manifest=$$(sample=$$(basename $* | cut -d. -f1 ); echo "$(*D)/$$sample.manifest"); \
			#if [ ! -s $$(echo "$$manifest") ] ; then \
			#	manifest=$$(sample=$$(basename $* | cut -d. -f1 ); echo "$$(dirname $(*D))/$$sample.manifest"); \
			#fi; \
			#if [ -s `echo "$$bed"` ] ; then \
			if [ -s $*.bed ] ; then \
				#echo "genes from BED '$$bed'" >> $*.test; \
				cut -f1,2,3 $*.bed > $@.manifest.bed ; \
			#elif [ -s `echo "$$manifest"` ] ; then \
			elif [ -s $*.manifest ] ; then \
				# manifest to bed ; \
				$(FATBAM_ManifestToBED) --input "$*.manifest" --output "$@.manifest.bed.tmp" --output_type "region" --type=PCR; \
				cut -f1,2,3 $@.manifest.bed.tmp > $@.manifest.bed ; \
				#rm $@.manifest.bed.tmp ; \
			fi; \
			if [ -s $@.manifest.bed ] ; then \
				echo "MANIFEST.bed to bedfile_genes_list : $bedfile_genes_list"; \
				$(BEDTOOLS)/bedtools intersect -wb -a $@.manifest.bed -b $(REFSEQ_GENES) | cut -f8 | sort -u > $@.manifest.bed.intersect; \
				sort -k5 $(REFSEQ_GENES) > $@.manifest.bed.refseq; \
				join -1 1 -2 5 $@.manifest.bed.intersect $@.manifest.bed.refseq -o 2.1,2.2,2.3,2.5 | sort -u -k1,2 | tr " " "\t" | awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t+\t"$$4}' > $$bedfile_genes_list; \
				#rm -f $*.manifest.bed.intersect $*.manifest.bed.refseq; \
			else \
				echo "#[ERROR] Generating GENES failed"; \
			fi;\
			rm -f $@.manifest*; \
		fi; \
		echo "bed_file list is : `echo $$bedfile_genes_list` "; \
		echo $$bedfile_genes_list | tr " " "\n" > $@ ; \
	else \
		echo "BEDFILE_GENES exists!!! "; \
	fi;


# Regnions COVERAGE METRICS
#############################

%.bam.metrics/metrics.regions_coverage: %.bam %.bam.bai %.list.genes %.design.bed %.bam.metrics/metrics.samtools.depthbed.coverage
	mkdir -p $(@D);
	touch $@;
	+if (( 1 )) ; then \
		if [ ! -z $*.list.genes ]; then \
			for one_bed in $$(cat $*.list.genes) $*.design.bed; \
			do \
				if [ -e $$one_bed ]; then \
					#bedfile_name=$$( basename $$one_bed | sed "s/\.genes$$//" ); \
					bedfile_name=$$( basename $$one_bed ); \
					$(NGSscripts)/genesCoverage.sh -f $*.bam -b $$one_bed -c "$(COVERAGE_CRITERIA)" --coverage-bases=$(@D)/$(*F).$$bedfile_name.depthbed --dp_fail=$(MINIMUM_DEPTH) --dp_warn=$(EXPECTED_DEPTH) --dp_threshold=$(DEPTH_COVERAGE_THRESHOLD) -n $(NB_BASES_AROUND) -t $(BEDTOOLS)/bedtools -s $(SAMTOOLS) --threads=$(THREADS) -o $(@D)/$(*F).$$bedfile_name; \
					echo "#[$$(date)] BAM Metrics on Regions coverage with $$bedfile_name BED file done" >> $@; \
				else \
					echo "#[$$(date)] BAM Metrics on Regions coverage with $$bedfile_name BED file FAILED" >> $@; \
				fi; \
			done; \
		else \
			echo "#[$$(date)] BAM Metrics on Regions coverage not generate because no BED list" >> $@; \
		fi; \
	else \
		echo "#[$$(date)] BAM Metrics on Regions coverage not done because BAM_METRICS=0" >> $@; \
	fi;
	if [ ! -e $@ ]; then echo "#[$$(date)] BAM Metrics on Regions coverage FAILED" > $@; fi;



# BAM CHECK rule (deprecated)
###############################

%.bam.metrics/metrics.bam_check: %.bam %.bam.bai
	# CHECK NUMBER of READS in BAM \
	# Create directory if not created
	mkdir -p $(@D);
	# Pattern
	echo $$(dirname $*)/$$(basename $$(dirname $*)) > $@.sample_pattern;
	# Header
	echo "# Check number of reads between $< and $$(cat $@.sample_pattern).R1.fastq.gz $$(cat $@.sample_pattern).R2.fastq.gz to ensure no reads removing (only mark reads)" >$(@D)/$(*F).bam_check;
	echo "# (possible new reads marked as secondary aligned)" >>$(@D)/$(*F).bam_check;
	# Check
	if [ -s $$(cat $@.sample_pattern).R1.fastq.gz ] && [ -s $$(cat $@.sample_pattern).R2.fastq.gz ]; then \
		echo "# Files $$(cat $@.sample_pattern).R1.fastq.gz and $$(cat $@.sample_pattern).R2.fastq.gz exists" >>$(@D)/$(*F).bam_check; \
		# Compare number of reads \
		echo "# $$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $<) reads for $<" >>$(@D)/$(*F).bam_check; \
		echo "# $$(zcat $$(cat $@.sample_pattern).R1.fastq.gz $$(cat $@.sample_pattern).R2.fastq.gz | grep ^@ -c ) reads for $$(cat $@.sample_pattern).R1.fastq.gz $$(cat $@.sample_pattern).R2.fastq.gz" >>$(@D)/$(*F).bam_check; \
		grep ".* reads for $<" $(@D)/$(*F).bam_check | cut -d" " -f2; \
		grep ".* reads for $$(cat $@.sample_pattern).R1.fastq.gz $$(cat $@.sample_pattern).R2.fastq.gz" $(@D)/$(*F).bam_check | cut -d" " -f2; \
		if [ "$$(grep ".* reads for $<" $(@D)/$(*F).bam_check | cut -d" " -f2)" != "$$(grep ".* reads for $$(cat $@.sample_pattern).R1.fastq.gz $$(cat $@.sample_pattern).R2.fastq.gz" $(@D)/$(*F).bam_check | cut -d" " -f2)" ]; then \
			# Error \
			echo "#[ERROR] KO Number of reads between $< and $$(cat $@.sample_pattern).R1.fastq.gz $$(cat $@.sample_pattern).R2.fastq.gz !!!" >>$(@D)/$(*F).bam_check; \
			echo "# BCFTOOLS STATS for $<" >>$(@D)/$(*F).bam_check;  \
			$(SAMTOOLS) stats $< | grep ^SN >>$(@D)/$(*F).bam_check; \
			$(SAMTOOLS) idxstats $< >>$(@D)/$(*F).bam_check; \
			#exit 1; \
		else \
			# OK \
			echo "# Number of reads OK between $< and $$(cat $@.sample_pattern).R1.fastq.gz $$(cat $@.sample_pattern).R2.fastq.gz"  >>$(@D)/$(*F).bam_check; \
		fi; \
	else \
		# Error \
		echo "#[ERROR] No $$(cat $@.sample_pattern).R1.fastq.gz or $$(cat $@.sample_pattern).R2.fastq.gz" >>$(@D)/$(*F).bam_check; \
	fi;
	cat $(@D)/$(*F).bam_check;
	echo "#[$$(date)] BAM Check done" > $@ ;
	#-rm -f $@.sample_pattern



# BAM CHECK rule for unaligned BAM (deprecated)
#################################################

%.bam.metrics/metrics.bam_check_fromUBAM: %.bam %.bam.bai
	# CHECK NUMBER of READS in BAM \
	# Create directory if not created
	mkdir -p $(@D);
	# Pattern
	echo $$(dirname $*)/$$(basename $$(dirname $*)) > $@.sample_pattern;
	# Header
	echo "# Check number of reads between $< and $$(cat $@.sample_pattern).unaligned.bam to ensure no reads removing (only mark reads)" >$(@D)/$(*F).bam_check;
	echo "# (possible new reads marked as secondary aligned)" >>$(@D)/$(*F).bam_check;
	# Check
	if [ -s $$(cat $@.sample_pattern).unaligned.bam ] && [ -s $$(cat $@.sample_pattern).unaligned.bam.bai ]; then \
		echo "# Files $$(cat $@.sample_pattern).unaligned.bam and $$(cat $@.sample_pattern).unaligned.bam.bai exists" >>$(@D)/$(*F).bam_check; \
		# Compare number of reads \
		#grep "Number of reads for /home1/IRC/DATA/DEV/RES/ALL/160923_M01656_0124_000000000-ATJUM/T_Horizon/T_Horizon.unaligned.bam" /home1/IRC/DATA/DEV/RES/ALL/160923_M01656_0124_000000000-ATJUM/T_Horizon/T_Horizon.bwamem.bam.metrics/T_Horizon.bwamem.bam_check | cut -d" " -f2 \
		echo "# $$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $<) reads for $<" >>$(@D)/$(*F).bam_check; \
		echo "# $$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $$(cat $@.sample_pattern).unaligned.bam) reads for $$(cat $@.sample_pattern).unaligned.bam" >>$(@D)/$(*F).bam_check; \
		grep ".* reads for $<" $(@D)/$(*F).bam_check | cut -d" " -f2; \
		grep ".* reads for $$(cat $@.sample_pattern).unaligned.bam" $(@D)/$(*F).bam_check | cut -d" " -f2; \
		if [ "$$(grep ".* reads for $<" $(@D)/$(*F).bam_check | cut -d" " -f2)" != "$$(grep ".* reads for $$(cat $@.sample_pattern).unaligned.bam" $(@D)/$(*F).bam_check | cut -d" " -f2)" ]; then \
			# Error \
			echo "#[ERROR] KO Number of reads between $< and $$(cat $@.sample_pattern).unaligned.bam !!!" >>$(@D)/$(*F).bam_check; \
			echo "# BCFTOOLS STATS for $<" >>$(@D)/$(*F).bam_check;  \
			$(SAMTOOLS) stats $< | grep ^SN >>$(@D)/$(*F).bam_check; \
			$(SAMTOOLS) idxstats $< >>$(@D)/$(*F).bam_check; \
			echo "# BCFTOOLS STATS for $$(cat $@.sample_pattern).unaligned.bam" >>$(@D)/$(*F).bam_check;  \
			$(SAMTOOLS) stats $$(cat $@.sample_pattern).unaligned.bam | grep ^SN >>$(@D)/$(*F).bam_check; \
			$(SAMTOOLS) idxstats $$(cat $@.sample_pattern).unaligned.bam >>$(@D)/$(*F).bam_check; \
			#exit 1; \
		else \
			# OK \
			echo "# Number of reads OK between $< and $$(cat $@.sample_pattern).unaligned.bam"  >>$(@D)/$(*F).bam_check; \
		fi; \
	else \
		# Error \
		echo "#[ERROR] No $$(cat $@.sample_pattern).unaligned.bam or $$(cat $@.sample_pattern).unaligned.bam.bai" >>$(@D)/$(*F).bam_check; \
	fi;
	cat $(@D)/$(*F).bam_check;
	echo "#[$$(date)] BAM Check done" > $@ ;
	-rm $@.sample_pattern


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# MAIN RULES '$(MK_RELEASE)' : basicaly to manage VCF, BAM, FASTQ, METRICS, MANIFEST, INTERVAL, BED... using SAMTOOLS, FASTQC GATK, PICARD, FATBAM, BWA, TABIX, IGV..."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )
