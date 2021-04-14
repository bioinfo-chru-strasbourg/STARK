############################
# Metrics Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.5.1b"
MK_DATE="27/09/2019"

# Release note
# 18/12/2015 - 0.9.2b : Force gzip metrics files
# 23/09/2016 - 0.9.3b : Add BAM Check metrics.bam_check. Change PICARD version
# 29/09/2016 - 0.9.4b : Add Amplicon coverage metrics.amplicon_coverage
# 29/09/2016 - 0.9.5b : Chenge metrics.genes rule
# 27/09/2019 - 0.9.5.1b: Change FATBAM to CAP tool, add HOWARD option




# OPTIONS
BAM_METRICS?=1
FULL_COVERAGE?=0
CAP_TMP_FOLDER?=$(TMP_FOLDER_TMP)
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
METRICS_FLAGS?=UNMAP,SECONDARY,QCFAIL,DUP

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
#echo $COVERAGE_CRITERIA | tr "," "\n" | sort -k1n -u | tail -n1
#SAMTOOLS_METRICS_MPILEUP_DEPTH_max_depth?=$(shell echo $(COVERAGE_CRITERIA) | tr " " "\n" | tr "," "\n" | sort -k1n -u | tail -n1)
SAMTOOLS_METRICS_MPILEUP_DEPTH_adjust_MQ?=0
SAMTOOLS_METRICS_MPILEUP_FLAGS?=$(METRICS_FLAGS)
SAMTOOLS_METRICS_MPILEUP_DEPTH_excl_flags?=$(shell if [ "$(SAMTOOLS_METRICS_MPILEUP_FLAGS)" == "" ]; then echo ""; else echo "--excl-flags $(SAMTOOLS_METRICS_MPILEUP_FLAGS)"; fi;)
SAMTOOLS_METRICS_MPILEUP_DEPTH_count_orphans?=--count-orphans
SAMTOOLS_METRICS_MPILEUP_PARAM?= $(SAMTOOLS_METRICS_MPILEUP_DEPTH_excl_flags) --min-BQ $(SAMTOOLS_METRICS_MPILEUP_DEPTH_base_q) --min-MQ $(SAMTOOLS_METRICS_MPILEUP_DEPTH_map_q) --max-depth $(SAMTOOLS_METRICS_MPILEUP_DEPTH_max_depth) --adjust-MQ $(SAMTOOLS_METRICS_MPILEUP_DEPTH_adjust_MQ) $(SAMTOOLS_METRICS_MPILEUP_DEPTH_count_orphans) $(shell if ! (( $(CLIP_OVERLAPPING_READS) )); then echo " -x "; fi )

# PICARD
PICARD_CollectHsMetrics_MINIMUM_MAPPING_QUALITY?=$(METRICS_MINIMUM_MAPPING_QUALITY)
PICARD_CollectHsMetrics_MINIMUM_BASE_QUALITY?=$(METRICS_MINIMUM_BASE_QUALITY)
PICARD_CollectHsMetrics_SAMPLE_SIZE?=10000
#PICARD_CollectHsMetrics_TMP_DIR=--TMP_DIR $(TMP_FOLDER_TMP)/PICARD_CollectHsMetrics_TMP_DIR
#PICARD_CollectHsMetrics_TMP_DIR=--TMP_DIR $(TMP_FOLDER_TMP)
PICARD_CollectHsMetrics_TMP_DIR=
#$(SAMTOOLS_METRICS_MPILEUP_DEPTH_max_depth)
#PICARD_CollectHsMetrics_PARAM?=MINIMUM_MAPPING_QUALITY=$(PICARD_CollectHsMetrics_MINIMUM_MAPPING_QUALITY) SAMPLE_SIZE=$(PICARD_CollectHsMetrics_SAMPLE_SIZE) MINIMUM_BASE_QUALITY=$(PICARD_CollectHsMetrics_MINIMUM_BASE_QUALITY) $(PICARD_CollectHsMetrics_TMP_DIR) CLIP_OVERLAPPING_READS=$(shell if (( $(CLIP_OVERLAPPING_READS) )); then echo "true"; else echo "false"; fi )
PICARD_CollectHsMetrics_PARAM?=-MINIMUM_MAPPING_QUALITY $(PICARD_CollectHsMetrics_MINIMUM_MAPPING_QUALITY) -SAMPLE_SIZE $(PICARD_CollectHsMetrics_SAMPLE_SIZE) -MINIMUM_BASE_QUALITY $(PICARD_CollectHsMetrics_MINIMUM_BASE_QUALITY) $(PICARD_CollectHsMetrics_TMP_DIR) -CLIP_OVERLAPPING_READS $(shell if (( $(CLIP_OVERLAPPING_READS) )); then echo "true"; else echo "false"; fi )

# CAP
CAP_METRICS_OPTIONS_CLIP_OVERLAPPING_READS?=$(shell if (( $(CLIP_OVERLAPPING_READS) )); then echo " --clip_overlapping_reads "; fi )
CAP_METRICS_OPTIONS_HSMETRICS_PARAMETERS?=--hsmetrics_parameters="MINIMUM_MAPPING_QUALITY=$(PICARD_CollectHsMetrics_MINIMUM_MAPPING_QUALITY);SAMPLE_SIZE=$(PICARD_CollectHsMetrics_SAMPLE_SIZE);MINIMUM_BASE_QUALITY=$(PICARD_CollectHsMetrics_MINIMUM_BASE_QUALITY)"
CAP_METRICS_OPTIONS?=$(CAP_METRICS_OPTIONS_CLIP_OVERLAPPING_READS) $(CAP_METRICS_OPTIONS_HSMETRICS_PARAMETERS)


GENESCOVERAGE_PRECISION?=2


################################
## BED / INTERVALS for Metrics #
################################


# BED from BAM
################

%.bam.bed: %.bam %.bam.bai
	#BAM.BED from BAM
	# samtools view P1335.bwamem.bam -b | /STARK/tools/bedtools/current/bin/bedtools genomecov -ibam stdin -bg | /STARK/tools/bedtools/current/bin/bedtools merge -i stdin
	#$(BEDTOOLS) bamtobed -i $< | $(BEDTOOLS) merge -i - | awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t+\t"$$1":"$$2"-"$$3}' > $@
	$(BEDTOOLS) bamtobed -i $< | $(BEDTOOLS) merge -i - | awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t"$$1":"$$2"-"$$3}' > $@



# DESIGN BED
##############
# BED for metrics (either a provided BED or generated from BAM)
# without header

%.design.bed: %.withoutheader.for_metrics_bed
	cp $< $@;



# BAM for metrics
###################

%.metrics.bam: %.bam %.bam.bai
	-ln -P $< $@;
	if [ ! -e $@ ]; then cp $< $@; fi; # if ln des not work



# BED for metrics
###################
# BED for metrics (either a provided BED or generated from BAM)
# Check if .metrics.region_clipped.bed and metrics.bed exists
# Else generate BED from BAM

%.for_metrics_bed: %.bam %.bam.bai %.bam.bed %.metrics.bed %.metrics.region_clipped.bed #%.bam.bed
	+if [ "`grep ^ -c $*.metrics.bed`" == "0" ] && [ "`grep ^ -c $*.metrics.region_clipped.bed`" == "0" ]; then \
		echo "# No BED file provided... BED from the BAM file"; \
		cp $*.bam.bed $@; \
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
		#cat $*.withoutheader.for_metrics_bed.3fields.bed | awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t+\t"$$1":"$$2"-"$$3}' >> $@; \
		cat $*.withoutheader.for_metrics_bed.3fields.bed | awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t"$$1":"$$2"-"$$3}' >> $@; \
		rm $*.withoutheader.for_metrics_bed.3fields.bed; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi



#############
## METRICS ##
#############



# BAM Validation
##################

BAM_VALIDATION_COMPRESSION?=4

%.validation.bam: %.bam %.bam.bai #%.list.genes %.design.bed
	# Create directory ;
	mkdir -p $(@D);
	# BAM Validation
	# samtools view -F 1284 F10.bwamem.bam
	$(SAMTOOLS) view $(SAMTOOLS_METRICS_VIEW_PARAM) -h $< -O BAM,level=$(BAM_VALIDATION_COMPRESSION) -@ $(THREADS) > $@ ;




# ALL METRICS
###############

#%.bam.metrics/metrics.depthbed

%.bam.metrics/metrics: %.bam.metrics/metrics.design %.bam.metrics/metrics.gatk %.bam.metrics/metrics.picard %.bam.metrics/metrics.samtools %.bam.metrics/metrics.regions_coverage %.bam.metrics/metrics.per_amplicon_coverage %.bam.metrics/metrics.post_alignment #%.bam.metrics/metrics.bam_check
#%.bam.metrics/metrics: %.bam.metrics/metrics.depthbed
	#Create directory
	mkdir -p $(@D)
	cat $^ > $@
	echo "#[INFO] BAM Metrics done" >> $@
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

#%.bam.metrics/metrics.markDuplicates: %.bam %.bam.bai
#	-cp -f $**.markduplicates.bam.metrics/* $(@D)/;
#	rm -rf $**.markduplicates.bam.metrics;
#	echo "#[INFO] MarkDuplicates metrics done. " > $@;


# UMIgroup metrics
##########################

#%.bam.metrics/metrics.UMIgroup: %.bam %.bam.bai
#	-cp -f $**.UMIgroup.bam.metrics/* $(@D)/;
#	rm -rf $**.UMIgroup.bam.metrics;
#	echo "#[INFO] UMIgroup metrics done. " > $@;


# Post align metrics
##########################

%.bam.metrics/metrics.post_alignment: %.bam %.bam.bai
	mkdir -p $(@D) ;
	-cp -f $*.*.bam.metrics/* $(@D)/;
	rm -rf $*.*.bam.metrics;
	echo "#[INFO] POST ALIGNMENT metrics done. " > $@;




# Amplicon metrics
####################
# From CAP

#%.bam.metrics/metrics.amplicon_coverage: %.bam %.bam.bai %.manifest %.genome
%.bam.metrics/metrics.per_amplicon_coverage: %.validation.bam %.validation.bam.bai %.manifest %.genome
	mkdir -p $(@D) ;
	+$(CAP) --function=coverage --env=$(CONFIG_TOOLS) --ref=$$(cat $*.genome) --bam=$< --output=$(@D)/$(*F).HsMetrics.per_amplicon_coverage.tmp --manifest=$*.manifest --threads=$(THREADS) $(CAP_METRICS_OPTIONS) --bedtools=$(BEDTOOLS) --samtools=$(SAMTOOLS) --picard=$(PICARD) --verbose --tmp=$(CAP_TMP_FOLDER) 1>$(@D)/$(*F).HsMetrics.per_amplicon_coverage.log 2>$(@D)/$(*F).HsMetrics.per_amplicon_coverage.err;
	awk -f $(STARK_FOLDER_BIN)/per_target_coverage_flag.awk -F"\t" -v EXPECTED_DEPTH=$(EXPECTED_DEPTH) -v MINIMUM_DEPTH=$(MINIMUM_DEPTH) $(@D)/$(*F).HsMetrics.per_amplicon_coverage.tmp > $(@D)/$(*F).HsMetrics.per_amplicon_coverage.flags; \
	rm -f $(@D)/$(*F).HsMetrics.per_amplicon_coverage.tmp; \
	#cat $(@D)/$(*F).amplicon_coverage.log $(@D)/$(*F).HsMetrics.per_amplicon_coverage.err;
	echo "#[INFO] BAM Amplicon Coverage Metrics done" > $@;


#%.validation.bam %.validation.bam.bai

# GATK metrics
################

GATKDOC_FLAGS= -rf BadCigar -allowPotentiallyMisencodedQuals
%.bam.metrics/metrics.gatk: %.validation.bam %.bam.bai %.genome %.for_metrics_bed %.3fields.for_metrics_bed
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
		echo "#[INFO] BAM GATK Metrics done" > $@; \
	else \
		echo "#[INFO] BAM GATK not done because BAM_METRICS=0" > $@; \
	fi;





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
			bed_subname="Design"; \
			[ "$$one_bed" != "$*.design.bed" ] && bed_subname="Panel."$$(basename $$one_bed); \
			# 4fields file \
			awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t"$$4}' $$one_bed > $(@D)/$(*F).$$(basename $$one_bed).4fields.tmp ; \
			# Clean bed with dict contig\
			grep -Po 'SN:([^\t]*)' $$(cat $*.dict) | cut -d: -f2 | sed "s/^/^/gi" | sed "s/$$/\t/gi" > $(@D)/$(*F).$$(basename $$one_bed).4fields.contig_from_dict ; \
			grep -f $(@D)/$(*F).$$(basename $$one_bed).4fields.contig_from_dict $(@D)/$(*F).$$(basename $$one_bed).4fields.tmp > $(@D)/$(*F).$$(basename $$one_bed).4fields ; \
			# BedToIntervalList \
			$(JAVA) $(JAVA_FLAGS_BY_SAMPLE) -jar $(PICARD) BedToIntervalList -I $(@D)/$(*F).$$(basename $$one_bed).4fields -O $(@D)/$(*F).$$(basename $$one_bed).interval -SD $$(cat $*.dict); \
			$(JAVA) $(JAVA_FLAGS_BY_SAMPLE) -jar $(PICARD) CollectHsMetrics -INPUT $*.validation.bam -OUTPUT $(@D)/$(*F).$$(basename $$one_bed).HsMetrics -R $$(cat $*.genome) -BAIT_INTERVALS $(@D)/$(*F).$$(basename $$one_bed).interval -TARGET_INTERVALS $(@D)/$(*F).$$(basename $$one_bed).interval -PER_TARGET_COVERAGE $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_target_coverage.tmp -PER_BASE_COVERAGE $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.tmp $(PICARD_CollectHsMetrics_PARAM) -VALIDATION_STRINGENCY SILENT 2>$(@D)/$(*F).$$(basename $$one_bed).HsMetrics.err ; \
			# If bed empty just touch \
			touch $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_target_coverage.tmp ; \
			# Flag HsMetrics per_target_coverage \
			awk -f $(STARK_FOLDER_BIN)/per_target_coverage_flag.awk -F"\t" -v EXPECTED_DEPTH=$(EXPECTED_DEPTH) -v MINIMUM_DEPTH=$(MINIMUM_DEPTH) $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_target_coverage.tmp > $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_target_coverage.flags; \
			# Flag HsMetrics per_target_coverage by FLAG \
			echo "#chrom	start	stop	target	mean	min	max	count" > $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.MISS.tsv; \
			cat $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.tmp | awk -F"\t" -v MISS=$(SEQUENCING_DEPTH) '($$4<MISS) {print $$1 "\t" $$2-1 "\t" $$2 "\t" $$3 "\t" $$4}' | $(BEDTOOLS) merge -c 4,5,5,5,5 -o distinct,mean,min,max,count >> $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.MISS.tsv 2>/dev/null; \
			echo "#chrom	start	stop	target	mean	min	max	count" > $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.FAIL.tsv; \
			cat $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.tmp | awk -F"\t" -v MISS=$(SEQUENCING_DEPTH) -v FAIL=$(MINIMUM_DEPTH) '($$4>=MISS && $$4<FAIL) {print $$1 "\t" $$2-1 "\t" $$2 "\t" $$3 "\t" $$4}' | $(BEDTOOLS) merge -c 4,5,5,5,5 -o distinct,mean,min,max,count >> $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.FAIL.tsv 2>/dev/null; \
			echo "#chrom	start	stop	target	mean	min	max	count" > $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.WARN.tsv; \
			cat $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.tmp | awk -F"\t" -v FAIL=$(MINIMUM_DEPTH) -v WARN=$(EXPECTED_DEPTH) '($$4>=FAIL && $$4<WARN) {print $$1 "\t" $$2-1 "\t" $$2 "\t" $$3 "\t" $$4}' | $(BEDTOOLS) merge -c 4,5,5,5,5 -o distinct,mean,min,max,count >> $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.WARN.tsv 2>/dev/null; \
			echo "#chrom	start	stop	target	mean	min	max	count" > $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.PASS.tsv; \
			cat $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.tmp | awk -F"\t" -v PASS=$(EXPECTED_DEPTH) '($$4>=PASS) {print $$1 "\t" $$2-1 "\t" $$2 "\t" $$3 "\t" $$4}' | $(BEDTOOLS) merge -c 4,5,5,5,5 -o distinct,mean,min,max,count >> $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.PASS.tsv 2>/dev/null; \
			# Flag HsMetrics per_target_coverage in one BED file \
			cat $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.tmp | awk -F"\t" -v MISS=$(SEQUENCING_DEPTH) -v FAIL=$(MINIMUM_DEPTH) -v WARN=$(EXPECTED_DEPTH) '($$4<MISS) {print $$1 "\t" $$2-1 "\t" $$2 "\t" $$3 "\t" $$4 "\t'$(MISS_COLOR_RGB)'\t+" }' | $(BEDTOOLS) merge -c 4,5,7,2,3,6 -o distinct,mean,first,first,last,distinct >> $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.bed.tmp 2>/dev/null; \
			cat $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.tmp | awk -F"\t" -v MISS=$(SEQUENCING_DEPTH) -v FAIL=$(MINIMUM_DEPTH) -v WARN=$(EXPECTED_DEPTH) '($$4>=MISS && $$4<FAIL) {print $$1 "\t" $$2-1 "\t" $$2 "\t" $$3 "\t" $$4 "\t'$(FAIL_COLOR_RGB)'\t+" }' | $(BEDTOOLS) merge -c 4,5,7,2,3,6 -o distinct,mean,first,first,last,distinct >> $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.bed.tmp 2>/dev/null; \
			cat $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.tmp | awk -F"\t" -v MISS=$(SEQUENCING_DEPTH) -v FAIL=$(MINIMUM_DEPTH) -v WARN=$(EXPECTED_DEPTH) '($$4>=FAIL && $$4<WARN) {print $$1 "\t" $$2-1 "\t" $$2 "\t" $$3 "\t" $$4 "\t'$(WARN_COLOR_RGB)'\t+" }' | $(BEDTOOLS) merge -c 4,5,7,2,3,6 -o distinct,mean,first,first,last,distinct >> $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.bed.tmp 2>/dev/null; \
			cat $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.tmp | awk -F"\t" -v MISS=$(SEQUENCING_DEPTH) -v FAIL=$(MINIMUM_DEPTH) -v WARN=$(EXPECTED_DEPTH) '($$4>=WARN) {print $$1 "\t" $$2-1 "\t" $$2 "\t" $$3 "\t" $$4 "\t'$(PASS_COLOR_RGB)'\t+" }' | $(BEDTOOLS) merge -c 4,5,7,2,3,6 -o distinct,mean,first,first,last,distinct >> $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.bed.tmp 2>/dev/null; \
			$(BEDTOOLS) sort -i $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.bed.tmp | cut -f1-9 > $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.bed; \
			# Copy Design bed \
			cp -p $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.bed $(@D)/$(*F).validation.flags.$$bed_subname.bed; \
			# HsMetrics per_target_coverage compression file file \
			$(GZ) -c $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.tmp > $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.gz; \
			rm $(@D)/$(*F).$$(basename $$one_bed).4fields* $(@D)/$(*F).$$(basename $$one_bed).interval $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_target_coverage.tmp $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.tmp rm -f $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.bed.tmp; \
		else \
			# BED empty \
			cp $*.empty.HsMetrics $(@D)/$(*F).$$(basename $$one_bed).HsMetrics; \
			echo "#[WARNING] BAM PICARD Metrics warning. Empty HsMetrics file generated. See '$(@D)/$(*F).$$(basename $$one_bed).HsMetrics.err'" >> $@; \
		fi; \
		if [ ! -s $(@D)/$(*F).$$(basename $$one_bed).HsMetrics ]; then \
			#cp $*.empty.HsMetrics $(@D)/$(*F).$$(basename $$one_bed).HsMetrics; \
			echo "#[ERROR] BAM PICARD Metrics failed. See '$(@D)/$(*F).$$(basename $$one_bed).HsMetrics.err'" >> $@; \
			# Exit to fail rule \
			exit_0; \
		else \
			echo "#[INFO] BAM PICARD Metrics done. HsMetrics file generated. See '$(@D)/$(*F).$$(basename $$one_bed).HsMetrics'" >> $@; \
		fi \
	done;
	echo "#[INFO] BAM PICARD Metrics done" >> $@


#echo "#chrom	start	stop	target	mean	min	max	count" > $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.MISS.bed; \
			cat $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.tmp | awk -F"\t" -v MISS=$(SEQUENCING_DEPTH) '($$4<MISS) {print $$1 "\t" $$2-1 "\t" $$2 "\t" $$3 "\t" $$4}' | $(BEDTOOLS) merge -c 4,5,5,5,5 -o distinct,mean,min,max,count >> $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.MISS.bed 2>/dev/null; \
			echo "#chrom	start	stop	target	mean	min	max	count" > $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.FAIL.bed; \
			cat $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.tmp | awk -F"\t" -v MISS=$(SEQUENCING_DEPTH) -v FAIL=$(MINIMUM_DEPTH) '($$4>=MISS && $$4<FAIL) {print $$1 "\t" $$2-1 "\t" $$2 "\t" $$3 "\t" $$4}' | $(BEDTOOLS) merge -c 4,5,5,5,5 -o distinct,mean,min,max,count >> $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.FAIL.bed 2>/dev/null; \
			echo "#chrom	start	stop	target	mean	min	max	count" > $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.WARN.bed; \
			cat $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.tmp | awk -F"\t" -v FAIL=$(MINIMUM_DEPTH) -v WARN=$(EXPECTED_DEPTH) '($$4>=FAIL && $$4<WARN) {print $$1 "\t" $$2-1 "\t" $$2 "\t" $$3 "\t" $$4}' | $(BEDTOOLS) merge -c 4,5,5,5,5 -o distinct,mean,min,max,count >> $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.WARN.bed 2>/dev/null; \
			echo "#chrom	start	stop	target	mean	min	max	count" > $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.PASS.bed; \
			cat $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.tmp | awk -F"\t" -v PASS=$(EXPECTED_DEPTH) '($$4>=PASS) {print $$1 "\t" $$2-1 "\t" $$2 "\t" $$3 "\t" $$4}' | $(BEDTOOLS) merge -c 4,5,5,5,5 -o distinct,mean,min,max,count >> $(@D)/$(*F).$$(basename $$one_bed).HsMetrics.per_base_coverage.PASS.bed 2>/dev/null; \


# SAMTOOLS metrics
####################


%.bam.metrics/metrics.samtools: %.bam.metrics/metrics.samtools.flagstat %.bam.metrics/metrics.samtools.idxstats %.bam.metrics/metrics.samtools.genomeCoverage %.bam.metrics/metrics.samtools.depth \
	%.bam.metrics/metrics.samtools.on.target %.bam.metrics/metrics.samtools.off.target %.bam.metrics/metrics.coverage #%.bam.metrics/metrics.samtools.genomeCoverageBed %.bam.metrics/metrics.samtools.depthbed.coverage
	# Create directory ;
	mkdir -p $(@D);
	cat $^ > $@
	echo "#[INFO] SAMTOOLS Metrics done" >> $@


%.bam.metrics/metrics.samtools.flagstat: %.bam %.bam.bai %.validation.bam %.validation.bam.bai
	# Create directory ;
	mkdir -p $(@D);
	> $@;
	$(SAMTOOLS) flagstat $*.bam > $(@D)/$(*F).flagstat ;
	$(SAMTOOLS) flagstat $*.validation.bam > $(@D)/$(*F).validation.flagstat ;
	if [ -s $(@D)/$(*F).flagstat ]; then \
		echo "#[INFO] SAMTOOLS flagstat done. See '$(@D)/$(*F).flagstat'. " >> $@; \
	else \
		echo "#[INFO] SAMTOOLS flagstat failed. " >> $@; \
	fi;


%.bam.metrics/metrics.samtools.idxstats: %.bam %.bam.bai %.validation.bam %.validation.bam.bai
	# Create directory ;
	mkdir -p $(@D);
	> $@;
	$(SAMTOOLS) idxstats $*.bam > $(@D)/$(*F).idxstats ;
	$(SAMTOOLS) idxstats $*.validation.bam > $(@D)/$(*F).validation.idxstats ;
	if [ -z $(@D)/$(*F).idxstats ]; then \
		echo "#[INFO] SAMTOOLS idxstats done. See '$(@D)/$(*F).idxstats'. " >> $@; \
	else \
		echo "#[INFO] SAMTOOLS idxstats failed. " >> $@; \
	fi;


%.bam.metrics/metrics.samtools.genomeCoverage: %.bam %.bam.bai %.validation.bam %.validation.bam.bai
	# Create directory ;
	mkdir -p $(@D);
	> $@;
	if (($(BAM_METRICS))); then \
		if (($(FULL_COVERAGE))); then \
			$(BEDTOOLS) genomecov -ibam $*.bam -dz > $(@D)/$(*F).genomeCoverage; \
			$(GZ) --best -f $(@D)/$(*F).genomeCoverage; \
			$(BEDTOOLS) genomecov -ibam $*.validation.bam -dz > $(@D)/$(*F).validation.genomeCoverage; \
			$(GZ) --best -f $(@D)/$(*F).validation.genomeCoverage; \
			if [ -s $(@D)/$(*F).genomeCoverage.gz ] && [ -s $(@D)/$(*F).validation.genomeCoverage.gz ]; then \
				echo "#[INFO] BEDTOOLS genomeCoverage done. See '$(@D)/$(*F).genomeCoverage.gz' and '$(@D)/$(*F).validation.genomeCoverage.gz'. " >> $@; \
			else \
				echo "#[INFO] SAMTOOLS genomeCoverage failed. " >> $@; \
			fi; \
		else \
			echo "#[INFO] SAMTOOLS genomeCoverage not done because FULL_COVERAGE=0. " >> $@; \
		fi; \
	else \
		echo "#[INFO] SAMTOOLS genomeCoverage not done because BAM_METRICS=0. " >> $@; \
	fi;


%.bam.metrics/metrics.samtools.depth: %.validation.bam %.bam.bai
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



%.bam.metrics/metrics.coverage: %.validation.bam %.validation.bam.bai %.list.genes %.design.bed %.bam.metrics/metrics.picard #%.bam.metrics/metrics.depthbed
	# Create directory ;
	mkdir -p $(@D);
	> $@;
	# For each BED GENES
	for one_bed in $$(cat $*.list.genes) $*.design.bed; do \
		# Test if BED exists \
		bedfile_name=$$( basename $$one_bed ); \
		if [ -e $$one_bed ]; then \
			echo -e "#Depth\tCoveredBases\tTotalBases\tPercent" > $(@D)/$(*F).$$(basename $$one_bed).coverage; \
			#$(UNGZ) -c $(@D)/$(*F).$$(basename $$one_bed).depthbed.gz | \
			$(UNGZ) -c $(@D)/$(*F).$$bedfile_name.HsMetrics.per_base_coverage.gz | awk 'NR!=1{print $$4"\t"$$3}' | \
			awk -v MDP=$$(echo $(COVERAGE_CRITERIA) | tr "," "\n" | sort -n | tail -n 1)  '{SUM++} { if ($$1>MDP) {DP[MDP]++} else {DP[$$1]++} } END { for (i=MDP; i>=0; i-=1) {print i" "DP[i]" SUM"SUM}}' | \
			sort -g -r | awk -v COVERAGE_CRITERIA=$(COVERAGE_CRITERIA) '{SUM+=$$2} {CUM[$$1]=SUM} {split(COVERAGE_CRITERIA,C,",")} END { for (j in C) {print C[j]"X\t"CUM[C[j]]"\t"SUM"\t"(CUM[C[j]]/SUM)} }' | \
			sort -g >> $(@D)/$(*F).$$bedfile_name.coverage; \
			echo "#[INFO] SAMTOOLS depthbed and coverage with '"$$bedfile_name"' done. See '$(@D)/$(*F).$$bedfile_name.depthbed' and  '$(@D)/$(*F).$$bedfile_name.coverage'. " >> $@; \
			if [ -s $(@D)/$(*F).$$bedfile_name.coverage ]; then \
				echo "#[INFO] SAMTOOLS coverage with '"$$bedfile_name"' done. See '$(@D)/$(*F).$$bedfile_name.coverage'. " >> $@; \
			else \
				echo "#[INFO] SAMTOOLS coverage with '"$$bedfile_name"' failed. " >> $@; \
			fi; \
		else \
			echo "#[INFO] SAMTOOLS depthbed and coverage with '"$$bedfile_name"' failed, because no '"$$bedfile_name"'." >> $@; \
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
			$(BEDTOOLS) intersect -abam $*.validation.bam -b $$one_bed -v | $(SAMTOOLS) view -c - > $(@D)/$(*F).$$(basename $$one_bed).off.target; \
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
				$(SAMTOOLS) view -b $*.validation.bam -L $$one_bed | $(SAMTOOLS) sort | $(BEDTOOLS) genomecov -ibam stdin -dz > $(@D)/$(*F).$$(basename $$one_bed).genomeCoverageBedbed; \
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

%.sequencing/metrics: %.sequencing/metrics.infos # %.sequencing/metrics.Q30 %.sequencing/metrics.fastqc %.sequencing/metrics.counts
	# create directory
	mkdir -p $(@D)
	cat $^ > $@
	-rm -f $^



# FastQC metrics
##################

# %.sequencing/metrics.fastqc: %.fastq.gz
# 	# create directory
# 	mkdir -p $(@D)
# 	# create link
# 	ln -s $< $(@D)/$(*F).fastq.gz
# 	# touch target
# 	touch $@;
# 	# FASTQC
# 	-if (($$(zcat $< | head -n 1 | wc -l))); then \
# 		#$(FASTQC) $< --outdir=$(@D) --casava --extract; \
# 		$(FASTQC) $(@D)/$(*F).fastq.gz --outdir=$(@D) --casava --extract --threads $(THREADS_BY_SAMPLE) ; \
# 		cp $(@D)/$(*F)_fastqc/fastqc_data.txt $(@D)/metrics.fastqc.txt; \
# 		echo "#[INFO] FASTQC done. See 'metrics.fastqc.txt' file." >> $@; \
# 	else \
# 		echo "#[ERROR] FASTQC can't be launched. No Reads in FASTQ file '$<'" >> $@; \
# 	fi;
# 	# create link
# 	-rm -f $(@D)/$(*F).fastq.gz



# Q30
########


# %.sequencing/metrics.Q30: %.sequencing
# 	cat $(@D)/*.fastp.json | python -c "import sys, json; print json.load(sys.stdin)['summary']['after_filtering']['q30_rate']" > $@.txt
# 	echo "#[INFO] Q30 calculation done. See 'metrics.Q30.txt' file." > $@;




# Reads Count
###############

# #%.fastqc/metrics.counts: %.fastq.gz #%.fastqc/metrics.fastqc
# %.sequencing/metrics.counts: %.fastq.gz #%.fastqc/metrics.fastqc
# 	# create directory
# 	mkdir -p $(@D)
# 	-if (($$(zcat $< | head -n 1 | wc -l))); then \
# 		# header \
# 		echo "total unique %unique maxRead count_maxRead %count_maxRead" > $@.tmp; \
# 		# Counts \
# 		if [ $$(zcat $< | head -n 1 | wc -l) ]; then \
# 			zcat $< | awk '{unique=0} ((NR-2)%4==0){read=$$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}' >> $@.tmp; \
# 		else \
# 			echo "0 0 - - - -" >> $@.tmp; \
# 		fi; \
# 		# transposition \
# 		awk '{ for (i=1; i<=NF; i++)  {a[NR,i] = $$i} } NF>p { p = NF } END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++){ str=str" "a[i,j]; } print str } }' $@.tmp | tr " " "\t" > $@.txt; \
# 		# nb bases \
# 		echo -e "count_bases\t"$$(zcat $< | paste - - - - | cut -f4 | wc -c) >> $@.txt; \
# 		# script \
# 		# echo "" > $*.fastqc/metrics.counts.txt; \
# 		echo "#[INFO] FASTQ Counts done. See 'metrics.counts.txt' file." > $@; \
# 	else \
# 		echo "#[ERROR] FASTQ Counts ERROR. No reads in '$<'..." > $@; \
# 	fi;
# 	-rm -f $@.tmp



# techno
########

%.sequencing/metrics.infos: %.manifest %.R2.fastq.gz
	# create directory
	mkdir -p $(@D)
	touch $@.txt
	echo -e $$((($$(zcat $*.R2.fastq.gz | head -n1 | wc -l))) && echo "mode\tPaired-End" || echo "mode\tSingle-End") >> $@.txt
	#echo -e $$((($$(grep "Amplicon Start" $*.manifest | grep -c "Upstream Probe Length"))) && echo "technology\tAmplicon" || echo "technology\tCapture") >> $@.txt
	echo -e $$((($$(grep -c "Upstream Probe Length\|ULSO Sequence" $*.manifest))) && echo "technology\tAmplicon" || echo "technology\tCapture") >> $@.txt
	echo "#[INFO] SEQUENCING INFOS done. See 'metrics.infos.txt' file." > $@;
	#-rm -f $@.*




# VCF METRICS
###############
# SNPEFF and BCFTOOLS metrics

# %.vcf.metrics/metrics: %.vcf.metrics/metrics.snpeff %.vcf.metrics/metrics.bcftools %.vcf.metrics/metrics.genes #%.vcf.metrics/metrics.info_field 
# 	cat $^ > $@



# SNPEFF Metrics
##################
# SNPEFF metrics through HOWARD

# %.vcf.metrics/metrics.snpeff: %.vcf
# 	mkdir -p $(@D);
# 	touch $@;
# 	if (($(METRICS_SNPEFF))); then \
# 		+$(HOWARD) --input=$< --output=$@.vcf --snpeff_stats=$@.html --annotation=null --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS)  --force; \
# 		echo "#[INFO] snpEff metrics done. See '$@.html' file." > $@; \
# 	else \
# 		echo "#[INFO] snpEff metrics NOT done." > $@; \
# 	fi;



# BCFTOOLS metrics
####################
# Stats from BCFTOOLS

# %.vcf.metrics/metrics.bcftools: %.vcf
# 	mkdir -p $(@D);
# 	touch $@;
# 	-$(BCFTOOLS) stats $< > $@.stats
# 	echo "#[INFO] BCFTOOLS metrics done. See '$@.stats' file." > $@;



# INFO stats
####################
# Stats of INFO field

# %.vcf.metrics/metrics.info_field: %.vcf
# 	mkdir -p $(@D);
# 	touch $@;
# 	grep -v "^#" $< | cut -f8 | tr ";" "\n" | sort | uniq -c | sed "s/^      / /gi" | tr "=" " " | awk '{print $$2"\t"$$3"\t"$$1} {a[$$2]+=$$1} END { for (key in a) { print "#\t" key "\t" a[key] } }' | sort > $@.stats
# 	echo "#[INFO] INFO field stats done. See '$@.stats' file." > $@;



# INFO stats
####################
# Stats of INFO field

# %.vcf.metrics/metrics.genesOLD: %.vcf.gz %.vcf.gz.tbi %.list.genes #%.design.bed
# 	mkdir -p $(@D);
# 	> $@;
# 	+for one_bed in $$(cat $*.list.genes) $*.design.bed; do \
# 		#cat "ONE_BED: "$$one_bed; \
# 		if [ -s $$one_bed ]; then \
# 			$(BCFTOOLS) view $*.vcf.gz -R $$one_bed | $(BCFTOOLS) norm --rm-dup exact > $@.$$(basename $$one_bed).vcf; \
# 		else \
# 			$(BCFTOOLS) view $*.vcf.gz | $(BCFTOOLS) norm --rm-dup exact > $@.$$(basename $$one_bed).vcf; \
# 		fi ; \
# 		$(HOWARD) --config=$(HOWARD_CONFIG) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION) --config_annotation=$(HOWARD_CONFIG_ANNOTATION) --pzfields="PZScore,PZFlag,PZComment,PZInfos" --translation=tab --fields="$(HOWARD_FIELDS)" --sort=$(HOWARD_SORT) --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)" --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS) --input=$@.$$(basename $$one_bed).vcf --output=$@.$$(basename $$one_bed).tsv --force; \
# 		$(HOWARD) --config=$(HOWARD_CONFIG) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION) --config_annotation=$(HOWARD_CONFIG_ANNOTATION) --pzfields="PZScore,PZFlag,PZComment,PZInfos" --translation=tab --fields="$(HOWARD_FIELDS_REPORT)" --sort=$(HOWARD_SORT) --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)" --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS) --input=$@.$$(basename $$one_bed).vcf --output=$@.$$(basename $$one_bed).report.tsv --stats=$@.$$(basename $$one_bed).report.info_field.stats --bcftools_stats=$@.$$(basename $$one_bed).report.bcftools.stats --force; \
# 		(($(METRICS_SNPEFF))) && $(HOWARD) --input=$@.$$(basename $$one_bed).vcf --output=$@.$$(basename $$one_bed).vcf.snpeff.vcf --snpeff_stats=$@.$$(basename $$one_bed).vcf.snpeff.html --annotation=null --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS)  --force ; \
# 		echo "#[INFO] VCF filtered by '$$one_bed' done. See '$@.$$(basename $$one_bed).vcf' file, '$@.$$(basename $$one_bed).tsv' file, and stats '$@.$$(basename $$one_bed).vcf.*'" >> $@; \
# 	done;


%.vcf.metrics/metrics: %.vcf.gz %.vcf.gz.tbi %.list.genes #%.design.bed
	mkdir -p $(@D);
	> $@;
	+for one_bed in $$(cat $*.list.genes) $*.design.bed; do \
		# "ONE_BED: "$$one_bed; \
		bed_subname="Design"; \
		[ "$$one_bed" != "$*.design.bed" ] && bed_subname="Panel."$$(basename $$one_bed); \
		if [ -s $$one_bed ] && [ "$$one_bed" != "$*.design.bed" ]; then \
			$(BCFTOOLS) view $*.vcf.gz -R $$one_bed | $(BCFTOOLS) norm --rm-dup exact > $@.$$bed_subname.vcf; \
		else \
			$(BCFTOOLS) view $*.vcf.gz | $(BCFTOOLS) norm --rm-dup exact > $@.$$bed_subname.vcf; \
		fi ; \
		# TSV \
		$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.$$bed_subname.vcf --output=$@.$$bed_subname.tsv --pzfields="PZScore,PZFlag,PZComment,PZInfos" --translation=TSV --fields="$(HOWARD_FIELDS)" --sort="$(HOWARD_SORT)" --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)" --stats=$@.$$bed_subname.info_field.stats.tsv --bcftools_stats=$@.$$bed_subname.bcftools.stats.tsv --force; \
		# TSV REPORT \
		$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.$$bed_subname.vcf --output=$@.$$bed_subname.report.tsv --pzfields="PZScore,PZFlag,PZComment,PZInfos" --translation=TSV --fields="$(HOWARD_FIELDS_REPORT)" --sort=$(HOWARD_SORT) --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)" --stats=$@.$$bed_subname.report.info_field.stats.tsv --bcftools_stats=$@.$$bed_subname.report.bcftools.stats.tsv --force; \
		# SNPEFF STATS \
		(($(METRICS_SNPEFF))) && $(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.$$bed_subname.vcf --output=$@.$$bed_subname.snpeff.vcf --snpeff_stats=$@.$$bed_subname.snpeff.html --annotation=null --force ; \
		if [ "$$(echo $* | rev | cut -d'.' -f 1 | rev)" == "final" ] || [ "$$(echo $* | rev | cut -d'.' -f 1 | rev)" == "full" ]; then \
			cp $@.$$bed_subname.tsv $*.$$bed_subname.tsv; \
			$(BGZIP) -c $@.$$bed_subname.vcf > $*.$$bed_subname.vcf.gz; \
		fi ; \
		echo "#[INFO] VCF filtered by '$$one_bed' done. See files '$@.$$bed_subname.vcf', '$@.$$bed_subname.tsv' and stats" >> $@; \
	done;
	



# List of Genes BED
#####################


%.bams.for_metrics_bed: $(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach ALIGNER,$(ALIGNERS),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ALIGNER).design.bed )) $(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach PIPELINE,$(PIPELINES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(call aligner,$(PIPELINE)).design.bed ))
	cat $^ | $(BEDTOOLS) sort -i - | $(BEDTOOLS) merge -i - | cut -f1,2,3 > $@
	echo "BAMS for_metrics_bed: $^"
	#cp $@ $*.test


%.list.genes: %.bed %.manifest %.bams.for_metrics_bed
	mkdir -p $(@D);
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
			echo "BEDFILE_GENES for $* generated from SAMPLE.bed, SAMPLE.manifest or from BAMs"; \
			bedfile_genes_list=`file=$$( basename $* | cut -d. -f1 ); echo "$$file.from_design.genes"`; \
			# Look for defined Design (bed or manifest) \
			if [ -s $*.bed ] ; then \
				#echo "genes from BED '$$bed'" >> $*.test; \
				cut -f1,2,3 $*.bed > $@.manifest.bed ; \
				bedfile_genes_list=`file=$$( basename $* | cut -d. -f1 ); echo "$$file.from_bed.genes"`; \
			#elif [ -s `echo "$$manifest"` ] ; then \
			elif [ -s $*.manifest ] ; then \
				# manifest to bed ; \
				$(CAP_ManifestToBED) --input "$*.manifest" --output "$@.manifest.bed.tmp" --output_type "region" --type=PCR; \
				cut -f1,2,3 $@.manifest.bed.tmp > $@.manifest.bed ; \
				bedfile_genes_list=`file=$$( basename $* | cut -d. -f1 ); echo "$$file.from_manifest.genes"`; \
				#rm $@.manifest.bed.tmp ; \
			elif [ -s $*.bams.for_metrics_bed ] ; then \
				cut -f1,2,3 $*.bams.for_metrics_bed > $@.manifest.bed ; \
				bedfile_genes_list=`file=$$( basename $* | cut -d. -f1 ); echo "$$file.from_alignments.genes"`; \
			fi; \
			# Generate BEDGENES \
			if [ -s $@.manifest.bed ] ; then \
				echo "MANIFEST.bed to bedfile_genes_list : $$bedfile_genes_list"; \
				$(BEDTOOLS) intersect -wb -a $@.manifest.bed -b $(REFSEQ_GENES) | cut -f7 | sort -u > $@.manifest.bed.intersect; \
				#sort -k5 $(REFSEQ_GENES) > $@.manifest.bed.refseq; \
				sort -k4 $(REFSEQ_GENES) > $@.manifest.bed.refseq; \
				#join -1 1 -2 5 $@.manifest.bed.intersect $@.manifest.bed.refseq -o 2.1,2.2,2.3,2.5 | sort -u -k1,2 | tr " " "\t" | awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t+\t"$$4}' > $$bedfile_genes_list; \
				join -1 1 -2 4 $@.manifest.bed.intersect $@.manifest.bed.refseq -o 2.1,2.2,2.3,2.4,2.5,2.6 | sort -u -k1,2 | tr " " "\t" | $(BEDTOOLS) sort | $(BEDTOOLS) merge -c 4,5,6 -o distinct,collapse,first | awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t0\t"$$6}' | $(STARK_BED_NORMALIZATION) > $$bedfile_genes_list; \
				#rm -f $*.manifest.bed.intersect $*.manifest.bed.refseq; \
			else \
				#echo "" > $$bedfile_genes_list ; \
				#echo "#[ERROR] Generating GENES failed. GENES file empty"; \
				echo "#[ERROR] Generating GENES failed."; \
			fi;\
			#bedfile_genes_list=$$(basename $$(echo $$bedfile_genes_list)); \
			rm -f $@.manifest*; \
		fi; \
		echo "bed_file list is : `echo $$bedfile_genes_list` "; \
		echo $$bedfile_genes_list | tr " " "\n" > $@ ; \
		#echo "List of genes list: "; cat $@ ; \
		#echo "List of genes: "; cat $$bedfile_genes_list ; \
	else \
		echo "BEDFILE_GENES exists!!! "; \
	fi;

# elif [ -s $*.manifest ] ; then \
# 				# manifest to bed ; \
# 				$(CAP_ManifestToBED) --input "$*.manifest" --output "$@.manifest.bed.tmp" --output_type "region" --type=PCR; \
# 				cut -f1,2,3 $@.manifest.bed.tmp > $(@D)/$(*F).from_manifest.bed ; \
# 				bedfile_genes_list=`file=$$( basename $* | cut -d. -f1 ); echo "$(@D)/$(*F).from_manifest.bed"`; \
# 				#rm $@.manifest.bed.tmp ; \
# 			elif [ -s $*.bams.for_metrics_bed ] ; then \
# 				cut -f1,2,3 $*.bams.for_metrics_bed > $(@D)/$(*F).from_alignments.bed ; \
# 				#bedfile_genes_list=`file=$$( basename $* | cut -d. -f1 ); echo "$$file.from_alignments.genes"`; \
# 				bedfile_genes_list=`file=$$( basename $* | cut -d. -f1 ); echo "$(@D)/$(*F).from_alignments.bed"`; \
# 			fi; \


# Regnions COVERAGE METRICS
#############################

%.bam.metrics/metrics.regions_coverage: %.validation.bam %.bam.bai %.list.genes %.design.bed %.bam.metrics/metrics.picard #%.bam.metrics/metrics.samtools.depthbed.coverage
	mkdir -p $(@D);
	touch $@;
	+if (( 1 )) ; then \
		if [ ! -z $*.list.genes ]; then \
			for one_bed in $$(cat $*.list.genes) $*.design.bed; \
			do \
				if [ -e $$one_bed ]; then \
					#bedfile_name=$$( basename $$one_bed | sed "s/\.genes$$//" ); \
					bedfile_name=$$( basename $$one_bed ); \
					#$(@D)/$(*F).$$bedfile_name.HsMetrics.per_base_coverage ; \
					$(UNGZ) -c $(@D)/$(*F).$$bedfile_name.HsMetrics.per_base_coverage.gz | awk 'NR!=1{print $$4"\t"$$3}' > $(@D)/$(*F).$$bedfile_name.HsMetrics.per_base_coverage.2cols ; \
					$(NGSscripts)/genesCoverage.sh -f $*.bam -b $$one_bed -c "$(COVERAGE_CRITERIA)" --coverage-bases=$(@D)/$(*F).$$bedfile_name.HsMetrics.per_base_coverage.2cols --dp_fail=$(MINIMUM_DEPTH) --dp_warn=$(EXPECTED_DEPTH) --dp_threshold=$(DEPTH_COVERAGE_THRESHOLD) --precision=$(GENESCOVERAGE_PRECISION) -n $(NB_BASES_AROUND) -t $(BEDTOOLS) -s $(SAMTOOLS) --threads=$(THREADS) -o $(@D)/$(*F).$$bedfile_name; \
					$(NGSscripts)/genesCoverage.sh -f $*.bam -b $$one_bed -c "1,$(MINIMUM_DEPTH),$(EXPECTED_DEPTH)" --coverage-bases=$(@D)/$(*F).$$bedfile_name.HsMetrics.per_base_coverage.2cols --dp_fail=$(MINIMUM_DEPTH) --dp_warn=$(EXPECTED_DEPTH) --dp_threshold=$(DEPTH_COVERAGE_THRESHOLD) --precision=$(GENESCOVERAGE_PRECISION) -n $(NB_BASES_AROUND) -t $(BEDTOOLS) -s $(SAMTOOLS) --threads=$(THREADS) -o $(@D)/$(*F).$$bedfile_name.report; \
					rm $(@D)/$(*F).$$bedfile_name.HsMetrics.per_base_coverage.2cols ; \
					echo "#[INFO] BAM Metrics on Regions coverage with $$bedfile_name BED file done" >> $@; \
				else \
					echo "#[INFO] BAM Metrics on Regions coverage with $$bedfile_name BED file FAILED" >> $@; \
				fi; \
			done; \
		else \
			echo "#[INFO] BAM Metrics on Regions coverage not generate because no BED list" >> $@; \
		fi; \
	else \
		echo "#[INFO] BAM Metrics on Regions coverage not done because BAM_METRICS=0" >> $@; \
	fi;
	if [ ! -e $@ ]; then echo "#[INFO] BAM Metrics on Regions coverage FAILED" > $@; fi;



# Global run metrics (Sam)
#############################
%.metrics: $(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach PIPELINE,$(PIPELINES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(call aligner,$(PIPELINE)).bam.metrics/metrics )) #$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).list.genes )
	# creates different files in the run directory:
	# <run>.reads.metrics
	# <run>.genes.metrics
	# <run>.design.metrics
	# <run>.amplicon.metrics if applicable
	# see python script for documentation
	#echo "reads.metrics test: ";
	#ls -l $$(dirname $@)/*/*.list.genes;
	#cat $$(dirname $@)/*/*.list.genes;
	$(PYTHON3) $(STARK_RUN_METRICS) --metricsFileList $$(echo $^ | tr " " ",") --outputPrefix $@. ;
	echo "#[INFO] All metrics files on Design and Panel(s), by targets and by genes, for global coverage, depth and coverage, are named $$(basename $@).*" > $@;
	#touch $@;



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# MAIN RULES '$(MK_RELEASE)' : basicaly to manage VCF, BAM, FASTQ, METRICS, MANIFEST, INTERVAL, BED... using SAMTOOLS, FASTQC GATK, PICARD, CAP, BWA, TABIX, IGV..."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )
