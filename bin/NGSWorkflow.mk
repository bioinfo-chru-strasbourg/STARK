#############################
# NGS workflow
# Release: 0.9.8.1
# Date: 13/04/2016
# Author: Antony Le Bechec
#############################


## INPUTS
###########

MK_PATH:=$(abspath $(lastword $(MAKEFILE_LIST)))
MK_DIR_PATH:=$(shell dirname $(MK_PATH))
#MK_DIR_PATH:=$(shell echo $(STARK_FOLDER_CONFIG))
PWD:=$(shell pwd -P)
#Define NGSEnv and config.mk file
NGSEnv?=/NGS
#CONFIG?=$(MK_DIR_PATH)/rules.mk/config.mk
CONFIG?=$(STARK_FOLDER_RULES)/config.mk
#FUNCTIONS?=$(MK_DIR_PATH)/functions.mk
FUNCTIONS?=$(STARK_FOLDER_BIN)/functions.mk
PARAM?=makefile.param
#RULES?=$(MK_DIR_PATH)/rules.mk/*.rules.mk
RULES?=$(STARK_FOLDER_RULES)/*.rules.mk
DATE:=$(shell date '+%Y%m%d-%H%M%S')
#ANALYSIS_DATE?=$(DATE)
ANALYSIS_REF?=$(DATE)
ANALYSIS_DATE?=$(ANALYSIS_REF)

# Functions
include $(FUNCTIONS)

# Parameters
include $(PARAM)

# Release
RELEASE_INFOS?=$(TMP_FOLDER_TMP)/$(ANALYSIS_DATE).release_infos # ANALYSIS_DATE
RELEASE_CMD := $(shell echo "" > $(RELEASE_INFOS) )

# pipelines
PIPELINES_INFOS?=$(TMP_FOLDER_TMP)/$(ANALYSIS_DATE).pipelines_infos # ANALYSIS_DATE
PIPELINES_CMD := $(shell echo "" > $(PIPELINES_INFOS) )
PIPELINES_COMMENT := "\#STEP_TYPE:STEP_NAME:STEP_DESCRIPTION:STEP_PARAMETERS"
PIPELINES_CMD := $(shell echo "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

# Load Config and all available rules
include $(CONFIG)

# Load rules
#include $(RULES)


#FORMAT RUNS_SAMPLES=RUN1:SAMPLE1 RUN1:SAMPLE2 RUN2:SAMPLEX
# Found in $(PARAM)
RUNS?=
RUNS_SAMPLES?=



#NUMBER CALCULATION
NB_RUN=$(shell echo $(RUNS) | wc -w)
NB_SAMPLE=$(shell echo $(RUNS_SAMPLES) | wc -w)
NB_ALIGNERS=$(shell echo $(ALIGNERS) | wc -w)
NB_CALLERS=$(shell echo $(CALLERS) | wc -w)
NB_ANNOTATORS=$(shell echo $(ANNOTATORS) | wc -w)
NB_WORKFLOWS=$(shell echo "$(NB_ALIGNERS) * $(NB_CALLERS) * $(NB_ANNOTATORS)" | bc)
NB_PIPELINES=$(shell echo $(PIPELINES) | wc -w)

#VARIABLES
SEP=:
DATE:=$(shell date '+%Y%m%d-%H%M%S')
VALIDATION?=0
REMOVE_INTERMEDIATE_SAM?=1 # Remove SAM file after BAM conversion


## FILES TO GENERATE and KEEP
###############################

#CLEAN :=

# NEEDED files/folder, i.e. Run's configuration and Demultiplexing OK
SAMPLESHEETS?=$(foreach RUN,$(RUNS),$(MISEQDIR)/$(RUN)/SampleSheet.csv )

#DemultiplexedBustardSummary=$(foreach RUN,$(RUNS),$(INPUTDIR)/$(RUN)/DemultiplexedBustardSummary.xml )
NEEDED= $(SAMPLESHEETS) #$(MAKEFILES) $(DemultiplexedBustardSummary)


SAMPLE=	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).SampleSheet.csv ) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).manifest ) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).bed ) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).manifest_name ) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).bed_name )

VCF=	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach ALIGNER,$(ALIGNERS),$(foreach CALLER,$(CALLERS),$(foreach ANNOTATOR,$(ANNOTATORS),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ALIGNER).$(CALLER).$(ANNOTATOR).vcf )))) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach ALIGNER,$(ALIGNERS),$(foreach CALLER,$(CALLERS),$(foreach ANNOTATOR,$(ANNOTATORS),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ALIGNER).$(CALLER).$(ANNOTATOR).vcf.gz )))) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach ALIGNER,$(ALIGNERS),$(foreach CALLER,$(CALLERS),$(foreach ANNOTATOR,$(ANNOTATORS),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ALIGNER).$(CALLER).$(ANNOTATOR).vcf.gz.tbi )))) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach ALIGNER,$(ALIGNERS),$(foreach CALLER,$(CALLERS),$(foreach ANNOTATOR,$(ANNOTATORS),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ALIGNER).$(CALLER).$(ANNOTATOR).tsv )))) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach ALIGNER,$(ALIGNERS),$(foreach CALLER,$(CALLERS),$(foreach ANNOTATOR,$(ANNOTATORS),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ALIGNER).$(CALLER).$(ANNOTATOR).vcf.metrics/metrics )))) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach PIPELINE,$(PIPELINES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(PIPELINE).vcf )) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach PIPELINE,$(PIPELINES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(PIPELINE).vcf.gz )) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach PIPELINE,$(PIPELINES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(PIPELINE).vcf.gz.tbi )) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach PIPELINE,$(PIPELINES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(PIPELINE).vcf.idx )) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach PIPELINE,$(PIPELINES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(PIPELINE).tsv )) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach PIPELINE,$(PIPELINES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(PIPELINE).vcf.metrics/metrics ))

BAM=	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach ALIGNER,$(ALIGNERS),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ALIGNER).bam )) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach ALIGNER,$(ALIGNERS),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ALIGNER).bam.bai )) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach ALIGNER,$(ALIGNERS),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ALIGNER).bam.metrics/metrics )) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach PIPELINE,$(PIPELINES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(call aligner,$(PIPELINE)).bam )) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach PIPELINE,$(PIPELINES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(call aligner,$(PIPELINE)).bam.bai )) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach PIPELINE,$(PIPELINES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(call aligner,$(PIPELINE)).bam.metrics/metrics ))


#$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach PIPELINE,$(PIPELINES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(call aligner,$(PIPELINE)).bam.bed ))
#$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach ALIGNER,$(ALIGNERS),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ALIGNER).bam.bed ))


UBAM=	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).unaligned.bam ) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).unaligned.bam.bai )



FASTQC_METRICS=$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).fastqc/metrics )

SEQUENCING_METRICS=$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).sequencing/metrics )


CRAM=	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).archive.cram ) \

FINAL=$(SAMPLE) $(BAM) $(VCF) $(CRAM) $(SEQUENCING_METRICS) #$(FASTQC_METRICS)
#$(UBAM)

REPORTS=$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).$(ANALYSIS_DATE).report )

#REPORT_FILES=	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).$(ANALYSIS_DATE).full.vcf ) \
#		$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).$(ANALYSIS_DATE).full.txt ) \
#		$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).$(ANALYSIS_DATE).final.vcf ) \
#		$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).$(ANALYSIS_DATE).final.txt ) \
#		$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).full.vcf ) \
#		$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).full.txt ) \
#		$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).final.vcf ) \
#		$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).final.txt ) \



REPORT_FILES=	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).full.vcf ) \
		$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).full.vcf.idx ) \
		$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).full.vcf.gz ) \
		$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).full.vcf.gz.tbi ) \
		$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).full.tsv ) \
		$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).final.vcf ) \
		$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).final.vcf.idx ) \
		$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).final.vcf.gz ) \
		$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).final.vcf.gz.tbi ) \
		$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).final.tsv ) \
		$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).final.vcf.metrics/metrics ) \


REPORT_FILES_EXT= final.vcf final.vcf.idx final.vcf.gz final.vcf.gz.tbi full.vcf full.vcf.idx full.vcf.gz full.vcf.gz.tbi final.tsv full.tsv $(ANALYSIS_DATE).final.vcf $(ANALYSIS_DATE).final.vcf.idx $(ANALYSIS_DATE).final.vcf.gz $(ANALYSIS_DATE).final.vcf.gz.tbi $(ANALYSIS_DATE).full.vcf $(ANALYSIS_DATE).full.vcf.idx $(ANALYSIS_DATE).full.vcf.gz $(ANALYSIS_DATE).full.vcf.gz.tbi $(ANALYSIS_DATE).final.tsv $(ANALYSIS_DATE).full.tsv


REPORT_FILES2=$(foreach EXT,$(REPORT_FILES_EXT), $(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).reports/$(call sample,$(RUN_SAMPLE)).$(EXT) ) )



ifdef FINAL_REPORT
	FINAL_REPORT_FILES=$(FINAL_REPORT) $(REPORT_FILES2)
endif


# Load rules
include $(RULES)

all: $(FINAL) $(FINAL_REPORT) $(FINAL_REPORT_FILES) $(FINAL_REPORT).samples.vcf.gz $(FINAL_REPORT).samples.tsv #$(FINAL_REPORT_FULL) $(FINAL_REPORT_FULL_VCF)
	# List of files to generate
	echo $^
	# CLEANING
	echo $(CLEAN)
	-rm -f $(CLEAN)
	#rm -rf $(OUTDIR)/*/*/*.for_metrics_bed
	#rm -rf $(OUTDIR)/*/*/*.*.manifest
	#echo $(CLEAN)$(SAMPLESHEETS)

# MAIN
########



$(RELEASE):
	#echo "RELEASE" > $@
	cat $(RELEASE_INFOS) > $@
	#head -n `echo $$(grep ^ $(RELEASE_INFOS) -c) / 2 | bc` $(RELEASE_INFOS) > $@
	#head -n `echo $$(grep ^ $(RELEASE_INFOS) -c) / 2 | bc` $(RELEASE_INFOS) > $@
	#echo "" >> $@

## List of BAMS
#%.bams.list: $(foreach ALIGNER,$(ALIGNERS),%.$(ALIGNER).bam)
%.bams.list: $(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach ALIGNER,$(ALIGNERS),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ALIGNER).bam )) $(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach PIPELINE,$(PIPELINES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(call aligner,$(PIPELINE)).bam ))
	ls $^ > $@
	echo "BAM_LIST2 $*.test: $^"
	#cp $@ $*.test

%.beds.list: $(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).bed )
	mkdir -p $(@D) ;
	ls $^ > $@;
	echo "BED LIST $*.test: $^";

%.SampleSheets.list: $(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).SampleSheet.csv )
	mkdir -p $(@D);
	ls $^ > $@;


#%.List of VCF
%.vcfs.list: $(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach ALIGNER,$(ALIGNERS),$(foreach CALLER,$(CALLERS),$(foreach ANNOTATOR,$(ANNOTATORS),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ALIGNER).$(CALLER).$(ANNOTATOR).vcf )))) $(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach PIPELINE,$(PIPELINES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(PIPELINE).vcf ))
	mkdir -p $(@D)
	ls $^ > $@

#%.List of VCF.GZ
%.vcfgzs.list: $(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach ALIGNER,$(ALIGNERS),$(foreach CALLER,$(CALLERS),$(foreach ANNOTATOR,$(ANNOTATORS),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ALIGNER).$(CALLER).$(ANNOTATOR).vcf.gz )))) $(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach PIPELINE,$(PIPELINES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(PIPELINE).vcf.gz ))
	mkdir -p $(@D)
	ls $^ > $@

# Files for Reports
#VARIANTS_FILES_FOR_REPORTS=$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ANALYSIS_DATE).variants.intersect.txt )

# REPORTS (intermediate?) to generate for the FINAL REPORT and FINAL REPORT "FULL"
#REPORTS=$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ANALYSIS_DATE).report )


#REPORTS_FULL=$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ANALYSIS_DATE).report.full )



## Create all reports for all SAMPLES
#$(FINAL_REPORT): %.$(ANALYSIS_DATE).report.summary $(REPORTS)
$(FINAL_REPORT): $(FINAL_REPORT).report.header $(REPORTS)
	@echo " " >> $@
	# Creating Report $@ with Reports $^
	-for report in $^; \
	do \
		echo "" >> $@; \
		cat $$report >> $@; \
		rm -f $$report; \
	done;
	# PDF creation
	@enscript -f $(FONT) --header="$(@F)||[`date '+%d/%m/%Y-%H:%M:%S'`]" -p$@.ps $@
	@ps2pdf $@.ps $@.pdf
	@rm -f $@.ps


#$(FINAL_REPORT_FULL): %.$(ANALYSIS_DATE).report.summary $(REPORTS_FULL)
#$(FINAL_REPORT_FULL): $(FINAL_REPORT_FULL).report.header $(REPORTS_FULL)
#	# Creating Report $@ with Reports $^
#	-for report in $^; \
#	do \
#		echo "" >> $@; \
#		cat $$report >> $@; \
#		rm $$report; \
#	done;
#	# PDF creation
#	@enscript -f $(FONT) --header="$(@F)||[`date '+%d/%m/%Y-%H:%M:%S'`]" -p$@.ps $@
#	@ps2pdf $@.ps $@.pdf
#	@rm $@.ps

# Warning !!! vcf-merge error occured !!!
#$(FINAL_REPORT_FULL_VCF): $(FINAL_REPORT_FULL_VCF).vcfs.list
#	-for VCF in `cat $<`; do \
#		if [ ! -e $$VCF.gz ]; then \
#			$(BGZIP) -f -c $$VCF > $$VCF.gz; \
#		fi; \
#		if [ ! -e $VCF.gz.tbi ]; then \
#			echo "# VCF tabix"; \
#			$(TABIX) -p vcf $$VCF.gz; \
#		fi; \
#	done;
#	$(VCFTOOLS)/vcf-merge $$(for VCF in $$(echo $$(cat $<)); do echo -n $$VCF'.gz '; done;) > $@;


# CLEAN

CLEAN=	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).*.*manifest ) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))*/*manifest ) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).*genes_from_manifest ) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).*manifest_by_samples.txt ) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).*metrics.bam ) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).*.metrics.bed ) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).*.transcripts ) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).*.intervals.bed ) \
	$(RELEASE_INFOS) \
	$(PIPELINES_INFOS) \
	#$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).*for_metrics_bed ) \
	#$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ANALYSIS_DATE)* ) \
	#$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(RELEASE_DATE) ) \
	#$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(RELEASE_REF) ) \
	$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).*.*genes ) \



clean: #$(CLEAN)
	echo $(CLEAN)
	-rm -f $(CLEAN)
	# Remove all symlink
	#find TEST -type l | xargs rm
	#rm -rf $(OUTDIR)/*/*/*.for_metrics_bed
	#rm -rf $(OUTDIR)/*/*/*.*.manifest
	#rm -rf $(OUTDIR)/*/*/*.*.genes
	#rm -rf $(OUTDIR)/*/*/$(ANALYSIS_DATE)
	#rm -rf $(OUTDIR)/*/*/$(RELEASE_DATE)
	#rm -rf $(OUTDIR)/*/*/$(RELEASE_REF)
