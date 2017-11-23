#############################
# NGS workflow validation
# Release: 0.9 
# Date: 31/03/2014
# Author: Antony Le Bechec
#############################


## INPUTS
###########

#Define NGSEnv and config.mk file
NGSEnv?=/NGS
CONFIG?=$(NGSEnv)/scripts/rules.mk/config.mk

#Functions
include $(NGSEnv)/scripts/functions.mk

#Parameters
PARAM?=makefile.param
include $(PARAM)

#FORMAT RUNS_SAMPLES=RUN1:SAMPLE1 RUN1:SAMPLE2 RUN2:SAMPLEX
RUNS_SAMPLES?=
INTERSEC?=2
NB_VARIANTS_TO_SHOW?=20
NB_VARIANTS_TO_SHOW_FULL?=10
FONT?=Courier5 #Times-Roman8

#PIPELINES: ALIGNERS CALLERS ANNOTATORS
ALIGNERS?=bwamem #bwasw #bwamem bwasw
CALLERS?=gatkUG gatkHC #mutect samtools 
ANNOTATORS?=trakxs #trakxs

#NUMBER CALCULATION
NB_SAMPLE=$(shell echo $(RUNS_SAMPLES) | wc -w)
NB_ALIGNERS=$(shell echo $(ALIGNERS) | wc -w)
NB_CALLERS=$(shell echo $(CALLERS) | wc -w)
NB_ANNOTATORS=$(shell echo $(ANNOTATORS) | wc -w)
NB_WORKFLOWS=$(shell echo "$(NB_ALIGNERS) * $(NB_CALLERS) * $(NB_CALLERS)" | bc)

#VARIABLES
SEP=:
DATE=$(shell date '+%Y%m%d-%H%M%S')

## FILES TO GENERATE and KEEP
###############################

# FINAL Files correspond to TAB, VCF, BAM, and METRICS
#RUN_SAMPLE_PATTERN=$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)) ))))
	

# FINAL REPORT and FINAL REPORT "FULL"
FINAL_REPORT?=$(OUTDIR)/analysis.V$(DATE).report
FINAL_REPORT_FULL?=$(FINAL_REPORT).full
FINAL_REPORT_VALIDATION?=$(OUTDIR)/analysis.V$(DATE).validation.report

# REPORTS (intermediate?) to generate for the FINAL REPORT and FINAL REPORT "FULL"
REPORTS=$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).report )
REPORTS_FULL=$(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).report.full )



all: $(FINAL_REPORT_VALIDATION)


# MAIN
########

#Load Config and all available rules
include $(CONFIG)
include $(NGSEnv)/scripts/rules.mk/*.rules.mk


$(FINAL_REPORT_VALIDATION): $(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)) )
	# Validation/check for each sample/VCF
	@echo "########################################### " >> $@
	@echo "###                                     ###" >> $@
	@echo "### Validation Report '$(DATE)' ###" >> $@
	@echo "###                                     ###" >> $@
	@echo "########################################### " >> $@
	#@echo "# Date: $(DATE) " >> $@
	#@echo "# " >> $@
	@echo " " >> $@
	@echo "# SUMMARY " >> $@
	@echo "########### " >> $@
	@echo "# $(NB_WORKFLOWS) WORKFLOWS:" >> $@
	@echo "#    $(NB_ALIGNERS) ALIGNERS     $(ALIGNERS) " >> $@
	@echo "#    $(NB_CALLERS) CALLERS      $(CALLERS) " >> $@
	@echo "#    $(NB_ANNOTATORS) ANNOTATORS   $(ANNOTATORS) " >> $@
	@echo "# REFERENCE GENOME: $(REF_GENOME) " >> $@ 
	@echo "# INTERSECTION THRESHOLD: $(INTERSEC) " >> $@
	@echo "# " >> $@
	@echo "# $(NB_SAMPLE) SAMPLES: $(RUNS_SAMPLES) " >> $@
	@echo " " >> $@
	-for pattern in $^; \
	do \
		echo "" >> $@; \
		echo "" >> $@; \
		echo "### Sample '"$$(basename $$pattern)"'" >> $@; \
		echo "########## " >> $@; \
		for vcf in $$(ls $$pattern/$$(basename $$pattern).*.vcf); \
		do \
			echo "" >> $@; \
			echo "# File '$$vcf' " >> $@; \
			$(TRAKXS)/check_samplemutation_vcf.pl --sample_file=$(TRAKXS)/sample_file_example.txt --sampleID=$$(basename $$vcf | awk -F"." '{print $$1}') --vcf_file=$$vcf --window=5 --show_not_found --comment=' ' >> $@; \
		done; \
	done;


