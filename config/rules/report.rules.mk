############################
# GATK Realignment Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.4b"
MK_DATE="02/10/2018"

# Release note
# 11/12/2015-0.9b: Create file
# 01/02/2016-0.9.1b: Debugging
# 13/04/2016-0.9.2b: update full.vcf, final.vcf
# 01/02/2018-0.9.3b: Add CALLING_QUALITY for rule generating full.vcf. Change pipeline name (remove sample name and vcf.gz extension)
# 02/10/2018-0.9.4b: Change Howard annotation, replace VCFTOOLS with BCFTOOLS, merge with multiallele not allowed

FONT?=Courier5 #Times-Roman8
NGSscripts?=$(NGS_SCRIPTS)
INTERSEC?=2
NB_VARIANTS_TO_SHOW?=20
NB_VARIANTS_TO_SHOW_FULL?=10
HOWARD?=$(NGSscripts)
#HOWARD_CONFIG?=$(HOWARD)/config.ini
HOWARD_CONFIG?="config.ini"
HOWARD_CONFIG_PRIORITIZATION?="config.prioritization.ini"
HOWARD_CONFIG_ANNOTATION?="config.annotation.ini"
HOWARD_ANNOTATION?=$(HOWARD)/VCFannotation.pl
ANNOTATION_TYPE_MINIMAL?="Symbol,location,outcome,hgvs"
ANNOVAR?=$(NGSscripts)
HOWARD_CALCULATION?=VAF,NOMEN,VAF_STATS,VARTYPE



# Analysis Summary
%.report.summary:
	mkdir -p $(@D)
	@echo "# SUMMARY " >> $@
	@echo "########### " >> $@
	@echo "# $(NB_PIPELINES) PIPELINES: $(PIPELINES)" >> $@
	@echo "# POST ALIGNMENT: $(POST_ALIGNMENT) " >> $@
	@echo "# OTHER OPTIONS: " >> $@
	@echo "#    VARIANT_RECALIBRATION: $(VARIANT_RECALIBRATION) " >> $@
	@echo "#    INTERVAL_PADDING: $(INTERVAL_PADDING) " >> $@
	@echo "#    PRIORITIZE_PIPELINES_LIST: $(PRIORITIZE_PIPELINES_LIST) " >> $@
	@echo "# REFERENCE GENOME: $(REF) " >> $@
	@echo "# INTERSECTION THRESHOLD: $(INTERSEC) " >> $@
	@echo "# " >> $@
	@echo " " >> $@


# Analysis Header
%.report.header: %.report.summary #%.config
	mkdir -p $(@D)
	@echo "######################################### " >> $@
	@echo "###                                   ###" >> $@
	@echo "### Analysis Report '$(ANALYSIS_REF)' ###" >> $@
	@echo "###                                   ###" >> $@
	@echo "######################################### " >> $@
	@echo " " >> $@
	@cat $*.report.summary >> $@
	@echo "# " >> $@
	@echo "# $(NB_RUN) RUNS: $(RUNS) " >> $@
	@echo "# $(NB_SAMPLE) SAMPLES: $(RUNS_SAMPLES) " >> $@
	@echo " " >> $@
	#@cat $*.config >> $@
	#@echo " " >> $@

%.config: %.report.header $(RELEASE)
	mkdir -p $(@D)
	#echo "CONFIG" > $@
	#cat $(RELEASE) > $@
	cat $^ > $@
	echo "" >> $@


## Report for a SAMPLE
%.$(ANALYSIS_DATE).report: $(FINAL) $(REPORT_FILES) %.$(ANALYSIS_DATE).config #$(REPORT_FILES) $(BAM) %.$(ANALYSIS_DATE).config
	#%.$(ANALYSIS_DATE).report.summary %.$(ANALYSIS_DATE).report.variants %.$(ANALYSIS_DATE).config
	@echo "######################################### " > $@
	@echo "### Sample Report '`echo $$(basename $$(dirname $(@D)))`/$(*F)' " >> $@
	@echo "### from Analysis '$(ANALYSIS_REF)' " >> $@
	@echo "######################################### " >> $@
	@echo " " >> $@

	# Report release
	-cp $*.$(ANALYSIS_DATE).config $*.$(ANALYSIS_DATE).config.txt

	# STARK REPORT
	echo "$(STARK_FOLDER_BIN)/stark_report.sh -f "`echo $$(basename $$(dirname $$(dirname $(@D))))`" -s "$(*F)" -e '$(ENV)' -i $$(echo $(PIPELINES) | tr " " ",") -d $(ANALYSIS_DATE) -r $(OUTDIR)"
	#$(STARK_FOLDER_BIN)/stark_report.sh -f "`echo $$(basename $$(dirname $$(dirname $(@D))))`" -s "$(*F)" -e "$(ENV)" -i $$(echo $(PIPELINES) | tr " " ",") -d $(ANALYSIS_DATE) -r $(OUTDIR)
	#$(STARK_FOLDER_BIN)/stark_report.sh -f "`echo $$(basename $$(dirname $$(dirname $(@D))))`" -s "$(*F)" -e "$(STARK_FOLDER_CONFIG)/config.app" -i $$(echo $(PIPELINES) | tr " " ",") -d $(ANALYSIS_DATE) -r $(OUTDIR)
	$(STARK_FOLDER_BIN)/stark_report.sh -f "`echo $$(basename $$(dirname $$(dirname $(@D))))`" -s "$(*F)" -e "$(CONFIG_TOOLS)" -i $$(echo $(PIPELINES) | tr " " ",") -d $(ANALYSIS_DATE) -r $(OUTDIR)
	#--env=$(CONFIG_TOOLS)


## list of vcf
%.$(ANALYSIS_DATE).final_variants_files_vcf_gz: %.$(ANALYSIS_DATE).vcfgzs.list $(VCF) %.$(ANALYSIS_DATE).config
	mkdir -p $(@D)
	# Include all vcf.gz of the sample in the list of all vcf.gz of the analysis
	cat $< | tr " " "\n" | grep "^$$(dirname $(@D))/" > $@;
	# TEST
	echo "FINALVARIANTSFILES: $@"
	cat $@


## list of vcf
%.final_variants_files_vcf_gz: %.$(ANALYSIS_DATE).final_variants_files_vcf_gz #$(VCF) %.$(ANALYSIS_DATE).config
	mkdir -p $(@D)
	# list all vcf.gz from this analysis, ordered
	cat $< > $@.full_vcf_gz_list.tmp;
	# Add all vcf.gz including from putative previous analyses
	ls $$(dirname $(@D))/*vcf.gz >> $@.full_vcf_gz_list.tmp;
	# Uniquifiy vcf.gz keeping order
	cat $@.full_vcf_gz_list.tmp | awk '{a[$$0]++} a[$$0]<2 {print $$0}' > $@;
	# TEST
	echo "FINALVARIANTSFILES.tmp: $@.full_vcf_gz_list.tmp"
	cat $@.full_vcf_gz_list.tmp
	echo "FINALVARIANTSFILES: $@"
	cat $@
	# cleaning
	rm -f $@.full_vcf_gz_list.tmp;



## MERGE OF GENERATED VCF
%.merge.unsorted.vcf: %.final_variants_files_vcf_gz %.genome
	# Generate pipeline name list
	cat $< | rev | cut -d/ -f1 | rev | sed s/\.vcf.gz//gi | cut -d. -f2- > $@.pipelines #| tr '\n' '\t' | sed 's/\t$$//'
	#echo "final_variants_files_vcf_gz:"; cat $<;
	#echo "pipelines:"; cat $@.pipelines;
	# Merge VCF, noramize and rehead with pipelines names
	$(BCFTOOLS) merge -l $< --force-samples -m none --info-rules - | $(BCFTOOLS) norm -m- -f `cat $*.genome` | $(BCFTOOLS) reheader -s $@.pipelines > $@;
	# Cleaning
	-rm -f $@.tmp* $@.pipelines


## FULLE VCF: ANNOTATION OF A MERGE FILE
%.full.vcf: %.merge.vcf %.transcripts
	# HOWARD annotation
	+$(HOWARD) --input=$< --output=$@.tmp --config=$(HOWARD_CONFIG) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION) --config_annotation=$(HOWARD_CONFIG_ANNOTATION) --annotation=$(HOWARD_ANNOTATION_REPORT) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS);
	# HOWARD calculation and prioritization
	+$(HOWARD) --input=$@.tmp --output=$@  --config=$(HOWARD_CONFIG) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION) --config_annotation=$(HOWARD_CONFIG_ANNOTATION) --calculation=$(HOWARD_CALCULATION_REPORT) --prioritization=$(HOWARD_PRIORITIZATION_REPORT) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --tmp=$(TMP_FOLDER_TMP)  --transcripts=$*.transcripts --force --multithreading --threads=$(THREADS) --env=$(CONFIG_TOOLS);
	# cleaning
	rm -rf $@.tmp

%.$(ANALYSIS_DATE).full.vcf: %.$(ANALYSIS_DATE).final_variants_files_vcf_gz %.full.vcf
	cat $< | rev | cut -d/ -f1 | rev | sed s/\.vcf.gz//gi | cut -d. -f2- > $@.pipelines
	$(BCFTOOLS) view -S $@.pipelines $*.full.vcf > $@
	-rm -f $@.tmp* $@.pipelines


## FINAL VCF  RULE
%.final.vcf: %.full.vcf
	-rm -f $<.tmp.*
	for S in $$(grep "^#CHROM" $< | cut -f10-); do \
		$(BCFTOOLS) view -U -s $$S $< | sed '/^#CHROM/s/'$$S'/'$$(echo $(@F) | cut -d\. -f1)'/' > $<.tmp.$$S; \
		$(BGZIP) -f $<.tmp.$$S; \
		$(TABIX) -f $<.tmp.$$S.gz; \
	done;
	#$(VCFTOOLS)/vcf-subset -c $$S $< | sed '/^#CHROM/s/'$$S'/'$$(echo $(@F) | cut -d\. -f1)'/' > $<.tmp.$$S; \
	$(BCFTOOLS) concat $<.tmp.*.gz -a -D > $@
	rm -f $<.tmp.*.gz*


## ALL samples VCF RULE
%.samples.vcf: $(REPORT_FILES)
	echo $(REPORT_FILES)
	echo $(REPORT_FILES) | tr " " "\n" | tr "\t" "\n" | grep "final.vcf.gz$$" > $@.vcf_list
	cat $@.vcf_list
	$(BCFTOOLS) merge -l $@.vcf_list --info-rules - > $@
	-rm -f $@.vcf_list
