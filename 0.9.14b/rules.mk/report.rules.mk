############################
# GATK Realignment Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.2b"
MK_DATE="13/04/2016"

# Release note
# 11/12/2015-0.9b: Create file
# 01/02/2016-0.9.1b: Debugging
# 13/04/2016-0.9.2b: update full.vcf, final.vcf

FONT?=Courier5 #Times-Roman8
NGSscripts?=$(NGS_SCRIPTS)
INTERSEC?=2
NB_VARIANTS_TO_SHOW?=20
NB_VARIANTS_TO_SHOW_FULL?=10
HOWARD?=$(NGSscripts)
HOWARD_CONFIG?=$(HOWARD)/config.ini
HOWARD_ANNOTATION?=$(HOWARD)/VCFannotation.pl
ANNOTATION_TYPE_MINIMAL?="Symbol,location,outcome,hgvs"
ANNOVAR?=$(NGSscripts)



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
%.$(ANALYSIS_DATE).report: $(REPORT_FILES) %.$(ANALYSIS_DATE).config
	#%.$(ANALYSIS_DATE).report.summary %.$(ANALYSIS_DATE).report.variants %.$(ANALYSIS_DATE).config
	@echo "######################################### " > $@
	@echo "### Sample Report '`echo $$(basename $$(dirname $(@D)))`/$(*F)' " >> $@
	@echo "### from Analysis '$(ANALYSIS_REF)' " >> $@
	@echo "######################################### " >> $@
	@echo " " >> $@
	
	# Report release
	#-enscript -f $(FONT) --color --header="`echo $$(basename $$(dirname $(@D)))`/$(*F)||[`date '+%d/%m/%Y-%H:%M:%S'`]" -p$(@F).config.ps $*.$(ANALYSIS_DATE).config
	#-ps2pdf $(@F).config.ps $(@F).config.pdf
	#-rm $*.reports/$(@F).config.ps
	-enscript -f $(FONT) --color --header="`echo $$(basename $$(dirname $(@D)))`/$(*F)||[`date '+%d/%m/%Y-%H:%M:%S'`]" -p$*.$(ANALYSIS_DATE).config.ps $*.$(ANALYSIS_DATE).config
	-ps2pdf $*.$(ANALYSIS_DATE).config.ps $*.$(ANALYSIS_DATE).config.pdf
	-rm $*.$(ANALYSIS_DATE).config.ps
	
	# STARK REPORT
	echo "$(STARK)/stark_report.sh -f "`echo $$(basename $$(dirname $$(dirname $(@D))))`" -s "$(*F)" -e $(ENV) -i $$(echo $(PIPELINES) | tr " " ",") -d $(ANALYSIS_DATE) -r $(OUTDIR)"
	$(STARK)/stark_report.sh -f "`echo $$(basename $$(dirname $$(dirname $(@D))))`" -s "$(*F)" -e $(ENV) -i $$(echo $(PIPELINES) | tr " " ",") -d $(ANALYSIS_DATE) -r $(OUTDIR)


## list of vcf
%.final_variants_files_vcf_gz: %.$(ANALYSIS_DATE).final_variants_files_vcf_gz #$(VCF) %.$(ANALYSIS_DATE).config
	mkdir -p $(@D)
	#echo $(VCF) | tr " " "\n" | grep  > $@.full_vcf_gz_list.tmp;
	# list all vcf.gz from this analysis, ordered
	cat $< > $@.full_vcf_gz_list.tmp;
	# Add all vcf.gz including from putative previous analyses
	#ls $(<D)/../*.vcf.gz >> $@.full_vcf_gz_list.tmp;
	ls $$(dirname $(@D))/*vcf.gz >> $@.full_vcf_gz_list.tmp;
	# Uniquifiy vcf.gz keeping order
	#cat $@.full_vcf_gz_list.tmp | uniq > $@;
	cat $@.full_vcf_gz_list.tmp | awk '{a[$$0]++} a[$$0]<2 {print $$0}' > $@;
	# TEST
	echo "FINALVARIANTSFILES.tmp: $@.full_vcf_gz_list.tmp"
	cat $@.full_vcf_gz_list.tmp
	echo "FINALVARIANTSFILES: $@"
	cat $@
	# cleaning
	rm -f $@.full_vcf_gz_list.tmp;
	

## list of vcf
%.$(ANALYSIS_DATE).final_variants_files_vcf_gz: $(VCF) %.$(ANALYSIS_DATE).config #%.final_variants_files_vcf_gz #$(VCF) 
	mkdir -p $(@D)
	touch $@
	# Include all vcf.gz of the sample in the list of all vcf.gz of the analysis
	-for file in $^; do if [[ $$file =~ $$(dirname $(@D))/.*\.gz$$ ]]; then echo $$file >> $@ ; fi; done;
	# TEST
	echo "FINALVARIANTSFILES: $@"
	cat $@


## MERGE OF GENERATED VCF
%.merge.vcf: %.final_variants_files_vcf_gz 
	cat $< | rev | cut -d/ -f1 | rev | sed s/\.vcf.gz/\.vcf/gi | tr '\n' '\t' | sed 's/\t$$//' > $@.pipelines
	$(BCFTOOLS) merge -l $< --force-samples | awk -F"\t" -v SN_PIPELINES="$$(cat $@.pipelines)" '$$1=="#CHROM"{print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"SN_PIPELINES} $$1!="#CHROM"{print}' > $@;
	rm -f $@.tmp  $@.pipelines

#%.$(ANALYSIS_DATE).merge.vcf: %.final_variants_files_vcf_gz 
#	cat $< | rev | cut -d/ -f1 | rev | sed s/\.vcf.gz/\.vcf/gi | tr '\n' '\t' | sed 's/\t$$//' > $@.pipelines
#	$(BCFTOOLS) merge -l $< --force-samples | awk -F"\t" -v SN_PIPELINES="$$(cat $@.pipelines)" '$$1=="#CHROM"{print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"SN_PIPELINES} $$1!="#CHROM"{print}' > $@;
#	rm -f $@.tmp  $@.pipelines
	

## MERGE OF ALL GENERATED VCF GZ including previous analyses
#%.merge.vcf: %.final_variants_files_vcf_gz 
#	ls $(<D)/../*.vcf.gz > $@.full_vcf_gz_list;
#	cat $@.full_vcf_gz_list | rev | cut -d/ -f1 | rev | sed s/\.vcf.gz/\.vcf/gi | tr '\n' '\t' | sed 's/\t$$//' > $@.pipelines
#	$(BCFTOOLS) merge -l $@.full_vcf_gz_list --force-samples | awk -F"\t" -v SN_PIPELINES="$$(cat $@.pipelines)" '$$1=="#CHROM"{print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"SN_PIPELINES} $$1!="#CHROM"{print}' > $@
#	rm $@.full_vcf_gz_list $@.pipelines
	
## FULLE VCF: ANNOTATION OF A MERGE FILE
%.full.vcf: %.merge.vcf
	+$(HOWARD) --input=$< --output=$@ --annotation=$(ANNOTATION_TYPE_MINIMAL) --filter=$(HOWARD_FILTER) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --snpeff_threads=$(THREADS_BY_SAMPLE)  --tmp=$(TMP_FOLDER_TMP) --calculation=FindByPipelines,GenotypeConcordance,VAF_STATS --force;


## FINAL VCF  RULE
%.final.vcf: %.full.vcf 
	-rm -f $<.tmp.*
	for S in $$(grep "^#CHROM" $< | cut -f10-); do \
		$(VCFTOOLS)/vcf-subset -c $$S $< | sed '/^#CHROM/s/'$$S'/'$$(echo $(@F) | cut -d\. -f1)'/' > $<.tmp.$$S; \
		$(BGZIP) -f $<.tmp.$$S; \
		$(TABIX) -f $<.tmp.$$S.gz; \
	done;
	$(BCFTOOLS) concat $<.tmp.*.gz -a -D > $@
	#+$(HOWARD) --input=$@.tmp --output=$@ --filter=$(HOWARD_FILTER) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --snpeff_threads=$(THREADS_BY_SAMPLE)  --tmp=$(TMP_FOLDER_TMP) --force;
	rm -f $<.tmp.*.gz*
	#$@.tmp
	

	










