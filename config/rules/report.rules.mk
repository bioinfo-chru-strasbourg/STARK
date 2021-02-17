############################
# GATK Realignment Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.4.1b"
MK_DATE="27/09/2019"

# Release note
# 11/12/2015-0.9b: Create file
# 01/02/2016-0.9.1b: Debugging
# 13/04/2016-0.9.2b: update full.vcf, final.vcf
# 01/02/2018-0.9.3b: Add CALLING_QUALITY for rule generating full.vcf. Change pipeline name (remove sample name and vcf.gz extension)
# 02/10/2018-0.9.4b: Change Howard annotation, replace VCFTOOLS with BCFTOOLS, merge with multiallele not allowed
# 27/09/2019-0.9.4.1b: Add HOWARD NOMEN field option

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
HOWARD_CALCULATION?=VAF,NOMEN,VAF_STATS,DP_STATS,VARTYPE
HOWARD_NOMEN_FIELDS?="hgvs"

REPORT_SECTIONS?=ALL


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
%.report.header: %.report.summary
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



%.config: %.report.header $(RELEASE)
	mkdir -p $(@D)
	cat $^ > $@
	echo "" >> $@



%.launch.json: %.config %.archive.cram %.manifest %.bed %.list.genes %.tag
	mkdir -p $(@D)
	> $@
	echo "{" >> $@;
	echo "\"sample\": [\"$(*F)\"]," >> $@;
	echo "\"reads\": [\"$(*F).archive.cram\"]," >> $@;
	echo "\"design\": [\"$(*F).$$([[ -s $*.manifest ]] && echo 'manifest' || echo 'bed')\"]," >> $@;
	if [ -s $*.list.genes ]; then echo "\"genes\": [\""$$(for g in $$(cat $*.list.genes); do basename $$g; done | tr '\n' '+' | sed s/+$$//gi)"\"]," >> $@; fi;
	if [ -s $*.transcripts ]; then echo "\"transcripts\": [\"$(*F).transcripts\"]," >> $@; fi;
	echo "\"sample_tag\": [\"$$(cat $*.tag | tr ' ' '!')\"]," >> $@;
	echo "\"application\": [\"$(APP)\"]" >> $@;
	echo "}" >> $@;



## Report for a SAMPLE
%.$(ANALYSIS_DATE).report: $(FINAL) $(REPORT_FILES) %.$(ANALYSIS_DATE).config #$(REPORT_FILES) $(BAM) %.$(ANALYSIS_DATE).config
	@echo "######################################### " > $@
	@echo "### Sample Report '`echo $$(basename $$(dirname $(@D)))`/$(*F)' " >> $@
	@echo "### from Analysis '$(ANALYSIS_REF)' " >> $@
	@echo "######################################### " >> $@
	@echo " " >> $@
	#
	# Report release
	-cp $*.$(ANALYSIS_DATE).config $*.$(ANALYSIS_DATE).config.txt
	#
	# STARK REPORT
	$(STARK_FOLDER_BIN)/STARK.report -f "`echo $$(basename $$(dirname $$(dirname $(@D))))`" -s "$(*F)" -e "$(ENV)" -i $$(echo $(PIPELINES) | tr " " ",") --sections="$(REPORT_SECTIONS)" -k $(ANALYSIS_DATE) -r $(OUTDIR) --verbose



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
%.merge.vcf: %.final_variants_files_vcf_gz %.genome
	# Generate pipeline name list
	cat $< | rev | cut -d/ -f1 | rev | sed s/\.vcf.gz//gi | cut -d. -f2- > $@.pipelines #| tr '\n' '\t' | sed 's/\t$$//'
	# Merge VCF, noramize and rehead with pipelines names
	$(BCFTOOLS) merge -l $< --force-samples -m none --info-rules - | $(BCFTOOLS) norm -m- -f $$(cat $*.genome) | $(BCFTOOLS) norm --rm-dup exact | $(BCFTOOLS) reheader -s $@.pipelines > $@;
	# Cleaning
	-rm -f $@.tmp* $@.pipelines

	

## FULL VCF: ANNOTATION OF A MERGE FILE
%.full.filtration.sorting.vcf: %.merge.vcf %.transcripts %.genome
	# HOWARD annotation
	+$(HOWARD) --input=$< --output=$@.tmp --config=$(HOWARD_CONFIG) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION) --config_annotation=$(HOWARD_CONFIG_ANNOTATION) --annotation=$(HOWARD_ANNOTATION_REPORT) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS) --norm=$$(cat $*.genome);
	# HOWARD calculation and prioritization
	+$(HOWARD) --input=$@.tmp --output=$@  --config=$(HOWARD_CONFIG) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION) --config_annotation=$(HOWARD_CONFIG_ANNOTATION) --calculation=$(HOWARD_CALCULATION_REPORT) --nomen_fields=$(HOWARD_NOMEN_FIELDS) --prioritization=$(HOWARD_PRIORITIZATION_REPORT) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --tmp=$(TMP_FOLDER_TMP)  --transcripts=$*.transcripts --force --multithreading --threads=$(THREADS) --env=$(CONFIG_TOOLS) --norm=$$(cat $*.genome);
	# cleaning
	rm -rf $@.tmp


## rehead full.vcf
%.$(ANALYSIS_DATE).full.vcf: %.$(ANALYSIS_DATE).final_variants_files_vcf_gz %.full.vcf
	cat $< | rev | cut -d/ -f1 | rev | sed s/\.vcf.gz//gi | cut -d. -f2- > $@.pipelines
	$(BCFTOOLS) view -S $@.pipelines $*.full.vcf > $@
	-rm -f $@.tmp* $@.pipelines


## FINAL VCF  RULE
%.final.sorting.vcf: %.full.vcf
	-rm -f $<.tmp.*
	for S in $$(grep "^#CHROM" $< | cut -f10-); do \
		$(BCFTOOLS) view -U -s $$S $< | sed '/^#CHROM/s/'$$S'/'$$(echo $(@F) | cut -d\. -f1)'/' > $<.tmp.$$S; \
		$(BGZIP) -f $<.tmp.$$S; \
		$(TABIX) -f $<.tmp.$$S.gz; \
	done;
	$(BCFTOOLS) concat $<.tmp.*.gz -a -D > $@
	rm -f $<.tmp.*.gz*


## ALL samples VCF RULE
%.samples.vcf: $(REPORT_FILES)
	echo $(REPORT_FILES)
	echo $(REPORT_FILES) | tr " " "\n" | tr "\t" "\n" | grep "final.vcf.gz$$" > $@.tmp.vcf_list
	$(BCFTOOLS) merge -l $@.tmp.vcf_list --info-rules - > $@.tmp.merged;
	+if [ "$(HOWARD_ANNOTATION_ANALYSIS)" != "" ]; then \
		$(HOWARD) --input=$@.tmp.merged --output=$@.tmp.annotated --config=$(HOWARD_CONFIG) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION) --config_annotation=$(HOWARD_CONFIG_ANNOTATION) --annotation=$(HOWARD_ANNOTATION_ANALYSIS) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS); \
	else \
		cp $@.tmp.merged $@.tmp.annotated; \
	fi;
	+if [ "$(HOWARD_CALCULATION_ANALYSIS)" != "" ]; then \
		$(HOWARD) --input=$@.tmp.annotated --output=$@.tmp.calculated --config=$(HOWARD_CONFIG) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION) --config_annotation=$(HOWARD_CONFIG_ANNOTATION) --calculation=$(HOWARD_CALCULATION_ANALYSIS) --nomen_fields=$(HOWARD_NOMEN_FIELDS) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS) --force; \
	else \
		cp $@.tmp.annotated $@.tmp.calculated; \
	fi;
	+if [ "$(HOWARD_PRIORITIZATION_ANALYSIS)" != "" ]; then \
		$(HOWARD) --input=$@.tmp.calculated --output=$@.tmp.prioritized --config=$(HOWARD_CONFIG) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION) --config_annotation=$(HOWARD_CONFIG_ANNOTATION) --prioritization=$(HOWARD_PRIORITIZATION_ANALYSIS) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS) --force; \
	else \
		cp $@.tmp.calculated $@.tmp.prioritized; \
	fi;
	cp $@.tmp.prioritized $@
	-if ((1)); then \
		echo "#DEV samples.vcf GENES filter:" ; \
		$(BGZIP) $@.tmp.prioritized -c > $@.tmp.prioritized.vcf.gz; \
		$(TABIX) $@.tmp.prioritized.vcf.gz; \
		for genes_file in $$(cat $$(cat $@.tmp.vcf_list | xargs dirname | sed s/reports$$/list.genes/) | cut -d. -f2- | sort -u); do \
			List_of_samples=$$(ls $(@D)/[^.]*/[^.]*.$$genes_file | xargs -l basename | cut -d. -f1 | tr "\n" "," | sed s/,$$//) ; \
			echo "samples: "$$List_of_samples; \
			echo "genes_file: "$$genes_file; \
			#echo "VCF files: "; \
			List_of_genes_files=$$(ls $$(for L in $$(echo $$List_of_samples | tr "," "\n"); do echo $(@D)/$$L/$$L.$$genes_file; done;)) ; \
			#echo $$(ls $(@D)/{$$(echo $$List_of_samples),}/{$$(echo $$List_of_samples),}*.$$genes_file) ; \
			echo "BED "; \
			#cat $$(ls $(@D)/{$$(ls $(@D)/[^.]*/[^.]*.$$genes_file | xargs -l basename | cut -d. -f1 | tr "\n" "," | sed s/,$$//)}/{$$(ls $(@D)/[^.]*/[^.]*.$$genes_file | xargs -l basename | cut -d. -f1 | tr "\n" "," | sed s/,$$//)}.$$genes_file) | $(BEDTOOLS) sort | $(BEDTOOLS) merge > $@.tmp.GENES.$$genes_file; \
			#cat $$(ls $(@D)/[^.]*/[^.]*.$$genes_file) | $(BEDTOOLS) sort | $(BEDTOOLS) merge > $@.tmp.GENES.$$genes_file; \
			cat $$(echo $$List_of_genes_files) | $(BEDTOOLS) sort | $(BEDTOOLS) merge > $@.tmp.GENES.$$genes_file; \
			echo "VCF "; \
			$(BCFTOOLS) view --samples $$List_of_samples $@.tmp.prioritized.vcf.gz -R $@.tmp.GENES.$$genes_file > $$(echo $@ | sed s/.vcf$$//).$$genes_file.vcf; \
			$(BGZIP) $$(echo $@ | sed s/.vcf$$//).$$genes_file.vcf ; \
			#$(TABIX) $$(echo $@ | sed s/.vcf$$//).$$genes_file.vcf.gz ; \
			$(HOWARD) --config=$(HOWARD_CONFIG) --config=$(HOWARD_CONFIG) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION) --config_annotation=$(HOWARD_CONFIG_ANNOTATION) --pzfields="PZScore,PZFlag,PZComment,PZInfos" --format=tab  --fields="$(HOWARD_FIELDS)" --sort=$(HOWARD_SORT) --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)"  --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS) --input=$$(echo $@ | sed s/.vcf$$//).$$genes_file.vcf.gz --output=$$(echo $@ | sed s/.vcf$$//).$$genes_file.tsv --force; \
		done; \
	fi;
	-rm -f $@.tmp.*



#
