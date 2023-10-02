############################
# GATK Realignment Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.5"
MK_DATE="27/09/2019"

# Release note
# 11/12/2015-0.9b: Create file
# 01/02/2016-0.9.1b: Debugging
# 13/04/2016-0.9.2b: update full.vcf, final.vcf
# 01/02/2018-0.9.3b: Add CALLING_QUALITY for rule generating full.vcf. Change pipeline name (remove sample name and vcf.gz extension)
# 02/10/2018-0.9.4b: Change Howard annotation, replace VCFTOOLS with BCFTOOLS, merge with multiallele not allowed
# 27/09/2019-0.9.4.1b: Add HOWARD NOMEN field option
# 06/02/2023-0.9.5: Add INFO_to_FORMAT, add threads on bcftools and bgzip

INTERSEC?=2
NB_VARIANTS_TO_SHOW?=20
NB_VARIANTS_TO_SHOW_FULL?=10
REPORT_VARIANTS_FULL?=0

REPORT_SECTIONS?=ALL


# Analysis Summary
%.report.summary:
	mkdir -p $(@D)
	@echo "" >> $@


# Analysis Header
%.report.header: %.report.summary
	mkdir -p $(@D)
	@echo "######################################### " >> $@
	@echo "###                                   ###" >> $@
	@echo "### Analysis Report '$(ANALYSIS_REF)' ###" >> $@
	@echo "###                                   ###" >> $@
	@echo "######################################### " >> $@



%.config: %.report.header $(RELEASE)
	mkdir -p $(@D)
	cat $*.report.header > $@
	cat $(RELEASE) >> $@
	echo "" >> $@


%.config.txt: %.config
	cp -p $< $@


%.analysis.json: %.archive.cram %.manifest %.bed %.list.genes %.tag
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
%.$(ANALYSIS_DATE).report: $(FINAL) $(REPORT_FILES) $(VCF_REPORT_FILES) %.$(ANALYSIS_DATE).config 
	@echo "######################################### " > $@
	@echo "### Sample Report '`echo $$(basename $$(dirname $(@D)))`/$(*F)' " >> $@
	@echo "### from Analysis '$(ANALYSIS_REF)' " >> $@
	@echo "######################################### " >> $@
	@echo " " >> $@
	# STARK REPORT
	rm -rf $*.report*
	if ! $(STARK_FOLDER_BIN)/STARK.report -f "`echo $$(basename $$(dirname $$(dirname $(@D))))`" -s "$(*F)" -e "$(ENV)" -i $$(echo $(PIPELINES) | tr " " ",") --sections="$(REPORT_SECTIONS)" -k $(ANALYSIS_DATE) -r $(OUTDIR) --tmp=$(TMP_FOLDER_TMP) --output=$*.report --verbose; then \
		$(STARK_FOLDER_BIN)/STARK.report -f "`echo $$(basename $$(dirname $$(dirname $(@D))))`" -s "$(*F)" -e "$(ENV)" -i $$(echo $(PIPELINES) | tr " " ",") --sections="$(REPORT_SECTIONS)" -k $(ANALYSIS_DATE) -r $(OUTDIR) --tmp=/tmp --output=$*.report --verbose; \
	fi;
	tar -C $(@D) -cvf $@.tar.gz $(*F).report.html $(*F).report.annex.html $(*F).report.html.folder/;
	


## list of vcf
%.$(ANALYSIS_DATE).final_variants_files_vcf_gz: %.$(ANALYSIS_DATE).vcfgzs.list $(VCF)
	mkdir -p $(@D)
	# Include all vcf.gz of the sample in the list of all vcf.gz of the analysis
	cat $< | tr " " "\n" | grep "^$$(dirname $(@D))/" > $@;
	# TEST
	echo "FINALVARIANTSFILES: $@"
	cat $@


## list of vcf
%.final_variants_files_vcf_gz: %.$(ANALYSIS_DATE).final_variants_files_vcf_gz
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
%.merge.vcf: %.final_variants_files_vcf_gz $(BAM)
	# Generate pipeline name list
	cat $< | rev | cut -d/ -f1 | rev | sed s/\.vcf.gz//gi | cut -d. -f2- > $@.pipelines
	# Merge VCF, normalize and rehead with pipelines names
	# bcftools norm -f <ref> <vcf> is not compatible with breakends. Do that treatment separately.
	# 1) merge | keep only breakends
	# 2) merge | exclude breakends | do the bcftools norm that crashes on breakends
	# 3) merge the two above | rest of normalization
	#
	# In case there is no SVTYPE, do the full command directly
	#
	$(BCFTOOLS) merge --threads=$(THREADS_BY_SAMPLE) -l $< --force-samples -m none --info-rules - > $@.merge_step0.vcf
	if (($(BCFTOOLS) head $@.merge_step0.vcf | grep ID=SVTYPE)); then \
	$(BCFTOOLS) view -i 'INFO/SVTYPE="BND"' $@.merge_step0.vcf > $@.bnd_only.tmp.vcf; \
	$(BGZIP) $@.bnd_only.tmp.vcf; \
	$(TABIX) $@.bnd_only.tmp.vcf.gz; \
	$(BCFTOOLS) view -e 'INFO/SVTYPE="BND"' $@.merge_step0.vcf | $(BCFTOOLS) norm --threads=$(THREADS_BY_SAMPLE) -m- -f $(GENOME) > $@.all_except_bnd.tmp.vcf; \
	$(BGZIP) $@.all_except_bnd.tmp.vcf; \
	$(TABIX) $@.all_except_bnd.tmp.vcf.gz; \
	$(BCFTOOLS) concat $@.bnd_only.tmp.vcf.gz $@.all_except_bnd.tmp.vcf.gz -a | $(BCFTOOLS) sort | $(BCFTOOLS) norm --threads=$(THREADS_BY_SAMPLE) --rm-dup exact | $(BCFTOOLS) +fill-tags -- -t AN,AC,AF,AC_Hemi,AC_Hom,AC_Het,ExcHet,HWE,MAF,NS | $(BCFTOOLS) reheader --threads=$(THREADS_BY_SAMPLE) -s $@.pipelines > $@.tmp.merged.vcf; \
	rm -f $@.merge_step0.vcf $@.bnd_only.tmp.vcf $@.all_except_bnd.tmp.vcf; \
	else \
	$(BCFTOOLS) merge --threads=$(THREADS_BY_SAMPLE) -l $< --force-samples -m none --info-rules - | $(BCFTOOLS) sort | $(BCFTOOLS) norm --threads=$(THREADS_BY_SAMPLE) --rm-dup exact | $(BCFTOOLS) +fill-tags -- -t AN,AC,AF,AC_Hemi,AC_Hom,AC_Het,ExcHet,HWE,MAF,NS | $(BCFTOOLS) reheader --threads=$(THREADS_BY_SAMPLE) -s $@.pipelines > $@.tmp.merged.vcf; \
	fi;
	# | $(BCFTOOLS) view --exclude 'FORMAT/GT="0/0"'
	# | $(BCFTOOLS) +setGT  -- -t . -n 0 
	# header file
	echo -e '##INFO=<ID=Validation_Depth_Flags,Number=.,Type=String,Description="Depth metrics flag from BAM validation">' > $@.tmp.annotate.Validation_Depth_Flags.hdr;
	echo -e '##INFO=<ID=Validation_Depth,Number=.,Type=Float,Description="Depth metrics from BAM validation">' > $@.tmp.annotate.Validation_Depth.hdr;
	# Add validation flags depth and alignments depth header
	grep "^##" $@.tmp.merged.vcf > $@.tmp.merged.reheaded.vcf;
	cat $@.tmp.annotate.Validation_Depth_Flags.hdr >> $@.tmp.merged.reheaded.vcf;
	cat $@.tmp.annotate.Validation_Depth.hdr >> $@.tmp.merged.reheaded.vcf;
	grep "^#CHROM" $@.tmp.merged.vcf >> $@.tmp.merged.reheaded.vcf;
	-grep "^#" -v $@.tmp.merged.vcf >> $@.tmp.merged.reheaded.vcf;
	# Add validation flags depth and alignments depth (if files exist)
	if ((1)); then \
		if (( $$(ls $$(dirname $$(dirname $@))/*bam.metrics/*.validation.flags.Design.bed | wc -l) )) || (( $$(ls $$(dirname $$(dirname $@))/*bam.metrics/*.design.bed.HsMetrics.per_base_coverage.gz | wc -l) )); then \
			if (( $$(ls $$(dirname $$(dirname $@))/*bam.metrics/*.validation.flags.Design.bed | wc -l) )); then \
				cat $$(dirname $$(dirname $@))/*bam.metrics/*.validation.flags.Design.bed | sort -k1,1 -k2,2 | $(BEDTOOLS) merge -c 10 -o distinct | $(BGZIP) --threads=$(THREADS_BY_SAMPLE) -c --index --index-name $@.tmp.flags.bed.gz.tbi > $@.tmp.flags.bed.gz; \
				$(BCFTOOLS) annotate --threads=$(THREADS_BY_SAMPLE) -a $@.tmp.flags.bed.gz -h $@.tmp.annotate.Validation_Depth_Flags.hdr -c CHROM,POS,TO,INFO/Validation_Depth_Flags -l Validation_Depth_Flags:unique $@.tmp.merged.reheaded.vcf > $@.tmp.merged.reheaded.flags.vcf; \
			else \
				cp $@.tmp.merged.reheaded.vcf $@.tmp.merged.reheaded.flags.vcf; \
			fi; \
			if (( $$(ls $$(dirname $$(dirname $@))/*bam.metrics/*.design.bed.HsMetrics.per_base_coverage.gz | wc -l) )); then \
				zcat $$(dirname $$(dirname $@))/*bam.metrics/*.design.bed.HsMetrics.per_base_coverage.gz | cut -f1,2,4 | grep ^chrom -v | sort -k1,1 -k2,2n | bgzip -c > $@.tmp.depth.tab.gz; \
				tabix -s1 -b2 -e2 $@.tmp.depth.tab.gz; \
				$(BCFTOOLS) annotate --threads=$(THREADS_BY_SAMPLE) -a $@.tmp.depth.tab.gz -h $@.tmp.annotate.Validation_Depth.hdr -c CHROM,POS,INFO/Validation_Depth -l Validation_Depth:avg $@.tmp.merged.reheaded.flags.vcf > $@; \
			else \
				cp $@.tmp.merged.reheaded.flags.vcf $@; \
			fi; \
		else \
			cp $@.tmp.merged.reheaded.vcf $@; \
		fi; \
	else \
		cp $@.tmp.merged.reheaded.vcf $@ ; \
	fi;
	# Cleaning
	rm -f $@.tmp* $@.pipelines



## FULL VCF: ANNOTATION OF A MERGE FILE
# %.full.vcf: %.merge.vcf %.transcripts
# 	cp $< $@.tmp0
# 	# Prevent comma in description in vcf header
# 	$(STARK_FOLDER_BIN)/fix_vcf_header.sh --input=$@.tmp0 --output=$@.tmp0 --threads=$(THREADS_BY_SAMPLE) --bcftools=$(BCFTOOLS) $(FIX_VCF_HEADER_REFORMAT_option);
# 	# HOWARD annotation
# 	+$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.tmp0 --output=$@.tmp1 --annotation=$(HOWARD_ANNOTATION_REPORT) --norm=$(GENOME);
# 	# HOWARD annotation dejavu (forced)
# 	# Prevent comma in description in vcf header
# 	$(STARK_FOLDER_BIN)/fix_vcf_header.sh --input=$@.tmp1 --output=$@.tmp1 --threads=$(THREADS_BY_SAMPLE) --bcftools=$(BCFTOOLS) $(FIX_VCF_HEADER_REFORMAT_option);
# 	+if [ "$(HOWARD_DEJAVU_ANNOTATION)" != "" ]; then \
# 		$(HOWARD) $(HOWARD_DEJAVU_CONFIG_OPTIONS) --input=$@.tmp1 --output=$@.tmp2 --annotation=$(HOWARD_DEJAVU_ANNOTATION) --norm=$(GENOME) --force; \
# 	else \
# 		mv $@.tmp1 $@.tmp2; \
# 	fi;
# 	# Prevent comma in description in vcf header
# 	$(STARK_FOLDER_BIN)/fix_vcf_header.sh --input=$@.tmp2 --output=$@.tmp2 --threads=$(THREADS_BY_SAMPLE) --bcftools=$(BCFTOOLS) $(FIX_VCF_HEADER_REFORMAT_option);
# 	# HOWARD calculation and prioritization
# 	+$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.tmp2 --output=$@.tmp3 --calculation=$(HOWARD_CALCULATION_REPORT) --nomen_fields=$(HOWARD_NOMEN_FIELDS) --prioritization=$(HOWARD_PRIORITIZATION_REPORT) --transcripts=$*.transcripts --force --norm=$(GENOME);
# 	# Prevent comma in description in vcf header
# 	$(STARK_FOLDER_BIN)/fix_vcf_header.sh --input=$@.tmp3 --output=$@.tmp3 --threads=$(THREADS_BY_SAMPLE) --bcftools=$(BCFTOOLS) $(FIX_VCF_HEADER_REFORMAT_option);
# 	+$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.tmp3 --output=$@ --prioritization=$(HOWARD_PRIORITIZATION_VARANK) --force --norm=$(GENOME);
# 	# cleaning
# 	rm -rf $@.tmp*

%.full.vcf: %.merge.vcf %.transcripts
	cp $< $@.tmp0
	# Prevent comma in description in vcf header
	$(STARK_FOLDER_BIN)/fix_vcf_header.sh --input=$@.tmp0 --output=$@ --threads=$(THREADS_BY_SAMPLE) --bcftools=$(BCFTOOLS) $(FIX_VCF_HEADER_REFORMAT_option);
	# cleaning
	rm -rf $@.tmp*


## rehead full.vcf
%.$(ANALYSIS_DATE).full.vcf: %.$(ANALYSIS_DATE).final_variants_files_vcf_gz %.full.vcf
	cat $< | rev | cut -d/ -f1 | rev | sed s/\.vcf.gz//gi | cut -d. -f2- > $@.pipelines
	$(BCFTOOLS) view --threads=$(THREADS_BY_SAMPLE) -S $@.pipelines $*.full.vcf > $@
	-rm -f $@.tmp* $@.pipelines


## FINAL VCF  RULE
%.final$(POST_CALLING_MERGING).vcf: %.full.vcf
	-rm -f $<.tmp.*
	for S in $$(grep "^#CHROM" $< | cut -f10-); do \
		$(BCFTOOLS) view --threads=$(THREADS_BY_SAMPLE) -U -s $$S $< | $(BCFTOOLS) view --threads=$(THREADS_BY_SAMPLE) -e 'FORMAT/GT="0/0"' | sed '/^#CHROM/s/'$$S'/'$$(echo $(@F) | cut -d\. -f1)'/' > $<.tmp.$$S; \
		$(BGZIP) --threads=$(THREADS_BY_SAMPLE) -f $<.tmp.$$S; \
		$(TABIX) -f $<.tmp.$$S.gz; \
	done;
	$(BCFTOOLS) concat $<.tmp.*.gz -a -D > $@.tmp;
	if [ "$(INFO_TO_FORMAT_ANNOTATIONS)" != "" ]; then \
		$(STARK_FOLDER_BIN)/INFO_to_FORMAT.sh --input=$@.tmp --output=$@ --annotations=$(INFO_TO_FORMAT_ANNOTATIONS) --bcftools=$(BCFTOOLS) --tabix=$(TABIX) --threads=$(THREADS_BY_SAMPLE); \
	else \
		mv $@.tmp $@; \
	fi;
	rm -f $<.tmp.*.gz* $@.tmp*


## ALL samples VCF RULE

# %.variants: $(VCF_REPORT_FILES) %.variants_full
# 	# List of final VCF files
# 	echo $^ | tr " " "\n" | tr "\t" "\n" | grep "final.vcf.gz$$" > $@.tmp.vcf_list
# 	$(BCFTOOLS) merge --threads=$(THREADS) -l $@.tmp.vcf_list -m none --force-samples | $(BCFTOOLS) filter --threads=$(THREADS) -S . -e 'GT=="0/0" | GT=="0|0"' > $@.tmp.merged;
# 	# Prevent comma in description in vcf header
# 	$(STARK_FOLDER_BIN)/fix_vcf_header.sh --input=$@.tmp.merged --output=$@.tmp.merged --threads=$(THREADS_BY_SAMPLE) --bcftools=$(BCFTOOLS) $(FIX_VCF_HEADER_REFORMAT_option);
# 	# Annotation
# 	+$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.tmp.merged --output=$@.tmp.annotated1 --annotation=$(HOWARD_ANNOTATION_ANALYSIS) --threads=$(THREADS);
# 	# Prevent comma in description in vcf header
# 	$(STARK_FOLDER_BIN)/fix_vcf_header.sh --input=$@.tmp.annotated1 --output=$@.tmp.annotated1 --threads=$(THREADS_BY_SAMPLE) --bcftools=$(BCFTOOLS) $(FIX_VCF_HEADER_REFORMAT_option);
# 	# Annotation dejavu (forced)
# 	+if [ "$(HOWARD_DEJAVU_ANNOTATION)" != "" ]; then \
# 		$(HOWARD) $(HOWARD_DEJAVU_CONFIG_OPTIONS) --input=$@.tmp.annotated1 --output=$@.tmp.annotated2 --annotation=$(HOWARD_DEJAVU_ANNOTATION) --norm=$$(cat $*.genome) --threads=$(THREADS) --force; \
# 	else \
# 		mv $@.tmp.annotated1 $@.tmp.annotated2; \
# 	fi;
# 	# Prevent comma in description in vcf header
# 	$(STARK_FOLDER_BIN)/fix_vcf_header.sh --input=$@.tmp.annotated2 --output=$@.tmp.annotated2 --threads=$(THREADS_BY_SAMPLE) --bcftools=$(BCFTOOLS) $(FIX_VCF_HEADER_REFORMAT_option);
# 	# Calculation and prioritization (forced)
# 	+$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.tmp.annotated2 --output=$@.tmp.calculated.prioritized --calculation=$(HOWARD_CALCULATION_ANALYSIS) --prioritization=$(HOWARD_PRIORITIZATION_ANALYSIS) --nomen_fields=$(HOWARD_NOMEN_FIELDS) --force --threads=$(THREADS);
# 	# Prevent comma in description in vcf header
# 	$(STARK_FOLDER_BIN)/fix_vcf_header.sh --input=$@.tmp.calculated.prioritized --output=$@.tmp.calculated.prioritized --threads=$(THREADS_BY_SAMPLE) --bcftools=$(BCFTOOLS) $(FIX_VCF_HEADER_REFORMAT_option);
# 	# Sort VCF
# 	mkdir -p $@.tmp.calculated.prioritized.SAMTOOLS_PREFIX
# 	$(BCFTOOLS) sort -T $@.tmp.calculated.prioritized.SAMTOOLS_PREFIX $@.tmp.calculated.prioritized > $@.tmp.calculated.prioritized.sorted
# 	rm -rf $@.tmp.calculated.prioritized.SAMTOOLS_PREFIX
# 	# Generate Design VCF
# 	$(BGZIP) --threads=$(THREADS) -c $@.tmp.calculated.prioritized.sorted > $@.tmp.calculated.prioritized.sorted.vcf.gz;
# 	$(TABIX) $@.tmp.calculated.prioritized.sorted.vcf.gz;
# 	cp $@.tmp.calculated.prioritized.sorted.vcf.gz $@.Design.vcf.gz;
# 	$(TABIX) $@.Design.vcf.gz; \
# 	# Generate Design TSV
# 	+$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.Design.vcf.gz --output=$@.Design.tsv --translation=TSV --fields="$(HOWARD_FIELDS)" --sort=$(HOWARD_SORT) --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)" --pzfields="PZScore,PZFlag,PZComment,PZInfos" --force --threads=$(THREADS);
# 	# Generate Panel(s) VCF and TSV from Design VCF ($@.Design.vcf.gz)
# 	+for genes_file in $$(cat $$(cat $@.tmp.vcf_list | xargs dirname | sed s/reports$$/list.genes/) | cut -d. -f2- | sort -u); do \
# 		# List of Samples with $$genes_file panel \
# 		echo "genes_file for variants is: "$$genes_file; \
# 		List_of_samples=$$(ls $(@D)/[^.]*/[^.]*.$$genes_file | xargs -l basename | cut -d. -f1 | sort -u | tr "\n" "," | sed s/,$$// 2>/dev/null) ; \
# 		echo "List_of_samples is: "$$List_of_samples; \
# 		if [ "$$List_of_samples" != "" ] || ((1)); then \
# 			# List of $$genes_file within Samples folders \
# 			List_of_genes_files=$$(ls $$(for L in $$(echo $$List_of_samples | tr "," "\n"); do echo $(@D)/$$L/$$L.$$genes_file; done;)) ; \
# 			# Merge all $$genes_file found into uniq BED file (but supposed to be the same) \
# 			cat $$(echo $$List_of_genes_files) | $(BEDTOOLS) sort | $(BEDTOOLS) merge > $@.tmp.GENES.$$genes_file; \
# 			# Generate VCF Panel from VCF Design (especially $@.tmp.calculated.prioritized.sorted.vcf.gz because tabix) with $$genes_file for List of Samples \
# 			if (( $$(grep ^ $@.tmp.GENES.$$genes_file -c) )); then \
# 				$(BCFTOOLS) view --threads=$(THREADS) --samples $$List_of_samples -U --force-samples $@.tmp.calculated.prioritized.sorted.vcf.gz -R $@.tmp.GENES.$$genes_file > $@.Panel.$$genes_file.vcf; \
# 			else \
# 				$(BCFTOOLS) view --threads=$(THREADS) --samples $$List_of_samples -U --force-samples $@.tmp.calculated.prioritized.sorted.vcf.gz > $@.Panel.$$genes_file.vcf; \
# 			fi; \
# 			# Compress VCF \
# 			$(BGZIP) --threads=$(THREADS) $@.Panel.$$genes_file.vcf; \
# 			$(TABIX) $@.Panel.$$genes_file.vcf.gz; \
# 			# Generate TSV Panel from VCF Panel compressed  \
# 			$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.Panel.$$genes_file.vcf.gz --output=$@.Panel.$$genes_file.tsv --translation=TSV --fields="$(HOWARD_FIELDS)" --sort=$(HOWARD_SORT) --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)" --pzfields="PZScore,PZFlag,PZComment,PZInfos" --force --threads=$(THREADS); \
# 		fi; \
# 	done;
# 	-rm -f $@.tmp.*
# 	echo "#[INFO] All variants files on Design and Panel(s) are named $$(basename $@).*" > $@;



# %.variants_full: $(VCF_REPORT_FILES)
# 	+if (($(REPORT_VARIANTS_FULL))); then \
# 		# List of full VCF files \
# 		echo $^ | tr " " "\n" | tr "\t" "\n" | grep "full.vcf.gz$$" > $@.tmp.vcf_list; \
# 		for f in $$(cat $@.tmp.vcf_list); do $(BCFTOOLS) view --threads=$(THREADS) -h $$f | grep "^#CHROM" | cut -f10- | sed  "s/^/$$(basename $$f | cut -d. -f1)./g;s/\t/\t$$(basename $$f | cut -d. -f1)./g";  done | tr "\t" "\n" > $@.tmp.vcf_list.pipelines; \
# 		$(BCFTOOLS) merge --threads=$(THREADS) -l $@.tmp.vcf_list -m none --force-samples | $(BCFTOOLS) filter --threads=$(THREADS) -S . -e 'GT=="0/0" | GT=="0|0"' | $(BCFTOOLS) reheader --threads=$(THREADS) -s $@.tmp.vcf_list.pipelines > $@.tmp.merged; \
# 		# Prevent comma in description in vcf header \
# 		$(STARK_FOLDER_BIN)/fix_vcf_header.sh --input=$@.tmp.merged --output=$@.tmp.merged --threads=$(THREADS_BY_SAMPLE) --bcftools=$(BCFTOOLS) $(FIX_VCF_HEADER_REFORMAT_option); \
# 		# Annotation; \
# 		$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.tmp.merged --output=$@.tmp.annotated0 --annotation=$(HOWARD_ANNOTATION_ANALYSIS) --norm=$$(cat $*.genome) --threads=$(THREADS); \
# 		# Prevent comma in description in vcf header \
# 		$(STARK_FOLDER_BIN)/fix_vcf_header.sh --input=$@.tmp.annotated0 --output=$@.tmp.annotated0 --threads=$(THREADS_BY_SAMPLE) --bcftools=$(BCFTOOLS) $(FIX_VCF_HEADER_REFORMAT_option); \
# 		# Annotation dejavu (forced); \
# 		if [ "$(HOWARD_DEJAVU_ANNOTATION)" != "" ]; then \
# 			$(HOWARD) $(HOWARD_DEJAVU_CONFIG_OPTIONS) --input=$@.tmp.annotated0 --output=$@.tmp.annotated1 --annotation=$(HOWARD_DEJAVU_ANNOTATION) --norm=$$(cat $*.genome) --threads=$(THREADS) --force; \
# 		else \
# 			mv @.tmp.annotated0 $@.tmp.annotated1; \
# 		fi; \
# 		# Prevent comma in description in vcf header \
# 		$(STARK_FOLDER_BIN)/fix_vcf_header.sh --input=$@.tmp.annotated1 --output=$@.tmp.annotated1 --threads=$(THREADS_BY_SAMPLE) --bcftools=$(BCFTOOLS) $(FIX_VCF_HEADER_REFORMAT_option); \
# 		# Calculation and prioritization (forced); \
# 		$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.tmp.annotated1 --output=$@.tmp.calculated.prioritized --calculation=$(HOWARD_CALCULATION_ANALYSIS) --prioritization=$(HOWARD_PRIORITIZATION_ANALYSIS) --nomen_fields=$(HOWARD_NOMEN_FIELDS) --force --threads=$(THREADS); \
# 		# Prevent comma in description in vcf header \
# 		$(STARK_FOLDER_BIN)/fix_vcf_header.sh --input=$@.tmp.calculated.prioritized --output=$@.tmp.calculated.prioritized --threads=$(THREADS_BY_SAMPLE) --bcftools=$(BCFTOOLS) $(FIX_VCF_HEADER_REFORMAT_option); \
# 		# Sort VCF; \
# 		mkdir -p $@.tmp.calculated.prioritized.SAMTOOLS_PREFIX; \
# 		$(BCFTOOLS) sort -T $@.tmp.calculated.prioritized.SAMTOOLS_PREFIX $@.tmp.calculated.prioritized > $@.tmp.calculated.prioritized.sorted; \
# 		rm -rf $@.tmp.calculated.prioritized.SAMTOOLS_PREFIX; \
# 		# Generate Design VCF; \
# 		$(BGZIP) --threads=$(THREADS) -c $@.tmp.calculated.prioritized.sorted > $@.tmp.calculated.prioritized.sorted.vcf.gz; \
# 		$(TABIX) $@.tmp.calculated.prioritized.sorted.vcf.gz; \
# 		cp $@.tmp.calculated.prioritized.sorted.vcf.gz $@.Design.vcf.gz; \
# 		$(TABIX) $@.Design.vcf.gz; \
# 		# Generate Design TSV; \
# 		$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.Design.vcf.gz --output=$@.Design.tsv --translation=TSV --fields="$(HOWARD_FIELDS)" --sort=$(HOWARD_SORT) --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)" --pzfields="PZScore,PZFlag,PZComment,PZInfos" --force --threads=$(THREADS); \
# 		# Generate Panel(s) VCF and TSV from Design VCF ($@.Design.vcf.gz); \
# 		for genes_file in $$(cat $$(cat $@.tmp.vcf_list | xargs dirname | sed s/reports$$/list.genes/) | cut -d. -f2- | sort -u); do \
# 			# List of Samples with $$genes_file panel \
# 			echo "genes_file for variants is: "$$genes_file; \
# 			List_of_samples=$$(ls $(@D)/[^.]*/[^.]*.$$genes_file | xargs -l basename | cut -d. -f1 | sort -u | tr "\n" "," | sed s/,$$// 2>/dev/null) ; \
# 			echo "List_of_samples is: "$$List_of_samples; \
# 			if [ "$$List_of_samples" != "" ] || ((1)); then \
# 				# List of $$genes_file within Samples folders \
# 				List_of_genes_files=$$(ls $$(for L in $$(echo $$List_of_samples | tr "," "\n"); do echo $(@D)/$$L/$$L.$$genes_file; done;)) ; \
# 				# Merge all $$genes_file found into uniq BED file (but supposed to be the same) \
# 				cat $$(echo $$List_of_genes_files) | $(BEDTOOLS) sort | $(BEDTOOLS) merge > $@.tmp.GENES.$$genes_file; \
# 				# Generate VCF Panel from VCF Design (especially $@.tmp.calculated.prioritized.sorted.vcf.gz because tabix) with $$genes_file for List of Samples \
# 				if (( $$(grep ^ $@.tmp.GENES.$$genes_file -c) )); then \
# 					$(BCFTOOLS) view --threads=$(THREADS)  --force-samples -U $@.tmp.calculated.prioritized.sorted.vcf.gz -R $@.tmp.GENES.$$genes_file > $@.Panel.$$genes_file.vcf; \
# 				else \
# 					$(BCFTOOLS) view --threads=$(THREADS)  --force-samples -U $@.tmp.calculated.prioritized.sorted.vcf.gz > $@.Panel.$$genes_file.vcf; \
# 				fi; \
# 				# Compress VCF \
# 				$(BGZIP) --threads=$(THREADS)  $@.Panel.$$genes_file.vcf; \
# 				$(TABIX) $@.Panel.$$genes_file.vcf.gz; \
# 				# Generate TSV Panel from VCF Panel compressed  \
# 				$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.Panel.$$genes_file.vcf.gz --output=$@.Panel.$$genes_file.tsv --translation=TSV --fields="$(HOWARD_FIELDS)" --sort=$(HOWARD_SORT) --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)" --pzfields="PZScore,PZFlag,PZComment,PZInfos" --force --threads=$(THREADS); \
# 			fi; \
# 		done; \
# 		rm -f $@.tmp.*; \
# 		echo "#[INFO] All variants files on Design and Panel(s) from full VCF are named $$(basename $@).*" > $@; \
# 	fi;



%.variants: $(VCF_REPORT_FILES) %.variants_full
	# List of final VCF files
	echo $^ | tr " " "\n" | tr "\t" "\n" | grep "final.vcf.gz$$" > $@.tmp.vcf_list
	$(BCFTOOLS) merge --threads=$(THREADS) -l $@.tmp.vcf_list -m none --force-samples | $(BCFTOOLS) filter --threads=$(THREADS) -S . -e 'GT=="0/0" | GT=="0|0"' > $@.tmp.merged;
	# Prevent comma in description in vcf header
	$(STARK_FOLDER_BIN)/fix_vcf_header.sh --input=$@.tmp.merged --output=$@.tmp.merged --threads=$(THREADS_BY_SAMPLE) --bcftools=$(BCFTOOLS) $(FIX_VCF_HEADER_REFORMAT_option);
	# Sort VCF
	mkdir -p $@.tmp.merged.SAMTOOLS_PREFIX
	$(BCFTOOLS) sort -T $@.tmp.merged.SAMTOOLS_PREFIX $@.tmp.merged > $@.tmp.calculated.prioritized.sorted
	rm -rf $@.tmp.merged.SAMTOOLS_PREFIX
	# Generate Design VCF
	$(BGZIP) --threads=$(THREADS) -c $@.tmp.calculated.prioritized.sorted > $@.tmp.calculated.prioritized.sorted.vcf.gz;
	$(TABIX) $@.tmp.calculated.prioritized.sorted.vcf.gz;
	cp $@.tmp.calculated.prioritized.sorted.vcf.gz $@.Design.vcf.gz;
	$(TABIX) $@.Design.vcf.gz; \
	# Generate Panel(s) VCF and TSV from Design VCF ($@.Design.vcf.gz)
	+for genes_file in $$(cat $$(cat $@.tmp.vcf_list | xargs dirname | sed s/reports$$/list.genes/) | cut -d. -f2- | sort -u); do \
		# List of Samples with $$genes_file panel \
		echo "genes_file for variants is: "$$genes_file; \
		List_of_samples=$$(ls $(@D)/[^.]*/[^.]*.$$genes_file | xargs -l basename | cut -d. -f1 | sort -u | tr "\n" "," | sed s/,$$// 2>/dev/null) ; \
		echo "List_of_samples is: "$$List_of_samples; \
		if [ "$$List_of_samples" != "" ] || ((1)); then \
			# List of $$genes_file within Samples folders \
			List_of_genes_files=$$(ls $$(for L in $$(echo $$List_of_samples | tr "," "\n"); do echo $(@D)/$$L/$$L.$$genes_file; done;)) ; \
			# Merge all $$genes_file found into uniq BED file (but supposed to be the same) \
			cat $$(echo $$List_of_genes_files) | $(BEDTOOLS) sort | $(BEDTOOLS) merge > $@.tmp.GENES.$$genes_file; \
			# Generate VCF Panel from VCF Design (especially $@.tmp.calculated.prioritized.sorted.vcf.gz because tabix) with $$genes_file for List of Samples \
			if (( $$(grep ^ $@.tmp.GENES.$$genes_file -c) )); then \
				$(BCFTOOLS) view --threads=$(THREADS) --samples $$List_of_samples -U --force-samples $@.tmp.calculated.prioritized.sorted.vcf.gz -R $@.tmp.GENES.$$genes_file > $@.Panel.$$genes_file.vcf; \
			else \
				$(BCFTOOLS) view --threads=$(THREADS) --samples $$List_of_samples -U --force-samples $@.tmp.calculated.prioritized.sorted.vcf.gz > $@.Panel.$$genes_file.vcf; \
			fi; \
			# Compress VCF \
			$(BGZIP) --threads=$(THREADS) $@.Panel.$$genes_file.vcf; \
			$(TABIX) $@.Panel.$$genes_file.vcf.gz; \
		fi; \
	done;
	-rm -f $@.tmp.*
	echo "#[INFO] All variants files on Design and Panel(s) are named $$(basename $@).*" > $@;



%.variants_full: $(VCF_REPORT_FILES)
	+if (($(REPORT_VARIANTS_FULL))); then \
		# List of full VCF files \
		echo $^ | tr " " "\n" | tr "\t" "\n" | grep "full.vcf.gz$$" > $@.tmp.vcf_list; \
		for f in $$(cat $@.tmp.vcf_list); do $(BCFTOOLS) view --threads=$(THREADS) -h $$f | grep "^#CHROM" | cut -f10- | sed  "s/^/$$(basename $$f | cut -d. -f1)./g;s/\t/\t$$(basename $$f | cut -d. -f1)./g";  done | tr "\t" "\n" > $@.tmp.vcf_list.pipelines; \
		$(BCFTOOLS) merge --threads=$(THREADS) -l $@.tmp.vcf_list -m none --force-samples | $(BCFTOOLS) filter --threads=$(THREADS) -S . -e 'GT=="0/0" | GT=="0|0"' | $(BCFTOOLS) reheader --threads=$(THREADS) -s $@.tmp.vcf_list.pipelines > $@.tmp.merged; \
		# Prevent comma in description in vcf header \
		$(STARK_FOLDER_BIN)/fix_vcf_header.sh --input=$@.tmp.merged --output=$@.tmp.merged --threads=$(THREADS_BY_SAMPLE) --bcftools=$(BCFTOOLS) $(FIX_VCF_HEADER_REFORMAT_option); \
		# Sort VCF; \
		mkdir -p $@.tmp.merged.SAMTOOLS_PREFIX; \
		$(BCFTOOLS) sort -T $@.tmp.merged.SAMTOOLS_PREFIX $@.tmp.merged > $@.tmp.calculated.prioritized.sorted; \
		rm -rf $@.tmp.merged.SAMTOOLS_PREFIX; \
		# Generate Design VCF; \
		$(BGZIP) --threads=$(THREADS) -c $@.tmp.calculated.prioritized.sorted > $@.tmp.calculated.prioritized.sorted.vcf.gz; \
		$(TABIX) $@.tmp.calculated.prioritized.sorted.vcf.gz; \
		cp $@.tmp.calculated.prioritized.sorted.vcf.gz $@.Design.vcf.gz; \
		$(TABIX) $@.Design.vcf.gz; \
		# Generate Panel(s) VCF and TSV from Design VCF ($@.Design.vcf.gz); \
		for genes_file in $$(cat $$(cat $@.tmp.vcf_list | xargs dirname | sed s/reports$$/list.genes/) | cut -d. -f2- | sort -u); do \
			# List of Samples with $$genes_file panel \
			echo "genes_file for variants is: "$$genes_file; \
			List_of_samples=$$(ls $(@D)/[^.]*/[^.]*.$$genes_file | xargs -l basename | cut -d. -f1 | sort -u | tr "\n" "," | sed s/,$$// 2>/dev/null) ; \
			echo "List_of_samples is: "$$List_of_samples; \
			if [ "$$List_of_samples" != "" ] || ((1)); then \
				# List of $$genes_file within Samples folders \
				List_of_genes_files=$$(ls $$(for L in $$(echo $$List_of_samples | tr "," "\n"); do echo $(@D)/$$L/$$L.$$genes_file; done;)) ; \
				# Merge all $$genes_file found into uniq BED file (but supposed to be the same) \
				cat $$(echo $$List_of_genes_files) | $(BEDTOOLS) sort | $(BEDTOOLS) merge > $@.tmp.GENES.$$genes_file; \
				# Generate VCF Panel from VCF Design (especially $@.tmp.calculated.prioritized.sorted.vcf.gz because tabix) with $$genes_file for List of Samples \
				if (( $$(grep ^ $@.tmp.GENES.$$genes_file -c) )); then \
					$(BCFTOOLS) view --threads=$(THREADS)  --force-samples -U $@.tmp.calculated.prioritized.sorted.vcf.gz -R $@.tmp.GENES.$$genes_file > $@.Panel.$$genes_file.vcf; \
				else \
					$(BCFTOOLS) view --threads=$(THREADS)  --force-samples -U $@.tmp.calculated.prioritized.sorted.vcf.gz > $@.Panel.$$genes_file.vcf; \
				fi; \
				# Compress VCF \
				$(BGZIP) --threads=$(THREADS)  $@.Panel.$$genes_file.vcf; \
				$(TABIX) $@.Panel.$$genes_file.vcf.gz; \
			fi; \
		done; \
		rm -f $@.tmp.*; \
		echo "#[INFO] All variants files on Design and Panel(s) from full VCF are named $$(basename $@).*" > $@; \
	fi;

