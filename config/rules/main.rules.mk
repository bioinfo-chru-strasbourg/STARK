############################
# Main Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.5.1b"
MK_DATE="27/09/2019"

# Release note
# 10/03/2015-0.9.4beta: change genome reference location, in the file %.genome
# 14/06/2015-0.9.4.1b: change interval from bed (grep to remove '@' lines)
# 25/11/2015-0.9.4.2b: Change VCF sorting, using vcftools instead of igvtools
# 25/01/2016-0.9.4.3b: Change VCF sorting vcftools parameter, adding -C chromosome ordering
# 18/02/2016-0.9.4.4b: change cp to cat form manifest copy
# 13/04/2016-0.9.4.5b: removing full.vcf rule
# 29/09/2016-0.9.4.6b: Update PICARD to releasr picard.jar
# 07/05/2018-0.9.4.7b: Add --force for translation vcf to txt
# 02/10/2018-0.9.5b: Change HOWARD translation, prioritization and hard filtering. Change Manifest/bed link generation
# 27/09/2019-0.9.5.1b: Change FATBAM to CAP tool

# TOOLS
IGVTOOLS?=$(NGSbin)/igvtools.jar
SAMTOOLS?=$(NGSbin)/samtools
JAVA?=java
TABIX?=$(NGSbin)/tabix
BWA?=$(NGSbin)/bwa
PICARDLIB?=$(NGSbin)/picard-tools
FASTQC?=$(NGSbin)/fastqc
BGZIP?=bgzip
GZ?=gzip

# HOWARD Prioritization
HOWARD_FILTER?="default"
HOWARD_CONFIG?="config.ini"
HOWARD_CONFIG_PRIORITIZATION?="config.prioritization.ini"
HOWARD_CONFIG_ANNOTATION?="config.annotation.ini"
#HOWARD_CONFIG_FILTER?=$(HOWARD_CONFIG_PRIORITIZATION)

# HOWARD Translation
HOWARD_ANNOTATION?="PZScore,PZFlag,PZComment,Symbol,hgvs,location,outcome,AlleleFrequency,AD,dbSNP,dbSNPNonFlagged,popfreq,database_CPSGEN,DP,AF,VF,GQ,Ensembl,TI,FC,GWASCatalog,COSMIC,LocalDB,1000genomesALL,1000genomesEUR,6500NHLBIALL,6500NHLBIEUR,PolyPhen2HumanVarPred,PolyPhen2HumanDivPred,MutationTasterPred,MutationAssessorPred,LTRPred,IARCTP53,SIFT,phastCons,PhyloP,SiPhy,FATHMM,LRT,GERP,PolyPhen2HumanVar,PolyPhen2HumanDiv,MutationTaster,MutationAssessor,TFBS,FilterComment,ALL"
#HOWARD_FIELDS?="PZScore,PZFlag,PZComment,Symbol,hgvs,location,outcome,AlleleFrequency,AD,dbSNP,dbSNPNonFlagged,popfreq,ALL"
HOWARD_FIELDS?="NOMEN,PZFlag,PZScore,PZComment,CNOMEN,PNOMEN,location,outcome,VAF_average,dbSNP,dbSNPNonFlagged,popfreq,ALL"
HOWARD_SORT?="PZFlag::DESC,PZScore:n:DESC"
HOWARD_SORT_BY?="PZFlag,PZScore"
HOWARD_ORDER_BY?="DESC,DESC"

# OPTIONS
REMOVE_INTERMEDIATE_SAM?=1
BED?=
PRIMER_BED?=
THREADS_SAMTOOLS?=$(THREADS_BY_SAMPLE)

BAM_COMPRESSION?=5
GZ?=gzip

## VCF files
##############

# VCF Compression
%.vcf.gz: %.vcf %.empty.vcf
	# If no VCF or empty file, create an empty VCF
	if [ ! -s $< ]; then cp $*.empty.vcf $<; fi;
	# VCF Sorting
	mkdir -p $@.SAMTOOLS_PREFIX
	$(BCFTOOLS) sort -T $@.SAMTOOLS_PREFIX $< > $@.sorted
	rm -rf $@.SAMTOOLS_PREFIX
	# VCF compression with BGZIP
	$(BGZIP) -f $@.sorted -c > $@
	rm -f $@.sorted
	# remove files
	-rm -f $*.empty.vcf*

# VCF Indexing with TABIX
%.vcf.gz.tbi: %.vcf.gz
	# Remove index if exists
	-rm -f $@
	# Indexing with TABIX
	$(TABIX) -f -p vcf $<

# VCF Indexing with IGVTOOLS
%.vcf.idx: %.vcf %.empty.vcf
	# If no VCF or empty file, create an empty VCF
	-if [ ! -s $< ]; then cp $*.empty.vcf $<; fi;
	# Indexing with IGVTOOLS
	-$(JAVA) -Xmx1g -jar $(IGVTOOLS) index $<
	# Empty index if fail
	if [ ! -e $@ ]; then touch $@; fi;
	# remove files
	-rm -f $*.empty.vcf*

# EMPTY VCF
%.empty.vcf:
	# Header
	echo "##fileformat=VCFv4.1" > $@
	# Head first line
	echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	"$$(echo $$(basename $@) | awk -F"." '{print $$1}') >> $@


# VCF SORTING
%.vcf: %.sorting.vcf
	grep "^#" $< > $@
	mkdir -p $<"_VCF_sort"
	grep -v "^#" $< | sort -k1,1V -k2,2n -T $<"_VCF_sort" >> $@
	rm -rf $<"_VCF_sort"


# VCF NORMALIZATION with BCFTOOLS
%.vcf: %.normalization.vcf %.genome
	$(BCFTOOLS) norm -m- -f `cat $*.genome` $< | $(BCFTOOLS) norm -d all > $@


# MERGE SNP and InDel VCF
%.vcf: %.SNP.vcf.gz %.InDel.vcf.gz %.SNP.vcf.gz.tbi %.InDel.vcf.gz.tbi
	# CONCAT
	$(BCFTOOLS) concat -a $*.SNP.vcf.gz $*.InDel.vcf.gz > $@.tmp
	# Sorting
	grep "^#" $@.tmp > $@
	mkdir -p $@.tmp"_VCF_sort"
	grep -v "^#" $@.tmp | sort -k1,1V -k2,2n -T $@.tmp"_VCF_sort" >> $@
	rm -rf $@.tmp*
	# Cleaning
	-rm -f $*.SNP.vcf.gz $*.InDel.vcf.gz
	-rm -f $*.SNP.vcf.gz.tbi $*.InDel.vcf.gz.tbi
	-rm -f $*.SNP.vcf.idx $*.InDel.vcf.idx
	-rm -f $*.idx


# HOWARD
#####################
# Translation, prioritiazzation, hard filtering


# VCF to tab delimiter
%.tsv: %.vcf
	# translation step
	$(HOWARD) --config=$(HOWARD_CONFIG) --config=$(HOWARD_CONFIG) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION) --config_annotation=$(HOWARD_CONFIG_ANNOTATION) --pzfields="PZScore,PZFlag,PZComment,PZInfos" --format=tab  --fields="$(HOWARD_FIELDS)" --sort=$(HOWARD_SORT) --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)"  --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS) --input=$< --output=$@ --force;
	# Touch
	if [ ! -e $@ ]; then touch $@; fi;
	# Cleaning

# Hard filtering
%.hard.tsv: %.vcf
	# translation step and hard filtering
	$(HOWARD) --config=$(HOWARD_CONFIG) --config=$(HOWARD_CONFIG) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION) --config_annotation=$(HOWARD_CONFIG_ANNOTATION) --prioritization=$(HOWARD_PRIORITIZATION) --pzfields="PZScore,PZFlag,PZComment,PZInfos" --format=tab  --fields=$(HOWARD_FIELDS) --sort=$(HOWARD_SORT) --sort_by=$(HOWARD_SORT_BY) --order_by=$(HOWARD_ORDER_BY)  --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS) --input=$< --output=$@ --hard --force;
	# Touch
	if [ ! -e $@ ]; then touch $@; fi;
	# Cleaning

# VCF to tab delimiter
%.txt: %.vcf
	# Translation step
	$(HOWARD) --config=$(HOWARD_CONFIG) --config=$(HOWARD_CONFIG) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION) --config_annotation=$(HOWARD_CONFIG_ANNOTATION) --prioritization=$(HOWARD_PRIORITIZATION) --pzfields="PZScore,PZFlag,PZComment,PZInfos" --format=tab  --fields=$(HOWARD_FIELDS) --sort=$(HOWARD_SORT) --sort_by=$(HOWARD_SORT_BY) --order_by=$(HOWARD_ORDER_BY)  --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS) --input=$< --output=$@ --force;
	# Touch
	if [ ! -e $@ ]; then touch $@; fi;
	# Cleaning


# Hard filtering
%.hard.txt: %.vcf
	# Translation and hard filtering
	$(HOWARD) --config=$(HOWARD_CONFIG) --config=$(HOWARD_CONFIG) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION) --config_annotation=$(HOWARD_CONFIG_ANNOTATION) --prioritization=$(HOWARD_PRIORITIZATION) --pzfields="PZScore,PZFlag,PZComment,PZInfos" --format=tab  --fields=$(HOWARD_FIELDS) --sort=$(HOWARD_SORT) --sort_by=$(HOWARD_SORT_BY) --order_by=$(HOWARD_ORDER_BY)  --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS) --input=$< --output=$@ --hard --force;
	# Touch
	if [ ! -e $@ ]; then touch $@; fi;
	# Cleaning


%.prioritized.vcf: %.vcf
	# Prioritization step
	$(HOWARD) --input=$< --output=$@ --config=$(HOWARD_CONFIG) --config=$(HOWARD_CONFIG) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION) --config_annotation=$(HOWARD_CONFIG_ANNOTATION) --filter=$(HOWARD_PRIORITIZATION)  --format=tab --annotation=$(HOWARD_ANNOTATION) --sort=$(HOWARD_SORT) --sort_by=$(HOWARD_SORT_BY) --env=$(CONFIG_TOOLS) --order_by=$(HOWARD_ORDER_BY)



## BAM/FASTQ Files
####################

# BAM Indexing
%.bam.bai: %.bam
	$(SAMTOOLS) index $<

# BAM from SAM
# sorting sam file by coordinate and output a bam file
%.bam : %.sam %.genome
	if ((! $$($(SAMTOOLS) view -H $< | grep "^@HD.*VN:.*SO:" | wc -l))) || (($$($(SAMTOOLS) view -H $< | grep "^@HD.*VN:.*SO:unsorted" | wc -l))) ; then \
		$(SAMTOOLS) sort $< -o $@ -O BAM -l 1 -T $<.SAMTOOLS_PREFIX -@ $(THREADS_SAMTOOLS); \
	else \
		$(SAMTOOLS) view -o $@ -b -1 -S -T `cat $*.genome` $< -@ $(THREADS_SAMTOOLS); \
	fi;
	-if [ $(REMOVE_INTERMEDIATE_SAM) -eq 1 ]; then rm -f $<; fi;
	#rm -f $<


# CRAM from SAM
# sorting sam file by coordinate and output a bam file
# BAM from CRAM
# sorting sam file by coordinate and output a bam file
%.cram: %.bam %.genome
	echo "test BAM to CRAM: $^"
	$(SAMTOOLS) view -o $@ -O CRAM -S -T `cat $*.genome` $*.bam -@ $(THREADS_SAMTOOLS);
	# test Empty output file
	# Remove intermediate SAM file


# BAM Sorting
#%.bam: %.uncompressed.bam
#	$(SAMTOOLS) sort $< -o $@ -l $(BAM_COMPRESSION) -@ $(THREADS_SAMTOOLS);

# BAM Sorting
%.bam: %.compress.bam
	$(SAMTOOLS) sort $< -o $@ -l $(BAM_COMPRESSION) -T $<.SAMTOOLS_PREFIX -@ $(THREADS_SAMTOOLS);


%.bam : %.sorting.bam
	if ((! $$($(SAMTOOLS) view -H $< | grep "^@HD.*VN:.*SO:" | wc -l))) || (($$($(SAMTOOLS) view -H $< | grep "^@HD.*VN:.*SO:unsorted" | wc -l))) ; then \
		$(SAMTOOLS) sort $< -o $@ -T $<.SAMTOOLS_PREFIX -@ $(THREADS_SAMTOOLS); \
	else \
		mv $< $@; \
	fi;
	rm -f $<


# SAM from FASTQ
%.sai : %.fastq.gz %.genome
	$(BWA) aln -f $@ `cat $*.genome` $<

# FASTQ Compression with GZIP
%.fastq.gz: %.fastq
	$(GZ) --fast $< -c > $@

# FASTQ(s) from BAM
%.R1.fastq %.R2.fastq: %.bam
	#$(JAVA) -jar $(PICARDLIB)/SamToFastq.jar INPUT=$< FASTQ=$*.R1.fastq SECOND_END_FASTQ=$*.R2.fastq
	$(JAVA) -jar $(PICARD) SamToFastq INPUT=$< FASTQ=$*.R1.fastq SECOND_END_FASTQ=$*.R2.fastq


# BAM reduction
GATKRR_FLAGS=
%.reduced.bam: %.bam %.bam.bai %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKRR_FLAGS) -T ReduceReads -R `cat $*.genome` -I $< -o $@



# SampleSheet Copy
%.SampleSheet.csv:
	-mkdir -p $(@D)
	-cp -p $(@D)/`echo $$(basename $(@D))`.SampleSheet.csv $@ || cp -p $(INPUTDIR)/`echo $$(basename $$(dirname $(@D)))`/SampleSheet.csv $@ || touch $@;



# BED / INTERVALS
###################


%.manifest: %.manifest_from_samplesheet
	-mkdir -p $(@D)
	# Create MANIFEST file
	# 1. test if main manifest file (SAMPLE.manifest) exists (use for Sample analysis)
	# 2. test if a manifest is found on SampleSheet (use for Run analysis)
	# 3. test if a manifest is defined in the APP (as an environment variable)
	#
	-if [ -e $(@D)/`echo $$(basename $(@D))`.manifest ]; then \
		echo "# MANIFEST for the sample '$(@D)/`echo $$(basename $(@D))`.manifest' exists" ; \
		if  [ "$(@D)/`echo $$(basename $(@D))`.manifest" != "$@" ]; then \
			echo "# MANIFEST for the sample '$(@D)/`echo $$(basename $(@D))`.manifest' exists and is different than the wanted MANIFEST file '$@'" ; \
			ln -s $(@D)/`echo $$(basename $(@D))`.manifest $@; \
			if [ ! -e $@ ]; then rm -f $@; cp $(@D)/`echo $$(basename $(@D))`.manifest $@; fi; \
		fi; \
	elif [ -s $*.manifest_from_samplesheet ]; then \
		echo "# MANIFEST for the SampleSheet found in the manifest" ; \
		cp $*.manifest_from_samplesheet $@; \
	elif [ -e "$(MANIFEST)" ] && [ "$(MANIFEST)" != "" ]; then \
		echo "# Input MANIFEST '$(MANIFEST)' exists" ; \
		cp $(MANIFEST) $@; \
	else \
		touch $@; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;
	# Clean


%.manifest_name: %.manifest_from_samplesheet_name
	-mkdir -p $(@D)
	# Create MANIFEST file
	# 1. test if main manifest name file (SAMPLE.manifest_name) exists (use for Sample analysis)
	# 2. test if a manifest name is found on SampleSheet (use for Run analysis)
	#
	-if [ -e $(@D)/`echo $$(basename $(@D))`.manifest_name ]; then \
		echo "# MANIFEST name for the sample '$(@D)/`echo $$(basename $(@D))`.manifest_name' exists" ; \
		if  [ "$(@D)/`echo $$(basename $(@D))`.manifest_name" != "$@" ]; then \
			echo "# MANIFEST name for the sample '$(@D)/`echo $$(basename $(@D))`.manifest_name' exists and is different than the wanted MANIFEST name file '$@'" ; \
			ln -s $(@D)/`echo $$(basename $(@D))`.manifest_name $@; \
			if [ ! -e $@ ]; then rm -f $@; cp $(@D)/`echo $$(basename $(@D))`.manifest_name $@; fi; \
		fi; \
	elif [ -e $(@D)/../`echo $$(basename $(@D))`.manifest_name ]; then \
		echo "# MANIFEST name for the sample '$(@D)/../`echo $$(basename $(@D))`.manifest_name' exists" ; \
		if  [ "$(@D)/../`echo $$(basename $(@D))`.manifest_name" != "$@" ]; then \
			echo "# MANIFEST name for the sample '$(@D)/../`echo $$(basename $(@D))`.manifest_name' exists and is different than the wanted MANIFEST name file '$@'" ; \
			ln -s $(@D)/../`echo $$(basename $(@D))`.manifest_name $@; \
			if [ ! -e $@ ]; then rm -f $@; cp $(@D)/../`echo $$(basename $(@D))`.manifest_name $@; fi; \
		fi; \
	elif [ -s $*.manifest_from_samplesheet_name ]; then \
		echo "# MANIFEST name for the SampleSheet found in the manifest" ; \
		cp $*.manifest_from_samplesheet_name $@; \
	elif [ -e "$(MANIFEST)" ] && [ "$(MANIFEST)" != "" ]; then \
		echo "# Input MANIFEST '$(MANIFEST)' exists" ; \
		echo -e "$(MANIFEST)\tfrom APP" > $@; \
	else \
		touch $@; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;
	# Clean



%.manifest_from_samplesheet: %.SampleSheet.csv #%.manifests_list.txt
	# 1. Create the sample Manifest
	# Find the index of Manifest field
	# Find the Manifest index (A,B...)
	# Find the Manifest (blabla.manifest)
	# copy the manifest into the sample folder
	# clean
	# Create directory if needed
	-mkdir -p $(@D)
	# SampleSheet empty
	if [ ! -e $< ]; then touch $<; fi;
	# Manifest List
	-if [ -s $(INPUTDIR)/`echo $$(basename $$(dirname $(@D)))`/manifests_list.txt ]; then \
		cp -p $(INPUTDIR)/`echo $$(basename $$(dirname $(@D)))`/manifests_list.txt $*.manifests_list.txt; \
	else \
		touch $*.manifests_list.txt; \
	fi;
	if [ ! -e $@ ]; then touch $*.manifests_list.txt; fi;
	#
	# Column number of sample manifest
	-grep -i ^Sample_ID $< | tr -d '\r\n' | sed -e 's/,/\n/g' | grep -i -n Manifest | cut -d \: -f 1 > $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I;
	-if [ ! -s $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I ]; then echo "1" > $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I; fi;
	# Manifest letter for the index
	-grep -e ^$$(basename $(@D)), $< | tail -n 1 | tr -d '\r\n' | cut -d \, -f `cat $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I` > $(@D)/$(*F).SAMPLE_MANIFEST_I;
	#
	# Manifest name for the sample
	-grep ^`cat $(@D)/$(*F).SAMPLE_MANIFEST_I`, $*.manifests_list.txt | cut -d \, -f 2 > $(@D)/$(*F).SAMPLE_MANIFEST;
	#
	# BED
	file=$$( echo "$$( dirname $(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST` )/$$( basename $(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST` ).bed" ); \
	echo "bed file is : $$file "; \
	if [ -s $$file ]; then \
		cat $$file > $*.bed; \
	else \
		echo "$$file don't exist !!!"; \
	fi;
	#
	# BED GENES
	file=$$( echo "$$( dirname $(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST` )/$$( basename $(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST` ).genes" ); \
	echo "genes file is : $$file "; \
	if [ -s $$file ]; then \
		cat $$file > $*.genes; \
	else \
		echo "$$file don't exist !!!"; \
	fi;
	#
	# TRANSCRIPTS
	file=$$( echo "$$( dirname $(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST` )/$$( basename $(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST` ).transcripts" ); \
	echo "transcripts file is : $$file "; \
	if [ -s $$file ]; then \
		cat $$file > $*.transcripts; \
	else \
		echo "$$file don't exist !!!"; \
	fi;
	#
	# Found manifest, or Default manifest is the first in the list
	if [ "`cat $(@D)/$(*F).SAMPLE_MANIFEST`" != "" ] && [ -e "$(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST`" ]; then \
		echo "# TEST Found Manifest '"`cat $(@D)/$(*F).SAMPLE_MANIFEST`"' in '$(MANIFEST_FOLDER)' for sample `echo $$(basename $$(dirname $(@D)))`/$(*F)"; \
		cat "$(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST`" > $@; \
	elif [ "`cut -d, -f2 $*.manifests_list.txt`" != "" ] && [ -e "$(MANIFEST_FOLDER)/`cut -d, -f2 $*.manifests_list.txt`" ]; then \
		echo "# TEST Manifest NOT found in '$(MANIFEST_FOLDER)'. Default Manifest used '"`cut -d, -f2 $*.manifests_list.txt`"' for sample `echo $$(basename $$(dirname $(@D)))`/$(*F)"; \
		cat "$(MANIFEST_FOLDER)/`cut -d, -f2 $*.manifests_list.txt`" > $@; \
	fi;
	#
	# IF root manifest exists
	if [ -s $(@D)/`echo $$(basename $(@D))`.manifest ] && [ "$(@D)/`echo $$(basename $(@D))`.manifest" != "$@" ]; then \
		cp -p $(@D)/`echo $$(basename $(@D))`.manifest $@; \
	fi;
	#
	# Empty manifest if failed!
	if [ ! -e $@ ]; then touch $@; fi;
	# Touch manifest to use time of SampleSheet
	if [ -e $< ]; then \
		touch $@ -r $<; \
	fi;
	for manifest in `cut -d"," -f1 $*.manifests_list.txt`; do grep ",$$manifest," $< | cut -d"," -f1 | sed "s/$$/,$$manifest/g" >> $*.manifest_by_samples.txt; echo "-----" >> $*.manifest_by_samples.txt; done; \
	# Clean
	-rm -f $*.manifests_list.txt
	# Remove intermediate files
	-rm -f $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I $(@D)/$(*F).SAMPLE_MANIFEST_I $(@D)/$(*F).SAMPLE_MANIFEST;


%.manifest_from_samplesheet_name: %.SampleSheet.csv
	# 1. Create the sample Manifest
	# Find the index of Manifest field
	# Find the Manifest index (A,B...)
	# Find the Manifest (blabla.manifest)
	# copy the manifest into the sample folder
	# clean
	# Create directory if needed
	-mkdir -p $(@D)
	# SampleSheet empty
	if [ ! -e $< ]; then touch $<; fi;
	# Manifest List
	-if [ -s $(INPUTDIR)/`echo $$(basename $$(dirname $(@D)))`/manifests_list.txt ]; then \
		cp -p $(INPUTDIR)/`echo $$(basename $$(dirname $(@D)))`/manifests_list.txt $*.manifests_list_name.txt; \
	else \
		touch $*.manifests_list_name.txt; \
	fi;
	if [ ! -e $@ ]; then touch $*.manifests_list_name.txt; fi;
	#
	# Column number of sample manifest
	-grep -i ^Sample_ID $< | tr -d '\r\n' | sed -e 's/,/\n/g' | grep -i -n Manifest | cut -d \: -f 1 > $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I_name;
	-if [ ! -s $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I_name ]; then echo "1" > $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I_name; fi;
	# Manifest letter for the index
	-grep -e ^$$(basename $(@D)), $< | tail -n 1 | tr -d '\r\n' | cut -d \, -f `cat $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I_name` > $(@D)/$(*F).SAMPLE_MANIFEST_I_name;
	#
	# Manifest name for the sample
	-grep ^`cat $(@D)/$(*F).SAMPLE_MANIFEST_I_name`, $*.manifests_list_name.txt | cut -d \, -f 2 > $(@D)/$(*F).SAMPLE_MANIFEST_name;
	#
	-if [ ! -s $(@D)/$(*F).SAMPLE_MANIFEST_name ]; then \
		> $@; \
	else \
		echo -e $$(cat $(@D)/$(*F).SAMPLE_MANIFEST_name)"\tfrom SampleSheet" > $@; \
	fi;
	#
	# Empty manifest if failed!
	if [ ! -e $@ ]; then touch $@; fi;
	#
	# Clean
	-rm -f $*.manifests_list_name.txt
	# Remove intermediate files
	-rm -f $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I_name $(@D)/$(*F).SAMPLE_MANIFEST_I_name $(@D)/$(*F).SAMPLE_MANIFEST_name;



# Interval from bed from manifest?
%.from_manifest.intervals: %.bed %.bam %.bam.bai %.bam.bed %.dict
	# manifest to interval (not needed?)
	#cat $<  | tr -d '\r' | sed -e "s/^M//" | awk -F"\t" '{print $$1":"$$2"-"$$3}' > $@
	# try to extract from the bam if exists, in order to not call in the whome genome
	-+if [ -s $< ]; then \
		cp $< $@.bed; \
	elif [ -s $*.bam ]; then \
		echo "# Generation of BED from BAM because bed/manifest is empty"; \
		cp $*.bam.bed $@.bed ; \
		if ((0)); then \
		if (($$($(SAMTOOLS) idxstats $*.bam | awk '{SUM+=$$3+$$4} END {print SUM}'))); then \
			rm -f $*.bam.genomeCoverageBed_for_intervals.mk $*.bam.genomeCoverageBed_for_intervals1.mk $*.bam.genomeCoverageBed_for_intervals2.mk $*.bam.genomeCoverageBed_for_intervals3.mk; \
			for chr in $$($(SAMTOOLS) idxstats $*.bam | grep -v "\*" | awk '{ if ($$3+$$4>0) print $$1 }'); do \
				echo "$*.bam.genomeCoverageBed_for_intervals.$$chr.bed: $*.bam" >> $*.bam.genomeCoverageBed_for_intervals1.mk; \
				echo "	$(SAMTOOLS) view $*.bam -b $$chr | $(BEDTOOLS) genomecov -ibam stdin -bg | $(BEDTOOLS) merge -i - > $*.bam.genomeCoverageBed_for_intervals.$$chr.bed " >> $*.bam.genomeCoverageBed_for_intervals1.mk; \
				echo -n " $*.bam.genomeCoverageBed_for_intervals.$$chr.bed" >> $*.bam.genomeCoverageBed_for_intervals2.mk; \
			done; \
			echo -n "$@.bed: " | cat - $*.bam.genomeCoverageBed_for_intervals2.mk > $*.bam.genomeCoverageBed_for_intervals3.mk; \
			echo ""  >> $*.bam.genomeCoverageBed_for_intervals3.mk; \
			echo "	cat $$^ " >> $*.bam.genomeCoverageBed_for_intervals3.mk; \
			echo "	cat $$^ > $@.bed " >> $*.bam.genomeCoverageBed_for_intervals3.mk; \
			cat $*.bam.genomeCoverageBed_for_intervals1.mk $*.bam.genomeCoverageBed_for_intervals3.mk >> $*.bam.genomeCoverageBed_for_intervals.mk; \
			make -i -f $*.bam.genomeCoverageBed_for_intervals.mk $@.bed -j $(THREADS) ; \
			rm $*.bam.genomeCoverageBed_for_intervals*; \
		fi; \
		fi; \
	else \
		echo "[ERROR] No intervals generated '$@'" ; \
	fi;
	#
	# INTERVAL WITH PICARD
	if [ -s $@.bed ]; then \
		echo "[INFO] Generate $@ from $@.bed with PICARD BedToIntervalList" ; \
		cut $@.bed -f1-3,5 > $@.bed.4fields ; \
		$(JAVA) -jar $(PICARD) BedToIntervalList I=$@.bed.4fields O=$@ SD=$$(cat $*.dict) ; \
		rm $@.bed.4fields ; \
	fi;
	# If error, try intervals with GREP/SED/AWK
	if [ ! -s $@ ] && [ -s $@.bed ]; then \
		echo "[INFO] Generate $@ from $@.bed with GREP/SED" ; \
		grep -v ^@ $@.bed  | tr -d '\r' | sed -e "s/^M//" | awk -F"\t" '{print $$1":"$$2"-"$$3}' | sed s/:0-/:1-/gi > $@; \
	fi;
	# touch
	if [ ! -e $@ ]; then touch $@; fi;
	# clean


# Interval from BED
%.bed.intervals: %.bed %.dict
	# BED to Intervals (not needed?)
	# INTERVAL WITH PICARD
	if [ -s $< ]; then \
		cut $< -f1-3,5 > $@.4fields ; \
		$(JAVA) -jar $(PICARD) BedToIntervalList I=$@.4fields O=$@ SD=$$(cat $*.dict) ; \
		rm $@.4fields ; \
	fi;
	# If error, try intervals with GREP/SED/AWK
	if [ ! -s $@ ] && [ -s $< ]; then \
		grep -v ^@ $<  | tr -d '\r' | sed -e "s/^M//" | awk -F"\t" '{print $$1":"$$2"-"$$3}' | sed s/:0-/:1-/gi > $@; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;

# BED primer file from a Manifest
%.primers.bed: %.manifest
	# Create BED file
	-if [ -e "$(PRIMER_BED)" ]; then \
		cp $(PRIMER_BED) $@; \
	elif [ -s $< ]; then \
		$(CAP_ManifestToBED) --input=$< --output=$@ --output_type=primer; \
	else \
		touch $@; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;

# BED file from a Manifest
%.bed: %.manifest
	# Create BED file
	# 1. test if main bed file (SAMPLE.bed) exists (use for Sample analysis)
	# 2. test if a bed can be generated from the manifest
	# 3. test if a manifest is defined in the APP (as an environment variable)
	#
	-rm -f $@
	-if [ -e $(@D)/`echo $$(basename $(@D))`.bed ]; then \
		echo "# BED for the sample '$(@D)/`echo $$(basename $(@D))`.bed' exists" ; \
		if  [ "$(@D)/`echo $$(basename $(@D))`.bed" != "$@" ]; then \
			echo "# BED for the sample '$(@D)/`echo $$(basename $(@D))`.bed' exists and is different than the wanted BED file '$@'" ; \
			ln -s $(@D)/`echo $$(basename $(@D))`.bed $@; \
			if [ ! -e $@ ]; then rm -f $@; cp $(@D)/`echo $$(basename $(@D))`.bed $@; fi; \
		fi; \
	elif [ -s $< ]; then \
		echo "# BED for the sample generated from the manifest '$<'" ; \
		rm -f $@.tmp $@.sorted.tmp; \
		$(CAP_ManifestToBED) --input=$< --output=$@.tmp --output_type=region_clipped; \
		$(BEDTOOLS) sort -i $@.tmp | $(BEDTOOLS) merge -c 4 -o collapse  | awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t+\t"$$4}' > $@; \
		rm -f $@.tmp; \
	elif [ -e "$(BED)" ] && [ "$(BED)" != "" ]; then \
		echo "# Input BED '$(BED)' exists" ; \
		cp $(BED) $@; \
	else \
		touch $@; \
	fi;
	#
	if [ ! -e $@ ]; then touch $@; fi;
	# Clean


# BED file from a Manifest
%.bed_name: %.manifest_name
	# Create BED file
	# 1. test if main bed file (SAMPLE.bed) exists (use for Sample analysis)
	# 2. test if a bed can be generated from the manifest
	# 3. test if a manifest is defined in the APP (as an environment variable)s
	#
	-if [ -e $(@D)/`echo $$(basename $(@D))`.bed_name ]; then \
		echo "# BED name for the sample '$(@D)/`echo $$(basename $(@D))`.bed_name' exists" ; \
		if  [ "$(@D)/`echo $$(basename $(@D))`.bed_name" != "$@" ]; then \
			echo "# BED name for the sample '$(@D)/`echo $$(basename $(@D))`.bed_name' exists and is different than the wanted BED name file '$@'" ; \
			ln -s $(@D)/`echo $$(basename $(@D))`.bed_name $@ || cp $(@D)/`echo $$(basename $(@D))`.bed_name $@; \
		fi; \
	elif [ -e $(@D)/../`echo $$(basename $(@D))`.bed_name ]; then \
		echo "# BED name for the sample '$(@D)/../`echo $$(basename $(@D))`.bed_name' exists" ; \
		if  [ "$(@D)/../`echo $$(basename $(@D))`.bed_name" != "$@" ]; then \
			echo "# BED name for the sample '$(@D)/../`echo $$(basename $(@D))`.bed_name' exists and is different than the wanted BED name file '$@'" ; \
			ln -s $(@D)/../`echo $$(basename $(@D))`.bed_name $@ || cp $(@D)/../`echo $$(basename $(@D))`.bed_name $@; \
		fi; \
	elif [ -s $< ]; then \
		echo "# BED name for the sample from the manifest '$<'" ; \
		echo -e "-\tfrom Manifest" > $@; \
	elif [ -e "$(BED)" ] && [ "$(BED)" != "" ]; then \
		echo "# Input BED '$(BED)' exists" ; \
		echo -e "$(BED)\tfrom APP" > $@; \
	else \
		touch $@; \
	fi;
	#
	if [ ! -e $@ ]; then touch $@; fi;
	# Clean


# BED file from a Manifest
%.transcripts: %.manifest
	# Create BED file
	# Merge bacause of PCR manifest defined only Amplicons
	-if [ -e "$(TRANSCRIPTS)" ] && [ "$(TRANSCRIPTS)" != "" ]; then \
		echo "# Input TRANSCRIPTS '$(TRANSCRIPTS)' exists" ; \
		cp $(TRANSCRIPTS) $@; \
	elif [ -e $(@D)/`echo $$(basename $(@D))`.transcripts ]; then \
		echo "# TRANSCRIPTS for the sample '$(@D)/`echo $$(basename $(@D))`.transcripts' exists" ; \
		if  [ "$(@D)/`echo $$(basename $(@D))`.transcripts" != "$@" ]; then \
			echo "# TRANSCRIPTS for the sample '$(@D)/`echo $$(basename $(@D))`.transcripts' exists and is different than the wanted TRANSCRIPTS file '$@'" ; \
			ln -s $(@D)/`echo $$(basename $(@D))`.transcripts $@ || cp $(@D)/`echo $$(basename $(@D))`.transcripts $@; \
		fi; \
	elif [ -e $(@D)/../`echo $$(basename $$(dirname $(@D)))`.transcripts ]; then \
		echo "# TRANSCRIPTS for the sample '$(@D)/../`echo $$(basename $$(dirname $(@D)))`.transcripts' exists" ; \
		if  [ "$(@D)/../`echo $$(basename $$(dirname $(@D)))`.transcripts" != "$@" ]; then \
			echo "# TRANSCRIPTS for the sample '$(@D)/../`echo $$(basename $$(dirname $(@D)))`.transcripts' exists and is different than the wanted TRANSCRIPTS file '$@'" ; \
			ln -s $(@D)/../`echo $$(basename $$(dirname $(@D)))`.transcripts $@ || cp $(@D)/../`echo $$(basename $$(dirname $(@D)))`.transcripts $@; \
		fi; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;
	# Clean


# BED clipped region file from a Manifest
%.region_clipped.bed: %.manifest %.bed
	# Create BED file
	# Merge bacause of PCR manifest defined only Amplicons
	-if [ -e "$(REGIONCLIPPEDBED)" ]; then \
		echo "# Input Region Clipped BED '$(REGIONCLIPPEDBED)' exists" ; \
		cp $(REGIONCLIPPEDBED) $@; \
	elif [ -e $(@D)/`echo $$(basename $(@D))`.region_clipped.bed ]; then \
		echo "# Region Clipped BED for the sample '$(@D)/`echo $$(basename $(@D))`.region_clipped.bed' exists" ; \
		if  [ "$(@D)/`echo $$(basename $(@D))`.region_clipped.bed" != "$@" ]; then \
			echo "# Region Clipped BED for the sample '$(@D)/`echo $$(basename $(@D))`.region_clipped.bed' exists and is different than the wanted BED file '$@'" ; \
			cp $(@D)/`echo $$(basename $(@D))`.region_clipped.bed $@; \
		fi; \
	elif [ -s $< ]; then \
		echo "# Region Clipped BED for the sample generated from the manifest '$<'" ; \
		rm -f $@.tmp $@.sorted.tmp; \
		$(CAP_ManifestToBED) --input=$< --output=$@.tmp --output_type=region_clipped; \
		$(BEDTOOLS) sort -i $@.tmp | $(BEDTOOLS) merge -c 4 -o collapse  | awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t+\t"$$4}' > $@; \
		rm -f $@.tmp; \
	else \
		echo "# Error creating region BED"; \
	fi;
	if [ "`grep ^ -c $@`" == "0" ]; then \
		echo "# Region Clipped BED generated is empty. BED will be used as Region Clipped BED" ; \
		cp $*.bed $@; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;




# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# MAIN RULES '$(MK_RELEASE)' : basicaly to manage VCF, BAM, FASTQ... using SAMTOOLS, FASTQC GATK, PICARD, CAP, BWA, TABIX, IGV, HOWARD..."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "POST_ALIGNMENT:sorting:BAM sorting"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "POST_ALIGNMENT:compress:BAM compression:BAM_COMPRESS='$(BAM_COMPRESSION)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


PIPELINES_COMMENT := "POST_CALLING:normalization:VCF Normalization."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "POST_ANNOTATION:normalization:VCF Normalization."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "POST_CALLING:sorting:VCF Sorting."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "POST_ANNOTATION:sorting:VCF Sorting."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
