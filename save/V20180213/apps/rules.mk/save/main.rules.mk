############################
# Main Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.4.7b"
MK_DATE="07/05/2018"

# Release note
# 10/03/2015-0.9.4beta: change genome reference location, in the file %.genome
# 14/06/2015-0.9.4.1b: change interval from bed (grep to remove '@' lines)
# 25/11/2015-0.9.4.2b: Change VCF sorting, using vcftools instead of igvtools
# 25/01/2016-0.9.4.3b: Change VCF sorting vcftools parameter, adding -C chromosome ordering
# 18/02/2016-0.9.4.4b: change cp to cat form manifest copy
# 13/04/2016-0.9.4.5b: removing full.vcf rule
# 29/09/2016-0.9.4.6b: Update PICARD to releasr picard.jar
# 07/05/2018-0.9.4.7b: Add --force for translation vcf to txt

# TOOLS
IGVTOOLS?=$(NGSbin)/igvtools.jar
SAMTOOLS?=$(NGSbin)/samtools
JAVA?=java
TABIX?=$(NGSbin)/tabix
BWA?=$(NGSbin)/bwa
PICARDLIB?=$(NGSbin)/picard-tools
FASTQC?=$(NGSbin)/fastqc
VCFTOOLS?=$(NGSbin)/vcftools
BGZIP?=bgzip
GZIP?=gzip

# HOWARD Prioritization
HOWARD_FILTER?="default"
HOWARD_CONFIG_FILTER?="config.filter.ini"

# HOWARD Translation
HOWARD_ANNOTATIONS?="PZScore,PZFlag,PZComment,Symbol,hgvs,location,outcome,AlleleFrequency,AD,dbSNP,dbSNPNonFlagged,popfreq,database_CPSGEN,DP,AF,VF,GQ,Ensembl,TI,FC,GWASCatalog,COSMIC,LocalDB,1000genomesALL,1000genomesEUR,6500NHLBIALL,6500NHLBIEUR,PolyPhen2HumanVarPred,PolyPhen2HumanDivPred,MutationTasterPred,MutationAssessorPred,LTRPred,IARCTP53,SIFT,phastCons,PhyloP,SiPhy,FATHMM,LRT,GERP,PolyPhen2HumanVar,PolyPhen2HumanDiv,MutationTaster,MutationAssessor,TFBS,FilterComment,ALL"
#HOWARD_FIELDS?="PZScore,PZFlag,PZComment,Symbol,hgvs,location,outcome,AlleleFrequency,AD,dbSNP,dbSNPNonFlagged,popfreq,ALL"
HOWARD_FIELDS?="NOMEN,PZFlag,PZScore,PZComment,CNOMEN,PNOMEN,location,outcome,VAF_average,dbSNP,dbSNPNonFlagged,popfreq,ALL"
HOWARD_SORT_BY?="PZFlag,PZScore"
HOWARD_ORDER_BY?="DESC,DESC"

# OPTIONS
REMOVE_INTERMEDIATE_SAM?=1
BED?=
PRIMER_BED?=
THREADS_SAMTOOLS?=$(THREADS_BY_SAMPLE)

BAM_COMPRESSION?=5

## VCF files
##############

# VCF Compression
%.vcf.gz: %.vcf %.empty.vcf
	# If no VCF or empty file, create an empty VCF
	if [ ! -s $< ]; then cp $*.empty.vcf $<; fi;
	# VCF Sorting
	#$(VCFTOOLS)/vcf-sort $^ > $<.sorted
	$(VCFTOOLS)/vcf-sort $< > $<.sorted
	# VCF compression with BGZIP
	$(BGZIP) -f $<.sorted -c > $@
	rm -f $<.sorted
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


# SORT VCF
%.vcf: %.unsorted.vcf
	mkdir -p $@"_BCFTOOLS_sort/"
	$(BCFTOOLS) sort -T $@"_BCFTOOLS_sort/" $< > $@
	rm -rf $@"_BCFTOOLS_sort/"


# NORMALIZE VCF with BCFTOOLS
%.norm.vcf: %.vcf %.genome 
	$(BCFTOOLS) norm -m- -f `cat $*.genome` $< > $@


# MERGE SNP and InDel VCF
%.unsorted.vcf: %.SNP.vcf %.InDel.vcf
	# CONCAT
	$(VCFTOOLS)/vcf-concat $^ > $@
	# Remove IDX
	#-rm $*.SNP.vcf $*.InDel.vcf
	-rm -f $*.SNP.vcf.idx $*.InDel.vcf.idx
	-rm -f $*.idx
	# SORTING
	#$(JAVA) -Xmx1g -jar $(IGVTOOLS) sort $@.unsorted.unsorted.vcf $@.unsorted.vcf
	# SORTING
	#$(JAVA) -Xmx1g -jar $(IGVTOOLS) sort $@.unsorted.vcf $@
	# CLEANING
	#-rm $@.unsorted.vcf


# VCF to tab delimiter
%.tsv: %.vcf
	# Prioritization, calculation and translation step
	$(HOWARD) --config=$(HOWARD_CONFIG) --config_filter=$(HOWARD_CONFIG_FILTER) --filter=$(HOWARD_FILTER) --pzfields="PZScore,PZFlag,PZComment,PZInfos" --format=tab  --fields="$(HOWARD_FIELDS)" --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)"  --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --input=$< --output=$@ --force;
	# Touch
	if [ ! -e $@ ]; then touch $@; fi;
	# Cleaning

# Hard filtering
%.hard.tsv: %.vcf
	# Prioritization, calculation and translation step and hard filtering
	$(HOWARD) --config=$(HOWARD_CONFIG) --config_filter=$(HOWARD_CONFIG_FILTER) --filter=$(HOWARD_FILTER) --pzfields="PZScore,PZFlag,PZComment,PZInfos" --format=tab  --fields="$(HOWARD_FIELDS)" --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)"  --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --input=$< --output=$@ --hard --force;
	# Touch
	if [ ! -e $@ ]; then touch $@; fi;
	# Cleaning
	



# VCF to tab delimiter
%.txt: %.vcf
	# Prioritization step
	#$(HOWARD_PRIORITIZATION) --config=$(HOWARD_CONFIG) --config_filter=$(HOWARD_CONFIG_FILTER) --filter=$(HOWARD_FILTER) --input_file=$< --output_file=$*.prioritized.vcf;
	# Translation step
	#$(HOWARD_TRANSLATION) --input_file=$*.prioritized.vcf --output_file=$@ --format=tab --annotation="$(HOWARD_ANNOTATIONS)" --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)"  ; # --noremove_filtered --output_file=$@ --output_format=tab
	#$(HOWARD) --config=$(HOWARD_CONFIG) --config_filter=$(HOWARD_CONFIG_FILTER) --filter=$(HOWARD_FILTER) --pzfields="PZScore,PZFlag,PZComment,PZInfos" --format=tab  --fields="$(HOWARD_ANNOTATIONS)" --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)"  --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --input=$< --output=$@;
	$(HOWARD) --config=$(HOWARD_CONFIG) --config_filter=$(HOWARD_CONFIG_FILTER) --filter=$(HOWARD_FILTER) --pzfields="PZScore,PZFlag,PZComment,PZInfos" --format=tab  --fields="$(HOWARD_FIELDS)" --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)"  --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --input=$< --output=$@ --force;
	# Touch
	if [ ! -e $@ ]; then touch $@; fi;
	# Cleaning
	-rm -f $*.prioritized.vcf



# Hard filtering
%.hard.txt: %.vcf
	# Prioritization and hard filtering
	#$(HOWARD) --input=$< --output=$@ --format=tab  --config=$(HOWARD_CONFIG) --config_filter=$(HOWARD_CONFIG_FILTER) --filter=$(HOWARD_FILTER)  --format=tab --annotation="$(HOWARD_ANNOTATIONS)" --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)"
	#$(HOWARD) --input=$< --output=$@ --format=tab  --config=$(HOWARD_CONFIG) --config_filter=$(HOWARD_CONFIG_FILTER) --filter=$(HOWARD_FILTER)  --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)" --fields="$(HOWARD_FIELDS)" --hard
	$(HOWARD) --config=$(HOWARD_CONFIG) --config_filter=$(HOWARD_CONFIG_FILTER) --filter=$(HOWARD_FILTER) --pzfields="PZScore,PZFlag,PZComment,PZInfos" --format=tab  --fields="$(HOWARD_FIELDS)" --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)"  --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --input=$< --output=$@ --hard --force;
	# Prioritization step
	#$(HOWARD_PRIORITIZATION) --config=$(HOWARD_CONFIG) --config_filter=$(HOWARD_CONFIG_FILTER) --filter=$(HOWARD_FILTER) --input_file=$< --output_file=$*.prioritized.hard.vcf --hard;
	# Translation step
	#$(HOWARD_TRANSLATION) --input_file=$*.prioritized.hard.vcf --output_file=$@ --format=tab --annotation="$(HOWARD_ANNOTATIONS)" --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)"  ; # --noremove_filtered --output_file=$@ --output_format=tab
	# Touch
	if [ ! -e $@ ]; then touch $@; fi;
	# Cleaning
	-rm -f $*.prioritized.hard.vcf




%.prioritized.vcf: %.vcf
	# Prioritization step
	#$(HOWARD_PRIORITIZATION) --config=$(HOWARD_CONFIG) --config_filter=$(HOWARD_CONFIG_FILTER) --filter=$(HOWARD_FILTER) --input_file=$< --output_file=$*.prioritized.hard.vcf --hard;
	$(HOWARD) --input=$< --output=$@ --config=$(HOWARD_CONFIG) --config_filter=$(HOWARD_CONFIG_FILTER) --filter=$(HOWARD_FILTER)  --format=tab --annotation="$(HOWARD_ANNOTATIONS)" --sort_by="$(HOWARD_SORT_BY)" --order_by="$(HOWARD_ORDER_BY)"



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
#%.bam: %.cram %.genome
#	$(SAMTOOLS) view -o $@ -O BAM -S -T `cat $*.genome` $< -@ $(THREADS_SAMTOOLS);
#	# test Empty output file
#	# Remove intermediate SAM file

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


# BAM Sorting
#%.bam : %.unsorted.bam
#	if ((! $$($(SAMTOOLS) view -H $< | grep "^@HD.*VN:.*SO:" | wc -l))); then \
#		$(SAMTOOLS) sort $< -o $@ -@ $(THREADS_SAMTOOLS); \
#		# CHECK NUMBER of READS in BAM \
#		#if [ "$$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS $<))" != "$$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $@)" ]; then echo "# ERROR in Number of reads between $< and $@ !!!"; exit 1; else echo "# Number of reads OK between $< and $@"; fi \
		# CHECK NUMBER of READS in BAM \
		#echo $$($(SAMTOOLS) view -c -@ $(THREADS_SAMTOOLS) -F 0x0100 $<)" READS for $<"; \
		#echo $$($(SAMTOOLS) view -c -@ $(THREADS_SAMTOOLS) -F 0x0100 $@)" READS for $@"; \
		#if [ "$$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $<)" != "$$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $@)" ]; then \
		#	echo "# ERROR in Number of reads between $< and $@ !!!"; \
		#	echo "# BCFTOOLS STATS for $<";  \
		#	$(SAMTOOLS) index $<; \
		#	$(SAMTOOLS) stats $< | grep SN; \
		#	$(SAMTOOLS) idxstats $<; \
		#	echo "# BCFTOOLS STATS for $@";  \
		#	$(SAMTOOLS) index $@; \
		#	$(SAMTOOLS) stats $@ | grep SN; \
		#	$(SAMTOOLS) idxstats $@; \
		#	exit 1; \
		#else \
		#	echo "# Number of reads OK between $< and $@"; \
		#fi; \
#	else \
#		mv $< $@; \
#	fi;
	
	#if [ "$$($(SAMTOOLS) view -c $<)" != "$$($(SAMTOOLS) view -c $@)"]; then echo "# ERROR in number of reads between $< and $@ !!!"; exit 0; fi
	# Remove unsorted file (if still exists...)
#	rm -f $<



# BAM Sorting
#%.bam : %.unsorted.sam
#	$(SAMTOOLS) sort $< $* -O bam -@ $(THREADS_SAMTOOLS)
#	# Remove unsorted file
#	rm $<

# SAM from FASTQ
%.sai : %.fastq.gz %.genome
	$(BWA) aln -f $@ `cat $*.genome` $<

# FASTQ Compression with GZIP
%.fastq.gz: %.fastq
	$(GZIP) $< -c > $@

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

#%.manifests_list.txt:
#	-mkdir -p $(@D)
#	-if [ -s $(INPUTDIR)/`echo $$(basename $$(dirname $(@D)))`/manifests_list.txt ]; then \
#		cp $(INPUTDIR)/`echo $$(basename $$(dirname $(@D)))`/manifests_list.txt $@; \
#	else \
#		touch $@; \
#	fi;
#	if [ ! -e $@ ]; then touch $@; fi;

%.manifest: %.manifest_from_samplesheet
	-mkdir -p $(@D)
	# Create MANIFEST file
	# 1. test if main manifest file (SAMPLE.manifest) exists (use for Sample analysis)
	# 2. test if a manifest is found on SampleSheet (use for Run analysis)
	# 3. test if a manifest is defined in the APP (as an environment variable)
	
	-if [ -e $(@D)/`echo $$(basename $(@D))`.manifest ]; then \
		echo "# MANIFEST for the sample '$(@D)/`echo $$(basename $(@D))`.manifest' exists" ; \
		if  [ "$(@D)/`echo $$(basename $(@D))`.manifest" != "$@" ]; then \
			echo "# MANIFEST for the sample '$(@D)/`echo $$(basename $(@D))`.manifest' exists and is different than the wanted MANIFEST file '$@'" ; \
			ln -s $(@D)/`echo $$(basename $(@D))`.manifest $@ || cp $(@D)/`echo $$(basename $(@D))`.manifest $@; \
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
	
	-if [ -e $(@D)/`echo $$(basename $(@D))`.manifest_name ]; then \
		echo "# MANIFEST name for the sample '$(@D)/`echo $$(basename $(@D))`.manifest_name' exists" ; \
		if  [ "$(@D)/`echo $$(basename $(@D))`.manifest_name" != "$@" ]; then \
			echo "# MANIFEST name for the sample '$(@D)/`echo $$(basename $(@D))`.manifest_name' exists and is different than the wanted MANIFEST name file '$@'" ; \
			ln -s $(@D)/`echo $$(basename $(@D))`.manifest_name $@ || cp $(@D)/`echo $$(basename $(@D))`.manifest_name $@; \
		fi; \
	elif [ -e $(@D)/../`echo $$(basename $(@D))`.manifest_name ]; then \
		echo "# MANIFEST name for the sample '$(@D)/../`echo $$(basename $(@D))`.manifest_name' exists" ; \
		if  [ "$(@D)/../`echo $$(basename $(@D))`.manifest_name" != "$@" ]; then \
			echo "# MANIFEST name for the sample '$(@D)/../`echo $$(basename $(@D))`.manifest_name' exists and is different than the wanted MANIFEST name file '$@'" ; \
			ln -s $(@D)/../`echo $$(basename $(@D))`.manifest_name $@ || cp $(@D)/../`echo $$(basename $(@D))`.manifest_name $@; \
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
	
	# Column number of sample manifest
	-grep -i ^Sample_ID $< | tr -d '\r\n' | sed -e 's/,/\n/g' | grep -i -n Manifest | cut -d \: -f 1 > $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I;
	-if [ ! -s $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I ]; then echo "1" > $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I; fi;
	# Manifest letter for the index
	-grep -e ^$$(basename $(@D)), $< | tail -n 1 | tr -d '\r\n' | cut -d \, -f `cat $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I` > $(@D)/$(*F).SAMPLE_MANIFEST_I;
	
	# Manifest name for the sample
	-grep ^`cat $(@D)/$(*F).SAMPLE_MANIFEST_I`, $*.manifests_list.txt | cut -d \, -f 2 > $(@D)/$(*F).SAMPLE_MANIFEST;
	
	# BED
	file=$$( echo "$$( dirname $(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST` )/$$( basename $(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST` ).bed" ); \
	echo "bed file is : $$file "; \
	if [ -s $$file ]; then \
		cat $$file > $*.bed; \
	else \
		echo "$$file don't exist !!!"; \
	fi;
	
	# BED GENES
	file=$$( echo "$$( dirname $(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST` )/$$( basename $(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST` ).genes" ); \
	echo "genes file is : $$file "; \
	if [ -s $$file ]; then \
		cat $$file > $*.genes; \
	else \
		echo "$$file don't exist !!!"; \
	fi;
	
	# TRANSCRIPTS
	file=$$( echo "$$( dirname $(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST` )/$$( basename $(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST` ).transcripts" ); \
	echo "transcripts file is : $$file "; \
	if [ -s $$file ]; then \
		cat $$file > $*.transcripts; \
	else \
		echo "$$file don't exist !!!"; \
	fi;
	
	# Found manifest, or Default manifest is the first in the list
	if [ "`cat $(@D)/$(*F).SAMPLE_MANIFEST`" != "" ] && [ -e "$(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST`" ]; then \
		echo "# TEST Found Manifest '"`cat $(@D)/$(*F).SAMPLE_MANIFEST`"' in '$(MANIFEST_FOLDER)' for sample `echo $$(basename $$(dirname $(@D)))`/$(*F)"; \
		cat "$(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST`" > $@; \
	elif [ "`cut -d, -f2 $*.manifests_list.txt`" != "" ] && [ -e "$(MANIFEST_FOLDER)/`cut -d, -f2 $*.manifests_list.txt`" ]; then \
		echo "# TEST Manifest NOT found in '$(MANIFEST_FOLDER)'. Default Manifest used '"`cut -d, -f2 $*.manifests_list.txt`"' for sample `echo $$(basename $$(dirname $(@D)))`/$(*F)"; \
		cat "$(MANIFEST_FOLDER)/`cut -d, -f2 $*.manifests_list.txt`" > $@; \
	fi;
	
	# IF root manifest exists
	if [ -s $(@D)/`echo $$(basename $(@D))`.manifest ] && [ "$(@D)/`echo $$(basename $(@D))`.manifest" != "$@" ]; then \
		cp -p $(@D)/`echo $$(basename $(@D))`.manifest $@; \
	fi;
	
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
	-if [ -s $(INPUTDIR)/`echo $$(basename $$(dirname $(@D)))`/manifests_list_name.txt ]; then \
		cp -p $(INPUTDIR)/`echo $$(basename $$(dirname $(@D)))`/manifests_list_name.txt $*.manifests_list_name.txt; \
	else \
		touch $*.manifests_list_name.txt; \
	fi;
	if [ ! -e $@ ]; then touch $*.manifests_list_name.txt; fi;
	
	# Column number of sample manifest
	-grep -i ^Sample_ID $< | tr -d '\r\n' | sed -e 's/,/\n/g' | grep -i -n Manifest | cut -d \: -f 1 > $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I_name;
	-if [ ! -s $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I_name ]; then echo "1" > $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I_name; fi;
	# Manifest letter for the index
	-grep -e ^$$(basename $(@D)), $< | tail -n 1 | tr -d '\r\n' | cut -d \, -f `cat $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I_name` > $(@D)/$(*F).SAMPLE_MANIFEST_I_name;
	
	# Manifest name for the sample
	-grep ^`cat $(@D)/$(*F).SAMPLE_MANIFEST_I_name`, $*.manifests_list_name.txt | cut -d \, -f 2 > $(@D)/$(*F).SAMPLE_MANIFEST_name;
	
	-if [ ! -s $(@D)/$(*F).SAMPLE_MANIFEST_name ]; then \
		> $@; \
	else \
		echo -e $$(cat $(@D)/$(*F).SAMPLE_MANIFEST_name)"\tfrom SampleSheet" > $@; \
	fi;
	
	# Empty manifest if failed!
	if [ ! -e $@ ]; then touch $@; fi;
	
	# Clean
	-rm -f $*.manifests_list_name.txt
	# Remove intermediate files
	-rm -f $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I_name $(@D)/$(*F).SAMPLE_MANIFEST_I_name $(@D)/$(*F).SAMPLE_MANIFEST_name;



# Interval from bed from manifest?
%.from_manifest.intervals: %.bed %.bam %.bam.bai %.dict
	# manifest to interval (not needed?)
	#cat $<  | tr -d '\r' | sed -e "s/^M//" | awk -F"\t" '{print $$1":"$$2"-"$$3}' > $@
	# try to extract from the bam if exists, in order to not call in the whome genome
	-+if [ ! -s $< ] && [ -s $*.bam ]; then \
		echo "# Generation of BED from BAM because bed/manifest is empty"; \
		#$(BEDTOOLS)/bamToBed -i $*.bam  | $(BEDTOOLS)/bedtools sort -i - | $(BEDTOOLS)/mergeBed -n -i - > $*.bed; \
		#$(BEDTOOLS)/genomeCoverageBed -ibam $*.bam -bg | $(BEDTOOLS)/mergeBed -n -i - > $*.bed; \
		if (($$($(SAMTOOLS) idxstats $*.bam | awk '{SUM+=$$3+$$4} END {print SUM}'))); then \
			rm -f $*.bam.genomeCoverageBed.mk $*.bam.genomeCoverageBed1.mk $*.bam.genomeCoverageBed2.mk $*.bam.genomeCoverageBed3.mk; \
			for chr in $$($(SAMTOOLS) idxstats $*.bam | grep -v "\*" | awk '{ if ($$3+$$4>0) print $$1 }'); do \
				echo "$*.bam.genomeCoverageBed.$$chr.bed: $*.bam" >> $*.bam.genomeCoverageBed1.mk; \
				echo "	$(SAMTOOLS) view $*.bam -b $$chr | $(BEDTOOLS)/genomeCoverageBed -ibam stdin -bg | $(BEDTOOLS)/mergeBed -n -i - > $*.bam.genomeCoverageBed.$$chr.bed " >> $*.bam.genomeCoverageBed1.mk; \
				echo -n " $*.bam.genomeCoverageBed.$$chr.bed" >> $*.bam.genomeCoverageBed2.mk; \
			done; \
			echo -n "$*.bed: " | cat - $*.bam.genomeCoverageBed2.mk > $*.bam.genomeCoverageBed3.mk; \
			echo ""  >> $*.bam.genomeCoverageBed3.mk; \
			#echo "	cat $$^ > $$@ " >> $*.bam.genomeCoverageBed3.mk; \
			echo "	cat $$^ > $*.bed " >> $*.bam.genomeCoverageBed3.mk; \
			echo "	-rm -f $$^ " >> $*.bam.genomeCoverageBed3.mk; \
			cat $*.bam.genomeCoverageBed1.mk $*.bam.genomeCoverageBed3.mk >> $*.bam.genomeCoverageBed.mk; \
			echo "TESTgenomeCoverageBed: "; \
			cat $*.bam.genomeCoverageBed.mk; \
			make -i -f $*.bam.genomeCoverageBed.mk $*.bed -j $(THREADS)  1>/dev/null 2>/dev/null; \
			rm $*.bam.genomeCoverageBed*.mk; \
		fi; \
	fi;
	# INTERVAL WITH PICARD
	if [ -s $< ]; then \
		$(JAVA) -jar $(PICARD) BedToIntervalList I=$< O=$@ SD=$$(cat $*.dict) ; \
	fi;
	# If error, try intervals with GREP/SED/AWK
	if [ ! -s $@ ] && [ -s $< ]; then \
		grep -v ^@ $<  | tr -d '\r' | sed -e "s/^M//" | awk -F"\t" '{print $$1":"$$2"-"$$3}' | sed s/:0-/:1-/gi > $@; \
	fi;
	# touch
	if [ ! -e $@ ]; then touch $@; fi;

# Interval from BED
%.bed.intervals: %.bed %.dict
	# BED to Intervals (not needed?)
	#cat $< | tr -d '\r' | sed -e "s/^M//" | awk -F"\t" '{print $$1":"$$2"-"$$3}' > $@
	# INTERVAL WITH PICARD
	if [ -s $< ]; then \
		$(JAVA) -jar $(PICARD) BedToIntervalList I=$< O=$@ SD=$$(cat $*.dict) ; \
	fi;
	# If error, try intervals with GREP/SED/AWK
	if [ ! -s $@ ] && [ -s $< ]; then \
		grep -v ^@ $<  | tr -d '\r' | sed -e "s/^M//" | awk -F"\t" '{print $$1":"$$2"-"$$3}' | sed s/:0-/:1-/gi > $@; \
	fi;
	#grep -v ^@ $< | tr -d '\r' | sed -e "s/^M//" | awk -F"\t" '{print $$1":"$$2"-"$$3}' > $@
	if [ ! -e $@ ]; then touch $@; fi;

# BED primer file from a Manifest
%.primers.bed: %.manifest
	# Create BED file
	-if [ -e "$(PRIMER_BED)" ]; then \
		cp $(PRIMER_BED) $@; \
	elif [ -s $< ]; then \
		$(FATBAM_ManifestToBED) --input=$< --output=$@ --output_type=primer; \
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
	
	
	-if [ -e $(@D)/`echo $$(basename $(@D))`.bed ]; then \
		echo "# BED for the sample '$(@D)/`echo $$(basename $(@D))`.bed' exists" ; \
		if  [ "$(@D)/`echo $$(basename $(@D))`.bed" != "$@" ]; then \
			echo "# BED for the sample '$(@D)/`echo $$(basename $(@D))`.bed' exists and is different than the wanted BED file '$@'" ; \
			ln -s $(@D)/`echo $$(basename $(@D))`.bed $@ || cp $(@D)/`echo $$(basename $(@D))`.bed $@; \
		fi; \
	elif [ -s $< ]; then \
		echo "# BED for the sample generated from the manifest '$<'" ; \
		rm -f $@.tmp $@.sorted.tmp; \
		$(FATBAM_ManifestToBED) --input=$< --output=$@.tmp --output_type=region; \
		$(BEDTOOLS)/sortBed -i $@.tmp > $@.sorted.tmp; \
		$(BEDTOOLS)/mergeBed -i $@.sorted.tmp -nms | awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t+\t"$$4}' > $@; \
		rm -f $@.tmp $@.sorted.tmp; \
	elif [ -e "$(BED)" ] && [ "$(BED)" != "" ]; then \
		echo "# Input BED '$(BED)' exists" ; \
		cp $(BED) $@; \
	else \
		touch $@; \
	fi;
	
	if [ ! -e $@ ]; then touch $@; fi;
	# Clean
	#-rm -f $*.manifest
	

# BED file from a Manifest
%.bed_name: %.manifest_name
	# Create BED file
	# 1. test if main bed file (SAMPLE.bed) exists (use for Sample analysis)
	# 2. test if a bed can be generated from the manifest
	# 3. test if a manifest is defined in the APP (as an environment variable)
	
	
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
	
	if [ ! -e $@ ]; then touch $@; fi;
	# Clean
	#-rm -f $*.manifest


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
	elif [ -e $(@D)/../`echo $$(basename $(@D))`.transcripts ]; then \
		echo "# TRANSCRIPTS for the sample '$(@D)/`echo $$(basename $(@D))`.transcripts' exists" ; \
		if  [ "$(@D)/../`echo $$(basename $(@D))`.transcripts" != "$@" ]; then \
			echo "# TRANSCRIPTS for the sample '$(@D)/../`echo $$(basename $(@D))`.transcripts' exists and is different than the wanted TRANSCRIPTS file '$@'" ; \
			ln -s $(@D)/../`echo $$(basename $(@D))`.transcripts $@ || cp $(@D)/../`echo $$(basename $(@D))`.transcripts $@; \
		fi; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;
	# Clean
	#-rm -f $*.manifest


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
		$(FATBAM_ManifestToBED) --input=$< --output=$@.tmp --output_type=region_clipped; \
		$(BEDTOOLS)/sortBed -i $@.tmp > $@.sorted.tmp; \
		$(BEDTOOLS)/mergeBed -i $@.sorted.tmp -nms  | awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t+\t"$$4}' > $@; \
		rm -f $@.tmp $@.sorted.tmp; \
	else \
		echo "# Error creating region BED"; \
	fi;
	if [ "`grep ^ -c $@`" == "0" ]; then \
		echo "# Region Clipped BED generated is empty. BED will be used as Region Clipped BED" ; \
		cp $*.bed $@; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;

	# Error
	#elif [ -s $*.bed ]; then
	#	cp $*.bed $@;


# REPORT
#%.full.vcf: %.bam %.bam.bai %.vcfs.list
#	# Find associated VCF
#	-for VCF in `cat $*.vcfs.list`; do \
#		if [ ! -e $VCF.gz ]; then \
#			$(BGZIP) -f -c $$VCF > $$VCF.gz; \
#		fi;
#		if [ ! -e $$VCF.gz.tbi ]; then
#			echo "# VCF tabix";
#			$(TABIX) -p -f vcf $$VCF.gz
#			$(TABIX) -p -f vcf $$VCF
#		#fi;
#		VCFS_GZ_OPTION=$$VCFS_GZ_OPTION" "$$VCF.gz
#		echo "# VCF end";
#	done;
#
#	#$NGS_FOLDER/bin/vcftools/vcf-isec -f --nfiles +1 $VCFS_GZ_OPTION > $OUTPUT;
#	$(VCFTOOLS)/vcf-merge $$VCFS_GZ_OPTION > $@;

#%.gatkUG.reduced.vcf: %.reduced.bam %.reduced.bam.bai %.reduced.from_manifest.intervals
#	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKUG_FLAGS) \
#		-T UnifiedGenotyper \
#		-R $(REF) \
#		-L $*.reduced.from_manifest.intervals \
#		-I $< \
#		-o $@
#	if [ ! -e $@ ]; then touch $@; fi;
#



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# MAIN RULES '$(MK_RELEASE)' : basicaly to manage VCF, BAM, FASTQ... using SAMTOOLS, FASTQC GATK, PICARD, FATBAM, BWA, TABIX, IGV..."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "POST_ALIGNMENT:sorting:BAM sorting"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "POST_ALIGNMENT:compress:BAM compression:BAM_COMPRESS='$(BAM_COMPRESSION)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


