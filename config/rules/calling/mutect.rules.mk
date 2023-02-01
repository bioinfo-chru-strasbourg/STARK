############################
# MUTECT Calling Rules
# Release: 0.9.1.1
# Date: 11/06/2021
# Author: Antony Le Bechec
############################




### FilterMutectCalls
# base_qual
# clustered_events
# contamination
# fragment
# germline
# haplotype
# map_qual
# multiallelic
# orientation
# panel_of_normals
# position
# slippage
# strand_bias
# weak_evidence


##########
# MUTECT #
##########

# Parameters
#MUTECT_INPUT_TYPE=tumor
MUTECT_INTERVAL_PADDING?=0

DPMIN_MUTECT?=30


%.MuTect$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	# Calling
	$(JAVA7) -jar $(MUTECT) \
		--analysis_type MuTect \
		--reference_sequence $$(cat $*.genome) \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "--intervals $*.design.bed.interval_list"; fi;) \
		--input_file:tumor $< \
		--vcf $@.tmp1 \
		-ip $(MUTECT_INTERVAL_PADDING);
	# Normalize
	grep "^##" $@.tmp1 > $@.tmp
	grep "^##" -v $@.tmp1 | grep "REJECT" -v | cut -f1-10 >> $@.tmp
	$(BCFTOOLS) view  -i 'FORMAT/DP>=$(DPMIN_MUTECT)' $@.tmp > $@.tmp2;
	# Empty
	if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp2; fi;
	# Copy
	if [ ! -e $@ ]; then cp $@.tmp2 $@; fi;
	# Clean
	rm -f $@.tmp* $@.idx





#################
# MUTECT2 GATK4 #
#################

#THREADS_GATK4?=$(THREADS_BY_CALLER)
THREADS_GATK4_MUTECT2?=$(THREADS_BY_CALLER)
MINPRUNING_GATK4_MUTECT2?=20
MAXREADS_GATK4_MUTECT2?=1000

DPMIN_MUTECT2?=30


GATK4_MUTECT2_FLAGS_SHARED?=--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter --max-reads-per-alignment-start $(MAXREADS_GATK4_MUTECT2) --dont-use-soft-clipped-bases true --min-pruning $(MINPRUNING_GATK4_MUTECT2) --callable-depth $(DPMIN_MUTECT2) --verbosity ERROR --native-pair-hmm-threads $(THREADS_GATK4_MUTECT2)

#--max-reads-per-alignment-start
#--dont-use-soft-clipped-bases

%.MuTect2$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.dict %.design.bed.interval_list #%.from_manifest.interval_list
	$(JAVA11) $(JAVA_FLAGS) -jar $(GATK4) Mutect2 $(GATK4_MUTECT2_FLAGS_SHARED) \
		-R $$(cat $*.genome) \
		-I $< \
		-tumor $$(basename $< | cut -d"." -f1) \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-O $@.tmp1;
	# Normalize
	grep "^##" $@.tmp1 | sed s/ID=TLOD,Number=A/ID=TLOD,Number=./gi > $@.tmp.vcf
	grep "^##" -v $@.tmp1 | cut -f1-10 >> $@.tmp.vcf
	# sort
	$(JAVA11) -jar $(PICARD) SortVcf -I $@.tmp.vcf -O $@.tmp2.vcf -SD $$(cat $*.dict);
	# DPMIN_MUTECT2
	#$(VCFUTILS) varFilter -d $(DPMIN_MUTECT2) $@.tmp2.vcf > $@.tmp3.vcf;
	$(BCFTOOLS) view  -i 'FORMAT/DP>=$(DPMIN_MUTECT2)' $@.tmp2.vcf > $@.tmp3.vcf;
	# Empty
	if [ ! -e $@.tmp3.vcf ]; then cp $*.empty.vcf $@.tmp3.vcf; fi;
	# Copy
	if [ ! -e $@ ]; then cp $@.tmp3.vcf $@; fi;
	# Clean
	rm -f $@.tmp* $@.idx





#########################
# MUTECT2 flagged GATK4 #
#########################

#THREADS_GATK4?=$(THREADS_BY_CALLER)
THREADS_GATK4_MUTECT2_FILTERED?=$(THREADS_BY_CALLER)
MINPRUNING_GATK4_MUTECT2_FILTERED?=20
MAXREADS_GATK4_MUTECT2_FILTERED?=1000

DPMIN_MUTECT2_FILTERED?=30
VAF_MUTECT2_FILTERED?=0.01
VAF_MUTECT2_FILTERED_HOM?=0.80


GATK4_MUTECT2_FILTERED_FLAGS_SHARED?=--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter --max-reads-per-alignment-start $(MAXREADS_GATK4_MUTECT2_FILTERED) --dont-use-soft-clipped-bases true --min-pruning $(MINPRUNING_GATK4_MUTECT2_FILTERED) --callable-depth $(DPMIN_MUTECT2_FILTERED) --verbosity ERROR --native-pair-hmm-threads $(THREADS_GATK4_MUTECT2_FILTERED)

%.MuTect2_filtered$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.dict %.design.bed.interval_list #%.from_manifest.interval_list
	# Calling by MuTect2
	$(JAVA11) $(JAVA_FLAGS) -jar $(GATK4) Mutect2 $(GATK4_MUTECT2_FILTERED_FLAGS_SHARED) \
		-R $$(cat $*.genome) \
		-I $< \
		-tumor $$(basename $< | cut -d"." -f1) \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-O $@.tmp.unfiltered.vcf;
	# Normalization
	grep "^##" $@.tmp.unfiltered.vcf | sed s/ID=TLOD,Number=A/ID=TLOD,Number=./gi > $@.tmp.unfiltered.TLOD.vcf
	grep "^##" -v $@.tmp.unfiltered.vcf | cut -f1-10 >> $@.tmp.unfiltered.TLOD.vcf
	# Sorting
	#$(JAVA11) -jar $(PICARD) SortVcf -I $@.tmp.unfiltered.TLOD.vcf -O $@.tmp.unfiltered.TLOD.sorted.vcf -SD $$(cat $*.dict);
	# HOWARD VAF calculation & Filtration by BCFTOOLS
	+if (($$($(BCFTOOLS) view -H $@.tmp.unfiltered.TLOD.vcf | wc -l ))); then \
		$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.tmp.unfiltered.TLOD.vcf --output=$@.tmp.unfiltered.TLOD.HOWARD.vcf --calculation="VAF"; \
		$(BCFTOOLS) view -i 'FORMAT/DP>=$(DPMIN_MUTECT2_FILTERED) && FORMAT/VAF>=$(VAF_MUTECT2_FILTERED) && FORMAT/VAF<$(VAF_MUTECT2_FILTERED_HOM)' $@.tmp.unfiltered.TLOD.HOWARD.vcf | $(BCFTOOLS) view -e 'GT="0/0"' -o $@.tmp.cleaned.vcf; \
	else \
		cp $@.tmp.unfiltered.TLOD.vcf $@.tmp.cleaned.vcf; \
	fi;
	# Sorting and contig
	$(JAVA11) -jar $(PICARD) SortVcf -I $@.tmp.cleaned.vcf -O $@ -SD $$(cat $*.dict);
	# Empty
	#if [ ! -e $@.tmp.FilterMutectCalls.TLOD.sorted.HOWARD.bcftools.vcf ]; then cp $*.empty.vcf $@.tmp.FilterMutectCalls.TLOD.sorted.HOWARD.bcftools.vcf; fi;
	# Copy
	#if [ ! -e $@ ]; then cp $@.tmp.FilterMutectCalls.TLOD.sorted.HOWARD.bcftools.vcf $@; fi;
	# Clean
	rm -f $@.tmp* $@.idx






###########################
# MUTECT2 stringent GATK4 #
###########################

#THREADS_GATK4?=$(THREADS_BY_CALLER)
THREADS_GATK4_MUTECT2_STRINGENT?=$(THREADS_BY_CALLER)
MINPRUNING_GATK4_MUTECT2_STRINGENT?=20
MAXREADS_GATK4_MUTECT2_STRINGENT?=1000

DPMIN_MUTECT2_STRINGENT?=30
VAF_MUTECT2_STRINGENT?=0.01
VAF_MUTECT2_STRINGENT_HOM?=0.80


GATK4_MUTECT2_STRINGENT_FLAGS_SHARED?=--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter --max-reads-per-alignment-start $(MAXREADS_GATK4_MUTECT2_STRINGENT) --dont-use-soft-clipped-bases true --min-pruning $(MINPRUNING_GATK4_MUTECT2_STRINGENT) --callable-depth $(DPMIN_MUTECT2_STRINGENT) --verbosity ERROR --native-pair-hmm-threads $(THREADS_GATK4_MUTECT2_STRINGENT)

#--max-reads-per-alignment-start
#--dont-use-soft-clipped-bases

%.MuTect2_stringent$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.dict %.design.bed.interval_list #%.from_manifest.interval_list
	# Calling by MuTect2
	$(JAVA11) $(JAVA_FLAGS) -jar $(GATK4) Mutect2 $(GATK4_MUTECT2_STRINGENT_FLAGS_SHARED) \
		-R $$(cat $*.genome) \
		-I $< \
		-tumor $$(basename $< | cut -d"." -f1) \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-O $@.tmp.unfiltered.vcf;
	# Filtration by MuTect2
	$(JAVA11) $(JAVA_FLAGS) -jar $(GATK4) FilterMutectCalls \
		-R $$(cat $*.genome) \
		-V $@.tmp.unfiltered.vcf \
		-O $@.tmp.FilterMutectCalls.vcf;
	# Normalization
	grep "^##" $@.tmp.FilterMutectCalls.vcf | sed s/ID=TLOD,Number=A/ID=TLOD,Number=./gi > $@.tmp.FilterMutectCalls.TLOD.vcf
	grep "^##" -v $@.tmp.FilterMutectCalls.vcf | cut -f1-10 >> $@.tmp.FilterMutectCalls.TLOD.vcf
	# Sorting
	$(JAVA11) -jar $(PICARD) SortVcf -I $@.tmp.FilterMutectCalls.TLOD.vcf -O $@.tmp.FilterMutectCalls.TLOD.sorted.vcf -SD $$(cat $*.dict);
	# Filtration by BCFTOOLS
	$(BCFTOOLS) view -i 'FORMAT/DP>=$(DPMIN_MUTECT2_STRINGENT) && FORMAT/AF>=$(VAF_MUTECT2_STRINGENT) && FORMAT/AF<$(VAF_MUTECT2_STRINGENT_HOM)' $@.tmp.FilterMutectCalls.TLOD.sorted.vcf | $(BCFTOOLS) view -f PASS  -e 'GT="0/0"' -o $@;
	# Clean
	rm -f $@.tmp* $@.idx






#################
# MUTECT2 GATK3 #
#################


#JAVA_FLAGS?= -Xmx8g
DFRAC_GATK3_MUTECT2?=1
MBQ_GATK3_MUTECT2?=17
VCFDBSNP_GATK3_MUTECT2?=$(VCFDBSNP)
GATK3_MUTECT2_FLAGS_SHARED?=--baq OFF --read_filter BadCigar --allow_potentially_misencoded_quality_scores --dontUseSoftClippedBases
MINPRUNING_GATK3_MUTECT2?=2
THREADS_GATK3_MUTECT2?=$(THREADS_BY_SAMPLE) #$(THREADS_GATK)
maxReadsInRegionPerSample_GATK3_MUTECT2?=1000
GATK3_MUTECT2_FLAGS= -nct $(THREADS_GATK3_MUTECT2) -stand_call_conf 10 -dfrac $(DFRAC_GATK3_MUTECT2) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample_GATK3_MUTECT2) --dbsnp $(VCFDBSNP_GATK3_MUTECT2) -mbq $(MBQ_GATK3_MUTECT2) -minPruning $(MINPRUNING_GATK3_MUTECT2)  $(GATK3_MUTECT2_FLAGS_SHARED)


%.GATK3_MuTect2$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.dict %.design.bed.interval_list #%.from_manifest.interval_list
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATK3_MUTECT2_FLAGS) \
	 	-T MuTect2 \
		-R $$(cat $*.genome) \
		-I:tumor $< \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-o $@.tmp1;
	# Clean header
	cat $@.tmp1 | sed s/ID=QSS,Number=A/ID=QSS,Number=./gi > $@.tmp.vcf
	# sort
	$(JAVA11) -jar $(PICARD) SortVcf -I $@.tmp.vcf -O $@.tmp2.vcf -SD $$(cat $*.dict);
	# DPMIN_MUTECT2
	#$(VCFUTILS) varFilter -d $(DPMIN_MUTECT2) $@.tmp2.vcf > $@.tmp3.vcf;
	$(BCFTOOLS) view  -i 'FORMAT/DP>=$(DPMIN_MUTECT2)' $@.tmp2.vcf > $@.tmp3.vcf;
	# Empty
	if [ ! -e $@.tmp3.vcf ]; then cp $*.empty.vcf $@.tmp3.vcf; fi;
	# Copy
	if [ ! -e $@ ]; then cp $@.tmp3.vcf $@; fi;
	# Clean
	rm -f $@.tmp* $@.idx





# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING: MUTECT to identify somatic variants and generate *.MuTect.vcf files. Parameteres MUTECT_INTERVAL_PADDING='$(MUTECT_INTERVAL_PADDING)' "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:MuTect:MuTect - by default, no filtering. MUTECT_INTERVAL_PADDING as parameter"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING: MUTECT2 with GATK4 to identify somatic variants and generate *.MuTect2.vcf files, with DP filtration. Parameters GATK4_MUTECT2_FLAGS_SHARED='$(GATK4_MUTECT2_FLAGS_SHARED)', DPMIN_MUTECT2='$(DPMIN_MUTECT2)' "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:MuTect2:MuTect2 with GATK4 - by default with DP filtration. "
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING: MUTECT2 with GATK4 to identify somatic variants and generate *.MuTect2_filtered.vcf files, with DP and VAF filtration, with FilterMutectCalls filteration by removing variant WT. Parameters GATK4_MUTECT2_FILTERED_FLAGS_SHARED='$(GATK4_MUTECT2_FILTERED_FLAGS_SHARED)', DPMIN_MUTECT2_FILTERED='$(DPMIN_MUTECT2_FILTERED)', VAF_MUTECT2_FILTERED='$(VAF_MUTECT2_FILTERED)', VAF_MUTECT2_FILTERED_HOM='$(VAF_MUTECT2_FILTERED_HOM)' "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:MuTect2_filtered:MuTect2 with GATK4 and FilterMutectCalls - by default with DP, VAF and WildType filtration. "
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING: MUTECT2 with GATK4 to identify somatic variants and generate *.MuTect2_stringent.vcf files, with DP and VAF filtration, with FilterMutectCalls filteration by removing variant WT and with quality not PASS. Parameters GATK4_MUTECT2_STRINGENT_FLAGS_SHARED='$(GATK4_MUTECT2_STRINGENT_FLAGS_SHARED)', DPMIN_MUTECT2_STRINGENT='$(DPMIN_MUTECT2_STRINGENT)', VAF_MUTECT2_STRINGENT='$(VAF_MUTECT2_STRINGENT)', VAF_MUTECT2_STRINGENT_HOM='$(VAF_MUTECT2_STRINGENT_HOM)' "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:MuTect2_stringent:MuTect2 with GATK4 and FilterMutectCalls - by default with DP, VAF, WildType and FilterMutectCalls filtration. "
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING: MUTECT2 with GATK3 to identify somatic variants and generate *.GATK3_MuTect2.vcf files. Parameters GATK3_MUTECT2_FLAGS='$(GATK3_MUTECT2_FLAGS)' "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:GATK3_MuTect2:MuTect2 with GATK3 - by default. "
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
