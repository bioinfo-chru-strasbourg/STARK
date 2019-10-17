############################
# MUTECT Calling Rules
# Release: 0.9
# Date: 17/10/2019
# Author: Antony Le Bechec
############################



##########
# MUTECT #
##########

# Parameters
#MUTECT_INPUT_TYPE=tumor
MUTECT_INTERVAL_PADDING?=0


%.MuTect$(POST_CALLING).vcf: %.bam %.bam.bai %.from_manifest.intervals %.empty.vcf %.genome
	# Calling
	$(JAVA7) -jar $(MUTECT) \
		--analysis_type MuTect \
		--reference_sequence $$(cat $*.genome) \
		$$(if [ "`grep ^ -c $*.from_manifest.intervals`" == "0" ]; then echo ""; else echo "--intervals $*.from_manifest.intervals"; fi;) \
		--input_file:tumor $< \
		--vcf $@.tmp1 \
		-ip $(MUTECT_INTERVAL_PADDING);
	# Normalize
	grep "^##" $@.tmp1 > $@.tmp
	grep "^##" -v $@.tmp1 | cut -f1-10 >> $@.tmp
	# Empty
	if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	# Copy
	if [ ! -e $@ ]; then cp $@.tmp $@; fi;
	# Clean
	rm -f $@.tmp* $@.idx





#################
# MUTECT2 GATK4 #
#################


GATK4_MUTECT2_FLAGS_SHARED?=--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter --verbosity ERROR


%.MuTect2$(POST_CALLING).vcf: %.bam %.bam.bai %.from_manifest.intervals %.empty.vcf %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK4) Mutect2 $(GATK4_MUTECT2_FLAGS_SHARED) \
		-R $$(cat $*.genome) \
		-I $< \
		-tumor $$(basename $< | cut -d"." -f1) \
		$$(if [ "`grep ^ -c $*.from_manifest.intervals`" == "0" ]; then echo ""; else echo "-L $*.from_manifest.intervals"; fi;) \
		-O $@.tmp1;
	# Normalize
	grep "^##" $@.tmp1 | sed s/ID=TLOD,Number=A/ID=TLOD,Number=./gi > $@.tmp
	grep "^##" -v $@.tmp1 | cut -f1-10 >> $@.tmp
	# Empty
	if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	# Copy
	if [ ! -e $@ ]; then cp $@.tmp $@; fi;
	# Clean
	rm -f $@.tmp* $@.idx





#################
# MUTECT2 GATK3 #
#################


JAVA_FLAGS?= -Xmx8g
DFRAC_GATK3_MUTECT2?=1
MBQ_GATK3_MUTECT2?=17
VCFDBSNP_GATK3_MUTECT2?=$(VCFDBSNP)
GATK3_MUTECT2_FLAGS_SHARED?=--baq OFF --read_filter BadCigar --allow_potentially_misencoded_quality_scores --dontUseSoftClippedBases
MINPRUNING_GATK3_MUTECT2?=2
THREADS_GATK3_MUTECT2?=$(THREADS_GATK)
maxReadsInRegionPerSample_GATK3_MUTECT2?=1000
GATK3_MUTECT2_FLAGS= -nct $(THREADS_GATK3_MUTECT2) -stand_call_conf 10 -dfrac $(DFRAC_GATK3_MUTECT2) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample_GATK3_MUTECT2) --dbsnp $(VCFDBSNP_GATK3_MUTECT2) -mbq $(MBQ_GATK3_MUTECT2) -minPruning $(MINPRUNING_GATK3_MUTECT2)  $(GATK3_MUTECT2_FLAGS_SHARED)


%.GATK3_MuTect2$(POST_CALLING).vcf: %.bam %.bam.bai %.from_manifest.intervals %.empty.vcf %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATK3_MUTECT2_FLAGS) \
	 	-T MuTect2 \
		-R $$(cat $*.genome) \
		-I:tumor $< \
		$$(if [ "`grep ^ -c $*.from_manifest.intervals`" == "0" ]; then echo ""; else echo "-L $*.from_manifest.intervals"; fi;) \
		-o $@.tmp1;
	# Clean header
	cat $@.tmp1 | sed s/ID=QSS,Number=A/ID=QSS,Number=./gi > $@.tmp
	# Empty
	if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	# Copy
	if [ ! -e $@ ]; then cp $@.tmp $@; fi;
	# Clean
	rm -f $@.tmp* $@.idx





# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING: MUTECT to identify variants and generate *.MuTect.vcf files. Parameteres MUTECT_INTERVAL_PADDING='$(MUTECT_INTERVAL_PADDING)' "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:MuTect:MuTect - by default, no filtering. MUTECT_INTERVAL_PADDING as parameter"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING: MUTECT2 with GATK4 to identify variants and generate *.MuTect2.vcf files. "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:MuTect2:MuTect2 with GATK4 - by default. "
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING: MUTECT2 with GATK3 to identify variants and generate *.GATK3_MuTect2.vcf files. Parameters GATK3_MUTECT2_FLAGS=$(GATK3_MUTECT2_FLAGS) "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:GATK3_MuTect2:MuTect2 with GATK3 - by default. "
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
