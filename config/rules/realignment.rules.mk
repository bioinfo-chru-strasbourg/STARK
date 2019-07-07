############################
# GATK Realignment Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.3beta"
MK_DATE="10/03/2015"

# Release note
# 10/03/2015: change genome reference location, in the file %.genome

# TOOLS
JAVA?=java
GATK?=$(NGSbin)/GenomeAnalysisTK.jar
SAMTOOLS?=$(NGSbin)/samtools

## INTERVALS
##############

# FLAGS and Options

THREADS_RTC?=$(THREADS_BY_SAMPLE)
GATKRealignerTargetCreatorFLAGS= -nt $(THREADS_RTC) #-nt 8
#GATKRealignerTargetCreatorOptions= -known $(KNOWN_ALLELES) -known $(VCFDBSNP)
#GATKRealignerTargetCreatorOptions= -known $(VCFDBSNP) -rf BadCigar -allowPotentiallyMisencodedQuals
GATKRealignerTargetCreatorOptions= -known $(VCFDBSNP) -allowPotentiallyMisencodedQuals
#-fixMisencodedQuals

# Create a .intervals file from the RealignerTargetCreator tool
# for IndelRealigner (detect suspicious small indels regions)
%.intervals: %.bam %.bam.bai %.genome %.from_manifest.intervals #$(KNOWN_ALLELES) $(KNOWN_ALLELES).idx %.genome
	if [ "`grep ^ -c $*.from_manifest.intervals`" == "0" ]; then \
		$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKRealignerTargetCreatorFLAGS) $(GATKRealignerTargetCreatorOptions) \
		-T RealignerTargetCreator \
		-R `cat $*.genome` \
		-I $< \
		-o $@; \
	else \
		$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKRealignerTargetCreatorFLAGS) $(GATKRealignerTargetCreatorOptions) \
		-T RealignerTargetCreator \
		-R `cat $*.genome` \
		-I $< \
		-o $@ \
		-L $*.from_manifest.intervals; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;


#%.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
#	mv $< $@

#%.bam: $.realignment.bam
#	mv $< $@


#%.TESTREALIGNMENT.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
%.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	# IF READS
	if (($($(SAMTOOLS) idxstats $< | awk '{SUM+=$$3+$$4} END {print SUM}'))); then \
		echo "$*.unmapped.bam: $*.realignment.bam" >> $*.realignment1.mk; \
		echo "	$(SAMTOOLS) view -b -f 12 $*.realignment.bam > $*.unmapped.bam;" >> $*.realignment1.mk; \
		echo -n "$*.unmapped.bam" >> $*.realignment2.mk; \
		for chr in $$($(SAMTOOLS) idxstats $< | grep -v "\*" | awk '{ if ($$3+$$4>0) print $$1 }'); do \
			#echo $$chr  >> $*.realignment.mk; \
			echo "$*.$$chr.bam: $*.realignment.bam" >> $*.realignment1.mk; \
			echo "	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $*.realignment.bam -o $*.$$chr.bam -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L $$chr" >> $*.realignment1.mk; \
			echo -n " $*.$$chr.bam" >> $*.realignment2.mk; \
		done; \
		echo -n "$@: " | cat - $*.realignment2.mk > $*.realignment3.mk; \
		echo ""  >> $*.realignment3.mk; \
		echo "	$(SAMTOOLS) merge -f $@ `cat $*.realignment2.mk` -@ $(THREADS_BY_SAMPLE)" >> $*.realignment3.mk; \
		cat $*.realignment1.mk $*.realignment3.mk >> $*.realignment.mk; \
		cat $*.realignment.mk; \
		+make -f $*.realignment.mk $@; \
	else \
		cp $< $@; \
	fi;
	# clean
	-rm -f $*.realignment.bam $*.realignment.bam.bai $*.realignment.* $*.realignment*.mk




# grep "^chr19:" /media/data2/RES/ALL/150306_M01658_0044_000000000-ABWBB/HL0052/HL0052.bwamem.unrecalibrated.unclipped.realignment.intervals | tr '\n' ',' | sed 's/,/ -L /g'
# -L chr >>>> -L $$(grep "^chr?:" $*.realignment.intervals | tr '\n' ',' | sed 's/,/ -L /g') chr

GATKIndelRealignerFLAGS=
#GATKIndelRealignerOptions= -known $(VCFDBSNP) --LODThresholdForCleaning 2.0 -compress 0 --maxReadsForRealignment 20000 --maxReadsForConsensuses 120 --maxReadsInMemory 2000000 --maxConsensuses 30 -model USE_READS -allowPotentiallyMisencodedQuals
GATKIndelRealignerOptions= -known $(VCFDBSNP) --LODThresholdForCleaning 2.0 -compress 0 --maxReadsForRealignment 50000 --maxReadsForConsensuses 120 --maxReadsInMemory 2000000 --maxConsensuses 30 -model USE_READS -allowPotentiallyMisencodedQuals -dfrac 1
# --maxConsensuses 100 --maxReadsForConsensuses 500 # Too much ressources needed!!!
#-fixMisencodedQuals

# Unmapped reads (0x04) and unmapped mate (0x08) > -f 0x012 or 12

%.unmapped.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(SAMTOOLS) view -b -f 12 $< > $@;

# Split IndelRealiger into Chromosomes chr1..22,X,Y,M

%.chr1.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr1

%.chr2.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr2

%.chr3.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr3

%.chr4.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr4

%.chr5.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr5

%.chr6.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr6

%.chr7.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr7

%.chr8.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr8

%.chr9.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr9

%.chr10.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr10

%.chr11.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr11

%.chr12.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr12

%.chr13.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr13

%.chr14.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr14

%.chr15.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr15

%.chr16.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr16

%.chr17.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr17

%.chr18.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr18

%.chr19.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr19

%.chr20.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr20

%.chr21.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr21

%.chr22.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chr22

%.chrX.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chrX

%.chrY.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chrY

%.chrM.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $< -o $@ -targetIntervals $*.realignment.intervals $(GATKIndelRealignerOptions) -L chrM



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# LOCAL REALIGNEMNT '$(MK_RELEASE)': Using GATK IndelRealigner and SAMTOOLS, BAM file is locally realigned. Options: GATKIndelRealignerFLAGS='$(GATKIndelRealignerFLAGS)', GATKIndelRealignerOptions='$(GATKIndelRealignerOptions)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "POST_ALIGNMENT:realignment:Local Realignment of reads in BAM"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
