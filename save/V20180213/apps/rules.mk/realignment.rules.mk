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
GATKRealignerTargetCreatorOptions= -known $(VCFDBSNP) -rf BadCigar -allowPotentiallyMisencodedQuals
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
	if (($(BAM_CHECK_STEPS))); then \
		if [ "$$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $<)" != "$$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $@)" ]; then \
			echo "# ERROR in Number of reads between $< and $@ !!!"; \
			echo "# BCFTOOLS STATS for $<.bam";  \
			$(SAMTOOLS) index $<.bam; \
			$(SAMTOOLS) stats $<.bam | grep SN; \
			$(SAMTOOLS) idxstats $<.bam; \
			echo "# BCFTOOLS STATS for $@";  \
			$(SAMTOOLS) index $@; \
			$(SAMTOOLS) stats $@ | grep SN; \
			$(SAMTOOLS) idxstats $@; \
			exit 1; \
		else \
			echo "# Number of reads OK between $< and $@"; \
		fi; \
	fi;
	#-rm $^
	#-rm $*.chr*.bai
	-rm $*.realignment.bam $*.realignment.bam.bai $*.realignment.*
	-rm -f $*.realignment*.mk
	



# Merge BAM
%.OLD.bam: %.chr1.bam %.chr2.bam %.chr3.bam %.chr4.bam \
	%.chr5.bam %.chr6.bam %.chr7.bam %.chr8.bam %.chr9.bam \
	%.chr10.bam %.chr11.bam %.chr12.bam %.chr13.bam %.chr14.bam %.chr15.bam \
	%.chr16.bam %.chr17.bam %.chr18.bam %.chr19.bam %.chr20.bam %.chr21.bam %.chr22.bam \
	%.chrX.bam %.chrY.bam %.chrM.bam %.unmapped.bam
	$(SAMTOOLS) merge -f $@ $^ -@ $(THREADS_BY_SAMPLE)
	# CHECK NUMBER of READS in BAM
	#echo $$($(SAMTOOLS) view -c -@ $(THREADS_SAMTOOLS) -F 0x0100 $*.realignment.bam)" READS for $*.realignment.bam";
	#echo $$($(SAMTOOLS) view -c -@ $(THREADS_SAMTOOLS) -F 0x0100 $@)" READS for $@";
	if (($(BAM_CHECK_STEPS))); then \
		if [ "$$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $*.realignment.bam)" != "$$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $@)" ]; then \
			echo "# ERROR in Number of reads between $.realignment.bam and $@ !!!"; \
			echo "# BCFTOOLS STATS for $*.realignment.bam";  \
			$(SAMTOOLS) index $*.realignment.bam; \
			$(SAMTOOLS) stats $*.realignment.bam | grep SN; \
			$(SAMTOOLS) idxstats $*.realignment.bam; \
			echo "# BCFTOOLS STATS for $@";  \
			$(SAMTOOLS) index $@; \
			$(SAMTOOLS) stats $@ | grep SN; \
			$(SAMTOOLS) idxstats $@; \
			exit 1; \
		else \
			echo "# Number of reads OK between $*.realignment.bam and $@"; \
		fi; \
	fi;
	-rm $^
	-rm $*.chr*.bai
	-rm $*.realignment.bam $*.realignment.bam.bai $*.realignment.*


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

