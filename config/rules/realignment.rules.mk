############################
# GATK Realignment Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.3beta"
MK_DATE="10/03/2015"

# Release note
# 10/03/2015: change genome reference location, in the file %.genome


## INTERVALS
##############

# FLAGS and Options

THREADS_RTC?=$(THREADS_BY_SAMPLE)
GATKRealignerTargetCreatorFLAGS= -nt $(THREADS_RTC) #-nt 8
#GATKRealignerTargetCreatorOptions= -known $(KNOWN_ALLELES) -known $(VCFDBSNP)
#GATKRealignerTargetCreatorOptions= -known $(VCFDBSNP) -rf BadCigar -allowPotentiallyMisencodedQuals
GATKRealignerTargetCreatorOptions= -known $(VCFDBSNP) -allowPotentiallyMisencodedQuals
GATKIndelRealignerFLAGS=
#GATKIndelRealignerOptions= -known $(VCFDBSNP) --LODThresholdForCleaning 2.0 -compress 0 --maxReadsForRealignment 20000 --maxReadsForConsensuses 120 --maxReadsInMemory 2000000 --maxConsensuses 30 -model USE_READS -allowPotentiallyMisencodedQuals
GATKIndelRealignerOptions= -known $(VCFDBSNP) --LODThresholdForCleaning 2.0 -compress 0 --maxReadsForRealignment 50000 --maxReadsForConsensuses 120 --maxReadsInMemory 2000000 --maxConsensuses 30 -model USE_READS -allowPotentiallyMisencodedQuals -dfrac 1

#-fixMisencodedQuals

# Create a .intervals file from the RealignerTargetCreator tool
# for IndelRealigner (detect suspicious small indels regions)
# %.intervals: %.bam %.bam.bai %.genome %.from_manifest.intervals #$(KNOWN_ALLELES) $(KNOWN_ALLELES).idx %.genome
# 	if [ "`grep ^ -c $*.from_manifest.intervals`" == "0" ]; then \
# 		$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKRealignerTargetCreatorFLAGS) $(GATKRealignerTargetCreatorOptions) \
# 		-T RealignerTargetCreator \
# 		-R `cat $*.genome` \
# 		-I $< \
# 		-o $@; \
# 	else \
# 		$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKRealignerTargetCreatorFLAGS) $(GATKRealignerTargetCreatorOptions) \
# 		-T RealignerTargetCreator \
# 		-R `cat $*.genome` \
# 		-I $< \
# 		-o $@ \
# 		-L $*.from_manifest.intervals; \
# 	fi;
# 	if [ ! -e $@ ]; then touch $@; fi;


#%.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
#	mv $< $@

#%.bam: $.realignment.bam
#	mv $< $@


#%.TESTREALIGNMENT.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
# (.bed, .list, .picard, .interval_list, or .intervals)
#%.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
#%.bam: %.realignment.bam %.realignment.bam.bai %.realignment.from_manifest.intervals %.genome
# %.bam.old: %.realignment.bam %.realignment.bam.bai %.realignment.design.bed %.genome
# 	# IF READS
# 	rm -f $*.realignment*.mk;
# 	+if (($$($(SAMTOOLS) idxstats $< | awk '{SUM+=$$3+$$4} END {print SUM}'))); then \
# 		echo "$*.unmapped.bam: $*.realignment.bam" >> $*.realignment1.mk; \
# 		echo "	$(SAMTOOLS) view -b -f 12 $*.realignment.bam > $*.unmapped.bam;" >> $*.realignment1.mk; \
# 		echo -n " $*.unmapped.bam " > $*.realignment2.mk; \
# 		for chr in $$($(SAMTOOLS) idxstats $< | grep -v "\*" | awk '{ if ($$3+$$4>0) print $$1 }'); do \
# 			#echo $$chr  >> $*.realignment.mk; \
# 			echo "$*.$$chr.bam: $*.realignment.bam" >> $*.realignment1.mk; \
# 			echo "	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $*.realignment.bam -o $*.$$chr.bam -targetIntervals $*.realignment.design.bed $(GATKIndelRealignerOptions) -L $$chr" >> $*.realignment1.mk; \
# 			echo -n " $*.$$chr.bam " >> $*.realignment2.mk; \
# 		done; \
# 		echo -n "$@: " | cat - $*.realignment2.mk > $*.realignment3.mk; \
# 		echo ""  >> $*.realignment3.mk; \
# 		echo "	$(SAMTOOLS) merge -f $@ $$(cat $*.realignment2.mk) -@ $(THREADS_BY_SAMPLE)" >> $*.realignment3.mk; \
# 		cat $*.realignment1.mk $*.realignment3.mk >> $*.realignment.mk; \
# 		cat $*.realignment.mk; \
# 		make -f $*.realignment.mk $@; \
# 	else \
# 		cp $< $@; \
# 	fi;
# 	# clean
# 	#-rm -f $*.realignment.bam $*.realignment.bam.bai $*.realignment.* $*.realignment*.mk
# 	-rm -f $*.realignment.bam $*.realignment.bam.bai $*.realignment*.mk



#%.TESTREALIGNMENT.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
# (.bed, .list, .picard, .interval_list, or .intervals)
#%.bam: %.realignment.bam %.realignment.bam.bai %.realignment.intervals %.genome
#%.bam: %.realignment.bam %.realignment.bam.bai %.realignment.from_manifest.intervals %.genome
%.bam: %.realignment.bam %.realignment.bam.bai %.realignment.design.bed %.genome
	# RealignerTargetCreator 
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKRealignerTargetCreatorFLAGS) $(GATKRealignerTargetCreatorOptions) \
			-T RealignerTargetCreator \
			-R $$(cat $*.genome) \
			-I $< \
			-o $*.realignment.bam.RealignerTargetCreator.intervals \
			$$(if (($$(grep ^ -c $*.realignment.design.bed))); then echo "-L $*.realignment.design.bed"; fi;);
	# IF READS
	rm -f $*.realignment*.mk;
	+if (($$($(SAMTOOLS) idxstats $< | awk '{SUM+=$$3+$$4} END {print SUM}'))); then \
		echo "$*.unmapped.bam: $*.realignment.bam" >> $*.realignment1.mk; \
		echo "	$(SAMTOOLS) view -b -f 12 $*.realignment.bam > $*.unmapped.bam;" >> $*.realignment1.mk; \
		echo -n " $*.unmapped.bam " > $*.realignment2.mk; \
		for chr in $$($(SAMTOOLS) idxstats $< | grep -v "\*" | awk '{ if ($$3+$$4>0) print $$1 }'); do \
			#echo $$chr  >> $*.realignment.mk; \
			echo "$*.$$chr.bam: $*.realignment.bam" >> $*.realignment1.mk; \
			echo "	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R $$(cat $*.genome) -I $*.realignment.bam -o $*.$$chr.bam -targetIntervals $*.realignment.bam.RealignerTargetCreator.intervals $(GATKIndelRealignerOptions) -L $$chr" >> $*.realignment1.mk; \
			echo -n " $*.$$chr.bam " >> $*.realignment2.mk; \
		done; \
		echo -n "$@: " | cat - $*.realignment2.mk > $*.realignment3.mk; \
		echo ""  >> $*.realignment3.mk; \
		echo "	$(SAMTOOLS) merge -f $@ $$(cat $*.realignment2.mk) -@ $(THREADS_BY_SAMPLE)" >> $*.realignment3.mk; \
		cat $*.realignment1.mk $*.realignment3.mk >> $*.realignment.mk; \
		cat $*.realignment.mk; \
		make -f $*.realignment.mk $@; \
	else \
		cp $< $@; \
	fi;
	# clean
	#-rm -f $*.realignment.bam $*.realignment.bam.bai $*.realignment.* $*.realignment*.mk $*.realignment.bam.RealignerTargetCreator.intervals
	-rm -f $*.realignment.bam $*.realignment.bam.bai $*.realignment*.mk






# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# LOCAL REALIGNEMNT '$(MK_RELEASE)': Using GATK IndelRealigner and SAMTOOLS, BAM file is locally realigned. Options: GATKIndelRealignerFLAGS='$(GATKIndelRealignerFLAGS)', GATKIndelRealignerOptions='$(GATKIndelRealignerOptions)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "POST_ALIGNMENT:realignment:Local Realignment of reads in BAM"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
