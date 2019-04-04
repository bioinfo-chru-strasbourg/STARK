############################
# MarkDuplicates Rules
# Author: PACHCHEK Sinthuja and SARAROLS Pauline
############################
# Release
MK_RELEASE="0.9.1beta"
MK_DATE="06/09/2016"

## Release note
# 25/07/2014: DEV into PROD. Creation V0.9. Splitting clipping. Option file from BED file. Checking Option file empty. Merge resulting SAM

# markDuplicates
# mark duplicates with PICARD Tools
#%.bam : %.markduplicates.bam
%.bam : %.markduplicates.bam %.markduplicates.bam.bai
	# Create Metrics Directory
	mkdir -p $(@D) ;
	# MarkDuplicates if BAM sorted
	if (( $$($(SAMTOOLS) view -H $< | grep "^@HD.*VN:.*SO:coordinate" | wc -l))); then \
		#$(JAVA) $(JAVA_FLAGS) -jar $(PICARDLIB)/MarkDuplicates.jar I=$*.unMD.bam O=$*.bam M=`echo $* | cut -d"." -f1,2`.markDuplicates.metrics.txt TMP_DIR=$(TMP_FOLDER_TMP) VALIDATION_STRINGENCY=SILENT || (echo "WARNING : MarkDuplicates.jar failed !"); \
		#$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) MarkDuplicates I=$< O=$@ M=$(@D)/$(*F).markDuplicates.metrics TMP_DIR=$(TMP_FOLDER_TMP) VALIDATION_STRINGENCY=SILENT 1>$(@D)/$(*F).markDuplicates.metrics.log 2>$(@D)/$(*F).markDuplicates.metrics.err  || (echo "WARNING : picard.jar MarkDuplicates failed !"); \
		mkdir -p $<.metrics; \
		$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) MarkDuplicates I=$< O=$@ M=$<.metrics/$(*F).markDuplicates.metrics TMP_DIR=$(TMP_FOLDER_TMP) VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=1 1>$<.metrics/$(*F).markDuplicates.metrics.log 2>$<.metrics/$(*F).markDuplicates.metrics.err  || (echo "WARNING : picard.jar MarkDuplicates failed !"); \
		test -s $*.bam || (echo "WARNING : $*.bam is empty after MarkDuplicates step !" ); \
		#rm $<; \
	else \
		echo "Markduplicate not created" ; \
	fi;
	-rm $<;

%.bam : %.markduplicatesSAMTOOLS.bam %.markduplicatesSAMTOOLS.bam.bai
	# Create Metrics Directory
	mkdir -p $(@D) ;
	# IF READS
	if (($($(SAMTOOLS) idxstats $< | awk '{SUM+=$$3+$$4} END {print SUM}'))); then \
		echo "$*.unmapped.bam: $*.markduplicates.bam" >> $*.markduplicates1.mk; \
		echo "	$(SAMTOOLS) view -b -f 12 $*.markduplicates.bam > $*.unmapped.bam;" >> $*.markduplicates1.mk; \
		echo -n "$*.unmapped.bam" >> $*.markduplicates2.mk; \
		for chr in $$($(SAMTOOLS) idxstats $< | grep -v "\*" | awk '{ if ($$3+$$4>0) print $$1 }'); do \
			#echo $$chr  >> $*.markduplicates.mk; \
			echo "$*.$$chr.bam: $*.markduplicates.bam" >> $*.markduplicates1.mk; \
			#echo "	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKIndelRealignerFLAGS) -T IndelRealigner -R `cat $*.genome` -I $*.markduplicates.bam -o $*.$$chr.bam -targetIntervals $*.markduplicates.intervals $(GATKIndelRealignerOptions) -L $$chr" >> $*.markduplicates1.mk; \
			echo " $(SAMTOOLS) view $*.markduplicates.bam $$chr | $(SAMTOOLS) rmdup - $*.$$chr.bam -S" >> $*.markduplicates1.mk; \
			echo -n " $*.$$chr.bam" >> $*.markduplicates2.mk; \
		done; \
		echo -n "$@: " | cat - $*.markduplicates2.mk > $*.markduplicates3.mk; \
		echo ""  >> $*.markduplicates3.mk; \
		echo "	$(SAMTOOLS) merge -f $@ `cat $*.markduplicates2.mk` -@ $(THREADS_BY_SAMPLE)" >> $*.markduplicates3.mk; \
		echo "	-rm -f `cat $*.markduplicates2.mk` " >> $*.markduplicates3.mk; \
		cat $*.markduplicates1.mk $*.markduplicates3.mk >> $*.markduplicates.mk; \
		cat $*.markduplicates.mk; \
		+make -f $*.markduplicates.mk $@; \
	else \
		cp $< $@; \
	fi;
	# clean
	-rm -f $<;

# MarkDuplicate Rule generic
#%.markduplicatesPICARD.bam : %.markduplicates.bam
#	mv $< $@;



# Default MarkDuplicates Rule
#%.markduplicates.bam: %.markduplicatesPICARD.bam:
#	mv $< $@;



#%.markduplicates.bam : %.unmarkduplicates.bam
#	mv $< $@
#
%.unMD.bam: %.unmarkduplicates.bam
	mv $< $@


#RELEASE_COMMENT := "\#\# MARK_DUPLICATES '$(MK_RELEASE)': Mark duplicated reads in BAM. default PICARD algorithm"
#RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# MARK_DUPLICATES: Mark duplicated reads in BAM with PICARD MarkDuplicates. Not parallelized"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# MARK_DUPLICATES SAMTOOLS: Mark duplicated reads in BAM with SAMTOOLS rmdup. Parallelized, but no mark duplicates for unaligned reads"

RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


#PIPELINES_COMMENT := "POST_ALIGNMENT:markduplicates:Mark duplicated reads in BAM. default PICARD"
#PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "POST_ALIGNMENT:markduplicates:Mark duplicated reads in BAM with PICARD MarkDuplicates"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "POST_ALIGNMENT:markduplicatesSAMTOOLS:Mark duplicated reads in BAM. With SAMTOOLS rmdup"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
