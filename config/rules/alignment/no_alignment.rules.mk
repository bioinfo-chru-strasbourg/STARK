############################
# NO ALIGNMENT Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE=1.0.0
MK_DATE=12/11/2026

# Release note
# 12/11/2026-1.0.0: Initial release with NO_ALIGNMENT rules


################
# NO_ALIGNMENT #
################

%.no_alignment.bam: %.R1.fastq.gz %.R2.fastq.gz %.bed
	# Generate dict
	# 4fields file
	awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t"$$4}' $*.bed > $@.tmp.4fields.tmp;
	# Clean bed with dict contig
	grep -Po 'SN:([^\t]*)' $(DICT) | cut -d: -f2 | sed "s/^/^/gi" | sed "s/$$/\t/gi" > $@.tmp.4fields.contig_from_dict;
	-grep -f $@.tmp.4fields.contig_from_dict $@.tmp.4fields.tmp > $@.tmp.4fields;
	# BedToIntervalList
	$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) BedToIntervalList -I $@.tmp.4fields -O $@.tmp.interval -SD $(DICT);
	# Alignment - Empty BAM
	grep "^@" $@.tmp.interval > $@.tmp.sam
	if (($$($(UNGZ) -c $*.R2.fastq.gz | head -n 1 | wc -l))); then \
		$(SAMTOOLS) import -1 $*.R1.fastq.gz -2 $*.R2.fastq.gz -T '*' -@ $(THREADS_SAMTOOLS) >> $@.tmp.sam; \
	else \
		$(SAMTOOLS) import -1 $*.R1.fastq.gz -T '*' -@ $(THREADS_SAMTOOLS) >> $@.tmp.sam; \
	fi;
	$(SAMTOOLS) view -O BAM,level=$(BAM_COMPRESSION) --no-PG -@ $(THREADS_SAMTOOLS) $@.tmp.sam > $@
	# Clean
	rm -rf $@.tmp*


# Documentation
$(file >>$(RELEASE_INFOS),## ALIGNER no_alignment '$(MK_RELEASE)': NO_ALIGNMENT generates a BAM file without performing any alignment. Parameters BAM_COMPRESSION='$(BAM_COMPRESSION)')
$(file >>$(PIPELINES_INFOS),ALIGNER:no_alignment:NO_ALIGNMENT - Generates BAM without alignment:BAM_COMPRESSION='$(BAM_COMPRESSION)')
