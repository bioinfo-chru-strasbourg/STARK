############################
# EMPTY ALIGNMENT Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE=1.0.0
MK_DATE=23/11/2026

# Release note
# 23/11/2026-1.0.0: Initial release with EMPTY_ALIGNMENT rules


###################
# EMPTY_ALIGNMENT #
###################

%.empty_alignment.bam: %.R1.fastq.gz %.R2.fastq.gz %.bed
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
	$(SAMTOOLS) view -O BAM,level=$(BAM_COMPRESSION) --no-PG -@ $(THREADS_SAMTOOLS) $@.tmp.sam > $@
	# Clean
	rm -rf $@.tmp*


# Documentation
$(file >>$(RELEASE_INFOS),## ALIGNER empty_alignment '$(MK_RELEASE)': EMPTY_ALIGNMENT generates an empty BAM file. Parameters BAM_COMPRESSION='$(BAM_COMPRESSION)')
$(file >>$(PIPELINES_INFOS),ALIGNER:empty_alignment:EMPTY_ALIGNMENT - Generates empty BAM:BAM_COMPRESSION='$(BAM_COMPRESSION)')