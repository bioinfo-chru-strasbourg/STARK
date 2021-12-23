############################
# No Pipeline Rules
# Release: 0.9.0.0
# Date: 29/09/2021
# Author: Antony Le Bechec
############################





################
# NO_ALIGNMENT #
################

%.no_alignment.bam: %.R1.fastq.gz %.R2.fastq.gz %.genome %.bed %.dict
	# Generate dict
	# 4fields file
	awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t"$$4}' $*.bed > $@.tmp.4fields.tmp;
	# Clean bed with dict contig
	grep -Po 'SN:([^\t]*)' $$(cat $*.dict) | cut -d: -f2 | sed "s/^/^/gi" | sed "s/$$/\t/gi" > $@.tmp.4fields.contig_from_dict;
	grep -f $@.tmp.4fields.contig_from_dict $@.tmp.4fields.tmp > $@.tmp.4fields;
	# BedToIntervalList
	#$(JAVA) $(JAVA_FLAGS_BY_SAMPLE) -jar $(PICARD) BedToIntervalList -I $@.tmp.4fields -O $@.tmp.interval -SD $$(cat $*.dict);
	$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) BedToIntervalList -I $@.tmp.4fields -O $@.tmp.interval -SD $$(cat $*.dict);
	# Alignment - Empty BAM
	#grep "^@" $@.tmp.interval | $(SAMTOOLS) view -b --no-PG > $@
	grep "^@" $@.tmp.interval > $@.tmp.sam
	if (($$(zcat $*.R2.fastq.gz | head -n 1 | wc -l))); then \
		$(SAMTOOLS) import -1 $*.R1.fastq.gz -2 $*.R2.fastq.gz -T '*' -@ $(THREADS_SAMTOOLS) >> $@.tmp.sam; \
	else \
		$(SAMTOOLS) import -1 $*.R1.fastq.gz -T '*' -@ $(THREADS_SAMTOOLS) >> $@.tmp.sam; \
	fi;
	$(SAMTOOLS) view -O BAM,level=$(BAM_COMPRESSION) --no-PG -@ $(THREADS_SAMTOOLS) $@.tmp.sam > $@
	# Clean
	rm -rf $@.tmp*
	


##############
# NO_CALLING #
##############

%.no_calling.vcf: %.bam %.bam.bai %.empty.vcf %.genome
	# Calling - Empty VCF
	cp $*.empty.vcf $@
	

#################
# NO_ANNOTATION #
#################

%.no_annotation.vcf: %.vcf
	# Annotation - Same VCF
	cp $< $@
	



# CONFIG/RELEASE

RELEASE_COMMENT := "\#\# ALIGNER: No alignemnt, generates an empty BAM "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "ALIGNER:no_alignment:No Alignment - generates an empty BAM"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


RELEASE_COMMENT := "\#\# CALLER: No calling, generates an empty VCF "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:no_calling:No Calling - generates an empty VCF"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


RELEASE_COMMENT := "\#\# ANNOTATOR: No annotation, no change in input VCF "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "ANNOTATOR:no_annotation:No Annotation - no change in VCF"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
