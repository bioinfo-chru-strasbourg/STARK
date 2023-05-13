############################
# Gencore Rules
# Authoring : Thomas LAVAUX
############################
# Release

MK_RELEASE=0.2.0.0
MK_DATE=12/09/2021

## Release note
# 0.1.0.0-12/12/2020 : DEV version
# 0.2.0.0-12/09/2021 : Pre-Prod version

%.bam: %.gencore.bam %.gencore.bam.bai %.genome %.gencore.design.bed 
	# gencore is a tool for fast and powerful deduplication for paired-end next-generation sequencing
	mkdir -p $<.metrics;
	cp -p $< $@;
	# Prepare BAM
	if [ "$$($(SAMTOOLS) view $< | head -n1)" != "$$($(SAMTOOLS) view $< | head -n1 | awk -v INPUT_FORMAT=BAM -v UMI_REFORMAT=1 -v UMI_TAG=1 -v UMI_SEP=_ -f $(FASTQ_CLEAN_HEADER))" ]; then \
		$(SAMTOOLS) view -h $< | awk -v INPUT_FORMAT=BAM -v UMI_REFORMAT=1 -v UMI_TAG=1 -v UMI_SEP=_ -f $(FASTQ_CLEAN_HEADER) | $(SAMTOOLS) view --bam --fast > $@.tmp.preformatted.bam; \
	else \
		ln -s $< $@.tmp.preformatted.bam; \
	fi;
	# gencore -i input.sorted.bam -o output.bam -r hg19.fasta -b test.bed
	$(GENCORE) -i $@.tmp.preformatted.bam -o $@.tmp.gencore.bam -r $$(cat $*.genome) -b $*.gencore.design.bed -j $<.metrics/$(*F).gencore.metrics.json -h $<.metrics/$(*F).gencore.metrics.html $(GENCORE_SUP_READS) $(GENCORE_SCORE_THREESHOLD) $(GENCORE_DIFF_THREESHOLD) $(GENCORE_RATIO_THREESHOLD) $(GENCORE_QUAL_THREESHOLD) $(GENCORE_QUAL_THREESHOLD);
	# reformat BAM TODO
	if [ "$$($(SAMTOOLS) view $@.tmp.gencore.bam | head -n1)" != "$$($(SAMTOOLS) view $@.tmp.gencore.bam | head -n1 | awk -v INPUT_FORMAT=BAM -v UMI_REFORMAT=1 -v UMI_TAG=1 -v UMI_SEP=- -f $(FASTQ_CLEAN_HEADER))" ]; then \
		$(SAMTOOLS) view -h $@.tmp.gencore.bam | awk -v INPUT_FORMAT=BAM -v UMI_REFORMAT=1 -v UMI_TAG=1 -v UMI_SEP=- -f $(FASTQ_CLEAN_HEADER) | $(SAMTOOLS) view --bam --fast > $@.tmp.reformatted.bam; \
	else \
		ln -s $@.tmp.gencore.bam $@.tmp.reformatted.bam; \
	fi;
	# Sort SAM by coordinate with Picard if gencore unsorted them during deduplication (Message in gencore log : WARNING: The output will be unordered)
	$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) SortSam -I $@.tmp.reformatted.bam -O $@ -SORT_ORDER coordinate;
	# Delete tempory files
	-rm $< $@.tmp.* ;




RELEASE_COMMENT := "\#\# GENCORE : gencore is a tool for fast and powerful deduplication for paired-end next-generation sequencing"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_ALIGNMENT:gencore:gencore is a tool for fast and powerful deduplication for paired-end next-generation sequencing"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )