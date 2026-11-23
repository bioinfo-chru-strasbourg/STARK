############################
# MarkDuplicates Rules
# Author: PACHCHEK Sinthuja and SARAROLS Pauline
############################
# Release
MK_RELEASE=1.0.0
MK_DATE=12/11/2026

## Release note
# 25/07/2014: DEV into PROD. Creation V0.9. Splitting clipping. Option file from BED file. Checking Option file empty. Merge resulting SAM
# 19/01/2022: Changing UMI options
# 12/11/2026: Add MARKDUPLICATES_FLAGS for MarkDuplicates configuration


##################
# MarkDuplicates #
##################

# Barcode tag
BARCODE_TAG?=
BARCODE_TAG_FLAG?=$(shell if [ "$(BARCODE_TAG)" != "" ]; then echo "BARCODE_TAG=$(BARCODE_TAG)"; fi;)
UMI_TYPE?=
UMI_TYPE_FLAG?=$(shell if [ "$(UMI_TYPE)" == "duplex" ]; then echo "DUPLEX_UMI=True"; fi;)

# Flags
MARKDUPLICATES_FLAGS?= \
	$(BARCODE_TAG_FLAG) \
	$(UMI_TYPE_FLAG) \
	$(PICARD_MARKDUP_OPTICAL_DEDUP)


# Rule markDuplicates
# mark duplicates with PICARD Tools
%.bam: %.markduplicates.bam %.markduplicates.bam.bai
	# Create Metrics Directory
	mkdir -p $(@D) ;
	# MarkDuplicates if BAM sorted
	if (( $$($(SAMTOOLS) view -H $< | grep "^@HD.*VN:.*SO:coordinate" | wc -l))); then \
		mkdir -p $<.metrics; \
		$(JAVA) $(JAVA_FLAGS) -jar $(PICARD) MarkDuplicates \
			I=$< \
			O=$@ \
			M=$<.metrics/$(*F).markDuplicates.metrics \
			$(MARKDUPLICATES_FLAGS) \
			VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=1 TMP_DIR=$(TMP_FOLDER_TMP) \
			1>$<.metrics/$(*F).markDuplicates.metrics.log 2>$<.metrics/$(*F).markDuplicates.metrics.err  || (echo "WARNING : picard.jar MarkDuplicates failed !"); \
		test -s $@ || (echo "WARNING : $@ is empty after MarkDuplicates step !" ); \
	else \
		echo "Markduplicate not created" ; \
	fi;
	-rm $<;


# Documentation
$(file >>$(RELEASE_INFOS),## POST_ALIGNMENT markduplicates '$(MK_RELEASE)': Mark duplicated reads in BAM with PICARD MarkDuplicates. Parameters MARKDUPLICATES_FLAGS='$(MARKDUPLICATES_FLAGS)')
$(file >>$(PIPELINES_INFOS),POST_ALIGNMENT:markduplicates:Mark duplicated reads in BAM with PICARD MarkDuplicates:MARKDUPLICATES_FLAGS='$(MARKDUPLICATES_FLAGS)')
