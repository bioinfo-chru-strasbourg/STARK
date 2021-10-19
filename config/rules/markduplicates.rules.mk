############################
# MarkDuplicates Rules
# Author: PACHCHEK Sinthuja and SARAROLS Pauline
############################
# Release
MK_RELEASE="0.9.1beta"
MK_DATE="06/09/2016"

## Release note
# 25/07/2014: DEV into PROD. Creation V0.9. Splitting clipping. Option file from BED file. Checking Option file empty. Merge resulting SAM


BARCODE_TAG?=

# markDuplicates
# mark duplicates with PICARD Tools
#%.bam : %.markduplicates.bam
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
		$$(if [ "$(BARCODE_TAG)" != "" ]; then echo "BARCODE_TAG=$(BARCODE_TAG)"; fi;) \
		$$(if [ "$(UMI_TYPE)" == "duplex" ]; then echo "DUPLEX_UMI=1"; fi;) \
		VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=1 TMP_DIR=$(TMP_FOLDER_TMP) \
		$(PICARD_MARKDUP_OPTICAL_DEDUP) \
		1>$<.metrics/$(*F).markDuplicates.metrics.log 2>$<.metrics/$(*F).markDuplicates.metrics.err  || (echo "WARNING : picard.jar MarkDuplicates failed !"); \
		test -s $@ || (echo "WARNING : $@ is empty after MarkDuplicates step !" ); \
		#rm $<; \
	else \
		echo "Markduplicate not created" ; \
	fi;
	-rm $<;



RELEASE_COMMENT := "\#\# MARK_DUPLICATES: Mark duplicated reads in BAM with PICARD MarkDuplicates. Not parallelized"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_ALIGNMENT:markduplicates:Mark duplicated reads in BAM with PICARD MarkDuplicates. Use BARCODE_TAG to specify tag"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

