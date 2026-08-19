############################
# STAR Aligner Rules
# Release: 1.0.0
# Date: 16/03/2026
# Author: Antony Le Béchec
############################

# Release note
# 1.0.0-16/03/2026: Create splitNCigar as a post alignment step

# Splits reads that contain Ns in their cigar string (e.g. spanning splicing events in RNAseq data).
# Identifies all N cigar elements and creates k+1 new reads (where k is the number of N cigar elements).
# The first read includes the bases that are to the left of the first N element, while the part of the read that is to the right of the N (including the Ns) is hard clipped and so on for the rest of the new reads.
# Used for post-processing RNA reads aligned against the full reference.
# Post alignment spécific for STAR: we need to use the bam with splitNcigar for SNV calling.
# However, splitNcigar is forbiden for fusion detection tools (Arriba and STARFusion) to work properly (see STARFusion.rules.mk and Arriba.rules.mk for explanation).
%.bam: %.splitncigar.bam %.splitncigar.bam.bai
	$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK4) SplitNCigarReads -R $(GENOME) -I $< -O $@



RELEASE_COMMENT := "\#\# BAM SPLITNCIGAR: GATK SplitNCigarReads is used to splits reads that contain Ns in their cigar string."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_ALIGNMENT:splitncigar:Splits reads that contain Ns in their cigar string."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )



