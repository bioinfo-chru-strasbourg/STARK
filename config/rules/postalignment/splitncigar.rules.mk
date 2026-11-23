############################
# STAR Aligner Rules
# Author: Antony Le Béchec
############################
# Release
MK_RELEASE=1.1.0
MK_DATE=12/11/2026

# Release note
# 1.0.0-16/03/2026: Create splitNCigar as a post alignment step
# 1.1.0-12/11/2026: Add GATK_SPLITNCIGAR_FLAGS for additional options, fix analysis_type for GATK


####################
# SplitNCigarReads #
####################

# Splits reads that contain Ns in their cigar string (e.g. spanning splicing events in RNAseq data).
# Identifies all N cigar elements and creates k+1 new reads (where k is the number of N cigar elements).
# The first read includes the bases that are to the left of the first N element, while the part of the read that is to the right of the N (including the Ns) is hard clipped and so on for the rest of the new reads.
# Used for post-processing RNA reads aligned against the full reference.
# Post alignment spécific for STAR: we need to use the bam with splitNcigar for SNV calling.
# However, splitNcigar is forbiden for fusion detection tools (Arriba and STARFusion) to work properly (see STARFusion.rules.mk and Arriba.rules.mk for explanation).

# Flags
GATK_SPLITNCIGAR_FLAGS?=


# Rule
%.bam: %.splitncigar.bam %.splitncigar.bam.bai
	# SplitNCigarReads with GATK
	$(JAVA) $(JAVA_FLAGS_GATK4_CALLING_STEP) -jar $(GATK) SplitNCigarReads -R $(GENOME) -I $< -O $@ $(GATK_SPLITNCIGAR_FLAGS)
	# Clean
	rm -rf $*.splitncigar.bam $*.splitncigar.bam.bai $@.bai


# Documentation
$(file >>$(RELEASE_INFOS),## POST_ALIGNMENT splitncigar '$(MK_RELEASE)': Splits reads that contain Ns in their cigar string. Parameters GATK_SPLITNCIGAR_FLAGS='$(GATK_SPLITNCIGAR_FLAGS)')
$(file >>$(PIPELINES_INFOS),POST_ALIGNMENT:splitncigar:Splits reads that contain Ns in their cigar string:GATK_SPLITNCIGAR_FLAGS='$(GATK_SPLITNCIGAR_FLAGS)')
