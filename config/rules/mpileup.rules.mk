############################
# VarScan Calling Rules
# Release: 0.9.5
# Date: 03/02/2023
# Author: Antony Le Bechec
############################

# Release notes:
# 0.9.1-24/04/2015: Add 'varscanlowfreqlowcovsnp' calling
# 0.9.2-05/10/2015: Add 'VarScan' calling, default VarScan parameters
# 0.9.3-07/10/2015: Add filters on VarScan vcf
# 0.9.3.1-10/11/2015: Add pipeline VarScanHEMATO
# 0.9.4-10/11/2015: Add Merge SNP and InDel
# 0.9.4.1-03/05/2016: Bug correction, if empty VCF after calling. Add Release infos
# 0.9.4.2-04/05/2016: Bug correction
# 0.9.4.3-17/05/2016: Bug correction
# 0.9.4.4-21/03/2019: Change parameters for mpileup by adding FATBAM Soft Clip to Q0
# 0.9.4.5-27/09/2019: Change FATBAM to CAP tool
# 0.9.4.6-25/05/2021: Change VARSCAN_FILTERS
# 0.9.4.7-28/06/2021: Homogenize VARSCAN thresholds 
# 0.9.4.8-29/07/2022: Homogenize VARSCAN rules path 
# 0.9.5-03/02/2023: Extract mpileup rule



####################
# SAMTOOLS mpileup
####################

# Parameters
MPILEUP_OPTIONS?= -E -d 10000000 -C 50 -Q 10 -q 1


%.bam.mpileup: %.bam %.bam.bai %.design.bed
	if (( $$(grep -v "^#" -c $*.design.bed) )); then \
		$(SAMTOOLS) view $< -L $*.design.bed -h | perl $(CAP_SOFTCLIPTOQ0) -v1 | $(SAMTOOLS) mpileup - -f $(GENOME) -l $*.design.bed $(MPILEUP_OPTIONS) > $@; \
	else \
		$(SAMTOOLS) view $< -h | perl $(CAP_SOFTCLIPTOQ0) -v1 | $(SAMTOOLS) mpileup - -f $(GENOME) $(MPILEUP_OPTIONS) > $@; \
	fi;



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# MPILEUP '$(MK_RELEASE)': MPILEUP tool to generate a mpileup file from a BAM file: MPILEUP_OPTIONS='$(MPILEUP_OPTIONS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )
