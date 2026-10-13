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
	# Clean and create split directory
	rm -rf $@.split && mkdir -p $@.split
	# Generate idxstats for the BAM file
	$(SAMTOOLS) idxstats $< | awk '$$3>0' > $@.split/idxstats.txt
	# Check if the design BED file has any non-comment lines and create the list of chromosomes accordingly
	if (( $$(grep -v -e "^#" -e "^$$" -c $*.design.bed) )); then \
		grep -v -e "^#" -e "^$$" $*.design.bed | cut -f1 | sort -u > $@.split/bedchroms.txt; \
		cut -f1 $@.split/idxstats.txt | grep -Fx -f $@.split/bedchroms.txt > $@.split/chroms.txt; \
		VIEW_L="-L $*.design.bed"; MP_L="-l $*.design.bed"; \
	else \
		cut -f1 $@.split/idxstats.txt > $@.split/chroms.txt; \
		VIEW_L=""; MP_L=""; \
	fi;
	# Order the chromosomes based on the idxstats information
	awk 'NR==FNR{k[$$1];next} $$1 in k{print $$3"\t"$$1}' $@.split/chroms.txt $@.split/idxstats.txt \
		| sort -k1,1nr | cut -f2 > $@.split/order.txt;
	# Generate the sub Makefile for parallel processing of each chromosome part
	{ \
	printf '.DELETE_ON_ERROR:\n'; \
	printf 'all: %s\n' "$$(sed 's|.*|$@.split/&.part|' $@.split/order.txt | tr '\n' ' ')"; \
	printf '$@.split/%%.part:\n'; \
	printf '\t$(SAMTOOLS) view $< "$$*" %s -h | perl $(CAP_SOFTCLIPTOQ0) -v1 | $(SAMTOOLS) mpileup - -f $(GENOME) %s $(MPILEUP_OPTIONS) > $$@\n' "$$VIEW_L" "$$MP_L"; \
	} > $@.split/sub.mk
	# Execute the sub Makefile to process all chromosome parts in parallel
	+$(MAKE) -f $@.split/sub.mk all
	# Concatenate all chromosome parts into the final mpileup file
	sed 's|.*|$@.split/&.part|' $@.split/chroms.txt | xargs -r cat > $@
	# Clean up the split directory
	rm -rf $@.split


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# MPILEUP '$(MK_RELEASE)': MPILEUP tool to generate a mpileup file from a BAM file: MPILEUP_OPTIONS='$(MPILEUP_OPTIONS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )
