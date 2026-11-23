############################
# Main Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE=1.0.0
MK_DATE=17/11/2026

# Release note
# 10/03/2015-0.9.4beta: change genome reference location, in the file %.genome
# 14/06/2015-0.9.4.1b: change interval from bed (grep to remove '@' lines)
# 25/11/2015-0.9.4.2b: Change VCF sorting, using vcftools instead of igvtools
# 25/01/2016-0.9.4.3b: Change VCF sorting vcftools parameter, adding -C chromosome ordering
# 18/02/2016-0.9.4.4b: change cp to cat form manifest copy
# 13/04/2016-0.9.4.5b: removing full.vcf rule
# 29/09/2016-0.9.4.6b: Update PICARD to releasr picard.jar
# 07/05/2018-0.9.4.7b: Add --force for translation vcf to txt
# 02/10/2018-0.9.5b: Change HOWARD translation, prioritization and hard filtering. Change Manifest/bed link generation
# 27/09/2019-0.9.5.1b: Change FATBAM to CAP tool
# 28/07/2022-0.9.5.2: Change SNP and InDel merge, 2 rules for POST_CALLING and Callers (such as VarScan)
# 25/05/2023-0.9.5.3: Cleaning
# 25/06/2025-0.9.6.0: STARK release 19 compatibility
# 12/10/2026-0.9.6.1: Add update of VCF sequence dictionary
# 17/11/2026-1.0.0: Extract specifics rules such as vcf, bam, design into separate files



# OPTIONS
REMOVE_INTERMEDIATE_SAM?=1


# Genome REF
##############

# Rule: Write the genome file in a file
%.genome:
	printf '%s\n' "$(GENOME)" > "$@";


## SampleSheet
################


# SampleSheet Copy
%.SampleSheet.csv:
	-mkdir -p $(@D);
	-cp -p $(@D)/`echo $$(basename $(@D))`.SampleSheet.csv $@ 2>/dev/null || cp -p $(INPUTDIR)/`echo $$(basename $$(dirname $(@D)))`/SampleSheet.csv $@ 2>/dev/null || touch $@;


