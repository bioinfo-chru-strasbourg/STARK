############################
# OutLyzer Calling Rules
# Release: 0.9
# Date: 23/10/2017
# Author: Antony Le Bechec
############################



#############
# OUTLYZER  #
#############



%.outLyzer$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.from_manifest.bed
	# create tmp Directory
	mkdir -p $@.outlyser_tmp
	# Generate VCF with OutLyser
	#$(OUTLYZER) calling -bed $*.from_manifest.bed -bam $< -ref $$(cat $*.genome) -output $@.outlyser_tmp/ -verbose 1
	$(PYTHON2) $(OUTLYZER) calling -pythonPath=$(PYTHON2) -samtools=$(SAMTOOLS) -bed $*.from_manifest.bed -bam $< -ref $$(cat $*.genome) -output $@.outlyser_tmp/ -verbose 1
	# Normalize OutLyzer output VCF
	cat $@.outlyser_tmp/*.vcf | awk -f $(STARK_FOLDER_BIN)/outlyzer_norm.awk > $@
	# Cleaning
	rm -rf $@.outlyser_tmp



RELEASE_CMD := $(shell echo "\#\# OutLyzer: OutLyzer caller to detect variant with low VAF, generate *.outLyzer.vcf files. VCF file is modified to include GT/Genotype information, depending on AF value - greater than 80 is 1/1, othewise 0/1." >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:outLyzer:OutLyzer - OutLyzer caller, GT info added:no parameters"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )