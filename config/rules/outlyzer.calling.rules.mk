############################
# OutLyzer Calling Rules
# Release: 0.9
# Date: 23/10/2017
# Author: Antony Le Bechec
############################


#############
# OUTLYZER  #
#############

DPMIN_OUTLYZER=30


THREADS_OUTLYZER2?=$(THREADS_BY_CALLER)

%.outLyzer_python2$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.dict %.design.bed
	# create tmp Directory
	mkdir -p $@.outlyser_tmp
	# Generate VCF with OutLyser
	if [ -s $*.design.bed ]; then \
		$(PYTHON2) $(OUTLYZER) calling -pythonPath=$(PYTHON2) -samtools=$(SAMTOOLS) -bed $*.design.bed -bam $< -ref $$(cat $*.genome) -output $@.outlyser_tmp/ -core $(THREADS_OUTLYZER2) -verbose 1 ; \
		# Normalize OutLyzer output VCF ; \
		cat $@.outlyser_tmp/*.vcf | awk -f $(STARK_FOLDER_BIN)/outlyzer_norm.awk > $@.tmp ; \
		# sort and contig ; \
		$(JAVA) -jar $(PICARD) SortVcf -I $@.tmp -O $@.tmp2 -SD $$(cat $*.dict) ; \
		$(VCFUTILS) varFilter -d $(DPMIN_OUTLYZER) $@.tmp2 > $@; \
	else \
		cp $*.empty.vcf $@ ; \
		echo "[ERROR] Empty BED. Empty OutLyzer VCF" ; \
	fi;
	# Cleaning
	rm -rf $@.outlyser_tmp $@.tmp* $@.idx


THREADS_OUTLYZER?=$(THREADS_BY_CALLER)

%.outLyzer$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.dict %.design.bed
	# create tmp Directory
	mkdir -p $@.outlyser_tmp
	# Generate VCF with OutLyser
	if [ -s $*.design.bed ]; then \
		$(PYTHON) $(OUTLYZER) calling -pythonPath=$(PYTHON) -samtools=$(SAMTOOLS) -bed $*.design.bed -bam $< -ref $$(cat $*.genome) -output $@.outlyser_tmp/ -core $(THREADS_OUTLYZER) -verbose 1; \
		# Normalize OutLyzer output VCF ; \
		cat $@.outlyser_tmp/*.vcf | awk -f $(STARK_FOLDER_BIN)/outlyzer_norm.awk > $@.tmp ; \
		# sort and contig ; \
		$(JAVA) -jar $(PICARD) SortVcf -I $@.tmp -O $@.tmp2 -SD $$(cat $*.dict) ; \
		$(VCFUTILS) varFilter -d $(DPMIN_OUTLYZER) $@.tmp2 > $@; \
	else \
		cp $*.empty.vcf $@ ; \
		echo "[ERROR] Empty BED. Empty OutLyzer VCF" ; \
	fi;
	# Cleaning
	rm -rf $@.outlyser_tmp $@.tmp* $@.idx


#$(JAVA) -jar $(PICARD) SortVcf I=$@.tmp O=$@ SEQUENCE_DICTIONARY=$$(cat $*.dict)

RELEASE_CMD := $(shell echo "\#\# OutLyzer: OutLyzer caller to detect variant with low VAF, generate *.outLyzer.vcf files. VCF file is modified to include GT/Genotype information, depending on AF value - greater than 80 is 1/1, othewise 0/1. DPMIN_OUTLYZER=$(DPMIN_OUTLYZER) " >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:outLyzer:OutLyzer - OutLyzer caller, GT info added:no parameters"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
