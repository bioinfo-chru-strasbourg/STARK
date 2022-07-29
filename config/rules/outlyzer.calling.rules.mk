############################
# OutLyzer Calling Rules
# Release: 0.9.1.0
# Date: 11/06/2021
# Author: Antony Le Bechec
############################


#############
# OUTLYZER  #
#############

DPMIN_OUTLYZER?=30
VAF_OUTLYZER?=0.01


THREADS_OUTLYZER2?=$(THREADS_BY_CALLER)

%.outLyzer_python2$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.dict %.design.bed
	# create tmp Directory
	mkdir -p $@.outlyser_tmp
	# Generate VCF with OutLyser
	if [ -s $*.design.bed ]; then \
		$(PYTHON2) $(OUTLYZER) calling -pythonPath=$(PYTHON2) -samtools=$(SAMTOOLS) -bed $*.design.bed -bam $< -ref $$(cat $*.genome) -output $@.outlyser_tmp/ -core $(THREADS_OUTLYZER2) -verbose 1 ; \
		# Normalize OutLyzer output VCF ; \
		cat $@.outlyser_tmp/*.vcf | awk -f $(STARK_FOLDER_BIN)/outlyzer_norm.awk > $@.tmp.vcf ; \
		# sort and contig ; \
		$(JAVA11) -jar $(PICARD) SortVcf -I $@.tmp.vcf -O $@.tmp2.vcf -SD $$(cat $*.dict) ; \
		$(BCFTOOLS) view  -i 'FORMAT/DP>=$(DPMIN_OUTLYZER)' $@.tmp2.vcf > $@; \
	else \
		cp $*.empty.vcf $@ ; \
		echo "[ERROR] Empty BED. Empty OutLyzer VCF" ; \
	fi;
	# Cleaning
	rm -rf $@.outlyser_tmp $@.tmp* $@.idx


# $(VCFUTILS) varFilter -d $(DPMIN_OUTLYZER) $@.tmp2.vcf > $@;
		


THREADS_OUTLYZER?=$(THREADS_BY_CALLER)

%.outLyzer$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.dict %.design.bed
	# create tmp Directory
	mkdir -p $@.outlyser_tmp
	# Generate VCF with OutLyser
	if [ -s $*.design.bed ]; then \
		$(PYTHON) $(OUTLYZER) calling -pythonPath=$(PYTHON) -samtools=$(SAMTOOLS) -bed $*.design.bed -bam $< -ref $$(cat $*.genome) -output $@.outlyser_tmp/ -core $(THREADS_OUTLYZER) -verbose 1; \
		# Normalize OutLyzer output VCF ; \
		cat $@.outlyser_tmp/*.vcf | awk -f $(STARK_FOLDER_BIN)/outlyzer_norm.awk > $@.tmp.vcf ; \
		# sort and contig ; \
		$(JAVA11) -jar $(PICARD) SortVcf -I $@.tmp.vcf -O $@.tmp2.vcf -SD $$(cat $*.dict) ; \
		$(BCFTOOLS) view  -i 'FORMAT/DP>=$(DPMIN_OUTLYZER)' $@.tmp2.vcf > $@; \
	else \
		cp $*.empty.vcf $@ ; \
		echo "[ERROR] Empty BED. Empty OutLyzer VCF" ; \
	fi;
	# Cleaning
	rm -rf $@.outlyser_tmp $@.tmp* $@.idx



THREADS_OUTLYZER?=$(THREADS_BY_CALLER)

%.outLyzer_filtered$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.dict %.design.bed
	# create tmp Directory
	mkdir -p $@.outlyser_tmp
	# Generate VCF with OutLyser
	if [ -s $*.design.bed ]; then \
		$(PYTHON) $(OUTLYZER) calling -pythonPath=$(PYTHON) -samtools=$(SAMTOOLS) -bed $*.design.bed -bam $< -ref $$(cat $*.genome) -output $@.outlyser_tmp/ -core $(THREADS_OUTLYZER) -verbose 1; \
		# Normalize OutLyzer output VCF ; \
		cat $@.outlyser_tmp/*.vcf | awk -f $(STARK_FOLDER_BIN)/outlyzer_norm.awk > $@.tmp.normalized.vcf ; \
		# sort and contig ; \
		$(JAVA11) -jar $(PICARD) SortVcf -I $@.tmp.normalized.vcf -O $@.tmp.normalized.sorted.vcf -SD $$(cat $*.dict) ; \
		# filtration by bcftools \
		$(BCFTOOLS) view -i 'FORMAT/DP>=$(DPMIN_OUTLYZER) && FORMAT/VAF>=$(VAF_OUTLYZER)' $@.tmp.normalized.sorted.vcf > $@; \
	else \
		cp $*.empty.vcf $@ ; \
		echo "[ERROR] Empty BED. Empty OutLyzer VCF" ; \
	fi;
	# Cleaning
	rm -rf $@.outlyser_tmp $@.tmp* $@.idx





	

#$(VCFUTILS) varFilter -d $(DPMIN_OUTLYZER) $@.tmp2.vcf > $@;

#$(JAVA11) -jar $(PICARD) SortVcf I=$@.tmp O=$@ SEQUENCE_DICTIONARY=$$(cat $*.dict)

RELEASE_CMD := $(shell echo "\#\# OutLyzer: OutLyzer caller to detect variant with low VAF, generate *.outLyzer.vcf files. VCF file is modified to include GT/Genotype information, depending on AF value - greater than 80 is 1/1, othewise 0/1. Filtered on DP. DPMIN_OUTLYZER='$(DPMIN_OUTLYZER)' " >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:outLyzer:OutLyzer - OutLyzer caller, GT info added, filter on DP:DP and VAF"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )



RELEASE_CMD := $(shell echo "\#\# OutLyzer: OutLyzer caller to detect variant with low VAF, generate *.outLyzer_filtered.vcf files. VCF file is modified to include GT/Genotype information, depending on AF value - greater than 80 is 1/1, othewise 0/1. Filtered on DP and VAF. DPMIN_OUTLYZER='$(DPMIN_OUTLYZER)', VAF_OUTLYZER='$(VAF_OUTLYZER)' " >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:outLyzer_filtered:OutLyzer - OutLyzer caller, GT info added, filter on DP and VAF:DP and VAF"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
