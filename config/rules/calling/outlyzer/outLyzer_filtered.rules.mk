############################
# OutLyzer Calling Rules
# Release: 0.9.2
# Date: 03/02/2023
# Author: Antony Le Bechec
############################


#####################
# OUTLYZER Filtered #
#####################

# Parameters
THREADS_OUTLYZER_FILTERED?=$(THREADS_BY_CALLER)
VAF_OUTLYZER_FILTERED?=0.01
DPMIN_OUTLYZER_FILTERED?=30

%.outLyzer_filtered$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.dict %.design.bed
	# create tmp Directory
	mkdir -p $@.outlyser_tmp
	# Generate VCF with OutLyser
	if [ -s $*.design.bed ]; then \
		$(PYTHON) $(OUTLYZER) calling -pythonPath=$(PYTHON) -samtools=$(SAMTOOLS) -bed $*.design.bed -bam $< -ref $$(cat $*.genome) -output $@.outlyser_tmp/ -core $(THREADS_OUTLYZER_FILTERED) -verbose 1; \
		# Normalize OutLyzer output VCF ; \
		cat $@.outlyser_tmp/*.vcf | awk -f $(STARK_FOLDER_BIN)/outlyzer_norm.awk > $@.tmp.normalized.vcf ; \
		# sort and contig ; \
		$(JAVA11) -jar $(PICARD) SortVcf -I $@.tmp.normalized.vcf -O $@.tmp.normalized.sorted.vcf -SD $$(cat $*.dict) ; \
		# filtration by bcftools \
		$(BCFTOOLS) view -i 'FORMAT/DP>=$(DPMIN_OUTLYZER_FILTERED) && FORMAT/VAF>=$(VAF_OUTLYZER_FILTERED)' $@.tmp.normalized.sorted.vcf > $@; \
	else \
		cp $*.empty.vcf $@ ; \
		echo "[ERROR] Empty BED. Empty OutLyzer VCF" ; \
	fi;
	# Cleaning
	rm -rf $@.outlyser_tmp $@.tmp* $@.idx



RELEASE_CMD := $(shell echo "\#\# OutLyzer: OutLyzer caller to detect variant with low VAF, generate *.outLyzer_filtered.vcf files. VCF file is modified to include GT/Genotype information, depending on AF value - greater than 80 is 1/1, othewise 0/1. Filtered on DP and VAF. DPMIN_OUTLYZER='$(DPMIN_OUTLYZER_FILTERED)', VAF_OUTLYZER='$(VAF_OUTLYZER_FILTERED)' " >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:outLyzer_filtered:OutLyzer - OutLyzer caller, GT info added, filter on DP and VAF:DP and VAF"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
