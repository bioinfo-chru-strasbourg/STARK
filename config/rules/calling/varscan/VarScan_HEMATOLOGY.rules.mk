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
# 0.9.5-03/02/2023: Extract VarScan rules HEMATOLOGY
# 0.9.6-12/05/2023: Remove VariantFiltration of GATK



######################
# VarScan HEMATOLOGY #
######################

VARSCAN_HEMATOLOGY_VAF=0.01
VARSCAN_HEMATOLOGY_ALT=5
VARSCAN_HEMATOLOGY_DP=50
VARSCAN_HEMATOLOGY_PVAL=1e-3
VARSCAN_HEMATOLOGY_BOTH_OPTIONS= --min-var-freq $(VARSCAN_HEMATOLOGY_VAF) --min-reads2 $(VARSCAN_HEMATOLOGY_ALT) --min-coverage $(VARSCAN_HEMATOLOGY_DP) --p-value $(VARSCAN_HEMATOLOGY_PVAL)
VARSCAN_HEMATOLOGY_SNP_OPTIONS= $(VARSCAN_HEMATOLOGY_BOTH_OPTIONS) --min-avg-qual 30
VARSCAN_HEMATOLOGY_INDEL_OPTIONS= $(VARSCAN_HEMATOLOGY_BOTH_OPTIONS) --min-avg-qual 10


%.VarScan_HEMATOLOGY$(POST_CALLING).SNP.vcf: %.bam.mpileup %.empty.vcf
	if [ -s $< ]; then \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.SNP.sample.txt; \
		cat $< | $(JAVA) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_HEMATOLOGY_SNP_OPTIONS) --vcf-sample-list $*.SNP.sample.txt > $@; \
		rm -f $*.SNP.sample.txt; \
	fi;
	if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	# Remove intermediate files
	rm -f $@.idx $@.unfiltered.vcf*


%.VarScan_HEMATOLOGY$(POST_CALLING).InDel.vcf: %.bam.mpileup %.empty.vcf
	if [ -s $< ]; then \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.InDel.sample.txt; \
		cat $< | $(JAVA) -jar $(VARSCAN) mpileup2indel --output-vcf 1 $(VARSCAN_HEMATOLOGY_INDEL_OPTIONS) --vcf-sample-list $*.InDel.sample.txt > $@; \
		rm -f $*.InDel.sample.txt; \
	fi;
	if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	# Remove intermediate files
	rm -f $@.idx $@.unfiltered.vcf*



# CONFIG/RELEASE
RELEASE_CMD := $(shell echo "\#\# VARSCAN HEMATOLOGY: identify variants and generate *.VarScan_HEMATOLOGY.vcf files with parameters SAMTOOLS mpileup and VARSCAN HEMATOLOGY application 'VARSCAN_HEMATOLOGY_SNP_OPTIONS: $(VARSCAN_HEMATOLOGY_SNP_OPTIONS)' 'VARSCAN_HEMATOLOGY_INDEL_OPTIONS: $(VARSCAN_HEMATOLOGY_INDEL_OPTIONS)' " >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:VarScan_HEMATOLOGY:SAMTOOLS/VARSCAN - design for HEMATOLOGY mutation, VAF\>$(VARSCAN_HEMATOLOGY_VAF) ALT\>$(VARSCAN_HEMATOLOGY_ALT) DP\>$(VARSCAN_HEMATOLOGY_DP) p-value\<$(VARSCAN_HEMATOLOGY_PVAL):VARSCAN_HEMATOLOGY_SNP_OPTIONS='$(VARSCAN_HEMATOLOGY_SNP_OPTIONS)', VARSCAN_HEMATOLOGY_INDEL_OPTIONS='$(VARSCAN_HEMATOLOGY_INDEL_OPTIONS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
