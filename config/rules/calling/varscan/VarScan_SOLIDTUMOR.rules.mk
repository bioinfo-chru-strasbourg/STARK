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
# 0.9.5-03/02/2023: Extract VarScan rules SOLIDTUMOR
# 0.9.6-12/05/2023: Remove VariantFiltration of GATK



######################
# VarScan SOLIDTUMOR #
######################

VARSCAN_SOLIDTUMOR_VAF=0.01
VARSCAN_SOLIDTUMOR_ALT=5
VARSCAN_SOLIDTUMOR_DP=50
VARSCAN_SOLIDTUMOR_PVAL=1e-3
VARSCAN_SOLIDTUMOR_BOTH_OPTIONS= --min-var-freq $(VARSCAN_SOLIDTUMOR_VAF) --min-reads2 $(VARSCAN_SOLIDTUMOR_ALT) --min-coverage $(VARSCAN_SOLIDTUMOR_DP) --p-value $(VARSCAN_SOLIDTUMOR_PVAL)
VARSCAN_SOLIDTUMOR_SNP_OPTIONS= $(VARSCAN_SOLIDTUMOR_BOTH_OPTIONS) --min-avg-qual 30
VARSCAN_SOLIDTUMOR_INDEL_OPTIONS= $(VARSCAN_SOLIDTUMOR_BOTH_OPTIONS) --min-avg-qual 10


%.VarScan_SOLIDTUMOR$(POST_CALLING).SNP.vcf: %.bam.mpileup %.empty.vcf %.genome
	if [ -s $< ]; then \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.SNP.sample.txt; \
		cat $< | $(JAVA) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_SOLIDTUMOR_SNP_OPTIONS) --vcf-sample-list $*.SNP.sample.txt > $@; \
		rm -f $*.SNP.sample.txt; \
	fi;
	if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	# Remove intermediate files
	rm -f $@.idx $@.unfiltered.vcf*


%.VarScan_SOLIDTUMOR$(POST_CALLING).InDel.vcf: %.bam.mpileup %.empty.vcf %.genome
	if [ -s $< ]; then \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.InDel.sample.txt; \
		cat $< | $(JAVA) -jar $(VARSCAN) mpileup2indel --output-vcf 1 $(VARSCAN_SOLIDTUMOR_INDEL_OPTIONS) --vcf-sample-list $*.InDel.sample.txt > $@; \
		rm -f $*.InDel.sample.txt; \
	fi;
	if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	# Remove intermediate files
	rm -f $@.idx $@.unfiltered.vcf*


# CONFIG/RELEASE
RELEASE_CMD := $(shell echo "\#\# VARSCAN SOLIDTUMOR: identify variants and generate *.VarScan_SOLIDTUMOR.vcf files with parameters SAMTOOLS mpileup and VARSCAN SOLIDTUMOR application 'VARSCAN_SOLIDTUMOR_SNP_OPTIONS: $(VARSCAN_SOLIDTUMOR_SNP_OPTIONS)' 'VARSCAN_SOLIDTUMOR_INDEL_OPTIONS: $(VARSCAN_SOLIDTUMOR_INDEL_OPTIONS)' " >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:VarScan_SOLIDTUMOR:SAMTOOLS/VARSCAN - design for SOLIDTUMOR mutation, VAF\>$(VARSCAN_SOLIDTUMOR_VAF) ALT\>$(VARSCAN_SOLIDTUMOR_ALT) DP\>$(VARSCAN_SOLIDTUMOR_DP) p-value\<$(VARSCAN_SOLIDTUMOR_PVAL):VARSCAN_SOLIDTUMOR_SNP_OPTIONS='$(VARSCAN_SOLIDTUMOR_SNP_OPTIONS)', VARSCAN_SOLIDTUMOR_INDEL_OPTIONS='$(VARSCAN_SOLIDTUMOR_INDEL_OPTIONS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
