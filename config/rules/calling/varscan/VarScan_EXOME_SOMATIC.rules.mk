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
# 0.9.5-03/02/2023: Extract VarScan rules EXOME SOMATIC



#########################
# VarScan EXOME_SOMATIC #
#########################

VARSCAN_EXOME_SOMATIC_VAF=0.01
VARSCAN_EXOME_SOMATIC_ALT=3
VARSCAN_EXOME_SOMATIC_DP=30
VARSCAN_EXOME_SOMATIC_PVAL=1e-3
VARSCAN_EXOME_SOMATIC_BOTH_OPTIONS= --min-var-freq $(VARSCAN_EXOME_SOMATIC_VAF) --min-reads2 $(VARSCAN_EXOME_SOMATIC_ALT) --min-coverage $(VARSCAN_EXOME_SOMATIC_DP) --p-value $(VARSCAN_EXOME_SOMATIC_PVAL)
VARSCAN_EXOME_SOMATIC_SNP_OPTIONS= $(VARSCAN_EXOME_SOMATIC_BOTH_OPTIONS) --min-avg-qual 30
VARSCAN_EXOME_SOMATIC_INDEL_OPTIONS= $(VARSCAN_EXOME_SOMATIC_BOTH_OPTIONS) --min-avg-qual 10
VARSCAN_EXOME_SOMATIC_FILTERS=--genotypeFilterExpression 'GQ < 20.0' --genotypeFilterName 'LowGQ'
VARSCAN_EXOME_SOMATIC_SNP_FILTERS=$(VARSCAN_EXOME_SOMATIC_FILTERS)
VARSCAN_EXOME_SOMATIC_INDEL_FILTERS=$(VARSCAN_EXOME_SOMATIC_FILTERS)

%.VarScan_EXOME_SOMATIC$(POST_CALLING).SNP.vcf: %.bam.mpileup %.empty.vcf %.genome #
	-if [ -s $< ]; then \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.SNP.sample.txt; \
		#timeout $(VARSCAN_TIMEOUT) cat $< | $(JAVA11) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_EXOME_SOMATIC_SNP_OPTIONS) --vcf-sample-list $*.SNP.sample.txt > $@.unfiltered.vcf; \
		cat $< | $(JAVA11) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_EXOME_SOMATIC_SNP_OPTIONS) --vcf-sample-list $*.SNP.sample.txt > $@.unfiltered.vcf; \
		rm -f $*.SNP.sample.txt; \
	fi;
	-if [ ! -s $@.unfiltered.vcf ]; then cp $*.empty.vcf $@.unfiltered.vcf; fi;
	-$(JAVA11) -jar $(IGVTOOLS) index $@.unfiltered.vcf;
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $@.unfiltered.vcf -o $@ $(VARSCAN_EXOME_SOMATIC_SNP_FILTERS);
	-if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -s $@ ]; then touch $@; fi;
	# Remove intermediate files
	-rm $@.unfiltered.vcf*
	# Remove IDX
	-rm $@.idx

%.VarScan_EXOME_SOMATIC$(POST_CALLING).InDel.vcf: %.bam.mpileup %.empty.vcf %.genome
	-if [ -s $< ]; then \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.InDel.sample.txt; \
		cat $< | $(JAVA11) -jar $(VARSCAN) mpileup2indel --output-vcf 1 $(VARSCAN_EXOME_SOMATIC_INDEL_OPTIONS) --vcf-sample-list $*.InDel.sample.txt > $@.unfiltered.vcf; \
		rm -f $*.InDel.sample.txt; \
	fi;
	-if [ ! -s $@.unfiltered.vcf ]; then cp $*.empty.vcf $@.unfiltered.vcf; fi;
	-$(JAVA11) -jar $(IGVTOOLS) index $@.unfiltered.vcf;
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $@.unfiltered.vcf -o $@ $(VARSCAN_EXOME_SOMATIC_INDEL_FILTERS);
	-if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -s $@ ]; then touch $@; fi;
	# Remove intermediate files
	-rm $@.unfiltered.vcf*
	# Remove IDX
	-rm $@.idx



# CONFIG/RELEASE
RELEASE_CMD := $(shell echo "\#\# VARSCAN EXOME_SOMATIC: identify variants and generate *.VarScan_EXOME_SOMATIC.vcf files with parameters SAMTOOLS mpileup and VARSCAN standard application 'VARSCAN_EXOME_SOMATIC_SNP_OPTIONS: $(VARSCAN_EXOME_SOMATIC_SNP_OPTIONS)' 'VARSCAN_EXOME_SOMATIC_INDEL_OPTIONS: $(VARSCAN_EXOME_SOMATIC_INDEL_OPTIONS)' and GATK VariantFiltration standard application 'VARSCAN_EXOME_SOMATIC_SNP_FILTERS: $(VARSCAN_EXOME_SOMATIC_SNP_FILTERS)' 'VARSCAN_EXOME_SOMATIC_INDEL_FILTERS: $(VARSCAN_EXOME_SOMATIC_INDEL_FILTERS)' " >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:VarScan_EXOME_SOMATIC:SAMTOOLS/VARSCAN - design for SOMATIC mutation on WES, VAF\>$(VARSCAN_EXOME_SOMATIC_VAF) ALT\>$(VARSCAN_EXOME_SOMATIC_ALT) DP\>$(VARSCAN_EXOME_SOMATIC_DP) p-value\<$(VARSCAN_EXOME_SOMATIC_PVAL):VARSCAN_EXOME_SOMATIC_SNP_OPTIONS='$(VARSCAN_EXOME_SOMATIC_SNP_OPTIONS)', VARSCAN_EXOME_SOMATIC_INDEL_OPTIONS='$(VARSCAN_EXOME_SOMATIC_INDEL_OPTIONS)', VARSCAN_EXOME_SOMATIC_SNP_FILTERS=$(VARSCAN_EXOME_SOMATIC_SNP_FILTERS), VARSCAN_EXOME_SOMATIC_INDEL_FILTERS=$(VARSCAN_EXOME_SOMATIC_INDEL_FILTERS)"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
