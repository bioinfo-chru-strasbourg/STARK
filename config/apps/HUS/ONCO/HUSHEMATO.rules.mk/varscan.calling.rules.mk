############################
# VarScan Calling Rules
# Release: 0.9.4.2
# Date: 04/05/2016
# Author: Antony Le Bechec
############################
MK_RELEASE="0.9.4.3b"
MK_DATE="17/05/2016"

# Release notes:
# 0.9.1-24/04/2015: Add 'varscanlowfreqlowcovsnp' calling
# 0.9.2-05/10/2015: Add 'VarScan' calling, default VarScan parameters
# 0.9.3-07/10/2015: Add filters on VarScan vcf
# 0.9.3.1-10/11/2015: Add pipeline VarScanHEMATO
# 0.9.4-10/11/2015: Add Merge SNP and InDel
# 0.9.4.1-03/05/2016: Bug correction, if empty VCF after calling. Add Release infos
# 0.9.4.2-04/05/2016: Bug correction
# 0.9.4.3-17/05/2016: Bug correction


# Parameters
# VARSCAN
VARSCAN?=$(NGSbin)/varscan.jar


####################
# SAMTOOLS mpileup
####################

MPILEUP_VARSCAN_HUSHEMATO_OPTIONS= -E -d 10000000 -L 10000000 -C 50 -Q 10 -q 1 --output-tags SP,DP,DP4,DV # -B  -d 100000 -L 100000 -C200 


%.bam.HUSHEMATO_mpileup: %.bam %.bam.bai %.genome
	#-$(SAMTOOLS) mpileup -f `cat $*.genome` $< $(MPILEUP_OPTIONS) -l $*.for_metrics.bed > $@
	-$(SAMTOOLS) mpileup -f `cat $*.genome` $< $(MPILEUP_VARSCAN_HUSHEMATO_OPTIONS) > $@



###################
# VarScan  HUSHEMATO
###################


VARSCAN_HUSHEMATO_BOTH_OPTIONS= --min-var-freq 0.03 --min-reads2 4 --min-coverage 50 --p-value 1e-1
VARSCAN_HUSHEMATO_SNP_OPTIONS= $(VARSCAN_HUSHEMATO_BOTH_OPTIONS) --min-avg-qual 30 
VARSCAN_HUSHEMATO_INDEL_OPTIONS= $(VARSCAN_HUSHEMATO_BOTH_OPTIONS) --min-avg-qual 10 
VARSCAN_HUSHEMATO_FILTERS=--genotypeFilterExpression 'GQ < 20.0' --genotypeFilterName 'LowGQ' --genotypeFilterExpression 'DP < 30' --genotypeFilterName 'LowCoverage' --genotypeFilterExpression 'DP < 10' --genotypeFilterName 'VeryLowCoverage'
VARSCAN_HUSHEMATO_SNP_FILTERS=$(VARSCAN_HUSHEMATO_FILTERS)
VARSCAN_HUSHEMATO_INDEL_FILTERS=$(VARSCAN_HUSHEMATO_FILTERS)
VARSCAN_HUSHEMATO_TIMEOUT=3600 # 1 = 1sec, 60=1min, 3600=1h


%.VarScan_HUSHEMATO.SNP.vcf: %.bam.HUSHEMATO_mpileup %.from_manifest.intervals %.empty.vcf %.genome #
	-if [ -s $< ]; then \
		echo "# GENERATION of '$@' start: "`date`; \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.SNP.sample.txt; \
		timeout $(VARSCAN_TIMEOUT) cat $< | $(JAVA) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_HUSHEMATO_SNP_OPTIONS) --vcf-sample-list $*.SNP.sample.txt > $@.unfiltered.vcf; \
		rm -f $*.SNP.sample.txt; \
		echo "# GENERATION of '$@' stop: "`date`; \
	fi;
	-if [ ! -s $@.unfiltered.vcf ]; then cp $*.empty.vcf $@.unfiltered.vcf; fi;
	echo "#NBVARIANT ($@.unfiltered.vcf after calling)"`grep -cv ^# $@.unfiltered.vcf`
	-$(JAVA) -jar $(IGVTOOLS) index $@.unfiltered.vcf;
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $@.unfiltered.vcf -o $@ $(VARSCAN_HUSHEMATO_SNP_FILTERS);
	echo "#NBVARIANT ($@ after filtering)"`grep -cv ^# $@`
	-if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -s $@ ]; then touch $@; fi;
	# Remove intermediate files
	-rm $@.unfiltered.vcf* 
	#-rm $*.empty.vcf
	# Remove IDX
	-rm $@.idx



%.VarScan_HUSHEMATO.InDel.vcf: %.bam.HUSHEMATO_mpileup %.from_manifest.intervals %.empty.vcf %.genome
	-if [ -s $< ]; then \
		echo "# GENERATION of '$@' start: "`date`; \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.InDel.sample.txt; \
		timeout $(VARSCAN_TIMEOUT) cat $< | $(JAVA) -jar $(VARSCAN) mpileup2indel --output-vcf 1 $(VARSCAN_HUSHEMATO_INDEL_OPTIONS) --vcf-sample-list $*.InDel.sample.txt > $@.unfiltered.vcf; \
		rm -f $*.InDel.sample.txt; \
		echo "# GENERATION of '$@' stop: "`date`; \
	fi;
	-if [ ! -s $@.unfiltered.vcf ]; then cp $*.empty.vcf $@.unfiltered.vcf; fi;
	echo "#NBVARIANT ($@.unfiltered.vcf after calling)"`grep -cv ^# $@.unfiltered.vcf`
	-$(JAVA) -jar $(IGVTOOLS) index $@.unfiltered.vcf;
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $@.unfiltered.vcf -o $@ $(VARSCAN_HUSHEMATO_INDEL_FILTERS);
	echo "#NBVARIANT ($@ after calling)"`grep -cv ^# $@`
	-if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -s $@ ]; then touch $@; fi;
	# Remove intermediate files
	-rm $@.unfiltered.vcf*
	#-rm $*.empty.vcf
	# Remove IDX
	-rm $@.idx





# CONFIG/RELEASE
RELEASE_CMD := $(shell echo "\#\# VARSCAN HUSHEMATO: identify variants and generate *.VarScan_HUSHEMATO.vcf files with parameters SAMTOOLS mpileup '$(MPILEUP_OPTIONS)' and VARSCAN HUSHEMATO application 'VARSCAN_HUSHEMATO_SNP_OPTIONS: $(VARSCAN_HUSHEMATO_SNP_OPTIONS)' 'VARSCAN_HUSHEMATO_SNP_OPTIONS: $(VARSCAN_HUSHEMATO_INDEL_OPTIONS)' and GATK VariantFiltration HUSHEMATO application 'VARSCAN_HUSHEMATO_SNP_FILTERS: $(VARSCAN_HUSHEMATO_SNP_FILTERS)' 'VARSCAN_HUSHEMATO_INDEL_FILTERS: $(VARSCAN_HUSHEMATO_INDEL_FILTERS)' " >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:VarScan_HUSHEMATO:SAMTOOLS/VARSCAN - design for HUSHEMATO mutation, VAF\>0.03 ALT\>4 DP\>50 p-value\<1e-1:MPILEUP_VARSCAN_HUSHEMATO_OPTIONS='$(MPILEUP_VARSCAN_HUSHEMATO_OPTIONS)', VARSCAN_HUSHEMATO_BOTH_OPTIONS=$(VARSCAN_HUSHEMATO_BOTH_OPTIONS)"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )





