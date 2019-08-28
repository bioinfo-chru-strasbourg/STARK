############################
# VarScan Calling Rules
# Author: Antony Le Bechec
############################
MK_RELEASE="0.9.4.5"
MK_DATE="27/08/2019"

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
# 0.9.4.5-27/08/2019: Change FATBAM to CAP


# Parameters
# VARSCAN
VARSCAN?=$(NGSbin)/varscan.jar
CAP?=/STARK/tools/cap/current/bin/CAP
CAP_SOFTCLIPTOQ0?=$(CAP)/CAP.SoftClipToQ0.pl


####################
# SAMTOOLS mpileup
####################

#MPILEUP_VARSCAN_HUSTUMSOL_OPTIONS= -E -d 10000000 -L 10000000 -C 50 -Q 10 -q 1 --output-tags SP,DP,DP4,DV # -B  -d 100000 -L 100000 -C200
MPILEUP_VARSCAN_HUSTUMSOL_OPTIONS= -E -d 10000000 -L 10000000 -C 50 -Q 10 -q 1 --output-tags SP,DP,DP4,DV,ADF,ADR,AD # -B  -d 100000 -L 100000 -C200


%.bam.HUSTUMSOL_mpileup: %.bam %.bam.bai %.genome
	#-$(SAMTOOLS) mpileup -f `cat $*.genome` $< $(MPILEUP_OPTIONS) -l $*.for_metrics.bed > $@
	$(SAMTOOLS) view $< -h | perl $(CAP_SOFTCLIPTOQ0) -v1 | $(SAMTOOLS) mpileup - -f `cat $*.genome` $(MPILEUP_VARSCAN_HUSTUMSOL_OPTIONS) > $@



#######################
# VarScan  HUSTUMSOL #
#######################


VARSCAN_HUSTUMSOL_BOTH_OPTIONS= --min-var-freq 0.01 --min-reads2 5 --min-coverage 30 --p-value 1e-1
VARSCAN_HUSTUMSOL_SNP_OPTIONS= $(VARSCAN_HUSTUMSOL_BOTH_OPTIONS) --min-avg-qual 30
VARSCAN_HUSTUMSOL_INDEL_OPTIONS= $(VARSCAN_HUSTUMSOL_BOTH_OPTIONS) --min-avg-qual 10
VARSCAN_HUSTUMSOL_FILTERS=--genotypeFilterExpression 'GQ < 200.0 && GQ >= 150' --genotypeFilterName 'LowGQ' --genotypeFilterExpression 'GQ < 150.0 && GQ >= 100' --genotypeFilterName 'VeryLowGQ'  --genotypeFilterExpression 'GQ < 100.0' --genotypeFilterName 'VeryVeryLowGQ'
VARSCAN_HUSTUMSOL_SNP_FILTERS=$(VARSCAN_HUSTUMSOL_FILTERS)
VARSCAN_HUSTUMSOL_INDEL_FILTERS=$(VARSCAN_HUSTUMSOL_FILTERS)
VARSCAN_HUSTUMSOL_TIMEOUT=3600 # 1 = 1sec, 60=1min, 3600=1h


%.VarScan_HUSTUMSOL.SNP.vcf: %.bam.HUSTUMSOL_mpileup %.from_manifest.intervals %.empty.vcf %.genome #
	-if [ -s $< ]; then \
		echo "# GENERATION of '$@' start: "`date`; \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.SNP.sample.txt; \
		timeout $(VARSCAN_TIMEOUT) cat $< | $(JAVA) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_HUSTUMSOL_SNP_OPTIONS) --vcf-sample-list $*.SNP.sample.txt > $@.unfiltered.vcf; \
		rm -f $*.SNP.sample.txt; \
		echo "# GENERATION of '$@' stop: "`date`; \
	fi;
	-if [ ! -s $@.unfiltered.vcf ]; then cp $*.empty.vcf $@.unfiltered.vcf; fi;
	echo "#NBVARIANT ($@.unfiltered.vcf after calling)"`grep -cv ^# $@.unfiltered.vcf`
	-$(JAVA) -jar $(IGVTOOLS) index $@.unfiltered.vcf;
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $@.unfiltered.vcf -o $@ $(VARSCAN_HUSTUMSOL_SNP_FILTERS);
	echo "#NBVARIANT ($@ after filtering)"`grep -cv ^# $@`
	-if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -s $@ ]; then touch $@; fi;
	# Remove intermediate files
	-rm $@.unfiltered.vcf*
	#-rm $*.empty.vcf
	# Remove IDX
	-rm $@.idx



%.VarScan_HUSTUMSOL.InDel.vcf: %.bam.HUSTUMSOL_mpileup %.from_manifest.intervals %.empty.vcf %.genome
	-if [ -s $< ]; then \
		echo "# GENERATION of '$@' start: "`date`; \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.InDel.sample.txt; \
		timeout $(VARSCAN_TIMEOUT) cat $< | $(JAVA) -jar $(VARSCAN) mpileup2indel --output-vcf 1 $(VARSCAN_HUSTUMSOL_INDEL_OPTIONS) --vcf-sample-list $*.InDel.sample.txt > $@.unfiltered.vcf; \
		rm -f $*.InDel.sample.txt; \
		echo "# GENERATION of '$@' stop: "`date`; \
	fi;
	-if [ ! -s $@.unfiltered.vcf ]; then cp $*.empty.vcf $@.unfiltered.vcf; fi;
	echo "#NBVARIANT ($@.unfiltered.vcf after calling)"`grep -cv ^# $@.unfiltered.vcf`
	-$(JAVA) -jar $(IGVTOOLS) index $@.unfiltered.vcf;
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $@.unfiltered.vcf -o $@ $(VARSCAN_HUSTUMSOL_INDEL_FILTERS);
	echo "#NBVARIANT ($@ after calling)"`grep -cv ^# $@`
	-if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -s $@ ]; then touch $@; fi;
	# Remove intermediate files
	-rm $@.unfiltered.vcf*
	#-rm $*.empty.vcf
	# Remove IDX
	-rm $@.idx


#######################
# VarScan CNV  HUSTUMSOL # NOT WORKING!!!!! # NEED normalized BAM !!!!
#######################

VARSCAN_CNV_NORMAL?=T_Horizon
VARSCAN_copyCaller?= --min-coverage=30 --amp-threshold=0.1 --del-threshold=0.1 --min-region-size=30

%.VarScanCNV_HUSTUMSOL.vcf: %.bam %.from_manifest.bed %.empty.vcf %.genome %.bams.list
	cp $*.empty.vcf $@
	echo $$(dirname $<)/../$(VARSCAN_CNV_NORMAL)/$(VARSCAN_CNV_NORMAL).bwamem.bam;
	#  | awk '{if($4 > 0) print $0}' | awk '{if($7 > 0) print $0}'  : TO AVOID a BUG
	-if [ "$(VARSCAN_CNV_NORMAL)" != "" ] && [ -s $$(dirname $<)/../$(VARSCAN_CNV_NORMAL)/$(VARSCAN_CNV_NORMAL).bwamem.bam  ]; then \
		$(SAMTOOLS) mpileup -f `cat $*.genome` $(MPILEUP_VARSCAN_OPTIONS) $$(dirname $<)/../$(VARSCAN_CNV_NORMAL)/$(VARSCAN_CNV_NORMAL).bwamem.bam $< | awk '{if($$4 > 0) print $$0}' | awk '{if($$7 > 0) print $$0}' | $(JAVA) -jar $(VARSCAN) copynumber $@.cnv.tmp --mpileup 1; \
		$(JAVA) -jar $(VARSCAN) copyCaller $@.cnv.tmp --output-file=$@.cnv $(VARSCAN_copyCaller) ; \
	else \
		echo "# No Normal BAM '$(VARSCAN_CNV_NORMAL)' for VarScan copynumber " ; \
	fi;
	# $(JAVA) -jar $(VARSCAN) copyCaller $@.cnv.tmp --output-file=$@.cnv --min-coverage=30 ;
	# --amp-threshold	Lower bound for log ratio to call amplification [0.25]
	#--del-threshold	Upper bound for log ratio to call deletion (provide as positive number) [0.25]
	#--min-region-size	Minimum size (in bases) for a region to be counted [10]
	#--recenter-up	Recenter data around an adjusted baseline > 0 [0]
	#--recenter-down	Recenter data around an adjusted baseline < 0 [0]
	# cleaning
	-rm $@.cnv.tmp




# CONFIG/RELEASE

RELEASE_CMD := $(shell echo "\#\# VARSCAN HUSTUMSOL: identify variants and generate *.VarScan_HUSTUMSOL.vcf files with parameters SAMTOOLS mpileup '$(MPILEUP_OPTIONS)' and VARSCAN HUSTUMSOL application 'VARSCAN_HUSTUMSOL_SNP_OPTIONS: $(VARSCAN_HUSTUMSOL_SNP_OPTIONS)' 'VARSCAN_HUSTUMSOL_SNP_OPTIONS: $(VARSCAN_HUSTUMSOL_INDEL_OPTIONS)' and GATK VariantFiltration HUSTUMSOL application 'VARSCAN_HUSTUMSOL_SNP_FILTERS: $(VARSCAN_HUSTUMSOL_SNP_FILTERS)' 'VARSCAN_HUSTUMSOL_INDEL_FILTERS: $(VARSCAN_HUSTUMSOL_INDEL_FILTERS)' " >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "CALLER:VarScan_HUSTUMSOL:SAMTOOLS/VARSCAN - design for HUSTUMSOL mutation, VAF\>0.01 ALT\>5 DP\>30 p-value\<1e-1:MPILEUP_VARSCAN_HUSTUMSOL_OPTIONS='$(MPILEUP_VARSCAN_HUSTUMSOL_OPTIONS)', VARSCAN_HUSTUMSOL_BOTH_OPTIONS=$(VARSCAN_HUSTUMSOL_BOTH_OPTIONS)"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
