############################
# VarScan Calling Rules
# Release: 0.9.4.6
# Date: 25/05/2021
# Author: Antony Le Bechec
############################
MK_RELEASE="0.9.4.6"
MK_DATE="25/05/2021"

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




####################
# SAMTOOLS mpileup
####################

MPILEUP_VARSCAN_OPTIONS= -E -d 10000000 -L 10000000 -C 50 -Q 10 -q 1 --output-tags SP,DP,DP4,DV,ADF,ADR,AD # -B  -d 100000 -L 100000 -C200


%.bam.mpileup: %.bam %.bam.bai %.genome %.design.bed
	#-$(SAMTOOLS) mpileup -f `cat $*.genome` $< $(MPILEUP_OPTIONS) -l $*.for_metrics.bed > $@
	#-$(SAMTOOLS) mpileup -f `cat $*.genome` $< $(MPILEUP_VARSCAN_OPTIONS) > $@
	#$(SAMTOOLS) view $< -h | perl $(CAP_SOFTCLIPTOQ0) -v1 | $(SAMTOOLS) mpileup - -f $$(cat $*.genome) $(MPILEUP_VARSCAN_OPTIONS) > $@
	if (( $$(grep -v "^#" -c $*.design.bed) )); then \
		$(SAMTOOLS) view $< -L $*.design.bed -h | perl $(CAP_SOFTCLIPTOQ0) -v1 | $(SAMTOOLS) mpileup - -f $$(cat $*.genome) -l $*.design.bed $(MPILEUP_VARSCAN_OPTIONS) > $@; \
	else \
		$(SAMTOOLS) view $< -h | perl $(CAP_SOFTCLIPTOQ0) -v1 | $(SAMTOOLS) mpileup - -f $$(cat $*.genome) $(MPILEUP_VARSCAN_OPTIONS) > $@; \
	fi;



####################
# VARSCAN Standard
####################

VARSCAN_BOTH_OPTIONS= --min-var-freq 0.01 --min-reads2 2 --min-coverage 30 --p-value 1e-1
VARSCAN_SNP_OPTIONS= $(VARSCAN_BOTH_OPTIONS) --min-avg-qual 30
VARSCAN_INDEL_OPTIONS= $(VARSCAN_BOTH_OPTIONS) --min-avg-qual 10
VARSCAN_FILTERS=--genotypeFilterExpression 'GQ < 20.0' --genotypeFilterName 'LowGQ'
#VARSCAN_FILTERS=--genotypeFilterExpression 'GQ < 200.0 && GQ >= 150' --genotypeFilterName 'LowGQ' --genotypeFilterExpression 'GQ < 150.0 && GQ >= 100' --genotypeFilterName 'VeryLowGQ'  --genotypeFilterExpression 'GQ < 100.0' --genotypeFilterName 'VeryVeryLowGQ'
VARSCAN_SNP_FILTERS=$(VARSCAN_FILTERS)
VARSCAN_INDEL_FILTERS=$(VARSCAN_FILTERS)
VARSCAN_TIMEOUT=3600 # 1 = 1sec, 60=1min, 3600=1h

%.VarScan.SNP$(POST_CALLING).vcf: %.bam.mpileup %.empty.vcf %.genome #
	-if [ -s $< ]; then \
		echo "# GENERATION of '$@' start: "`date`; \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.SNP.sample.txt; \
		#timeout $(VARSCAN_TIMEOUT) cat $< | $(JAVA) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_SNP_OPTIONS) --vcf-sample-list $*.SNP.sample.txt > $@.unfiltered.vcf; \
		cat $< | $(JAVA) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_SNP_OPTIONS) --vcf-sample-list $*.SNP.sample.txt > $@.unfiltered.vcf; \
		rm -f $*.SNP.sample.txt; \
		echo "# GENERATION of '$@' stop: "`date`; \
	fi;
	-if [ ! -s $@.unfiltered.vcf ]; then cp $*.empty.vcf $@.unfiltered.vcf; fi;
	echo "#NBVARIANT ($@.unfiltered.vcf after calling)"`grep -cv ^# $@.unfiltered.vcf`
	-$(JAVA) -jar $(IGVTOOLS) index $@.unfiltered.vcf;
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $@.unfiltered.vcf -o $@ $(VARSCAN_SNP_FILTERS);
	echo "#NBVARIANT ($@ after calling)"`grep -cv ^# $@`
	-if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -s $@ ]; then touch $@; fi;
	# Remove intermediate files
	-rm $@.unfiltered.vcf*
	#-rm $*.empty.vcf
	# Remove IDX
	-rm $@.idx


%.VarScan.InDel$(POST_CALLING).vcf: %.bam.mpileup %.empty.vcf %.genome
	-if [ -s $< ]; then \
		echo "# GENERATION of '$@' start: "`date`; \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.InDel.sample.txt; \
		#timeout $(VARSCAN_TIMEOUT) cat $< | $(JAVA) -jar $(VARSCAN) mpileup2indel --output-vcf 1 $(VARSCAN_INDEL_OPTIONS) --vcf-sample-list $*.InDel.sample.txt > $@.unfiltered.vcf; \
		cat $< | $(JAVA) -jar $(VARSCAN) mpileup2indel --output-vcf 1 $(VARSCAN_INDEL_OPTIONS) --vcf-sample-list $*.InDel.sample.txt > $@.unfiltered.vcf; \
		rm -f $*.InDel.sample.txt; \
		echo "# GENERATION of '$@' stop: "`date`; \
	fi;
	-if [ ! -s $@.unfiltered.vcf ]; then cp $*.empty.vcf $@.unfiltered.vcf; fi;
	echo "#NBVARIANT ($@.unfiltered.vcf after calling)"`grep -cv ^# $@.unfiltered.vcf`
	-$(JAVA) -jar $(IGVTOOLS) index $@.unfiltered.vcf;
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $@.unfiltered.vcf -o $@ $(VARSCAN_INDEL_FILTERS);
	echo "#NBVARIANT ($@ after calling)"`grep -cv ^# $@`
	-if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -s $@ ]; then touch $@; fi;
	# Remove intermediate files
	-rm $@.unfiltered.vcf*
	#-rm $*.empty.vcf
	# Remove IDX
	-rm $@.idx


#######################
# VarScan  SOMATIC #
#######################

VARSCAN_SOMATIC_VAF=0.01
VARSCAN_SOMATIC_ALT=4
VARSCAN_SOMATIC_DP=50
VARSCAN_SOMATIC_PVAL=1e-1

VARSCAN_SOMATIC_BOTH_OPTIONS= --min-var-freq $(VARSCAN_SOMATIC_VAF) --min-reads2 $(VARSCAN_SOMATIC_ALT) --min-coverage $(VARSCAN_SOMATIC_DP) --p-value $(VARSCAN_SOMATIC_PVAL)
VARSCAN_SOMATIC_SNP_OPTIONS= $(VARSCAN_SOMATIC_BOTH_OPTIONS) --min-avg-qual 30
VARSCAN_SOMATIC_INDEL_OPTIONS= $(VARSCAN_SOMATIC_BOTH_OPTIONS) --min-avg-qual 10
VARSCAN_SOMATIC_FILTERS=--genotypeFilterExpression 'GQ < 20.0' --genotypeFilterName 'LowGQ'
#--genotypeFilterExpression "GQ < 200.0 && GQ >= 150" --genotypeFilterName "LowGQ" --genotypeFilterExpression "GQ < 150.0 && GQ >= 100" --genotypeFilterName "VeryLowGQ"  --genotypeFilterExpression "GQ < 100.0" --genotypeFilterName "VeryVeryLowGQ" --filterExpression "DP < 30" --filterName "LowCoverage" --filterExpression "DP < 10" --filterName "VeryLowCoverage"
VARSCAN_SOMATIC_SNP_FILTERS=$(VARSCAN_SOMATIC_FILTERS)
VARSCAN_SOMATIC_INDEL_FILTERS=$(VARSCAN_SOMATIC_FILTERS)
VARSCAN_SOMATIC_TIMEOUT=3600 # 1 = 1sec, 60=1min, 3600=1h


%.VarScan_SOMATIC$(POST_CALLING).SNP.vcf: %.bam.mpileup %.empty.vcf %.genome #
	-if [ -s $< ]; then \
		echo "# GENERATION of '$@' start: "`date`; \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.SNP.sample.txt; \
		#timeout $(VARSCAN_TIMEOUT) cat $< | $(JAVA) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_SOMATIC_SNP_OPTIONS) --vcf-sample-list $*.SNP.sample.txt > $@.unfiltered.vcf; \
		cat $< | $(JAVA) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_SOMATIC_SNP_OPTIONS) --vcf-sample-list $*.SNP.sample.txt > $@.unfiltered.vcf; \
		rm -f $*.SNP.sample.txt; \
		echo "# GENERATION of '$@' stop: "`date`; \
	fi;
	-if [ ! -s $@.unfiltered.vcf ]; then cp $*.empty.vcf $@.unfiltered.vcf; fi;
	echo "#NBVARIANT ($@.unfiltered.vcf after calling)"`grep -cv ^# $@.unfiltered.vcf`
	-$(JAVA) -jar $(IGVTOOLS) index $@.unfiltered.vcf;
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $@.unfiltered.vcf -o $@ $(VARSCAN_SOMATIC_SNP_FILTERS);
	echo "#NBVARIANT ($@ after filtering)"`grep -cv ^# $@`
	-if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -s $@ ]; then touch $@; fi;
	# Remove intermediate files
	-rm $@.unfiltered.vcf*
	#-rm $*.empty.vcf
	# Remove IDX
	-rm $@.idx


%.VarScan_SOMATIC$(POST_CALLING).InDel.vcf: %.bam.mpileup %.empty.vcf %.genome
	-if [ -s $< ]; then \
		echo "# GENERATION of '$@' start: "`date`; \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.InDel.sample.txt; \
		#timeout $(VARSCAN_TIMEOUT) cat $< | $(JAVA) -jar $(VARSCAN) mpileup2indel --output-vcf 1 $(VARSCAN_SOMATIC_INDEL_OPTIONS) --vcf-sample-list $*.InDel.sample.txt > $@.unfiltered.vcf; \
		cat $< | $(JAVA) -jar $(VARSCAN) mpileup2indel --output-vcf 1 $(VARSCAN_SOMATIC_INDEL_OPTIONS) --vcf-sample-list $*.InDel.sample.txt > $@.unfiltered.vcf; \
		rm -f $*.InDel.sample.txt; \
		echo "# GENERATION of '$@' stop: "`date`; \
	fi;
	-if [ ! -s $@.unfiltered.vcf ]; then cp $*.empty.vcf $@.unfiltered.vcf; fi;
	echo "#NBVARIANT ($@.unfiltered.vcf after calling)"`grep -cv ^# $@.unfiltered.vcf`
	-$(JAVA) -jar $(IGVTOOLS) index $@.unfiltered.vcf;
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $@.unfiltered.vcf -o $@ $(VARSCAN_SOMATIC_INDEL_FILTERS);
	echo "#NBVARIANT ($@ after calling)"`grep -cv ^# $@`
	-if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -s $@ ]; then touch $@; fi;
	# Remove intermediate files
	-rm $@.unfiltered.vcf*
	#-rm $*.empty.vcf
	# Remove IDX
	-rm $@.idx


#######################
# VarScan EXOME_SOMATIC #
#######################

VARSCAN_EXOME_SOMATIC_VAF=0.05
VARSCAN_EXOME_SOMATIC_ALT=5
VARSCAN_EXOME_SOMATIC_DP=50
VARSCAN_EXOME_SOMATIC_PVAL=1e-1

VARSCAN_EXOME_SOMATIC_BOTH_OPTIONS= --min-var-freq $(VARSCAN_EXOME_SOMATIC_VAF) --min-reads2 $(VARSCAN_EXOME_SOMATIC_ALT) --min-coverage $(VARSCAN_EXOME_SOMATIC_DP) --p-value $(VARSCAN_EXOME_SOMATIC_PVAL)
VARSCAN_EXOME_SOMATIC_SNP_OPTIONS= $(VARSCAN_EXOME_SOMATIC_BOTH_OPTIONS) --min-avg-qual 30
VARSCAN_EXOME_SOMATIC_INDEL_OPTIONS= $(VARSCAN_EXOME_SOMATIC_BOTH_OPTIONS) --min-avg-qual 10
VARSCAN_EXOME_SOMATIC_FILTERS=--genotypeFilterExpression 'GQ < 20.0' --genotypeFilterName 'LowGQ'
#--genotypeFilterExpression "GQ < 200.0 && GQ >= 150" --genotypeFilterName "LowGQ" --genotypeFilterExpression "GQ < 150.0 && GQ >= 100" --genotypeFilterName "VeryLowGQ"  --genotypeFilterExpression "GQ < 100.0" --genotypeFilterName "VeryVeryLowGQ" --filterExpression "DP < 30" --filterName "LowCoverage" --filterExpression "DP < 10" --filterName "VeryLowCoverage"
VARSCAN_EXOME_SOMATIC_SNP_FILTERS=$(VARSCAN_EXOME_SOMATIC_FILTERS)
VARSCAN_EXOME_SOMATIC_INDEL_FILTERS=$(VARSCAN_EXOME_SOMATIC_FILTERS)
VARSCAN_EXOME_SOMATIC_TIMEOUT=36000 # 1 = 1sec, 60=1min, 3600=1h


%.VarScan_EXOME_SOMATIC.SNP$(POST_CALLING).vcf: %.bam.mpileup %.empty.vcf %.genome #
	-if [ -s $< ]; then \
		echo "# GENERATION of '$@' start: "`date`; \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.SNP.sample.txt; \
		#timeout $(VARSCAN_TIMEOUT) cat $< | $(JAVA) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_EXOME_SOMATIC_SNP_OPTIONS) --vcf-sample-list $*.SNP.sample.txt > $@.unfiltered.vcf; \
		cat $< | $(JAVA) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_EXOME_SOMATIC_SNP_OPTIONS) --vcf-sample-list $*.SNP.sample.txt > $@.unfiltered.vcf; \
		rm -f $*.SNP.sample.txt; \
		echo "# GENERATION of '$@' stop: "`date`; \
	fi;
	-if [ ! -s $@.unfiltered.vcf ]; then cp $*.empty.vcf $@.unfiltered.vcf; fi;
	echo "#NBVARIANT ($@.unfiltered.vcf after calling)"`grep -cv ^# $@.unfiltered.vcf`
	-$(JAVA) -jar $(IGVTOOLS) index $@.unfiltered.vcf;
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $@.unfiltered.vcf -o $@ $(VARSCAN_EXOME_SOMATIC_SNP_FILTERS);
	echo "#NBVARIANT ($@ after filtering)"`grep -cv ^# $@`
	-if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -s $@ ]; then touch $@; fi;
	# Remove intermediate files
	-rm $@.unfiltered.vcf*
	#-rm $*.empty.vcf
	# Remove IDX
	-rm $@.idx


%.VarScan_EXOME_SOMATIC.InDel$(POST_CALLING).vcf: %.bam.mpileup %.empty.vcf %.genome
	-if [ -s $< ]; then \
		echo "# GENERATION of '$@' start: "`date`; \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.InDel.sample.txt; \
		#timeout $(VARSCAN_TIMEOUT) cat $< | $(JAVA) -jar $(VARSCAN) mpileup2indel --output-vcf 1 $(VARSCAN_EXOME_SOMATIC_INDEL_OPTIONS) --vcf-sample-list $*.InDel.sample.txt > $@.unfiltered.vcf; \
		cat $< | $(JAVA) -jar $(VARSCAN) mpileup2indel --output-vcf 1 $(VARSCAN_EXOME_SOMATIC_INDEL_OPTIONS) --vcf-sample-list $*.InDel.sample.txt > $@.unfiltered.vcf; \
		rm -f $*.InDel.sample.txt; \
		echo "# GENERATION of '$@' stop: "`date`; \
	fi;
	-if [ ! -s $@.unfiltered.vcf ]; then cp $*.empty.vcf $@.unfiltered.vcf; fi;
	echo "#NBVARIANT ($@.unfiltered.vcf after calling)"`grep -cv ^# $@.unfiltered.vcf`
	-$(JAVA) -jar $(IGVTOOLS) index $@.unfiltered.vcf;
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $@.unfiltered.vcf -o $@ $(VARSCAN_EXOME_SOMATIC_INDEL_FILTERS);
	echo "#NBVARIANT ($@ after calling)"`grep -cv ^# $@`
	-if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -s $@ ]; then touch $@; fi;
	# Remove intermediate files
	-rm $@.unfiltered.vcf*
	#-rm $*.empty.vcf
	# Remove IDX
	-rm $@.idx




#######################
# VarScan  HEMATOLOGY #
#######################

VARSCAN_HEMATOLOGY_VAF=0.03
VARSCAN_HEMATOLOGY_ALT=4
VARSCAN_HEMATOLOGY_DP=50
VARSCAN_HEMATOLOGY_PVAL=1e-1

VARSCAN_HEMATOLOGY_BOTH_OPTIONS= --min-var-freq $(VARSCAN_HEMATOLOGY_VAF) --min-reads2 $(VARSCAN_HEMATOLOGY_ALT) --min-coverage $(VARSCAN_HEMATOLOGY_DP) --p-value $(VARSCAN_HEMATOLOGY_PVAL)
VARSCAN_HEMATOLOGY_SNP_OPTIONS= $(VARSCAN_HEMATOLOGY_BOTH_OPTIONS) --min-avg-qual 30
VARSCAN_HEMATOLOGY_INDEL_OPTIONS= $(VARSCAN_HEMATOLOGY_BOTH_OPTIONS) --min-avg-qual 10
VARSCAN_HEMATOLOGY_FILTERS=--genotypeFilterExpression 'GQ < 20.0' --genotypeFilterName 'LowGQ'
#--genotypeFilterExpression "GQ < 200.0 && GQ >= 150" --genotypeFilterName "LowGQ" --genotypeFilterExpression "GQ < 150.0 && GQ >= 100" --genotypeFilterName "VeryLowGQ"  --genotypeFilterExpression "GQ < 100.0" --genotypeFilterName "VeryVeryLowGQ" --filterExpression "DP < 30" --filterName "LowCoverage" --filterExpression "DP < 10" --filterName "VeryLowCoverage"
VARSCAN_HEMATOLOGY_SNP_FILTERS=$(VARSCAN_HEMATOLOGY_FILTERS)
VARSCAN_HEMATOLOGY_INDEL_FILTERS=$(VARSCAN_HEMATOLOGY_FILTERS)
VARSCAN_HEMATOLOGY_TIMEOUT=3600 # 1 = 1sec, 60=1min, 3600=1h


%.VarScan_HEMATOLOGY.SNP$(POST_CALLING).vcf: %.bam.mpileup %.empty.vcf %.genome
	-if [ -s $< ]; then \
		echo "# GENERATION of '$@' start: "`date`; \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.SNP.sample.txt; \
		#timeout $(VARSCAN_TIMEOUT) cat $< | $(JAVA) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_HEMATOLOGY_SNP_OPTIONS) --vcf-sample-list $*.SNP.sample.txt > $@.unfiltered.vcf; \
		cat $< | $(JAVA) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_HEMATOLOGY_SNP_OPTIONS) --vcf-sample-list $*.SNP.sample.txt > $@.unfiltered.vcf; \
		rm -f $*.SNP.sample.txt; \
		echo "# GENERATION of '$@' stop: "`date`; \
	fi;
	-if [ ! -s $@.unfiltered.vcf ]; then cp $*.empty.vcf $@.unfiltered.vcf; fi;
	echo "#NBVARIANT ($@.unfiltered.vcf after calling)"`grep -cv ^# $@.unfiltered.vcf`
	-$(JAVA) -jar $(IGVTOOLS) index $@.unfiltered.vcf;
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $@.unfiltered.vcf -o $@ $(VARSCAN_HEMATOLOGY_SNP_FILTERS);
	echo "#NBVARIANT ($@ after filtering)"`grep -cv ^# $@`
	-if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -s $@ ]; then touch $@; fi;
	# Remove intermediate files
	-rm $@.unfiltered.vcf*
	#-rm $*.empty.vcf
	# Remove IDX
	-rm $@.idx


%.VarScan_HEMATOLOGY.InDel$(POST_CALLING).vcf: %.bam.mpileup %.empty.vcf %.genome
	-if [ -s $< ]; then \
		echo "# GENERATION of '$@' start: "`date`; \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.InDel.sample.txt; \
		#timeout $(VARSCAN_TIMEOUT) cat $< | $(JAVA) -jar $(VARSCAN) mpileup2indel --output-vcf 1 $(VARSCAN_HEMATOLOGY_INDEL_OPTIONS) --vcf-sample-list $*.InDel.sample.txt > $@.unfiltered.vcf; \
		cat $< | $(JAVA) -jar $(VARSCAN) mpileup2indel --output-vcf 1 $(VARSCAN_HEMATOLOGY_INDEL_OPTIONS) --vcf-sample-list $*.InDel.sample.txt > $@.unfiltered.vcf; \
		rm -f $*.InDel.sample.txt; \
		echo "# GENERATION of '$@' stop: "`date`; \
	fi;
	-if [ ! -s $@.unfiltered.vcf ]; then cp $*.empty.vcf $@.unfiltered.vcf; fi;
	echo "#NBVARIANT ($@.unfiltered.vcf after calling)"`grep -cv ^# $@.unfiltered.vcf`
	-$(JAVA) -jar $(IGVTOOLS) index $@.unfiltered.vcf;
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $@.unfiltered.vcf -o $@ $(VARSCAN_HEMATOLOGY_INDEL_FILTERS);
	echo "#NBVARIANT ($@ after calling)"`grep -cv ^# $@`
	-if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -s $@ ]; then touch $@; fi;
	# Remove intermediate files
	-rm $@.unfiltered.vcf*
	#-rm $*.empty.vcf
	# Remove IDX
	-rm $@.idx


#######################
# VarScan  SOLIDTUMOR #
#######################

VARSCAN_SOLIDTUMOR_VAF=0.01
VARSCAN_SOLIDTUMOR_ALT=5
VARSCAN_SOLIDTUMOR_DP=30
VARSCAN_SOLIDTUMOR_PVAL=1e-1

VARSCAN_SOLIDTUMOR_BOTH_OPTIONS= --min-var-freq $(VARSCAN_SOLIDTUMOR_VAF) --min-reads2 $(VARSCAN_SOLIDTUMOR_ALT) --min-coverage $(VARSCAN_SOLIDTUMOR_DP) --p-value $(VARSCAN_SOLIDTUMOR_PVAL)
VARSCAN_SOLIDTUMOR_SNP_OPTIONS= $(VARSCAN_SOLIDTUMOR_BOTH_OPTIONS) --min-avg-qual 30
VARSCAN_SOLIDTUMOR_INDEL_OPTIONS= $(VARSCAN_SOLIDTUMOR_BOTH_OPTIONS) --min-avg-qual 10
VARSCAN_SOLIDTUMOR_FILTERS=--genotypeFilterExpression 'GQ < 20.0' --genotypeFilterName 'LowGQ'
#VARSCAN_SOLIDTUMOR_FILTERS=--genotypeFilterExpression 'GQ < 200.0 && GQ >= 150' --genotypeFilterName 'LowGQ' --genotypeFilterExpression 'GQ < 150.0 && GQ >= 100' --genotypeFilterName 'VeryLowGQ'  --genotypeFilterExpression 'GQ < 100.0' --genotypeFilterName 'VeryVeryLowGQ'
VARSCAN_SOLIDTUMOR_SNP_FILTERS=$(VARSCAN_SOLIDTUMOR_FILTERS)
VARSCAN_SOLIDTUMOR_INDEL_FILTERS=$(VARSCAN_SOLIDTUMOR_FILTERS)
VARSCAN_SOLIDTUMOR_TIMEOUT=3600 # 1 = 1sec, 60=1min, 3600=1h


%.VarScan_SOLIDTUMOR.SNP$(POST_CALLING).vcf: %.bam.mpileup %.empty.vcf %.genome
	-if [ -s $< ]; then \
		echo "# GENERATION of '$@' start: "`date`; \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.SNP.sample.txt; \
		#timeout $(VARSCAN_TIMEOUT) cat $< | $(JAVA) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_SOLIDTUMOR_SNP_OPTIONS) --vcf-sample-list $*.SNP.sample.txt > $@.unfiltered.vcf; \
		cat $< | $(JAVA) -jar $(VARSCAN) mpileup2snp --output-vcf 1 $(VARSCAN_SOLIDTUMOR_SNP_OPTIONS) --vcf-sample-list $*.SNP.sample.txt > $@.unfiltered.vcf; \
		rm -f $*.SNP.sample.txt; \
		echo "# GENERATION of '$@' stop: "`date`; \
	fi;
	-if [ ! -s $@.unfiltered.vcf ]; then cp $*.empty.vcf $@.unfiltered.vcf; fi;
	echo "#NBVARIANT ($@.unfiltered.vcf after calling)"`grep -cv ^# $@.unfiltered.vcf`
	-$(JAVA) -jar $(IGVTOOLS) index $@.unfiltered.vcf;
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $@.unfiltered.vcf -o $@ $(VARSCAN_SOLIDTUMOR_SNP_FILTERS);
	echo "#NBVARIANT ($@ after filtering)"`grep -cv ^# $@`
	-if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -s $@ ]; then touch $@; fi;
	# Remove intermediate files
	-rm $@.unfiltered.vcf*
	#-rm $*.empty.vcf
	# Remove IDX
	-rm $@.idx


%.VarScan_SOLIDTUMOR.InDel$(POST_CALLING).vcf: %.bam.mpileup %.empty.vcf %.genome
	-if [ -s $< ]; then \
		echo "# GENERATION of '$@' start: "`date`; \
		sample=`basename $* | cut -d"." -f1`; \
		echo $$sample > $*.InDel.sample.txt; \
		#timeout $(VARSCAN_TIMEOUT) cat $< | $(JAVA) -jar $(VARSCAN) mpileup2indel --output-vcf 1 $(VARSCAN_SOLIDTUMOR_INDEL_OPTIONS) --vcf-sample-list $*.InDel.sample.txt > $@.unfiltered.vcf; \
		cat $< | $(JAVA) -jar $(VARSCAN) mpileup2indel --output-vcf 1 $(VARSCAN_SOLIDTUMOR_INDEL_OPTIONS) --vcf-sample-list $*.InDel.sample.txt > $@.unfiltered.vcf; \
		rm -f $*.InDel.sample.txt; \
		echo "# GENERATION of '$@' stop: "`date`; \
	fi;
	-if [ ! -s $@.unfiltered.vcf ]; then cp $*.empty.vcf $@.unfiltered.vcf; fi;
	echo "#NBVARIANT ($@.unfiltered.vcf after calling)"`grep -cv ^# $@.unfiltered.vcf`
	-$(JAVA) -jar $(IGVTOOLS) index $@.unfiltered.vcf;
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $@.unfiltered.vcf -o $@ $(VARSCAN_SOLIDTUMOR_INDEL_FILTERS);
	echo "#NBVARIANT ($@ after calling)"`grep -cv ^# $@`
	-if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -s $@ ]; then touch $@; fi;
	# Remove intermediate files
	-rm $@.unfiltered.vcf*
	#-rm $*.empty.vcf
	# Remove IDX
	-rm $@.idx





# CONFIG/RELEASE
#RELEASE_CMD := $(shell echo "\#\# VARSCAN: identify variants and generate *.varscan*.vcf files with parameters SAMTOOLS mpileup '$(MPILEUP_OPTIONS)'  and VARSCAN '$(VARSCAN_OPTIONS)' " >> $(RELEASE_INFOS) )
RELEASE_COMMENT := "\#\# CALLING VARSCAN '$(MK_RELEASE)': VARSCAN tool identify variants from mpileup file, generated from aligned BAM xith parameters: MPILEUP_OPTIONS='$(MPILEUP_OPTIONS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_CMD := $(shell echo "\#\# VARSCAN: identify variants and generate *.VarScan.vcf files with parameters SAMTOOLS mpileup '$(MPILEUP_OPTIONS)' and VARSCAN standard application 'VARSCAN_SNP_OPTIONS: $(VARSCAN_SNP_OPTIONS)' 'VARSCAN_SNP_OPTIONS: $(VARSCAN_INDEL_OPTIONS)' and GATK VariantFiltration standard application 'VARSCAN_SNP_FILTERS: $(VARSCAN_SNP_FILTERS)' 'VARSCAN_INDEL_FILTERS: $(VARSCAN_INDEL_FILTERS)' " >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\# VARSCAN SOMATIC: identify variants and generate *.VarScan_SOMATIC.vcf files with parameters SAMTOOLS mpileup '$(MPILEUP_OPTIONS)' and VARSCAN standard application 'VARSCAN_SNP_OPTIONS: $(VARSCAN_SOMATIC_SNP_OPTIONS)' 'VARSCAN_SNP_OPTIONS: $(VARSCAN_SOMATIC_INDEL_OPTIONS)' and GATK VariantFiltration standard application 'VARSCAN_SNP_FILTERS: $(VARSCAN_SOMATIC_SNP_FILTERS)' 'VARSCAN_INDEL_FILTERS: $(VARSCAN_SOMATIC_INDEL_FILTERS)' " >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\# VARSCAN EXOME_SOMATIC: identify variants and generate *.VarScan_EXOME_SOMATIC.vcf files with parameters SAMTOOLS mpileup '$(MPILEUP_OPTIONS)' and VARSCAN standard application 'VARSCAN_SNP_OPTIONS: $(VARSCAN_EXOME_SOMATIC_SNP_OPTIONS)' 'VARSCAN_SNP_OPTIONS: $(VARSCAN_EXOME_SOMATIC_INDEL_OPTIONS)' and GATK VariantFiltration standard application 'VARSCAN_SNP_FILTERS: $(VARSCAN_EXOME_SOMATIC_SNP_FILTERS)' 'VARSCAN_INDEL_FILTERS: $(VARSCAN_EXOME_SOMATIC_INDEL_FILTERS)' " >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\# VARSCAN HEMATOLOGY: identify variants and generate *.VarScan_HEMATOLOGY.vcf files with parameters SAMTOOLS mpileup '$(MPILEUP_OPTIONS)' and VARSCAN HEMATOLOGY application 'VARSCAN_HEMATOLOGY_SNP_OPTIONS: $(VARSCAN_HEMATOLOGY_SNP_OPTIONS)' 'VARSCAN_HEMATOLOGY_SNP_OPTIONS: $(VARSCAN_HEMATOLOGY_INDEL_OPTIONS)' and GATK VariantFiltration HEMATOLOGY application 'VARSCAN_HEMATOLOGY_SNP_FILTERS: $(VARSCAN_HEMATOLOGY_SNP_FILTERS)' 'VARSCAN_HEMATOLOGY_INDEL_FILTERS: $(VARSCAN_HEMATOLOGY_INDEL_FILTERS)' " >> $(RELEASE_INFOS) )
RELEASE_CMD := $(shell echo "\#\# VARSCAN SOLIDTUMOR: identify variants and generate *.VarScan_SOLIDTUMOR.vcf files with parameters SAMTOOLS mpileup '$(MPILEUP_OPTIONS)' and VARSCAN SOLIDTUMOR application 'VARSCAN_SOLIDTUMOR_SNP_OPTIONS: $(VARSCAN_SOLIDTUMOR_SNP_OPTIONS)' 'VARSCAN_SOLIDTUMOR_SNP_OPTIONS: $(VARSCAN_SOLIDTUMOR_INDEL_OPTIONS)' and GATK VariantFiltration SOLIDTUMOR application 'VARSCAN_SOLIDTUMOR_SNP_FILTERS: $(VARSCAN_SOLIDTUMOR_SNP_FILTERS)' 'VARSCAN_SOLIDTUMOR_INDEL_FILTERS: $(VARSCAN_SOLIDTUMOR_INDEL_FILTERS)' " >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "CALLER:VarScan:SAMTOOLS/VARSCAN - by default:MPILEUP_OPTIONS='$(MPILEUP_OPTIONS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:VarScan_SOMATIC:SAMTOOLS/VARSCAN - design for SOMATIC mutation, VAF\>$(VARSCAN_SOMATIC_VAF) ALT\>$(VARSCAN_SOMATIC_ALT) DP\>$(VARSCAN_SOMATIC_DP) p-value\<$(VARSCAN_SOMATIC_PVAL):MPILEUP_OPTIONS='$(MPILEUP_OPTIONS)', VARSCAN_SOMATIC_BOTH_OPTIONS=$(VARSCAN_SOMATIC_BOTH_OPTIONS)"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:VarScan_EXOME_SOMATIC:SAMTOOLS/VARSCAN - design for SOMATIC mutation on WES, VAF\>$(VARSCAN_EXOME_SOMATIC_VAF) ALT\>$(VARSCAN_EXOME_SOMATIC_ALT) DP\>$(VARSCAN_EXOME_SOMATIC_DP) p-value\<$(VARSCAN_EXOME_SOMATIC_PVAL):MPILEUP_OPTIONS='$(MPILEUP_OPTIONS)', VARSCAN_EXOME_SOMATIC_BOTH_OPTIONS=$(VARSCAN_EXOME_SOMATIC_BOTH_OPTIONS)"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:VarScan_SOLIDTUMOR:SAMTOOLS/VARSCAN - design for SOLIDTUMOR mutation, VAF\>$(VARSCAN_SOLIDTUMOR_VAF) ALT\>$(VARSCAN_SOLIDTUMOR_ALT) DP\>$(VARSCAN_SOLIDTUMOR_DP) p-value\<$(VARSCAN_SOLIDTUMOR_PVAL):MPILEUP_OPTIONS='$(MPILEUP_OPTIONS)', VARSCAN_SOLIDTUMOR_BOTH_OPTIONS=$(VARSCAN_SOLIDTUMOR_BOTH_OPTIONS)"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:VarScan_HEMATOLOGY:SAMTOOLS/VARSCAN - design for HEMATOLOGY mutation, VAF\>$(VARSCAN_HEMATOLOGY_VAF) ALT\>$(VARSCAN_HEMATOLOGY_ALT) DP\>$(VARSCAN_HEMATOLOGY_DP) p-value\<$(VARSCAN_HEMATOLOGY_PVAL):MPILEUP_OPTIONS='$(MPILEUP_OPTIONS)', VARSCAN_HEMATOLOGY_BOTH_OPTIONS=$(VARSCAN_HEMATOLOGY_BOTH_OPTIONS)"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
