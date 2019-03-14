############################
# SAMTOOLS Calling Rules
# Release: 0.9.1 
# Date: 21/09/2016
# Author: Antony Le Bechec
############################

#SAMTOOLS=$(NGSbin)/samtools


####################
# SAMTOOLS mpileup
####################

MPILEUP_SAMTOOLS_HUSHEMATO_OPTIONS= -E -d 10000000 -L 10000000 -C 50 -Q 10 -q 1 --output-tags SP,DP,DP4,DV  # -B  -d 100000 -L 100000 -C200 #MPILEUP_OPTIONS= -B  -d 100000 -o 0 -L 100000



###############################
# SAMTOOLS/BCFTOOLS HUSHEMATO #
###############################

# Filter
SAMTOOLS_HUSHEMATO_BCFTOOLS_FILTERS= -i 'INFO/DP>30'
SAMTOOLS_HUSHEMATO_FILTERS=--genotypeFilterExpression 'GQ < 99.0 && GQ >= 90.0' --genotypeFilterName 'LowGQ' --genotypeFilterExpression 'GQ < 90.0 && GQ >= 50' --genotypeFilterName 'VeryLowGQ'  --genotypeFilterExpression 'GQ < 50.0' --genotypeFilterName 'VeryVeryLowGQ' 


%.samtools_HUSHEMATO.vcf: %.bam %.bam.bai %.empty.vcf %.genome %.bam.bed %.from_manifest.bed
	# Calling	
	$(SAMTOOLS) mpileup -uf `cat $*.genome` $< $(MPILEUP_SAMTOOLS_HUSHEMATO_OPTIONS) --positions $*.from_manifest.bed | $(BCFTOOLS) call -mv -Ov -f GQ -f GP | $(BCFTOOLS) filter $(SAMTOOLS_HUSHEMATO_BCFTOOLS_FILTERS) > $@.unfiltered.vcf
	# Filtering
	-if [ ! -s $@.unfiltered.vcf ]; then cp $*.empty.vcf $@.unfiltered.vcf; fi;
	echo "#NBVARIANT ($@.unfiltered.vcf after calling)"`grep -cv ^# $@.unfiltered.vcf`
	-$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantFiltration -R `cat $*.genome` --variant $@.unfiltered.vcf -o $@ $(SAMTOOLS_HUSHEMATO_FILTERS);
	echo "#NBVARIANT ($@ after filtering)"`grep -cv ^# $@`
	# Fix
	-if [ ! -s $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -s $@ ]; then touch $@; fi;
	# Remove intermediate files
	-rm $@.unfiltered.*vcf* $@.idx




# CONFIG/RELEASE
RELEASE_CMD := $(shell echo "\#\# CALLING: SAMTOOLS HUSHEMATO ($(SAMTOOLS_VERSION)), BCFTOOLS ($(BCFTOOLS_VERSION)) and GATK ($(GATK_VERSION)) VariantFiltration ($(SAMTOOLS_HUSHEMATO_FILTERS)) to identify variants and generate *.samtools_HUSHEMATO.vcf files " >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "CALLER:samtools_HUSHEMATO:SAMTOOLS/BCLTOOLS - filtering variant designed for HUSHEMATO depth\<30"
# :MPILEUP_SAMTOOLS_OPTIONS='$(MPILEUP_SAMTOOLS_OPTIONS)',SAMTOOLS_FILTERS='$(SAMTOOLS_FILTERS)'
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

