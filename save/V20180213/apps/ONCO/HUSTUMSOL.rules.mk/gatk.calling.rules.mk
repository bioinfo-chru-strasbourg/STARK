############################
# GATK Calling Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.3.5b"
MK_DATE="04/05/2016"

# Release note
# 0.9.1beta-10/03/2015: change genome reference location, in the file %.genome
# 0.9.1.1beta-30/04/2015: VariantFiltration correction
# 0.9.2b-23/11/2015: adding gatkUG_CPSGEN_MASTR pipeline
# 0.9.3b-25/01/2016: adding gatkUG_DIAG and gatkHC_DIAG pipelines
# 0.9.3.1b-08/02/2016: adding gatkUG_DIAGGEN pipeline
# 0.9.3.2b-03/03/2016: adding gatkUG_DIAG_PARIS pipeline
# 0.9.3.3b-11/04/2016: change gatkHC pipeline to fit with standard analysis
# 0.9.3.4b-03/05/2016: Modification of gatkUG_HUSHEMATO rule
# 0.9.3.5b-04/05/2016: Rewrite rules and update release information

# TOOLS
GATK?=$(NGSbin)/GenomeAnalysisTK.jar
GATKIG?=$(NGSbin)/IndelGenotyper.jar
VCFTOOLS?=$(NGSbin)/vcftools

# OPTIONS
JAVA_FLAGS?= -Xmx8g
#DCOV=1000000
DFRAC_HUSTUMSOL=1
MBQ_HC_HUSTUMSOL=17
MBQ_UG_HUSTUMSOL=17
#STAND_EMIT_CONF=20
#STAND_CALL_CONF=20



#########################
# GATK UnifiedGenotyper #
#########################

# Call SNPs and indels on a per-locus basis
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php

# This tool uses a Bayesian genotype likelihood model to estimate simultaneously the most likely genotypes and allele frequency in a population of N samples, emitting a genotype for each sample. The system can either emit just the variant sites or complete genotypes (which includes homozygous reference calls) satisfying some phred-scaled confidence value.

#
# -minIndelFrac,--min_indel_fraction_per_sample --min_indel_fraction_per_sample--
# Minimum fraction of all reads at a locus that must contain an indel (of any allele) for
# that sample to contribute to the indel count for alleles
# -minIndelCnt,--min_indel_count_for_genotyping --min_indel_count_for_genotyping--
# Minimum number of consensus indels required to trigger genotyping run
#-deletions,--max_deletion_fraction --max_deletion_fraction--
# Maximum fraction of reads with deletions spanning this locus for it to be callable
# --min_base_quality_score / -mbq (default HC10/UG17): Minimum base quality required to consider a base for calling



####################
# GATKUG_HUSTUMSOL #
####################

# GATKUG_HUSTUMSOL Flags
THREADS_GATK?=$(THREADS_BY_SAMPLE)
THREADS_GATKUG_HUSTUMSOL?=$(THREADS_GATK)
#GATKUG_HUSTUMSOL_THREADS=4
GATKUG_HUSTUMSOL_FLAGS= -nct $(THREADS_GATKUG_HUSTUMSOL) -glm BOTH \
		-minIndelFrac 0.01 \
		-minIndelCnt 2 \
		-deletions 0.01 \
		-baq OFF \
		-stand_call_conf 10 -dfrac $(DFRAC_HUSTUMSOL) --dbsnp $(VCFDBSNP) -mbq $(MBQ_UG_HUSTUMSOL) -rf BadCigar -dt NONE \
		-allowPotentiallyMisencodedQuals
# Minimum coverage for a variant called by GATKUG_HUSTUMSOL
DPMIN_HUSTUMSOL=30

%.gatkUG_HUSTUMSOL.unfiltered.unrecalibrated.vcf: %.bam %.bam.bai %.from_manifest.intervals %.empty.vcf %.genome
	
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKUG_HUSTUMSOL_FLAGS) \
		-T UnifiedGenotyper \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.from_manifest.intervals`" == "0" ]; then echo ""; else echo "-L $*.from_manifest.intervals"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 		# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_HUSTUMSOL) $@.tmp > $@ 		# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 		# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx



####################
# gatkHC_HUSTUMSOL #
####################


MINPRUNING_SOLIDTUMOR?=20
THREADS_GATKHC_HUSTUMSOL?=$(THREADS_GATK)
maxReadsInRegionPerSample=8000
GATKHC_HUSTUMSOL_FLAGS= -nct $(THREADS_GATKHC_HUSTUMSOL) -baq OFF -stand_call_conf 10 -dfrac $(DFRAC_HUSTUMSOL) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample) --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC_HUSTUMSOL) -rf BadCigar -minPruning $(MINPRUNING_SOLIDTUMOR) -allowPotentiallyMisencodedQuals

%.gatkHC_HUSTUMSOL.unfiltered.unrecalibrated.vcf: %.bam %.bam.bai %.from_manifest.intervals %.empty.vcf %.genome
	
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKHC_HUSTUMSOL_FLAGS) \
		-T HaplotypeCaller \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.from_manifest.intervals`" == "0" ]; then echo ""; else echo "-L $*.from_manifest.intervals"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@;
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -e $@ ]; then touch $@; fi;
	-rm -f $@.idx



RELEASE_COMMENT := "\#\# CALLING GATKUG HUSTUMSOL identify variants and generate *.gatkUG_HUSTUMSOL.vcf files with parameters: GATKUG_HUSTUMSOL_FLAGS='$(GATKUG_HUSTUMSOL_FLAGS)', DPMIN_HUSTUMSOL='$(DPMIN_HUSTUMSOL)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKHC HUSTUMSOL identify variants and generate *.gatkHC_HUSTUMSOL.vcf files from *.gatkHC_SOLDITUMOR.vcf"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkUG_HUSTUMSOL:GATK Unified Genotyper - designed for HUSTUMSOL variant discovery:GATKUG_SOLIDTUMOR_FLAGS='$(GATKUG_HUSTUMSOL_FLAGS)', DPMIN_SOLIDTUMOR='$(DPMIN_HUSTUMSOL)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkHC_HUSTUMSOL:GATK Haplotype Caller - designed for HUSTUMSOL variant discovery:GATKHC_SOLIDTUMOR_FLAGS='$(GATKHC_HUSTUMSOL_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )





