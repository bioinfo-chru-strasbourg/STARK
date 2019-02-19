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
DCOV=1000000
MBQ_HC_HUSDIAGGEN=17
MBQ_UG_HUSDIAGGEN=17
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


##################
# gatkUG_HUSDIAGGEN #
##################

# GATKUG Flags
GATKUG_THREADS_HUSDIAGGEN?=1 #$(THREADS_GATK)
DCOV_HUSDIAGGEN=1000000
DFRAC_HUSDIAGGEN=1
MBQ_UG_HUSDIAGGEN=17
STAND_EMIT_CONF_HUSDIAGGEN=10
STAND_CALL_CONF_HUSDIAGGEN=10
minIndelFrac_HUSDIAGGEN=0.05
minIndelCnt_HUSDIAGGEN=5
deletions_HUSDIAGGEN=0.05
GATKUG_HUSDIAGGEN_FLAGS= -nct $(GATKUG_THREADS_HUSDIAGGEN) -glm BOTH \
		-minIndelFrac $(minIndelFrac_HUSDIAGGEN)  \
		-minIndelCnt $(minIndelCnt_HUSDIAGGEN) \
		-deletions $(deletions_HUSDIAGGEN) \
		-baq OFF \
		-stand_call_conf $(STAND_CALL_CONF_HUSDIAGGEN) -dfrac $(DFRAC_HUSDIAGGEN) --dbsnp $(VCFDBSNP) -mbq $(MBQ_UG_HUSDIAGGEN) -rf BadCigar -dt NONE
# Minimum coverage for a variant called by gatkUG
DPMIN_HUSDIAGGEN=4

%.gatkUG_HUSDIAGGEN.unfiltered.unrecalibrated.vcf: %.bam %.bam.bai %.from_manifest.intervals %.empty.vcf %.genome
	
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKUG_HUSDIAGGEN_FLAGS) \
		-T UnifiedGenotyper \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.from_manifest.intervals`" == "0" ]; then echo ""; else echo "-L $*.from_manifest.intervals"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 			# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_HUSDIAGGEN) $@.tmp > $@ 	# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 			# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx



#####################
# gatkHC_HUSDIAGGEN #
#####################

MINPRUNING_HUSDIAGGEN?=20
THREADS_GATKHC_HUSDIAGGEN?=$(THREADS_GATK)
maxReadsInRegionPerSample=8000
GATKHC_HUSDIAGGEN_FLAGS= -nct $(THREADS_GATKHC_HUSDIAGGEN) -baq OFF -stand_call_conf 10 -dfrac $(DFRAC) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample) --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC_HUSDIAGGEN) -rf BadCigar -minPruning $(MINPRUNING_HUSDIAGGEN) -allowPotentiallyMisencodedQuals

%.gatkHC_HUSDIAGGEN.unfiltered.unrecalibrated.vcf: %.bam %.bam.bai %.from_manifest.intervals %.empty.vcf %.genome
	
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKHC_HUSDIAGGEN_FLAGS) \
		-T HaplotypeCaller \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.from_manifest.intervals`" == "0" ]; then echo ""; else echo "-L $*.from_manifest.intervals"; fi;) \
		-I $< \
		-o $@;
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -e $@ ]; then touch $@; fi;
	-rm -f $@.idx





RELEASE_COMMENT := "\#\# CALLING GATKUG HUSDIAGGEN identify variants and generate *.gatkUG_HUSDIAGGEN.vcf files with parameters: GATKUG_HUSDIAGGEN_FLAGS='$(GATKUG_HUSDIAGGEN_FLAGS)', DPMIN_HUSDIAGGEN='$(DPMIN_HUSDIAGGEN)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKHC HUSDIAGGEN identify variants and generate *.gatkHC_HUSDIAGGEN.vcf files with parameters: GATKHC_HUSDIAGGEN_FLAGS='$(GATKHC_HUSDIAGGEN_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "CALLER:gatkUG_HUSDIAGGEN:GATK Unified Genotyper - designed for HUSDIAGGEN variant  discovery:GATKUG_HUSDIAGGEN_FLAGS='$(GATKUG_HUSDIAGGEN_FLAGS)', DPMIN_HUSDIAGGEN='$(DPMIN_HUSDIAGGEN)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkHC_HUSDIAGGEN:GATK Haplotype Caller - designed for HUSDIAGGEN variant analysis:GATKHC_HUSDIAGGEN_FLAGS='$(GATKHC_HUSDIAGGEN_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )



