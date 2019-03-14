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
DFRAC_HUSHEMATO=1
MBQ_HC_HUSHEMATO=17
MBQ_UG_HUSHEMATO=17
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


#####################
# GATKUG_HUSHEMATO #
#####################

# GZTKUG_HUSHEMATO Flags
THREADS_GATK?=$(THREADS_BY_SAMPLE)
THREADS_GATKUG_HUSHEMATO?=$(THREADS_GATK)
#GZTKUG_HUSHEMATO_THREADS=4
GATKUG_HUSHEMATO_FLAGS= -nct $(THREADS_GATKUG_HUSHEMATO) -glm BOTH \
		-minIndelFrac 0.01 \
		-minIndelCnt 2 \
		-deletions 0.01 \
		-baq OFF \
		-stand_call_conf 10 -dfrac $(DFRAC_HUSHEMATO) --dbsnp $(VCFDBSNP) -mbq $(MBQ_UG_HUSHEMATO) -rf BadCigar -dt NONE \
		-allowPotentiallyMisencodedQuals
# Minimum coverage for a variant called by GATKUG_HUSHEMATO
DPMIN_HUSHEMATO=4

%.gatkUG_HUSHEMATO.unfiltered.unrecalibrated.vcf: %.bam %.bam.bai %.from_manifest.intervals %.empty.vcf %.genome
	
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKUG_HUSHEMATO_FLAGS) \
		-T UnifiedGenotyper \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.from_manifest.intervals`" == "0" ]; then echo ""; else echo "-L $*.from_manifest.intervals"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 		# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_HUSHEMATO) $@.tmp > $@ 		# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 		# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx


#####################
# GATKUG_HUSHEMATOHALOPLEX #
#####################

# GZTKUG_HUSHEMATO Flags
THREADS_GATK?=$(THREADS_BY_SAMPLE)
THREADS_GATKUG_HUSHEMATOHALOPLEX?=1
#$(THREADS_GATK)
#GZTKUG_HUSHEMATO_THREADS=4
GATKUG_HUSHEMATOHALOPLEX_FLAGS= -nct $(THREADS_GATKUG_HUSHEMATOHALOPLEX) -glm BOTH \
		-minIndelFrac 0.01 \
		-minIndelCnt 2 \
		-deletions 0.01 \
		-baq OFF \
		-stand_call_conf 10 -dfrac $(DFRAC_HUSHEMATO) --dbsnp $(VCFDBSNP) -mbq $(MBQ_UG_HUSHEMATO) -rf BadCigar -dt NONE \
		-allowPotentiallyMisencodedQuals
# Minimum coverage for a variant called by GATKUG_HUSHEMATO
DPMIN_HUSHEMATOHALOPLEX=1

%.gatkUG_HUSHEMATOHALOPLEX.unfiltered.unrecalibrated.vcf: %.bam %.bam.bai %.from_manifest.intervals %.empty.vcf %.genome
	
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKUG_HUSHEMATOHALOPLEX_FLAGS) \
		-T UnifiedGenotyper \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.from_manifest.intervals`" == "0" ]; then echo ""; else echo "-L $*.from_manifest.intervals"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 		# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_HUSHEMATOHALOPLEX) $@.tmp > $@ 		# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 		# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx


#####################
# gatkHC_HUSHEMATO #
#####################

MINPRUNING_HUSHEMATO?=20
THREADS_GATKHC_HUSHEMATO?=$(THREADS_GATK)
maxReadsInRegionPerSample=8000
GATKHC_HUSHEMATO_FLAGS= -nct $(THREADS_GATKHC_HUSHEMATO) -baq OFF -stand_call_conf 10 -dfrac $(DFRAC_HUSHEMATO) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample) --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC_HUSHEMATO) -rf BadCigar -minPruning $(MINPRUNING_HUSHEMATO) -allowPotentiallyMisencodedQuals

%.gatkHC_HUSHEMATO.unfiltered.unrecalibrated.vcf: %.bam %.bam.bai %.from_manifest.intervals %.empty.vcf %.genome
	
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKHC_HUSHEMATO_FLAGS) \
		-T HaplotypeCaller \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.from_manifest.intervals`" == "0" ]; then echo ""; else echo "-L $*.from_manifest.intervals"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@;
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -e $@ ]; then touch $@; fi;
	-rm -f $@.idx




RELEASE_COMMENT := "\#\# CALLING GATKUG HUSHEMATO identify variants and generate *.gatkUG_HUSHEMATO.vcf files with parameters: GATKUG_HUSHEMATO_FLAGS='$(GATKUG_HUSHEMATO_FLAGS)', DPMIN_HUSHEMATO='$(DPMIN_HUSHEMATO)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKHC HUSHEMATO identify variants and generate *.gatkHC_HUSHEMATO.vcf files from *.gatkHC_HUSHEMATO.vcf"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "CALLER:gatkUG_HUSHEMATO:GATK Unified Genotyper - designed for HUSHEMATO variant discovery:GATKUG_HUSHEMATO_FLAGS='$(GATKUG_HUSHEMATO_FLAGS)', DPMIN_HUSHEMATO='$(DPMIN_HUSHEMATO)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkHC_HUSHEMATO:GATK Haplotype Caller - designed for HUSHEMATO variant discovery:GATKHC_HUSHEMATO_FLAGS='$(GATKHC_HUSHEMATO_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )





