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
DFRAC=1
MBQ_HC=10
MBQ_UG=17
STAND_EMIT_CONF=20
STAND_CALL_CONF=20
GATKHC_FLAGS_DIAG__SHARED=--baq OFF --read_filter BadCigar --allow_potentially_misencoded_quality_scores --dontUseSoftClippedBases



#########################
# GATK UnifiedGenotyper #
#########################

# Call SNPs and indels on a per-locus basis
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php

# This tool uses a Bayesian genotype likelihood model to estimate simultaneously the most likely genotypes and allele frequency in a population of N samples, emitting a genotype for each sample. The system can either emit just the variant sites or complete genotypes (which includes homozygous reference calls) satisfying some phred-scaled confidence value.

#
# -minIndelFrac,--min_indel_fraction_per_sample
# Minimum fraction of all reads at a locus that must contain an indel (of any allele) for
# that sample to contribute to the indel count for alleles
# -minIndelCnt,--min_indel_count_for_genotyping
# Minimum number of consensus indels required to trigger genotyping run
#-deletions,--max_deletion_fraction
# Maximum fraction of reads with deletions spanning this locus for it to be callable
# --min_base_quality_score / -mbq (default HC10/UG17): Minimum base quality required to consider a base for calling


###############
# gatkUG_DIAG_IP #
###############

# GATKUG Flags
GATKUG_THREADS_DIAG_IP?=1 #$(THREADS_GATK)
DFRAC_UG_DIAG_IP=1
MBQ_UG_UG_DIAG_IP=17
STAND_EMIT_CONF_UG_DIAG_IP=30
STAND_CALL_CONF_UG_DIAG_IP=30
minIndelFrac_UG_DIAG_IP=0.1
minIndelCnt_UG_DIAG_IP=4
deletions_UG_DIAG_IP=1
DPMIN_UG_DIAG_IP=4
GATKUG_DIAG_IP_FLAGS= -nct $(GATKUG_THREADS_DIAG_IP) -glm BOTH \
		-minIndelFrac $(minIndelFrac_UG_DIAG_IP)  \
		-minIndelCnt $(minIndelCnt_UG_DIAG_IP) \
		-deletions $(deletions_UG_DIAG_IP) \
		-baq OFF \
		-stand_call_conf $(STAND_CALL_CONF_UG_DIAG_IP) -dfrac $(DFRAC_UG_DIAG_IP) --dbsnp $(VCFDBSNP) -mbq $(MBQ_UG_UG_DIAG_IP) -rf BadCigar -dt NONE
# Minimum coverage for a variant called by gatkUG

%.gatkUG_DIAG_IP.unfiltered.unrecalibrated.vcf: %.bam %.bam.bai %.from_manifest.intervals %.empty.vcf %.genome

	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKUG_DIAG_IP_FLAGS) \
		-T UnifiedGenotyper \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.from_manifest.intervals`" == "0" ]; then echo ""; else echo "-L $*.from_manifest.intervals"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 			# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_UG_DIAG_IP) $@.tmp > $@ 	# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 			# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx



###############
# gatkUG_DIAG_MOSAIC #
###############

# GATKUG Flags
GATKUG_THREADS_DIAG_MOSAIC?=1 #$(THREADS_GATK)
DFRAC_UG_DIAG_MOSAIC=1
MBQ_UG_DIAG_MOSAIC=17
STAND_EMIT_CONF_UG_DIAG_MOSAIC=30
STAND_CALL_CONF_UG_DIAG_MOSAIC=30
minIndelFrac_UG_DIAG_MOSAIC=0.00
minIndelCnt_UG_DIAG_MOSAIC=1
deletions_UG_DIAG_MOSAIC=1
# Minimum coverage for a variant called by gatkUG
DPMIN_UG_DIAG_MOSAIC=4

GATKUG_DIAG_MOSAIC_FLAGS= -nct $(GATKUG_THREADS_DIAG_MOSAIC) -glm BOTH \
		-minIndelFrac $(minIndelFrac_UG_DIAG_MOSAIC)  \
		-minIndelCnt $(minIndelCnt_UG_DIAG_MOSAIC) \
		-deletions $(deletions_UG_DIAG_MOSAIC) \
		-baq OFF \
		-stand_call_conf $(STAND_CALL_CONF_UG_DIAG_MOSAIC) -dfrac $(DFRAC_UG_DIAG_MOSAIC) --dbsnp $(VCFDBSNP) -mbq $(MBQ_UG_DIAG_MOSAIC) -rf BadCigar -dt NONE
# Minimum coverage for a variant called by gatkUG

%.gatkUG_DIAG_MOSAIC.unfiltered.unrecalibrated.vcf: %.bam %.bam.bai %.from_manifest.intervals %.empty.vcf %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKUG_DIAG_MOSAIC_FLAGS) \
		-T UnifiedGenotyper \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.from_manifest.intervals`" == "0" ]; then echo ""; else echo "-L $*.from_manifest.intervals"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 			# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_UG_DIAG_MOSAIC) $@.tmp > $@ 	# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 			# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx




########################
# GATK HaplotypeCaller #
########################

# Call germline SNPs and indels via local re-assembly of haplotypes
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php

# The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other. It also makes the HaplotypeCaller much better at calling indels than position-based callers like UnifiedGenotyper.
#
#In the so-called GVCF mode used for scalable variant calling in DNA sequence data, HaplotypeCaller runs per-sample to generate an intermediate genomic gVCF (gVCF), which can then be used for joint genotyping of multiple samples in a very efficient way, which enables rapid incremental processing of samples as they roll off the sequencer, as well as scaling to very large cohort sizes (e.g. the 92K exomes of ExAC).
#
#In addition, HaplotypeCaller is able to handle non-diploid organisms as well as pooled experiment data. Note however that the #algorithms used to calculate variant likelihoods is not well suited to extreme allele frequencies (relative to ploidy) so its #use is not recommended for somatic (cancer) variant discovery. For that purpose, use MuTect2 instead.
#
# Finally, HaplotypeCaller is also able to correctly handle the splice junctions that make RNAseq a challenge for most variant callers.



##################
# gatkHC_DIAG_IP #
##################

MINPRUNING_HC_DIAG_IP?=4
THREADS_GATKHC_DIAG_IP?=$(THREADS_GATK)
maxReadsInRegionPerSample_HC_DIAG_IP=1000
MBQ_HC_DIAG_IP=$(MBQ_HC)
STAND_EMIT_CONF_HC_DIAG_IP=30
STAND_CALL_CONF_HC_DIAG_IP=30
DFRAC_HC_DIAG_IP=1
DPMIN_HC_DIAG_IP=4


GATKHC_DIAG_IP_FLAGS= -nct $(THREADS_GATKHC_DIAG_IP) -stand_call_conf $(STAND_CALL_CONF_HC_DIAG_IP) -dfrac $(DFRAC_HC_DIAG_IP) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample_HC_DIAG_IP) --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC_DIAG_IP) -minPruning $(MINPRUNING_HC_DIAG_IP) $(GATKHC_FLAGS_DIAG__SHARED)

%.gatkHC_DIAG_IP.unfiltered.unrecalibrated.vcf: %.bam %.bam.bai %.from_manifest.intervals %.empty.vcf %.genome
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKHC_DIAG_IP_FLAGS) \
		-T HaplotypeCaller \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.from_manifest.intervals`" == "0" ]; then echo ""; else echo "-L $*.from_manifest.intervals"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 			# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_HC_DIAG_IP) $@.tmp > $@ 	# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 			# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx



######################
# gatkHC_DIAG_MOSAIC #
######################

JAVA_FLAGS_HC_DIAG_MOSAIC?=-Xmx8g
MINPRUNING_HC_DIAG_MOSAIC?=1
THREADS_GATKHC_DIAG_MOSAIC?=$(THREADS_GATK)
maxReadsInRegionPerSample_HC_DIAG_MOSAIC=1000
MBQ_HC_DIAG_MOSAIC=$(MBQ_HC)
DPMIN_HC_DIAG_MOSAIC=4
STAND_EMIT_CONF_HC_DIAG_MOSAIC=30
STAND_CALL_CONF_HC_DIAG_MOSAIC=30

GATKHC_DIAG_MOSAIC_FLAGS= -nct $(THREADS_GATKHC_DIAG_MOSAIC) -stand_call_conf $(STAND_CALL_CONF_HC_DIAG_MOSAIC) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample_HC_DIAG_MOSAIC) --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC_DIAG_MOSAIC) -minPruning $(MINPRUNING_HC_DIAG_MOSAIC) $(GATKHC_FLAGS_DIAG__SHARED)

%.gatkHC_DIAG_MOSAIC.unfiltered.unrecalibrated.vcf: %.bam %.bam.bai %.from_manifest.intervals %.empty.vcf %.genome
	$(JAVA) $(JAVA_FLAGS_HC_DIAG_MOSAIC) -jar $(GATK) $(GATKHC_DIAG_MOSAIC_FLAGS) \
		-T HaplotypeCaller \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.from_manifest.intervals`" == "0" ]; then echo ""; else echo "-L $*.from_manifest.intervals"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 			# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_HC_DIAG_MOSAIC) $@.tmp > $@ 	# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 			# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx




RELEASE_COMMENT := "\#\# CALLING GATKUG DIAG_IP identify variants and generate *.gatkUG_DIAG_IP.vcf files with parameters: GATKUG_DIAG_IP_FLAGS='$(GATKUG_DIAG_IP_FLAGS)', DPMIN_UG_DIAG_IP='$(DPMIN_UG_DIAG_IP)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKUG DIAG_MOSAIC identify variants and generate *.gatkUG_DIAG_MOSAIC.vcf files with parameters: GATKUG_DIAG_MOSAIC_FLAGS='$(GATKUG_DIAG_MOSAIC_FLAGS)', DPMIN_UG_DIAG_MOSAIC='$(DPMIN_UG_DIAG_MOSAIC)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKHC DIAG_IP identify variants and generate *.gatkHC_DIAG_IP.vcf files with parameters: GATKHC_DIAG_IP_FLAGS='$(GATKHC_DIAG_IP_FLAGS)', DPMIN_HC_DIAG_IP='$(DPMIN_HC_DIAG_IP)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKHC DIAG_MOSAIC identify variants and generate *.gatkHC_DIAG_MOSAIC.vcf files with parameters: GATKHC_DIAG_MOSAIC_FLAGS='$(GATKHC_DIAG_MOSAIC_FLAGS)', DPMIN_HC_DIAG_MOSAIC='$(DPMIN_HC_DIAG_MOSAIC)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )



PIPELINES_COMMENT := "CALLER:gatkUG_DIAG_IP:GATK Unified Genotyper - designed for DIAG_IP discovery:GATKUG_DIAG_IP_FLAGS='$(GATKUG_DIAG_IP_FLAGS)', DPMIN_UG_DIAG_IP='$(DPMIN_UG_DIAG_IP)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkUG_DIAG_MOSAIC:GATK Unified Genotyper - designed for DIAG_MOSAIC variant discovery:GATKUG_DIAG_MOSAIC_FLAGS='$(GATKUG_DIAG_MOSAIC_FLAGS)', DPMIN_UG_DIAG_MOSAIC='$(DPMIN_UG_DIAG_MOSAIC)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkHC_DIAG_IP:GATK Haplotype Caller - designed for DIAG_IP analysis:GATKHC_DIAG_IP_FLAGS='$(GATKHC_DIAG_IP_FLAGS)', DPMIN_HC_DIAG_IP='$(DPMIN_HC_DIAG_IP)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkHC_DIAG_MOSAIC:GATK Haplotype Caller - designed for DIAG_MOSAIC discovery:GATKHC_DIAG_MOSAIC_FLAGS='$(GATKHC_DIAG_MOSAIC_FLAGS)', DPMIN_HC_DIAG_MOSAIC='$(DPMIN_HC_DIAG_MOSAIC)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
