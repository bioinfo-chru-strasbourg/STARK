############################
# GATK Calling Rules
# Release: 0.9.5
# Date: 31/10/2025
# Author: Antony Le Bechec
############################


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
# 0.9.3.6b-10/11/2017: adding gatkUG_ONCOGENET pipeline
# 0.9.3.7b-22/03/2019: Add --dontUseSoftClippedBases for GATKHC
# 0.9.3.8-29/07/2022: Remove --dontUseSoftClippedBases for GATKHC
# 0.9.3.9-29/07/2022: Add --dontUseSoftClippedBases for GATKHC, add GATKUG_LONG_INDELS and GATKHC_LONG_INDELS
# 0.9.4-03/02/2023: Extract gatkHC
# 0.9.5-31/10/2025: Change GATKHC to GATK3HC



# OPTIONS
GATK3HC_FLAGS_SHARED=--baq OFF --read_filter BadCigar --allow_potentially_misencoded_quality_scores --dontUseSoftClippedBases


########################
# GATK HaplotypeCaller #
########################

# Call germline SNPs and indels via local re-assembly of haplotypes
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php

# The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other. It also makes the HaplotypeCaller much better at calling indels than position-based callers like UnifiedGenotyper.
# In the so-called GVCF mode used for scalable variant calling in DNA sequence data, HaplotypeCaller runs per-sample to generate an intermediate genomic gVCF (gVCF), which can then be used for joint genotyping of multiple samples in a very efficient way, which enables rapid incremental processing of samples as they roll off the sequencer, as well as scaling to very large cohort sizes (e.g. the 92K exomes of ExAC).
# In addition, HaplotypeCaller is able to handle non-diploid organisms as well as pooled experiment data. Note however that the #algorithms used to calculate variant likelihoods is not well suited to extreme allele frequencies (relative to ploidy) so its #use is not recommended for somatic (cancer) variant discovery. For that purpose, use MuTect2 instead.
# Finally, HaplotypeCaller is also able to correctly handle the splice junctions that make RNAseq a challenge for most variant callers.


###########
# gatk3HC #
###########

# GATK3HC Flags
GATK3HC_FLAGS=$(GATK3HC_FLAGS_SHARED) \
	--num_cpu_threads_per_data_thread $(THREADS_BY_CALLER) \
	--dbsnp $(VCFDBSNP) \
	--standard_min_confidence_threshold_for_calling 10 \
	--downsample_to_fraction 1 \
	--maxReadsInRegionPerSample 250 \
	--min_base_quality_score 17 \
	--minPruning 4


%.gatk3HC$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.design.bed.interval_list #%.from_manifest.interval_list
	$(JAVA8) $(JAVA_FLAGS) -XX:ParallelGCThreads=$(THREADS_BY_CALLER) -jar $(GATK3) $(GATK3HC_FLAGS) \
		-T HaplotypeCaller \
		-R $(GENOME) \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@;
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -e $@ ]; then touch $@; fi;
	-rm -f $@.idx


RELEASE_COMMENT := "\#\# CALLING GATK Haplotype Caller '$(MK_RELEASE)': GATK Haplotype Caller tool identify variants from aligned BAM with shared parameters: GATK='$(GATK3)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATK3HC identify variants and generate *.gatk3HC.vcf files with parameters: GATK3HC_FLAGS='$(GATK3HC_FLAGS)', INTERVAL_PADDING='$(INTERVAL_PADDING)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:gatk3HC:GATK Haplotype Caller - by default:GATK3HC_FLAGS='$(GATK3HC_FLAGS)', INTERVAL_PADDING='$(INTERVAL_PADDING)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
