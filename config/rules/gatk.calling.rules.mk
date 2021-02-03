############################
# GATK Calling Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.3.7b"
MK_DATE="22/03/2019"

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


# TOOLS
GATK?=$(NGSbin)/GenomeAnalysisTK.jar
GATKIG?=$(NGSbin)/IndelGenotyper.jar

# OPTIONS
JAVA_FLAGS?= -Xmx8g
DCOV=1000
DFRAC=1
MBQ_HC=17
MBQ_UG=17
STAND_EMIT_CONF=20
STAND_CALL_CONF=20
GATKHC_FLAGS_SHARED=--baq OFF --read_filter BadCigar --allow_potentially_misencoded_quality_scores --dontUseSoftClippedBases



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

##########
# gatkUG #
##########

# GATKUG Flags
THREADS_GATK?=$(THREADS_BY_SAMPLE)
THREADS_GATKUG?=$(THREADS_GATK)
INTERVAL_PADDING?=0
#GATKUG_THREADS=4
GATKUG_FLAGS= -nct $(THREADS_GATKUG) -glm BOTH \
		-minIndelFrac 0.01 \
		-minIndelCnt 2 \
		-deletions 0.01 \
		-baq OFF \
		-stand_call_conf 10 -dfrac $(DFRAC) --dbsnp $(VCFDBSNP) -mbq $(MBQ_UG) -rf BadCigar -dt NONE \
		-allowPotentiallyMisencodedQuals
# Minimum coverage for a variant called by gatkUG
DPMIN=1

%.gatkUG$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	#
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKUG_FLAGS) \
		-T UnifiedGenotyper \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 		# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN) $@.tmp > $@ 		# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 		# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx

####################
# gatkUG_GERMLINE #
####################

# GATKUG Flags
GATKUG_THREADS_GERMLINE?=1 #$(THREADS_GATK)
DFRAC_UG_GERMLINE=1
MBQ_UG_UG_GERMLINE=17
STAND_EMIT_CONF_UG_GERMLINE=30
STAND_CALL_CONF_UG_GERMLINE=30
minIndelFrac_UG_GERMLINE=0.1
minIndelCnt_UG_GERMLINE=4
deletions_UG_GERMLINE=1
DPMIN_UG_GERMLINE=4
GATKUG_GERMLINE_FLAGS= -nct $(GATKUG_THREADS_GERMLINE) -glm BOTH \
		-minIndelFrac $(minIndelFrac_UG_GERMLINE)  \
		-minIndelCnt $(minIndelCnt_UG_GERMLINE) \
		-deletions $(deletions_UG_GERMLINE) \
		-baq OFF \
		-stand_call_conf $(STAND_CALL_CONF_UG_GERMLINE) -dfrac $(DFRAC_UG_GERMLINE) --dbsnp $(VCFDBSNP) -mbq $(MBQ_UG_UG_GERMLINE) -rf BadCigar -dt NONE
# Minimum coverage for a variant called by gatkUG

%.gatkUG_GERMLINE$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	#
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKUG_GERMLINE_FLAGS) \
		-T UnifiedGenotyper \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 			# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_UG_GERMLINE) $@.tmp > $@ 	# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 			# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx



####################
# gatkUG_EXOME #
####################

# GATKUG Flags
GATKUG_THREADS_EXOME?=1 #$(THREADS_GATK)
DFRAC_UG_EXOME=1
MBQ_UG_UG_EXOME=17
STAND_EMIT_CONF_UG_EXOME=30
STAND_CALL_CONF_UG_EXOME=30
minIndelFrac_UG_EXOME=0.1
minIndelCnt_UG_EXOME=4
deletions_UG_EXOME=1
DPMIN_UG_EXOME=4
GATKUG_EXOME_FLAGS= -nct $(GATKUG_THREADS_EXOME) -glm BOTH \
		-minIndelFrac $(minIndelFrac_UG_EXOME)  \
		-minIndelCnt $(minIndelCnt_UG_EXOME) \
		-deletions $(deletions_UG_EXOME) \
		-baq OFF \
		-stand_call_conf $(STAND_CALL_CONF_UG_EXOME) -dfrac $(DFRAC_UG_EXOME) --dbsnp $(VCFDBSNP) -mbq $(MBQ_UG_UG_EXOME) -rf BadCigar -dt NONE
# Minimum coverage for a variant called by gatkUG

%.gatkUG_EXOME$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	#
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKUG_EXOME_FLAGS) \
		-T UnifiedGenotyper \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 			# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_UG_EXOME) $@.tmp > $@ 	# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 			# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx



####################
# gatkUG_EXOME_SOMATIC #
####################

# GATKUG Flags
GATKUG_THREADS_EXOME_SOMATIC?=1 #$(THREADS_GATK)
DFRAC_UG_EXOME_SOMATIC=1
MBQ_UG_UG_EXOME_SOMATIC=17
STAND_EMIT_CONF_UG_EXOME_SOMATIC=30
STAND_CALL_CONF_UG_EXOME_SOMATIC=30
minIndelFrac_UG_EXOME_SOMATIC=0.01
minIndelCnt_UG_EXOME_SOMATIC=2
deletions_UG_EXOME_SOMATIC=0.01
DPMIN_UG_EXOME_SOMATIC=50
GATKUG_EXOME_SOMATIC_FLAGS= -nct $(GATKUG_THREADS_EXOME_SOMATIC) -glm BOTH \
		-minIndelFrac $(minIndelFrac_UG_EXOME_SOMATIC)  \
		-minIndelCnt $(minIndelCnt_UG_EXOME_SOMATIC) \
		-deletions $(deletions_UG_EXOME_SOMATIC) \
		-baq OFF \
		-stand_call_conf $(STAND_CALL_CONF_UG_EXOME_SOMATIC) -dfrac $(DFRAC_UG_EXOME_SOMATIC) --dbsnp $(VCFDBSNP) -mbq $(MBQ_UG_UG_EXOME_SOMATIC) -rf BadCigar -dt NONE
# Minimum coverage for a variant called by gatkUG

%.gatkUG_EXOME_SOMATIC$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	#
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKUG_EXOME_SOMATIC_FLAGS) \
		-T UnifiedGenotyper \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 			# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_UG_EXOME_SOMATIC) $@.tmp > $@ 	# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 			# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx


####################
# gatkUG_GENOME #
####################

# GATKUG Flags
GATKUG_THREADS_GENOME?=1 #$(THREADS_GATK)
DFRAC_UG_GENOME=1
MBQ_UG_UG_GENOME=17
STAND_EMIT_CONF_UG_GENOME=30
STAND_CALL_CONF_UG_GENOME=30
minIndelFrac_UG_GENOME=0.1
minIndelCnt_UG_GENOME=4
deletions_UG_GENOME=1
DPMIN_UG_GENOME=4
GATKUG_GENOME_FLAGS= -nct $(GATKUG_THREADS_GENOME) -glm BOTH \
		-minIndelFrac $(minIndelFrac_UG_GENOME)  \
		-minIndelCnt $(minIndelCnt_UG_GENOME) \
		-deletions $(deletions_UG_GENOME) \
		-baq OFF \
		-stand_call_conf $(STAND_CALL_CONF_UG_GENOME) -dfrac $(DFRAC_UG_GENOME) --dbsnp $(VCFDBSNP) -mbq $(MBQ_UG_UG_GENOME) -rf BadCigar -dt NONE
# Minimum coverage for a variant called by gatkUG

%.gatkUG_GENOME$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	#
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKUG_GENOME_FLAGS) \
		-T UnifiedGenotyper \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 			# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_UG_GENOME) $@.tmp > $@ 	# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 			# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx


#####################
# GATKUG_HEMATOLOGY #
#####################

# GZTKUG_HEMATOLOGY Flags
THREADS_GATK?=$(THREADS_BY_SAMPLE)
THREADS_GATKUG_HEMATOLOGY?=$(THREADS_GATK)
#GZTKUG_HEMATOLOGY_THREADS=4
GATKUG_HEMATOLOGY_FLAGS= -nct $(THREADS_GATKUG_HEMATOLOGY) -glm BOTH \
		-minIndelFrac 0.01 \
		-minIndelCnt 2 \
		-deletions 0.01 \
		-baq OFF \
		-stand_call_conf 10 -dfrac $(DFRAC) --dbsnp $(VCFDBSNP) -mbq $(MBQ_UG) -rf BadCigar -dt NONE \
		-allowPotentiallyMisencodedQuals
# Minimum coverage for a variant called by GATKUG_HEMATOLOGY
DPMIN_HEMATOLOGY=1

%.gatkUG_HEMATOLOGY$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	#
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKUG_HEMATOLOGY_FLAGS) \
		-T UnifiedGenotyper \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 		# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_HEMATOLOGY) $@.tmp > $@ 		# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 		# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx


##########
# GATKUG_SOLIDTUMOR #
##########

# GATKUG_SOLIDTUMOR Flags
THREADS_GATK?=$(THREADS_BY_SAMPLE)
THREADS_GATKUG_SOLIDTUMOR?=$(THREADS_GATK)
#GATKUG_SOLIDTUMOR_THREADS=4
GATKUG_SOLIDTUMOR_FLAGS= -nct $(THREADS_GATKUG_SOLIDTUMOR) -glm BOTH \
		-minIndelFrac 0.01 \
		-minIndelCnt 2 \
		-deletions 0.01 \
		-baq OFF \
		-stand_call_conf 10 -dfrac $(DFRAC) --dbsnp $(VCFDBSNP) -mbq $(MBQ_UG) -rf BadCigar -dt NONE \
		-allowPotentiallyMisencodedQuals
# Minimum coverage for a variant called by GATKUG_SOLIDTUMOR
DPMIN_SOLIDTUMOR=1

%.gatkUG_SOLIDTUMOR$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	#
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKUG_SOLIDTUMOR_FLAGS) \
		-T UnifiedGenotyper \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 		# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_SOLIDTUMOR) $@.tmp > $@ 		# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 		# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx



##################
# GATKUG_SOMATIC #
##################

# GATKUG_SOMATIC Flags
THREADS_GATK?=$(THREADS_BY_SAMPLE)
THREADS_GATKUG_SOMATIC?=$(THREADS_GATK)
#GATKUG_SOMATIC_THREADS=4
GATKUG_SOMATIC_FLAGS= -nct $(THREADS_GATKUG_SOMATIC) -glm BOTH \
		-minIndelFrac 0.01 \
		-minIndelCnt 2 \
		-deletions 0.01 \
		-baq OFF \
		-stand_call_conf 10 -dfrac $(DFRAC) --dbsnp $(VCFDBSNP) -mbq $(MBQ_UG) -rf BadCigar -dt NONE \
		-allowPotentiallyMisencodedQuals
# Minimum coverage for a variant called by GATKUG_SOMATIC
DPMIN_SOMATIC=1

%.gatkUG_SOMATIC$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	#
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKUG_SOMATIC_FLAGS) \
		-T UnifiedGenotyper \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 		# in case of no vcf creation, to not kill the pipeline
	#-$(VCFUTILS) varFilter -d $(DPMIN_SOMATIC) $@.tmp > $@ 		# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	$(BCFTOOLS) view -i "FORMAT/DP>$(DPMIN_SOMATIC)" $@.tmp > $@
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 		# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx





####################
# gatkUG_ONCOGENET #
####################

# GATKUG Flags
GATKUG_THREADS_ONCOGENET?=$(THREADS_GATK)
DCOV_ONCOGENET=1000000
DFRAC_ONCOGENET=1
MBQ_UG_ONCOGENET=17
STAND_EMIT_CONF_ONCOGENET=10
STAND_CALL_CONF_ONCOGENET=10
minIndelFrac_ONCOGENET=0.05
minIndelCnt_ONCOGENET=5
deletions_ONCOGENET=0.05
GATKUG_ONCOGENET_FLAGS= -nct $(GATKUG_THREADS_ONCOGENET) -glm BOTH \
		-minIndelFrac $(minIndelFrac_ONCOGENET)  \
		-minIndelCnt $(minIndelCnt_ONCOGENET) \
		-deletions $(deletions_ONCOGENET) \
		-baq OFF \
		-stand_call_conf $(STAND_CALL_CONF_ONCOGENET) -dfrac $(DFRAC_ONCOGENET) --dbsnp $(VCFDBSNP) -mbq $(MBQ_UG_ONCOGENET) -rf BadCigar -dt NONE
# Minimum coverage for a variant called by gatkUG
DPMIN_ONCOGENET=4

%.gatkUG_ONCOGENET$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	#
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKUG_ONCOGENET_FLAGS) \
		-T UnifiedGenotyper \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 			# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_ONCOGENET) $@.tmp > $@ 	# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
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


##########
# gatkHC #
##########

MINPRUNING?=4
THREADS_GATKHC?=$(THREADS_GATK)
maxReadsInRegionPerSample=250
GATKHC_FLAGS= -nct $(THREADS_GATKHC) -stand_call_conf 10 -dfrac $(DFRAC) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample) --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC) -minPruning $(MINPRUNING) $(GATKHC_FLAGS_SHARED)

%.gatkHC$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	#
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKHC_FLAGS) \
		-T HaplotypeCaller \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@;
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -e $@ ]; then touch $@; fi;
	-rm -f $@.idx


###################
# gatkHC_GERMLINE #
###################

MINPRUNING_HC_GERMLINE?=4
THREADS_GATKHC_GERMLINE?=$(THREADS_GATK)
maxReadsInRegionPerSample_HC_GERMLINE=1000
MBQ_HC_GERMLINE=$(MBQ_HC)
STAND_EMIT_CONF_HC_GERMLINE=30
STAND_CALL_CONF_HC_GERMLINE=30
DFRAC_HC_GERMLINE=1
DPMIN_HC_GERMLINE=4


GATKHC_GERMLINE_FLAGS= -nct $(THREADS_GATKHC_GERMLINE) -stand_call_conf $(STAND_CALL_CONF_HC_GERMLINE) -dfrac $(DFRAC_HC_GERMLINE) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample_HC_GERMLINE) --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC_GERMLINE) -minPruning $(MINPRUNING_HC_GERMLINE)  $(GATKHC_FLAGS_SHARED)

%.gatkHC_GERMLINE$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKHC_GERMLINE_FLAGS) \
		-T HaplotypeCaller \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 			# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_HC_GERMLINE) $@.tmp > $@ 	# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 			# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx


###################
# gatkHC_EXOME #
###################

MINPRUNING_HC_EXOME?=4
THREADS_GATKHC_EXOME?=$(THREADS_GATK)
maxReadsInRegionPerSample_HC_EXOME=1000
MBQ_HC_EXOME=$(MBQ_HC)
STAND_EMIT_CONF_HC_EXOME=30
STAND_CALL_CONF_HC_EXOME=30
DFRAC_HC_EXOME=1
DPMIN_HC_EXOME=4


GATKHC_EXOME_FLAGS= -nct $(THREADS_GATKHC_EXOME) -stand_call_conf $(STAND_CALL_CONF_HC_EXOME) -dfrac $(DFRAC_HC_EXOME) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample_HC_EXOME) --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC_EXOME) -minPruning $(MINPRUNING_HC_EXOME)  $(GATKHC_FLAGS_SHARED)

%.gatkHC_EXOME$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKHC_EXOME_FLAGS) \
		-T HaplotypeCaller \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 			# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_HC_EXOME) $@.tmp > $@ 	# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 			# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx


########################
# gatkHC_EXOME_SOMATIC #
########################

MINPRUNING_HC_EXOME_SOMATIC?=4
THREADS_GATKHC_EXOME_SOMATIC?=$(THREADS_GATK)
maxReadsInRegionPerSample_HC_EXOME_SOMATIC=1000
MBQ_HC_EXOME_SOMATIC=$(MBQ_HC)
STAND_EMIT_CONF_HC_EXOME_SOMATIC=30
STAND_CALL_CONF_HC_EXOME_SOMATIC=30
DFRAC_HC_EXOME_SOMATIC=1
DPMIN_HC_EXOME_SOMATIC=50


GATKHC_EXOME_SOMATIC_FLAGS= -nct $(THREADS_GATKHC_EXOME_SOMATIC) -stand_call_conf $(STAND_CALL_CONF_HC_EXOME_SOMATIC) -dfrac $(DFRAC_HC_EXOME_SOMATIC) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample_HC_EXOME_SOMATIC) --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC_EXOME_SOMATIC) -minPruning $(MINPRUNING_HC_EXOME_SOMATIC)  $(GATKHC_FLAGS_SHARED)

%.gatkHC_EXOME_SOMATIC$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKHC_EXOME_SOMATIC_FLAGS) \
		-T HaplotypeCaller \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 			# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_HC_EXOME_SOMATIC) $@.tmp > $@ 	# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 			# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx



##################
# gatkHC_GENOME #
#################

MINPRUNING_HC_GENOME?=2
THREADS_GATKHC_GENOME?=$(THREADS_GATK)
maxReadsInRegionPerSample_HC_GENOME=1000
MBQ_HC_GENOME=$(MBQ_HC)
STAND_EMIT_CONF_HC_GENOME=30
STAND_CALL_CONF_HC_GENOME=30
DFRAC_HC_GENOME=1
DPMIN_HC_GENOME=4


GATKHC_GENOME_FLAGS= -nct $(THREADS_GATKHC_GENOME) -stand_call_conf $(STAND_CALL_CONF_HC_GENOME) -dfrac $(DFRAC_HC_GENOME) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample_HC_GENOME) --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC_GENOME) -minPruning $(MINPRUNING_HC_GENOME)  $(GATKHC_FLAGS_SHARED)

%.gatkHC_GENOME$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	#
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKHC_GENOME_FLAGS) \
		-T HaplotypeCaller \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@.tmp;
	-if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp; fi;
	-if [ ! -e $@.tmp ]; then touch $@.tmp; fi; 			# in case of no vcf creation, to not kill the pipeline
	-$(VCFUTILS) varFilter -d $(DPMIN_HC_GENOME) $@.tmp > $@ 	# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp $@; fi; 			# in case of error in previous line
	-rm -f $@.tmp $@.tmp.idx $@.idx



#####################
# gatkHC_HEMATOLOGY #
#####################

MINPRUNING_HEMATOLOGY?=20
THREADS_GATKHC_HEMATOLOGY?=$(THREADS_GATK)
maxReadsInRegionPerSample_HEMATOLOGY=1000
GATKHC_HEMATOLOGY_FLAGS= -nct $(THREADS_GATKHC_HEMATOLOGY) -stand_call_conf 10 -dfrac $(DFRAC) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample_HEMATOLOGY) --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC) -minPruning $(MINPRUNING_HEMATOLOGY)  $(GATKHC_FLAGS_SHARED)

%.gatkHC_HEMATOLOGY$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	#
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKHC_HEMATOLOGY_FLAGS) \
		-T HaplotypeCaller \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@;
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -e $@ ]; then touch $@; fi;
	-rm -f $@.idx



#####################
# gatkHC_SOLIDTUMOR #
#####################

MINPRUNING_SOLIDTUMOR?=20
THREADS_GATKHC_SOLIDTUMOR?=$(THREADS_GATK)
maxReadsInRegionPerSample_SOLIDTUMOR=1000
GATKHC_SOLIDTUMOR_FLAGS= -nct $(THREADS_GATKHC_SOLIDTUMOR) -stand_call_conf 10 -dfrac $(DFRAC) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample_SOLIDTUMOR) --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC) -minPruning $(MINPRUNING_SOLIDTUMOR)  $(GATKHC_FLAGS_SHARED)

%.gatkHC_SOLIDTUMOR$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	#
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKHC_SOLIDTUMOR_FLAGS) \
		-T HaplotypeCaller \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@;
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -e $@ ]; then touch $@; fi;
	-rm -f $@.idx




##################
# gatkHC_SOMATIC #
##################

MINPRUNING_SOMATIC?=20
THREADS_GATKHC_SOMATIC?=$(THREADS_GATK)
maxReadsInRegionPerSample_SOMATIC=1000
GATKHC_SOMATIC_FLAGS= -nct $(THREADS_GATKHC_SOMATIC) -stand_call_conf 10 -dfrac $(DFRAC) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample_SOMATIC) --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC) -minPruning $(MINPRUNING_SOMATIC)  $(GATKHC_FLAGS_SHARED)

%.gatkHC_SOMATIC$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	#
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKHC_SOMATIC_FLAGS) \
		-T HaplotypeCaller \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@;
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -e $@ ]; then touch $@; fi;
	-rm -f $@.idx





####################
# gatkHC_ONCOGENET #
####################

MBQ_HC_ONCOGENET=17
MINPRUNING_ONCOGENET?=20
THREADS_GATKHC_ONCOGENET?=$(THREADS_GATK)
maxReadsInRegionPerSample=8000
GATKHC_ONCOGENET_FLAGS= -nct $(THREADS_GATKHC_ONCOGENET) -stand_call_conf 10 -dfrac $(DFRAC_ONCOGENET) --maxReadsInRegionPerSample $(maxReadsInRegionPerSample) --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC_ONCOGENET) -minPruning $(MINPRUNING_ONCOGENET)  $(GATKHC_FLAGS_SHARED)

%.gatkHC_ONCOGENET$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	#
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKHC_ONCOGENET_FLAGS) \
		-T HaplotypeCaller \
		-R `cat $*.genome` \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-o $@;
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -e $@ ]; then touch $@; fi;
	-rm -f $@.idx





##########
# gatk4HC #
##########

MINPRUNING?=4
THREADS_GATK4HC?=$(THREADS_GATK)
maxReadsInRegionPerSample=250
GATK4HC_FLAGS= --dbsnp $(VCFDBSNP) -mbq $(MBQ_HC) --min-pruning $(MINPRUNING) --native-pair-hmm-threads $(THREADS_GATKHC_ONCOGENET) $(GATK4HC_FLAGS_SHARED)

%.gatk4HC$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	#
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK4) HaplotypeCaller \
		$(GATK4HC_FLAGS) \
		-R $$(cat $*.genome) \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
		-I $< \
		-ip $(INTERVAL_PADDING) \
		-O $@;
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -e $@ ]; then touch $@; fi;
	-rm -f $@.idx




####################
# IndelGenotyperV2 #
####################

GATKIG_FLAGS= -minFraction 0.01 -minCnt 1 -B:dbsnp,vcf $(VCFDBSNP) -rf BadCigar
#-mnr 10000000

%.gatkIG.vcf: %.bam %.bam.bai %.genome %.design.bed.interval_list #%.from_manifest.interval_list
	if [ ! "$(GATKIG)" == ""] && [ -e $(GATKIG) ]; \
	then \
		$(JAVA) $(JAVA_FLAGS) -jar $(GATKIG) $(GATKIG_FLAGS) \
			-T IndelGenotyperV2 \
			-R `cat $*.genome` \
			$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "-L $*.design.bed.interval_list"; fi;) \
			-I $< \
			-o $@ ; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;
	-rm $@.idx




RELEASE_COMMENT := "\#\# CALLING GATK '$(MK_RELEASE)': GATK tool identify variants from aligned BAM with shared parameters: GATK='$(GATK)', filterExpression='$(filterExpression)', genotypeFilterExpression='$(genotypeFilterExpression)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


RELEASE_COMMENT := "\#\# CALLING GATKUG identify variants and generate *.gatkUG.vcf files with parameters: GATKUG_FLAGS='$(GATKUG_FLAGS)', DPMIN='$(DPMIN)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKUG GERMLINE identify variants and generate *.gatkUG_GERMLINE.vcf files with parameters: GATKUG_GERMLINE_FLAGS='$(GATKUG_GERMLINE_FLAGS)', DPMIN_UG_GERMLINE='$(DPMIN_UG_GERMLINE)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKUG EXOME identify variants and generate *.gatkUG_EXOME.vcf files with parameters: GATKUG_EXOME_FLAGS='$(GATKUG_EXOME_FLAGS)', DPMIN_UG_EXOME='$(DPMIN_UG_EXOME)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKUG EXOME_SOMATIC identify variants and generate *.gatkUG_EXOME_SOMATIC.vcf files with parameters: GATKUG_EXOME_SOMATIC_FLAGS='$(GATKUG_EXOME_SOMATIC_FLAGS)', DPMIN_UG_EXOME_SOMATIC='$(DPMIN_UG_EXOME_SOMATIC)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKUG GENOME identify variants and generate *.gatkUG_GENOME.vcf files with parameters: GATKUG_GENOME_FLAGS='$(GATKUG_GENOME_FLAGS)', DPMIN_UG_GENOME='$(DPMIN_UG_GENOME)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKUG HEMATOLOGY identify variants and generate *.gatkUG_HEMATOLOGY.vcf files with parameters: GATKUG_HEMATOLOGY_FLAGS='$(GATKUG_HEMATOLOGY_FLAGS)', DPMIN_HEMATOLOGY='$(DPMIN_HEMATOLOGY)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKUG SOLIDTUMOR identify variants and generate *.gatkUG_SOLIDTUMOR.vcf files with parameters: GATKUG_SOLIDTUMOR_FLAGS='$(GATKUG_SOLIDTUMOR_FLAGS)', DPMIN_SOLIDTUMOR='$(DPMIN_SOLIDTUMOR)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKUG SOMATIC identify variants and generate *.gatkUG_SOMATIC.vcf files with parameters: GATKUG_SOMATIC_FLAGS='$(GATKUG_SOMATIC_FLAGS)', DPMIN_SOMATIC='$(DPMIN_SOMATIC)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKUG ONCOGENET identify variants and generate *.gatkUG_ONCOGENET.vcf files with parameters: GATKUG_ONCOGENET_FLAGS='$(GATKUG_ONCOGENET_FLAGS)', DPMIN_SOLIDTUMOR='$(DPMIN_ONCOGENET)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


RELEASE_COMMENT := "\#\# CALLING GATKHC identify variants and generate *.gatkHC.vcf files with parameters: GATKHC_FLAGS='$(GATKHC_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKHC GERMLINE identify variants and generate *.gatkHC_GERMLINE.vcf files with parameters: GATKHC_GERMLINE_FLAGS='$(GATKHC_GERMLINE_FLAGS)', DPMIN_HC_GERMLINE='$(DPMIN_HC_GERMLINE)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKHC EXOME identify variants and generate *.gatkHC_EXOME.vcf files with parameters: GATKHC_EXOME_FLAGS='$(GATKHC_EXOME_FLAGS)', DPMIN_HC_EXOME='$(DPMIN_HC_EXOME)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKHC EXOME_SOMATIC identify variants and generate *.gatkHC_EXOME_SOMATIC.vcf files with parameters: GATKHC_EXOME_SOMATIC_FLAGS='$(GATKHC_EXOME_SOMATIC_FLAGS)', DPMIN_HC_EXOME_SOMATIC='$(DPMIN_HC_EXOME_SOMATIC)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKHC GENOME identify variants and generate *.gatkHC_GENOME.vcf files with parameters: GATKHC_GENOME_FLAGS='$(GATKHC_GENOME_FLAGS)', DPMIN_HC_GENOME='$(DPMIN_HC_GENOME)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKHC HEMATOLOGY identify variants and generate *.gatkHC_HEMATOLOGY.vcf files with parameters: GATKHC_HEMATOLOGY_FLAGS='$(GATKHC_HEMATOLOGY_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKHC SOLIDTUMOR identify variants and generate *.gatkHC_SOLIDTUMOR.vcf files with parameters: GATKHC_SOLIDTUMOR_FLAGS='$(GATKHC_SOLIDTUMOR_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING GATKHC ONCOGENET identify variants and generate *.gatkHC_ONCOGENET.vcf files with parameters: GATKHC_ONCOGENET_FLAGS='$(GATKHC_ONCOGENET_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


RELEASE_COMMENT := "\#\# CALLING GATK4HC identify variants and generate *.gatk4HC.vcf files with parameters: GATK4HC_FLAGS='$(GATK4HC_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "CALLER:gatkUG:GATK Unified Genotyper - by default:GATKUG_FLAGS='$(GATKUG_FLAGS)', DPMIN='$(DPMIN)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkUG_GERMLINE:GATK Unified Genotyper - designed for GERMLINE discovery:GATKUG_GERMLINE_FLAGS='$(GATKUG_GERMLINE_FLAGS)', DPMIN_UG_GERMLINE='$(DPMIN_UG_GERMLINE)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkUG_EXOME:GATK Unified Genotyper - designed for EXOME discovery:GATKUG_EXOME_FLAGS='$(GATKUG_EXOME_FLAGS)', DPMIN_UG_EXOME='$(DPMIN_UG_EXOME)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkUG_EXOME_SOMATIC:GATK Unified Genotyper - designed for EXOME_SOMATIC discovery:GATKUG_EXOME_SOMATIC_FLAGS='$(GATKUG_EXOME_SOMATIC_FLAGS)', DPMIN_UG_EXOME_SOMATIC='$(DPMIN_UG_EXOME_SOMATIC)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkUG_GENOME:GATK Unified Genotyper - designed for GENOME discovery:GATKUG_GENOME_FLAGS='$(GATKUG_GENOME_FLAGS)', DPMIN_UG_GENOME='$(DPMIN_UG_GENOME)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkUG_HEMATOLOGY:GATK Unified Genotyper - designed for HEMATOLOGY variant discovery:GATKUG_HEMATOLOGY_FLAGS='$(GATKUG_HEMATOLOGY_FLAGS)', DPMIN_HEMATOLOGY='$(DPMIN_HEMATOLOGY)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkUG_SOLIDTUMOR:GATK Unified Genotyper - designed for SOLIDTUMOR variant discovery:GATKUG_SOLIDTUMOR_FLAGS='$(GATKUG_SOLIDTUMOR_FLAGS)', DPMIN_SOLIDTUMOR='$(DPMIN_SOLIDTUMOR)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkUG_SOMATIC:GATK Unified Genotyper - designed for SOMATIC variant discovery:GATKUG_SOMATIC_FLAGS='$(GATKUG_SOMATIC_FLAGS)', DPMIN_SOMATIC='$(DPMIN_SOMATIC)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkUG_ONCOGENET:GATK Unified Genotyper - designed for ONCOGENET variant discovery:GATKUG_ONCOGENET_FLAGS='$(GATKUG_ONCOGENET_FLAGS)', DPMIN_SOLIDTUMOR='$(DPMIN_ONCOGENET)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


PIPELINES_COMMENT := "CALLER:gatkHC:GATK Haplotype Caller - by default:GATKHC_FLAGS='$(GATKHC_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkHC_GERMLINE:GATK Haplotype Caller - designed for GERMLINE discovery:GATKHC_GERMLINE_FLAGS='$(GATKHC_GERMLINE_FLAGS)', DPMIN_HC_GERMLINE='$(DPMIN_HC_GERMLINE)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkHC_EXOME:GATK Haplotype Caller - designed for EXOME analysis:GATKHC_EXOME_FLAGS='$(GATKHC_EXOME_FLAGS)', DPMIN_HC_EXOME='$(DPMIN_HC_EXOME)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkHC_EXOME_SOMATIC:GATK Haplotype Caller - designed for EXOME_SOMATIC analysis:GATKHC_EXOME_SOMATIC_FLAGS='$(GATKHC_EXOME_SOMATIC_FLAGS)', DPMIN_HC_EXOME_SOMATIC='$(DPMIN_HC_EXOME_SOMATIC)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkHC_GENOME:GATK Haplotype Caller - designed for GENOME analysis:GATKHC_GENOME_FLAGS='$(GATKHC_GENOME_FLAGS)', DPMIN_HC_GENOME='$(DPMIN_HC_GENOME)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkHC_HEMATOLOGY:GATK Haplotype Caller - designed for HEMATOLOGY variant discovery:GATKHC_HEMATOLOGY_FLAGS='$(GATKHC_HEMATOLOGY_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkHC_SOLIDTUMOR:GATK Haplotype Caller - designed for SOLIDTUMOR variant discovery:GATKHC_SOLIDTUMOR_FLAGS='$(GATKHC_SOLIDTUMOR_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkHC_SOMATIC:GATK Haplotype Caller - designed for SOMATIC variant discovery:GATKHC_SOMATIC_FLAGS='$(GATKHC_SOMATIC_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkHC_ONCOGENET:GATK Haplotype Caller - designed for ONCOGENET variant discovery:GATKHC_ONCOGENET_FLAGS='$(GATKHC_ONCOGENET_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


PIPELINES_COMMENT := "CALLER:gatk4HC:GATK4 Haplotype Caller - by default:GATK4HC_FLAGS='$(GATK4HC_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
