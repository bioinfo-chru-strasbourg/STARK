############################
# GATK Calling Rules
# Release: 1.0.0
# Date: 15/03/2026
# Author: Antony Le Bechec
############################


# Release note
# 1.0.0-15/03/2026: Create DeepVariant rule


###############
# DeepVariant #
###############


###############
# deepvariant #
###############

# Deepvariant image
DEEPVARIANT_DOCKER?=google/deepvariant:1.10.0

# Deepvariant haplotype contigs: Heterozygous variants in these contigs will be re-genotyped as the most likely
DEEPVARIANT_HAPLOID_CONTIGS?=

# Deepvariant PAR BED: BED file with PAR regions to be excluded from genotype adjustment in haploid contigs (chrX and chrY for human genome)
DEEPVARIANT_PAR_BED?=

# DeepVariant Flags
# Available and useful options:
# --model_type=WES 						# (Mandatory) Replace this string with exactly one of the following [WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA]**
# --haploid_contigs="chrX,chrY"			# Heterozygous variants in these contigs will be re-genotyped as the most likely of reference or homozygous alternates. For a sample with karyotype XY, it should be set to "chrX,chrY" for GRCh38 and "X,Y" for GRCh37. For a sample with karyotype XX, this should not be used.
# --par_regions_bed						# If --haploid_contigs is set, then this can be used to provide PAR regions to be excluded from genotype adjustment.
# --vcf_stats_report=false				# Creates VCF statistics report in html file. Default is false.
# --disable_small_model=false			# Disables the small model from make_examples stage. Default is false.
# --num_shards=$(nproc) 				# This will use all your cores to run make_examples. Feel free to change.**
# --logging_dir=/output/logs 			# This saves the log output for each stage separately.

DEEPVARIANT_FLAGS?=--model_type=WES \
	--vcf_stats_report=false \
	--disable_small_model=false \
	$(shell if [ ! -z "$(DEEPVARIANT_HAPLOID_CONTIGS)" ]; then echo ' --haploid_contigs="$(DEEPVARIANT_HAPLOID_CONTIGS)" '; fi) \
	$(shell if [ ! -z "$(DEEPVARIANT_PAR_BED)" ] && [ -e "$(DEEPVARIANT_PAR_BED)" ]; then echo ' --par_regions_bed="$(DEEPVARIANT_PAR_BED)" '; fi) \
	--dry_run=false

# DeepVariant depth filter (DP) for variants in the output VCF file (e.g., DP>=4)
DEEPVARIANT_DPMIN?=4

# DeepVariant rule
%.deepvariant$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.design.bed
	$(DOCKER_RUN) --rm --name "deepvariant-$(@F)-$$(date +%s)" $(DEEPVARIANT_DOCKER) \
		/opt/deepvariant/bin/run_deepvariant \
		$(DEEPVARIANT_FLAGS) \
		--ref=$(GENOME) \
		--reads=$< \
		--output_vcf=$@.tmp.vcf \
		--num_shards=$(THREADS_BY_CALLER) \
		--logging_dir=$@.logs \
		--regions="$*.design.bed"
	# Filter out missing genotypes (GT=./.) and keep only hom and het variants with a called genotype (GT=0/1, 1/1, etc.)
	# Filter on DP (read depth) if specified
	$(BCFTOOLS) view $@.tmp.vcf -g ^miss --include  '(GT="het" || GT="hom") && GT!="0/0" && GT!="0|0" && FORMAT/DP >= $(DEEPVARIANT_DPMIN)' --threads=$(THREADS_BY_CALLER) > $@
	# Empty if no file and clean up temporary files
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -e $@ ]; then touch $@; fi;
	-rm -rf $@.idx $@.tmp* $@.logs


RELEASE_COMMENT := "\#\# CALLING DeepVariant Caller '$(MK_RELEASE)': DeepVariant tool identify variants from aligned BAM with shared parameters: DEEPVARIANT='$(DEEPVARIANT_DOCKER)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING DeepVariant identify variants and generate *.deepvariant.vcf files with parameters: DEEPVARIANT_FLAGS='$(DEEPVARIANT_FLAGS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:deepvariant:DeepVariant - by default:DEEPVARIANT_FLAGS='$(DEEPVARIANT_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
