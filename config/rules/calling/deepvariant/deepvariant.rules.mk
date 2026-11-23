############################
# DeepVariant Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE=1.0.0
MK_DATE=19/08/2026

# Release note
# 1.0.0-19/08/2026: Create DeepVariant rule


###############
# deepvariant #
###############

# DeepVariant Flags
# Available and useful options:
# --model_type=WES 						# (Mandatory) Replace this string with exactly one of the following [WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA]**
# --haploid_contigs="chrX,chrY"			# Heterozygous variants in these contigs will be re-genotyped as the most likely of reference or homozygous alternates. For a sample with karyotype XY, it should be set to "chrX,chrY" for GRCh38 and "X,Y" for GRCh37. For a sample with karyotype XX, this should not be used.
# --par_regions_bed						# If --haploid_contigs is set, then this can be used to provide PAR regions to be excluded from genotype adjustment.
# --vcf_stats_report=false				# Creates VCF statistics report in html file. Default is false.
# --disable_small_model=false			# Disables the small model from make_examples stage. Default is false.
# --num_shards=$(nproc) 				# This will use all your cores to run make_examples. Feel free to change.**
# --logging_dir=/output/logs 			# This saves the log output for each stage separately.

# Deepvariant image
DEEPVARIANT_DOCKER?=google/deepvariant:1.10.0

# Deepvariant haplotype contigs: Heterozygous variants in these contigs will be re-genotyped as the most likely
DEEPVARIANT_HAPLOID_CONTIGS?=
DEEPVARIANT_HAPLOID_CONTIGS_FLAGS?=$(shell if [ ! -z "$(DEEPVARIANT_HAPLOID_CONTIGS)" ]; then echo ' --haploid_contigs="$(DEEPVARIANT_HAPLOID_CONTIGS)" '; fi)

# Deepvariant PAR BED: BED file with PAR regions to be excluded from genotype adjustment in haploid contigs (chrX and chrY for human genome)
DEEPVARIANT_PAR_BED?=
DEEPVARIANT_PAR_BED_FLAGS?=$(shell if [ ! -z "$(DEEPVARIANT_PAR_BED)" ] && [ -e "$(DEEPVARIANT_PAR_BED)" ]; then echo ' --par_regions_bed="$(DEEPVARIANT_PAR_BED)" '; fi)

# Flags
# This variable can be defined directly on Apps
DEEPVARIANT_FLAGS?= \
	$(DEEPVARIANT_HAPLOID_CONTIGS_FLAGS) \
	$(DEEPVARIANT_PAR_BED_FLAGS) \
	--model_type=WES \
	--vcf_stats_report=false \
	--disable_small_model=false \
	--dry_run=false

# Filters Flags
# This variable can be defined directly on Apps
DEEPVARIANT_FILTERS_FLAGS?= \
	--exclude-genotypes '^miss' \
	--exclude '' \
	--include '(GT="het" || GT="hom") && GT!="0/0" && GT!="0|0" && FORMAT/DP>=4' \
	--filter '' \
	--sort


# Rule
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
	# Filters
	bash $(STARK_FOLDER_BIN)/filter_vcf.sh --input $@.tmp.vcf --output $@ --threads $(THREADS_BY_CALLER) --bcftools $(BCFTOOLS) $(DEEPVARIANT_FILTERS_FLAGS)
	#$(BCFTOOLS) view $@.tmp.vcf -g ^miss --include  '(GT="het" || GT="hom") && GT!="0/0" && GT!="0|0" && FORMAT/DP >= $(DEEPVARIANT_DPMIN)' --threads=$(THREADS_BY_CALLER) > $@
	# Empty if no file and clean up temporary files
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	-if [ ! -e $@ ]; then touch $@; fi;
	-rm -rf $@.idx $@.tmp* $@.logs


# Documentation
$(file >>$(RELEASE_INFOS),## CALLER deepvariant '$(MK_RELEASE)': DeepVariant caller tool identify variants from aligned BAM, generates *.deepvariant.vcf files. Parameters DEEPVARIANT_DOCKER='$(DEEPVARIANT_DOCKER)', DEEPVARIANT_HAPLOID_CONTIGS='$(DEEPVARIANT_HAPLOID_CONTIGS)', DEEPVARIANT_PAR_BED='$(DEEPVARIANT_PAR_BED)', DEEPVARIANT_FLAGS='$(DEEPVARIANT_FLAGS)', DEEPVARIANT_FILTERS_FLAGS='$(DEEPVARIANT_FILTERS_FLAGS)')
$(file >>$(PIPELINES_INFOS),CALLER:deepvariant:DeepVariant - by default:DEEPVARIANT_DOCKER='$(DEEPVARIANT_DOCKER)', DEEPVARIANT_HAPLOID_CONTIGS='$(DEEPVARIANT_HAPLOID_CONTIGS)', DEEPVARIANT_PAR_BED='$(DEEPVARIANT_PAR_BED)', DEEPVARIANT_FLAGS='$(DEEPVARIANT_FLAGS)', DEEPVARIANT_FILTERS_FLAGS='$(DEEPVARIANT_FILTERS_FLAGS)')
