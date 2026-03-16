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

# DeepVariant Flags
#DEEPVARIANT_FLAGS?=--model_type=WES
DEEPVARIANT_FLAGS?=--model_type=WES \
	--vcf_stats_report=false \
	--disable_small_model=false \
	--haploid_contigs="chrX,chrY" \
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
		--par_regions_bed="$*.design.bed" \
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
