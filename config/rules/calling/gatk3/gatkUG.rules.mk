############################
# GATK Calling Rules
# Release: 0.9.4
# Date: 03/02/2023
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
# 0.9.5.0-31/10/2025: Parallelised GATKUG


GATKUG_FLAGS_SHARED?=

##########
# gatkUG #
##########

# GATKUG DPMIN filter
DPMIN_GATKUG=4

# GATKUG Flags
GATKUG_FLAGS=$(GATKUG_FLAGS_SHARED) \
	$(VCFDBSNP_WITH_GATK) \
	--genotype_likelihoods_model BOTH \
	--read_filter BadCigar \
	--baq OFF \
	--downsampling_type NONE \
	--min_indel_fraction_per_sample 0.01 \
	--min_indel_count_for_genotyping 2 \
	--max_deletion_fraction 0.01 \
	--standard_min_confidence_threshold_for_calling 10 \
	--downsample_to_fraction 1 \
	--min_base_quality_score 17
		

#%.gatkUG$(POST_CALLING).vcf: %.softclippedtoq0.bam %.softclippedtoq0.bam.bai %.empty.vcf %.design.bed  #%.bam %.bam.bai %.empty.vcf %.design.bed
%.gatkUG$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.design.bed
	rm -f $@.for_gatkUG*.mk;
	+if (($$($(SAMTOOLS) idxstats $< | awk '{SUM+=$$3+$$4} END {print SUM}'))); then \
		for chr in $$($(SAMTOOLS) idxstats $< | grep -v "\*" | awk '{ if ($$3+$$4>0) print $$1 }'); do \
			grep -P "^$$chr\t" $*.design.bed > $@.for_gatkUG.design.bed.$$chr.bed; \
			if [ -s $@.for_gatkUG.design.bed.$$chr.bed ]; then \
				echo "$@.for_gatkUG.$$chr.vcf.gz: $<" >> $@.for_gatkUG.generation_vcf.mk; \
				echo "	$(JAVA8) $(JAVA_FLAGS) -jar $(GATK3) $(GATKUG_FLAGS) -T UnifiedGenotyper -R $(GENOME) $$(if [ "`grep ^ -c $@.for_gatkUG.design.bed.$$chr.bed`" == "0" ]; then echo ""; else echo "-L $@.for_gatkUG.design.bed.$$chr.bed"; fi;) -I $< -ip $(INTERVAL_PADDING) -o $@.for_gatkUG.$$chr.vcf.gz " >> $@.for_gatkUG.generation_vcf.mk; \
				echo -n " $@.for_gatkUG.$$chr.vcf.gz " >> $@.for_gatkUG.list_of_vcfs.mk; \
			else \
				echo "#[INFO] No reads on chromosome $$chr for $<:"; \
				continue; \
			fi; \
		done; \
		echo -n "$@.tmp.vcf: " | cat - $@.for_gatkUG.list_of_vcfs.mk > $@.for_gatkUG.final_vcf.mk; \
		echo ""  >> $@.for_gatkUG.final_vcf.mk; \
		if [ -s $@.for_gatkUG.list_of_vcfs.mk ]; then \
			echo "	$(BCFTOOLS) concat $$(cat $@.for_gatkUG.list_of_vcfs.mk) -a -d all --threads $(THREADS_BY_CALLER) | $(BCFTOOLS) norm -m -any -o $@.tmp.vcf --threads $(THREADS_BY_CALLER)" >> $@.for_gatkUG.final_vcf.mk; \
		else \
			echo "	touch $@.tmp.vcf " >> $@.for_gatkUG.final_vcf.mk; \
		fi; \
		cat $@.for_gatkUG.generation_vcf.mk $@.for_gatkUG.final_vcf.mk >> $@.for_gatkUG.mk; \
		cat $@.for_gatkUG.mk; \
		make -f $@.for_gatkUG.mk $@.tmp.vcf; \
		rm -rf $@.for_gatkUG*; \
	else \
		touch $@.tmp.vcf; \
	fi;
	-if [ ! -e $@.tmp.vcf ]; then cp $*.empty.vcf $@.tmp.vcf; fi;
	-if [ ! -e $@.tmp.vcf ]; then touch $@.tmp.vcf; fi; 						# in case of no vcf creation, to not kill the pipeline
	$(BCFTOOLS) view -i "FORMAT/DP>$(DPMIN_GATKUG)" $@.tmp.vcf --threads $(THREADS_BY_CALLER) > $@ 	# filter on DP, cause UG is too relax, espacially because the clipping can generate few errors
	-if [ ! -e $@ ]; then cp $@.tmp.vcf $@; fi; 							# in case of error in previous line
	-rm -f $@.tmp.vcf $@.tmp.vcf.idx $@.idx


RELEASE_COMMENT := "\#\# CALLING GATKUG identify variants and generate *.gatkUG.vcf files with parameters: GATKUG_FLAGS='$(GATKUG_FLAGS)', INTERVAL_PADDING='$(INTERVAL_PADDING)', DPMIN_UG='$(DPMIN_UG)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:gatkUG:GATK Unified Genotyper - designed for PARALLEL discovery:GATKUG_FLAGS='$(GATKUG_FLAGS)', INTERVAL_PADDING='$(INTERVAL_PADDING)', DPMIN_UG='$(DPMIN_UG)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
