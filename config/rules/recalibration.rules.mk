############################
# GATK Recalibration Rules
# Release: 0.9.1
# Date: 19/01/2015
# Author: Antony Le Bechec
############################




# BAM RECALIBRATION

%.bam.grp: %.bam %.bam.bai %.from_manifest.interval_list %.genome
	# Generate BaseRecalibrator grp file for recalibration
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T BaseRecalibrator -R $$(cat $*.genome) -knownSites $(VCFDBSNP) -I $< -o $@ -L $*.from_manifest.interval_list -nct $(THREADS_BY_SAMPLE) -U -compress 0

%.bam: %.recalibration.bam %.genome %.recalibration.bam.bai %.recalibration.bam.grp #%.recalibration.bam.bai %.from_manifest.intervals %.recalibration.from_manifest.intervals
	# Recalibrate BAM with BaseRecalibrator grp file
	#$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T PrintReads -R $$(cat $*.genome) -I $< -BQSR $*.recalibration.bam.grp -o $@ -L $*.recalibration.from_manifest.intervals -nct $(THREADS_BY_SAMPLE) -U
	$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T PrintReads -R $$(cat $*.genome) -I $< -BQSR $*.recalibration.bam.grp -o $@ -nct $(THREADS_BY_SAMPLE) -U -EOQ
	# clean
	-rm -f $<;




# VCF VariantRecalibration (optionnal)
# all the GATK commands are preceded of "&&" in order to stop the process if the previous command failed, and if one of the GATK command failed, we move %.unfiltered.unrecalibrated.vcf in %.unfiltered.vcf ( || symbol )
#####################

# DBs options for variant recalibration on SNPs
#recalibration_DBs_SNPs_options?=-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $(HAPMAP) -resource:omni,known=false,training=true,truth=true,prior=12.0 $(OMNI) -resource:1000G,known=false,training=true,truth=false,prior=10.0 $(PHASE1_1000G) -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $(VCFDBSNP137VCF)
#recalibration_DBs_indels_options?=-resource:mills,known=false,training=true,truth=true,prior=12.0 $(VCFMILLS1000G) -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $(VCFDBSNP137VCF)
#recalibration_DBs_SNPs_options?=-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $(VCFDBSNP137VCF)
#recalibration_DBs_indels_options?=-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $(VCFDBSNP137VCF)
recalibration_DBs_SNPs_options?=-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $(VCFDBSNP_RECALIBRATION)
#recalibration_DBs_indels_options?=-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $(VCFDBSNP_RECALIBRATION)
recalibration_DBs_indels_options?=$(VCFDBSNP_RECALIBRATION)
#VCFDBSNP_RECALIBRATION

# options for variant recalibration on SNPs
recalibration_SNPs_options=$(recalibration_DBs_SNPs_options) -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an DP -minNumBad 1000 --maxGaussians 4 -mode SNP
# options for variant recalibration on indels
recalibration_indels_options=$(recalibration_DBs_indels_options) --maxGaussians 4 -minNumBad 1000  -an DP -an FS -an ReadPosRankSum -an MQRankSum -mode INDEL

#%.unfiltered.vcf: %.unfiltered.unrecalibrated.vcf %.unfiltered.unrecalibrated.vcf.idx %.empty.vcf %.genome
%.vcf: %.recalibration.vcf %.recalibration.vcf.idx %.empty.vcf %.genome
	#if (($(VARIANT_RECALIBRATION))); then
	if ((0)); then \
		$(BCFTOOLS) view --types=indels $< > $*.indels.vcf; \
		$(BCFTOOLS) view --exclude-types=indels $< > $*.snps.vcf; \
		$(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantRecalibrator -R $$(cat $*.genome) -input $*.snps.vcf -recalFile $*.snps.recal -tranchesFile $*.snps.tranches -ip $(INTERVAL_PADDING) $(recalibration_SNPs_options) \
		&& $(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T VariantRecalibrator -R $$(cat $*.genome) -input $*.indels.vcf -recalFile $*.indels.recal -tranchesFile $*.indels.tranches -ip $(INTERVAL_PADDING) $(recalibration_indels_options) \
		&& $(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T ApplyRecalibration -R $$(cat $*.genome) -input $*.snps.vcf -tranchesFile $*.snps.tranches -recalFile $*.snps.recal -ip $(INTERVAL_PADDING) --ts_filter_level 99.5 -mode "SNP" -o $*.recalibrated.snps.vcf \
		&& $(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T ApplyRecalibration -R $$(cat $*.genome) -input $*.indels.vcf -tranchesFile $*.indels.tranches -recalFile $*.indels.recal -ip $(INTERVAL_PADDING) --ts_filter_level 99.0 	-mode "INDEL" -o $*.recalibrated.indels.vcf \
		&& $(JAVA) $(JAVA_FLAGS) -jar $(GATK) -T CombineVariants -R $$(cat $*.genome) --variant $*.recalibrated.snps.vcf --variant $*.recalibrated.indels.vcf -o $@ \
		|| mv $< $@ && echo "Variant Recalibration failed, so no Variant Recalibration is done on your data"; \
		rm -f $*.indels.vcf $*.snps.vcf $*.indels.vcf.idx $*.snps.vcf.idx $*.recalibrated.snps.vcf $*.recalibrated.indels.vcf $*.snps.recal $*.snps.tranches $*.indels.recal $*.indels.tranches; \
	else \
		mv $< $@; \
		echo "Variant Recalibration not available yet"; \
	fi;

# VariantRecalibrator GATK4
# /STARK/tools/java/current/bin/java  -Xmx7g  -Dsnappy.disable=true -Dsamjdk.try_use_intel_deflater=false   -Dorg.xerial.snappy.tempdir=/STARK/data/STARKData/V20210401b/output/tmp -Djava.io.tmpdir=/STARK/data/STARKData/V20210401b/output/tmp -jar $GATK4 VariantRecalibrator 
# -R $REF
# -V $VCF.snps.vcf.gz
# --resource:dbsnp,known=true,training=true,truth=true,prior=2.0 $DBSNP
# -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR
# -mode SNP
# -O output.recal
# --tranches-file output.tranches
# --rscript-file output.plots.R



RELEASE_COMMENT := "\#\# BAM RECALIBRATION: GATK BaseRecalibrator and PrintReads are used to recalibrate BAM files."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_ALIGNMENT:recalibration:BaseRecalibrator of reads in BAM. Warning: step BAM destructive, i.e. remove reads"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


RELEASE_COMMENT := "\#\# VCF RECALIBRATION: GATK BaseRecalibrator and PrintReads are used to recalibrate BAM files."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "POST_CALLING:recalibration:VariantRecalibrator of VCF."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

#PIPELINES_COMMENT := "POST_ANNOTATION:recalibration:VariantRecalibrator of VCF."
#PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
