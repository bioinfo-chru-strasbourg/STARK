############################
# Metrics Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.5b"
MK_DATE="17/03/2019"

# Release note
# 18/12/2015 - 0.9.2b : Force gzip metrics files
# 23/09/2016 - 0.9.3b : Add BAM Check metrics.bam_check. Change PICARD version
# 29/09/2016 - 0.9.4b : Add Amplicon coverage metrics.amplicon_coverage
# 29/09/2016 - 0.9.5b : Chenge metrics.genes rule


# TOOLS
SAMTOOLS?=$(NGSbin)/samtools
JAVA?=java
PICARDLIB?=$(NGSbin)/picard-tools
FASTQC?=$(NGSbin)/fastqc
NGSscripts?=$(NGS_SCRIPTS)

# OPTIONS
BED?=
PRIMER_BED?=
BAM_METRICS?=1

FATBAM_TMP_FOLDER?=$(TMP_FOLDER_TMP)

METRICS_SNPEFF?=0

GZ?=gzip

DP_FAIL?=30
DP_WARN?=100
DP_THRESHOLD?=1


# BED / INTERVALS for Metrics
###############################


# INTERVAL from BED

%.bam.bed: %.bam %.bam.bai
	#BAM.BED from BAM
	-+if ((1)); then \
	if (($$($(SAMTOOLS) idxstats $< | awk '{SUM+=$$3+$$4} END {print SUM}'))); then \
		rm -f $<.genomeCoverageBed.mk $<.genomeCoverageBed1.mk $<.genomeCoverageBed2.mk $<.genomeCoverageBed3.mk; \
		for chr in $$($(SAMTOOLS) idxstats $< | grep -v "\*" | awk '{ if ($$3+$$4>0) print $$1 }'); do \
			echo "$<.genomeCoverageBed.$$chr.bed: $<" >> $<.genomeCoverageBed1.mk; \
			echo "	$(SAMTOOLS) view $< -b $$chr | $(BEDTOOLS)/genomeCoverageBed -ibam stdin -bg | $(BEDTOOLS)/mergeBed -n -i - > $<.genomeCoverageBed.$$chr.bed " >> $<.genomeCoverageBed1.mk; \
			echo -n " $<.genomeCoverageBed.$$chr.bed" >> $<.genomeCoverageBed2.mk; \
		done; \
		echo -n "$@: " | cat - $<.genomeCoverageBed2.mk > $<.genomeCoverageBed3.mk; \
		echo ""  >> $<.genomeCoverageBed3.mk; \
		echo "	cat $$^ > $@ " >> $<.genomeCoverageBed3.mk; \
		echo "	-rm -f $$^ " >> $<.genomeCoverageBed3.mk; \
		#echo "	-rm -f \$^ " >> $<.genomeCoverageBed3.mk; \
		cat $<.genomeCoverageBed1.mk $<.genomeCoverageBed3.mk >> $<.genomeCoverageBed.mk; \
		make -j $(THREADS) -i -f $<.genomeCoverageBed.mk $@ 1>/dev/null 2>/dev/null; \
		rm $<.genomeCoverageBed*.mk; \
	fi; \
	fi;


# BED for metrics
%.metrics.bam: %.bam
	-ln -P $< $@;
	if [ ! -e $@ ]; then cp $< $@; fi; # if ln des not work

%.for_metrics_bed: %.bam %.bam.bai %.metrics.from_manifest.intervals %.metrics.bed %.metrics.region_clipped.bed %.bam.bed
	echo "# Create Metrics BED... "
	if [ "`grep ^ -c $*.metrics.bed`" == "0" ]; then \
		echo "# No BED file provided... BED from the BAM file"; \
		cp $*.bam.bed $@; \
	else \
		echo "# BED file provided..."; \
		echo "# Extract BAM Header"; \
		$(SAMTOOLS) view -H $< > $@; \
		if [ -e $*.metrics.region_clipped.bed ] && [ "`grep ^ -c $*.metrics.region_clipped.bed`" != "0" ]; then \
			echo "# BED file from clipped BED"; \
			grep -v ^@ $*.metrics.region_clipped.bed >> $@; \
		else \
			if [ -e $*.metrics.bed ]; then \
				echo "# BED FILE from BED"; \
				grep -v ^@ $*.metrics.bed >> $@; \
			fi; \
		fi; \
	fi; \
	if [ ! -e $@ ]; then touch $@; fi


# BED without header for metrics

%.withoutheader.for_metrics_bed: %.for_metrics_bed
	if [ -s $< ]; then \
		grep -v ^@ $< > $@; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi

# BED 3 fields for Picards

%.3fields.for_metrics_bed: %.for_metrics_bed #%.withoutheader.for_metrics_bed
	if [ -s $< ]; then \
		grep ^@ $< > $@; \
		grep -v ^@ $< > $*.withoutheader.for_metrics_bed.3fields.bed; \
		cat $*.withoutheader.for_metrics_bed.3fields.bed | awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t+\t"$$1":"$$2"-"$$3}' >> $@; \
		rm $*.withoutheader.for_metrics_bed.3fields.bed; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi


## METRICS
############

# ALL METRICS

%.bam.metrics/metrics: %.bam.metrics/metrics.gatk %.bam.metrics/metrics.picard %.bam.metrics/metrics.samtools %.bam.metrics/metrics.genes %.bam.metrics/metrics.amplicon_coverage %.bam.metrics/metrics.bam_check
	#Create directory
	mkdir -p $(@D)
	cat $^ > $@
	echo "#[$$(date)] BAM Metrics done" >> $@
	-rm -f $*.for_metrics_bed


# FATBAM Metrics
%.bam.metrics/metrics.amplicon_coverage: %.bam %.bam.bai %.manifest %.genome
	mkdir -p $(@D) ;
	-+$(FATBAM_COVERAGE) --env=$(CONFIG_TOOLS) --ref=`cat $*.genome`  --bam=$< --output=$(@D)/$(*F).amplicon_coverage --manifest=$*.manifest --multithreading --threads=$(THREADS) -v --tmp=$(FATBAM_TMP_FOLDER);
	echo "#[$$(date)] BAM Amplicon Coverage Metrics done" > $@;

# GATK METRICS

GATKDOC_FLAGS= -rf BadCigar -allowPotentiallyMisencodedQuals
%.bam.metrics/metrics.gatk: %.bam %.bam.bai %.genome %.for_metrics_bed %.3fields.for_metrics_bed
	# TODO: speed up ! Too loog for exome/genome...
	# use split algorithm with makefile and $(SAMTOOLS) view $< -b $$chr | $(BEDTOOLS)/genomeCoverageBed -ibam stdin -bg
	# see rule %.bam.bed
	#if (($(BAM_METRICS))); then \
	#if (($(BAM_METRICS)) && (0)); then \
	#if (($(BAM_METRICS))) && ((0)); then \
	if (($(BAM_METRICS))) && ((1)); then \
		#%.withoutheader.for_metrics_bed ; \
		#Create directory ; \
		mkdir -p $(@D) ; \
		grep -v ^@ $*.for_metrics_bed > $*.withoutheader.for_metrics_bed.gatk.bed ; \
		# GATK DepthOfCoverage needs BED without HEADER!!! ; \
		if [ ! -e $(@D)/$(*F) ]; then \
			$(JAVA) $(JAVA_FLAGS) -jar $(GATK) $(GATKDOC_FLAGS) \
				-T DepthOfCoverage \
				-R `cat $*.genome` \
				-o $(@D)/$(*F) \
				-I $< \
				-L $*.withoutheader.for_metrics_bed.gatk.bed; \
		fi; \
		rm $*.withoutheader.for_metrics_bed.gatk.bed; \
		echo "#[$$(date)] BAM GATK Metrics done" > $@; \
	else \
		#Create directory ; \
		mkdir -p $(@D) ; \
		echo "#[$$(date)] BAM GATK not done because BAM_METRICS=0" > $@; \
	fi;

# PICARD Metrics (DO NOT WORK, problem with bait file)

%.empty.HsMetrics:
	echo "BAIT_SET	GENOME_SIZE	BAIT_TERRITORY	TARGET_TERRITORY	BAIT_DESIGN_EFFICIENCY	TOTAL_READS	PF_READS	PF_UNIQUE_READS	PCT_PF_READS	PCT_PF_UQ_READS	PF_UQ_READS_ALIGNED	PCT_PF_UQ_READS_ALIGNED	PF_UQ_BASES_ALIGNED	ON_BAIT_BASES	NEAR_BAIT_BASES	OFF_BAIT_BASES	ON_TARGET_BASES	PCT_SELECTED_BASES	PCT_OFF_BAIT	ON_BAIT_VS_SELECTED	MEAN_BAIT_COVERAGE	MEAN_TARGET_COVERAGE	PCT_USABLE_BASES_ON_BAIT	PCT_USABLE_BASES_ON_TARGET	FOLD_ENRICHMENT	ZERO_CVG_TARGETS_PCT	FOLD_80_BASE_PENALTY	PCT_TARGET_BASES_2X	PCT_TARGET_BASES_10X	PCT_TARGET_BASES_20X	PCT_TARGET_BASES_30X	PCT_TARGET_BASES_40X	PCT_TARGET_BASES_50X	PCT_TARGET_BASES_100X	HS_LIBRARY_SIZE	HS_PENALTY_10X	HS_PENALTY_20X	HS_PENALTY_30X	HS_PENALTY_40X	HS_PENALTY_50X	HS_PENALTY_100X	AT_DROPOUT	GC_DROPOUT	SAMPLE	LIBRARY	READ_GROUP" > $@
	echo $$(echo $$(basename $@) | awk -F"." '{print $$1}')"	?	?	?	?	0	0	0	?	?	0	?	0	0	0	0	0	?	?	?	0	?	?	?	?	1	?	0	0	0	0	0	0	0	?	0	0	0	0	0	0	0	0" >> $@


%.bam.metrics/metrics.picard: %.bam %.bam.bai %.for_metrics_bed %.3fields.for_metrics_bed %.empty.HsMetrics %.genome
	#%.withoutheader.for_metrics_bed
	#JAVA_MEMORY_BY_SAMPLE ???
	# METRIC_ACCUMULATION_LEVEL=ALL_READS VALIDATION_STRINGENCY=LENIENT
	#Create directory
	mkdir -p $(@D)
	touch $@
	# Picard Metrics needs BED with Header and 3+3 fields!!!
	-if [ -s $*.3fields.for_metrics_bed ]; then \
		#$(JAVA) $(JAVA_FLAGS) -jar $(PICARDLIB)/CalculateHsMetrics.jar INPUT=$< OUTPUT=$(@D)/$(*F).HsMetrics BAIT_INTERVALS=$*.3fields.for_metrics_bed TARGET_INTERVALS=$*.3fields.for_metrics_bed VALIDATION_STRINGENCY=LENIENT 2>$(@D)/$(*F).HsMetrics.err; \
		#$(JAVA) $(JAVA_FLAGS_BY_SAMPLE) -jar $(PICARD) CollectHsMetrics INPUT=$< OUTPUT=$(@D)/$(*F).HsMetrics R=`cat $*.genome` BAIT_INTERVALS=$*.3fields.for_metrics_bed TARGET_INTERVALS=$*.3fields.for_metrics_bed PER_TARGET_COVERAGE=$(@D)/$(*F).HsMetrics.per_target_coverage 2>$(@D)/$(*F).HsMetrics.err VALIDATION_STRINGENCY=LENIENT; \
		$(JAVA) $(JAVA_FLAGS_BY_SAMPLE) -jar $(PICARD) CollectHsMetrics INPUT=$< OUTPUT=$(@D)/$(*F).HsMetrics R=`cat $*.genome` BAIT_INTERVALS=$*.for_metrics_bed TARGET_INTERVALS=$*.for_metrics_bed PER_TARGET_COVERAGE=$(@D)/$(*F).HsMetrics.per_target_coverage 2>$(@D)/$(*F).HsMetrics.err VALIDATION_STRINGENCY=LENIENT; \
	fi;
	# VALIDATION_STRINGENCY=LENIENT
	if [ ! -s $(@D)/$(*F).HsMetrics ]; then cp $*.empty.HsMetrics $(@D)/$(*F).HsMetrics; echo "#[ERROR] BAM PICARD Metrics failed. Empty HsMetrics file generated. See '$(@D)/$(*F).HsMetrics.err'" >> $@; fi
	echo "#[$$(date)] BAM PICARD Metrics done" >> $@


# SAMTOOLS METRICS
%.bam.metrics/metrics.samtools: %.bam %.bam.bai %.genome %.for_metrics_bed %.withoutheader.for_metrics_bed %.3fields.for_metrics_bed
	#Create directory ;
	mkdir -p $(@D);
	# Mandatory metrics
	# SAMTOOLS depthbed. mandatory for base coverage
	if [ -s $*.withoutheader.for_metrics_bed ]; then \
		#$(SAMTOOLS) depth -b $*.withoutheader.for_metrics_bed $< > $(@D)/$(*F).depthbed; \
		# CLEAN BAM \
		$(SAMTOOLS) view -F 1024 -F 4 -q 10 -h $<  -1 -@ $(THREADS) > $<.cleaned.bam ; \
		# DEPTHBED ON \
		$(SAMTOOLS) depth -b $*.withoutheader.for_metrics_bed $<.cleaned.bam > $(@D)/$(*F).depthbed; \
		# DEPTHBED OFF \
		#$(BEDTOOLS)/intersectBed -abam $<.cleaned.bam -b $*.withoutheader.for_metrics_bed -v | $(SAMTOOLS) depth - > $(@D)/$(*F).off.depthbed; \
		# ON NBReads \
		$(SAMTOOLS) view -c $<.cleaned.bam -L $*.withoutheader.for_metrics_bed > $(@D)/$(*F).on.nbreads; \
		# OFF NBReads \
		$(BEDTOOLS)/intersectBed -abam $<.cleaned.bam -b $*.withoutheader.for_metrics_bed -v | $(SAMTOOLS) view -c - > $(@D)/$(*F).off.nbreads; \
		rm -f $<.cleaned.bam; \
		#$(SAMTOOLS) view -F 1024 -q 10 -h $< | $(SAMTOOLS) depth -b $*.withoutheader.for_metrics_bed - > $(@D)/$(*F).depthbed; \
		#$(SAMTOOLS) view -F 1024 -q 10 -h $< -1 -@ $(THREADS) | $(BEDTOOLS)/intersectBed -abam - -b $*.withoutheader.for_metrics_bed -v | $(SAMTOOLS) depth - > $(@D)/$(*F).off.depthbed; \
		echo "#[$$(date)] COVERAGE with SAMTOOLS with BED file '$*.withoutheader.for_metrics_bed', without duplicates reads and reads with quality mapping < 10" >> $@; \
	else \
		echo "#[$$(date)] COVERAGE with SAMTOOLS with a BED file not done due to a lack of region defined in BED '$*.withoutheader.for_metrics_bed'" >> $@; \
	fi;
	if (($(BAM_METRICS))); then \
		#%.withoutheader.for_metrics_bed ; \
		#Create directory ; \
		#mkdir -p $(@D) ; \
		# SAMTOOLS Metrics ; \
		$(SAMTOOLS) flagstat $< > $(@D)/$(*F).flagstat ; \
		$(SAMTOOLS) idxstats $< > $(@D)/$(*F).idxstats ; \
		# Coverage with SAMTOOLS and BETTOOLS ; \
		if [ -s $*.withoutheader.for_metrics_bed ]; then \
			#$(SAMTOOLS) depth -b $*.withoutheader.for_metrics_bed $< > $(@D)/$(*F).depthbed; \
			$(SAMTOOLS) view -b $< -L $*.withoutheader.for_metrics_bed | $(BEDTOOLS)/genomeCoverageBed -ibam stdin -g `cat $*.genome` -dz > $(@D)/$(*F).genomeCoverageBedbed; \
		else \
			echo "#[$$(date)] COVERAGE with SAMTOOLS and BEDTOOLS not done due to a lack of region defined in BED '$*.withoutheader.for_metrics_bed'" >> $@; \
		fi; \
		if [ "$(FULL_COVERAGE)" == "1" ]; then \
			$(SAMTOOLS) depth $< > $(@D)/$(*F).depth; \
			$(GZ) --best -f $(@D)/$(*F).depth; \
			$(BEDTOOLS)/genomeCoverageBed -ibam $< -g `cat $*.genome` -dz > $(@D)/$(*F).genomeCoverageBed; \
			$(GZ) --best -f $(@D)/$(*F).genomeCoverageBed; \
		else \
			echo "#[$$(date)] FULL COVERAGE with SAMTOOLS and BEDTOOLS not done due to option disabled." >> $@; \
		fi; \
		echo "#[$$(date)] BAM SAMTOOLS/BEDTOOLS Metrics done" >> $@ ; \
	else \
		mkdir -p $(@D) ; \
		# SAMTOOLS Metrics ; \
		$(SAMTOOLS) flagstat $< > $(@D)/$(*F).flagstat ; \
		$(SAMTOOLS) idxstats $< > $(@D)/$(*F).idxstats ; \
		echo "#[$$(date)] BAM BEDTOOLS not done because BAM_METRICS=0" >> $@ ; \
	fi;


# FASTQC METRICS

%.fastqc/metrics: %.fastqc/metrics.fastqc %.fastqc/metrics.counts
	cat $< > $@


%.fastqc/metrics.fastqc: %.fastq.gz
	# create directory
	mkdir -p $(@D)
	# create link
	ln -s $< $(@D)/$(*F).unaligned.fastq.gz
	# touch target
	touch $@;
	# FASTQC
	-if (($$(zcat $< | head -n 1 | wc -l))); then \
		#$(FASTQC) $< --outdir=$(@D) --casava --extract; \
		$(FASTQC) $(@D)/$(*F).unaligned.fastq.gz --outdir=$(@D) --casava --extract --threads $(THREADS_BY_SAMPLE) ; \
		cp $(@D)/$(*F).unaligned_fastqc/fastqc_data.txt $*.fastqc/metrics.fastqc.txt; \
		echo "#FASTQC done. See 'metrics.fastqc.txt' file." >> $@; \
	else \
		echo "#FASTQC can't be launched. No Reads in the unaligned BAM file '$<'" >> $@; \
	fi;
	# create link
	-rm -f $(@D)/$(*F).unaligned.fastq.gz

%.fastqc/metrics.counts: %.fastq.gz #%.fastqc/metrics.fastqc
	# create directory
	mkdir -p $(@D)
	-if (($$(zcat $< | head -n 1 | wc -l))); then \
		# header \
		echo "total unique %unique maxRead count_maxRead %count_maxRead" > $*.counts.tmp; \
		# Counts \
		if [ $$(zcat $< | head -n 1 | wc -l) ]; then \
			zcat $< | awk '{unique=0} ((NR-2)%4==0){read=$$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}' >> $*.counts.tmp; \
		else \
			echo "0 0 - - - -" >> $*.counts.tmp; \
		fi; \
		# transposition \
		awk '{ for (i=1; i<=NF; i++)  {a[NR,i] = $$i} } NF>p { p = NF } END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++){ str=str" "a[i,j]; } print str } }' $*.counts.tmp | tr " " "\t" > $*.fastqc/metrics.counts.txt; \
		# Q30 \
		# script \
		# echo "" > $*.fastqc/metrics.counts.txt; \
		echo "#FASTQ Counts done. See 'metrics.counts.txt' file." > $@; \
	else \
		echo "#FASTQ Counts ERROR. No reads in '$<'..." > $@; \
	fi;
	-rm -f $*.counts.tmp

%.from_unaligned_bam.fastqc/metrics: %.unaligned.bam
	mkdir -p $(@D)
	touch $@;
	-if [ "$$($(SAMTOOLS) view $< | grep ^ -c)" == "0" ]; then \
		echo "#FASTQGC can't be launched. No Reads in the unaligned BAM file '$<'" >> $@; \
	else \
		#$(FASTQC) $< --outdir=$(@D) --casava --extract; \
		$(FASTQC) $< --outdir=$(@D) --casava --extract --threads $(THREADS) ; \
		echo "#FASTQGC done" >> $@; \
	fi;


## VCF METRICS
%.vcf.metrics/metrics: %.vcf.metrics/metrics.snpeff %.vcf.metrics/metrics.bcftools
	cat $< > $@

%.vcf.metrics/metrics.snpeff: %.vcf
	mkdir -p $(@D);
	touch $@;
	if (($(METRICS_SNPEFF))); then \
		+$(HOWARD) --input=$< --output=$@.vcf --snpeff_stats=$@.html --annotation=null --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS)  --force; \
		echo "#snpEff metrics done. See '$@.html' file." >> $@; \
	else \
		echo "# snpEff metrics NOT done." >> $@; \
	fi;


%.vcf.metrics/metrics.bcftools: %.vcf
	mkdir -p $(@D);
	touch $@;
	-$(BCFTOOLS) stats $< > $@.stats
	echo "#BCFTOOLS metrics done. See '$@.stats' file." >> $@;


#

# GENES COVERAGE METRICS

%.bam.metrics/metrics.genes: %.bam %.bam.bai
	# join -1 1 -2 5 <(sort test) <(sort -k5 /home1/TOOLS/db/RefSeq.hg19.bed) -o 2.1,2.2,2.3,2.5 | tr " " "\t"
	# grep  -P $GENES_LIST /home1/TOOLS/db/RefSeq.hg19.be
	mkdir -p $(@D);
	touch $@;
	+if (($(BAM_METRICS))) || [ -s $(BEDFILE_GENES) ] ; then \
		mkdir -p $(@D); \
		touch $@; \
		if [ -s `file=$$( echo $* | cut -d. -f1 ); echo "$$file.genes"` ] ; then \
			echo "BEDFILE_GENES from SAMPLE.genes"; \
			bedfile_genes_list=`file=$$( echo $* | cut -d. -f1 ); echo "$$file.genes"`; \
		elif [ ! -z $(BEDFILE_GENES) ]; then \
			echo "BEDFILE_GENES from BEDFILE_GENES variable "; \
			bedfile_genes_list=""; \
			for bedfile_genes in $$(echo $(BEDFILE_GENES)); \
			do \
				bedfile_name=$$( basename $$bedfile_genes | sed "s/\.genes$$//" ); \
				cp $$bedfile_genes $(@D)/$(*F)_$${bedfile_name}.bed; \
				bedfile_genes_list="$$bedfile_genes_list $(@D)/$(*F)_$${bedfile_name}.bed"; \
			done; \
		else \
			echo "BEDFILE_GENES generated from SAMPLE.manifest"; \
			bedfile_genes_list=`file=$$( echo $* | cut -d. -f1 ); echo "$$file.genes_from_manifest"`; \
			# ici on va faire une intersection avec le manifest et ce sera notre bed par defaut sil nexiste pas; \
			bed=`file=$$( echo $* | cut -d. -f1 ); echo "$$file.bed"`; \
			manifest=`file=$$( echo $* | cut -d. -f1 ); echo "$$file.manifest"`; \
			if [ -s `echo "$$bed"` ] ; then \
				#echo "genes from BED '$$bed'" >> $*.test; \
				cut -f1,2,3 $${bed} > $${manifest}.bed ; \
			elif [ -s `echo "$$manifest"` ] ; then \
				#echo "genes from BED '$$manifest'" >> $*.test; \
				# manifest to bed ; \
				$(FATBAM_ManifestToBED) --input "$$manifest" --output "$${manifest}.bed_tmp" --output_type "region" | $(FATBAM_ManifestToBED) --input "$$manifest" --output "$${manifest}.bed_tmp" --output_type "region" --type=PCR; \
				cut -f1,2,3 $${manifest}.bed_tmp > $${manifest}.bed ; \
				rm $${manifest}.bed_tmp ; \
			fi; \
			#echo "final genes file '$${manifest}.bed'" >> $*.test; \
			if [ -s `echo "$${manifest}.bed"` ] ; then \
				echo "MANIFEST.bed to bedfile_genes_list : $bedfile_genes_list"; \
				$(BEDTOOLS)/bedtools intersect -wb -a $${manifest}.bed -b $(REFSEQ_GENES) | cut -f8 | sort -u > $${manifest}.bed.intersect; \
				sort -k5 $(REFSEQ_GENES) > $${manifest}.bed.refseq; \
				join -1 1 -2 5 $${manifest}.bed.intersect $${manifest}.bed.refseq -o 2.1,2.2,2.3,2.5 | sort -u -k1,2 | tr " " "\t" | awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t+\t"$$4}' > $$bedfile_genes_list; \
				rm -f $${manifest}.bed.intersect $${manifest}.bed.refseq; \
			else \
				echo "#[$$(date)] BAM Metrics on Genes coverage not generate, No manifest found, ." >> $@; \
			fi;\
		fi; \
		mkdir -p $(@D); \
		touch $@; \
		echo "bed_file list is : `echo $$bedfile_genes_list` "; \
		#cat `echo $$bedfile_genes_list`; \
		if [ ! -z $$bedfile_genes_list ]; then \
			for bedfile_genes in $$(echo $$bedfile_genes_list); \
			do \
				if [ -e $$bedfile_genes ]; then \
					bedfile_name=$$( basename $$bedfile_genes | sed "s/\.genes$$//" ); \
					$(NGSscripts)/genesCoverage.sh -f $*.bam -b $$bedfile_genes -c "$(COVERAGE_CRITERIA)" --dp_fail=$(DP_FAIL) --dp_warn=$(DP_WARN) --dp_threshold=$(DP_THRESHOLD) -n $(NB_BASES_AROUND) -t $(BEDTOOLS)/bedtools -s $(SAMTOOLS) --threads=$(THREADS) -o $(@D)/$$bedfile_name; \
					echo "#[$$(date)] BAM Metrics on Genes coverage with $$bedfile_name bedfile done" >> $@; \
				else \
					echo "#[$$(date)] BAM Metrics on Genes coverage with $$bedfile_name bedfile FAILED" >> $@; \
				fi; \
			done; \
		else \
			echo "#[$$(date)] BAM Metrics on Genes coverage not generate because bedfile_genes variable is empty" >> $@; \
		fi; \
	else \
		echo "#[$$(date)] BAM Metrics on Genes coverage not done because BAM_METRICS=0" >> $@; \
	fi;
	if [ ! -e $@ ]; then echo "#[$$(date)] BAM Metrics on Genes coverage FAILED" > $@; fi;


	#$(NGSscripts)/genesCoverage.sh -f $*.bam -b $$bedfile_genes -c "$(COVERAGE_CRITERIA)" -n $(NB_BASES_AROUND) -t $(BEDTOOLS) -u $(BEDTOOLS2) -s $(SAMTOOLS) --threads=$(THREADS) -o $(@D)/$$bedfile_name; \
	echo "join -1 1 -2 5 $${manifest}.bed.intersect $${manifest}.bed.refseq -o 2.1,2.2,2.3,2.5 | sort -u -k1,2 | tr " " "\t" | awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t+\t"$$4}' > $$bedfile_genes_list;" ; \


%.bam.metrics/metrics.genes_OLD: %.bam %.bam.bai
	mkdir -p $(@D);
	touch $@;
	if (($(BAM_METRICS))); then \
		mkdir -p $(@D); \
		touch $@; \
		if [ -s `file=$$( echo $* | cut -d. -f1 ); echo "$$file.genes"` ] ; then \
			bedfile_genes_list=`file=$$( echo $* | cut -d. -f1 ); echo "$$file.genes"`; \
		elif [ ! -z $(BEDFILE_GENES) ]; then \
			bedfile_genes_list=""; \
			for bedfile_genes in $$(echo $(BEDFILE_GENES)); \
			do \
				bedfile_name=$$( basename $$bedfile_genes | sed "s/\.genes$$//" ); \
				cp $$bedfile_genes $(@D)/$(*F)_$${bedfile_name}.bed; \
				bedfile_genes_list="$$bedfile_genes_list $(@D)/$(*F)_$${bedfile_name}.bed"; \
			done; \
		else \
			bedfile_genes_list=`file=$$( echo $* | cut -d. -f1 ); echo "$$file.genes_from_manifest"`; \
			# ici on va faire une intersection avec le manifest et ce sera notre bed par defaut sil nexiste pas; \
			manifest=`file=$$( echo $* | cut -d. -f1 ); echo "$$file.manifest"`; \
			if [ -s `echo "$$manifest"` ] ; then \
				# manifest to bed ; \
				$(FATBAM_ManifestToBED) --input "$$manifest" --output "$${manifest}.bed_tmp" --output_type "region" | $(FATBAM_ManifestToBED) --input "$$manifest" --output "$${manifest}.bed_tmp" --output_type "region" --type=PCR; \
				cut -f1,2,3 $${manifest}.bed_tmp > $${manifest}.bed ; \
				rm $${manifest}.bed_tmp ; \
				# recover genes list from manifest in refseq file ; \
				genes_list=`$(BEDTOOLS)/bedtools intersect -wb -a $${manifest}.bed -b $(REFSEQ_GENES) | cut -f8 | sort -u `; \
				# if strand information : genes_list=`$(BEDTOOLS)/bedtools intersect -s -wb -a $${manifest}.bed -b $(REFSEQ_GENES) | cut -f8 | sort -u `; \
				touch $${bedfile_genes_list}.tmp; \
				# recover exon coordinates for all genes in refseq file and create bed file for genes coverage report ; \
				for gene in $$(echo $$genes_list); \
				do \
					awk -v pattern="$$gene" '$$5 == pattern {print $$0}' $(REFSEQ_GENES) >> $${bedfile_genes_list}.tmp; \
				done; \
				touch $${bedfile_genes_list}; \
				cut -f1,2,3,5 $${bedfile_genes_list}.tmp > $${bedfile_genes_list}; \
				rm $${bedfile_genes_list}.tmp ; \
			else \
				echo "No manifest found, genes coverage not generate." >> $@; \
			fi;\
		fi; \
		mkdir -p $(@D); \
		touch $@; \
		echo "bed_file list is : `echo $$bedfile_genes_list` "; \
		if [ ! -z $$bedfile_genes_list ]; then \
			for bedfile_genes in $$(echo $$bedfile_genes_list); \
			do \
				if [ -e $$bedfile_genes ]; then \
					bedfile_name=$$( basename $$bedfile_genes | sed "s/\.genes$$//" ); \
					$(NGSscripts)/genesCoverage.sh -f $*.bam -b $$bedfile_genes -c "$(COVERAGE_CRITERIA)" -n $(NB_BASES_AROUND) -t $(BEDTOOLS) -u $(BEDTOOLS2) -s $(SAMTOOLS) -o $(@D)/$$bedfile_name; \
					echo "#[$$(date)] Metrics on genes coverage with $$bedfile_name bedfile done" >> $@; \
				else \
					echo "#[$$(date)] Metrics on genes coverage with $$bedfile_name bedfile FAILED" >> $@; \
				fi; \
			done; \
		else \
			echo "# Genes coverage not generate because bedfile_genes variable is empty" >> $@; \
		fi; \
	else \
		echo "# Genes coverage not done because BAM_METRICS=0" >> $@; \
	fi;
	if [ ! -e $@ ]; then echo "#[$$(date)] Metrics on genes coverage FAILED" > $@; fi;



%.bam.metrics/metrics.bam_check: %.bam %.bam.bai #%.R1.fastq.gz %.R2.fastq.gz
	# CHECK NUMBER of READS in BAM \
	# Create directory if not created
	mkdir -p $(@D);
	# Pattern
	echo $$(dirname $*)/$$(basename $$(dirname $*)) > $@.sample_pattern;
	# Header
	echo "# Check number of reads between $< and $$(cat $@.sample_pattern).R1.fastq.gz $$(cat $@.sample_pattern).R2.fastq.gz to ensure no reads removing (only mark reads)" >$(@D)/$(*F).bam_check;
	echo "# (possible new reads marked as secondary aligned)" >>$(@D)/$(*F).bam_check;
	# Check
	if [ -s $$(cat $@.sample_pattern).R1.fastq.gz ] && [ -s $$(cat $@.sample_pattern).R2.fastq.gz ]; then \
		echo "# Files $$(cat $@.sample_pattern).R1.fastq.gz and $$(cat $@.sample_pattern).R2.fastq.gz exists" >>$(@D)/$(*F).bam_check; \
		# Compare number of reads \
		echo "# $$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $<) reads for $<" >>$(@D)/$(*F).bam_check; \
		echo "# $$(zcat $$(cat $@.sample_pattern).R1.fastq.gz $$(cat $@.sample_pattern).R2.fastq.gz | grep ^@ -c ) reads for $$(cat $@.sample_pattern).R1.fastq.gz $$(cat $@.sample_pattern).R2.fastq.gz" >>$(@D)/$(*F).bam_check; \
		grep ".* reads for $<" $(@D)/$(*F).bam_check | cut -d" " -f2; \
		grep ".* reads for $$(cat $@.sample_pattern).R1.fastq.gz $$(cat $@.sample_pattern).R2.fastq.gz" $(@D)/$(*F).bam_check | cut -d" " -f2; \
		if [ "$$(grep ".* reads for $<" $(@D)/$(*F).bam_check | cut -d" " -f2)" != "$$(grep ".* reads for $$(cat $@.sample_pattern).R1.fastq.gz $$(cat $@.sample_pattern).R2.fastq.gz" $(@D)/$(*F).bam_check | cut -d" " -f2)" ]; then \
			# Error \
			echo "#[ERROR] KO Number of reads between $< and $$(cat $@.sample_pattern).R1.fastq.gz $$(cat $@.sample_pattern).R2.fastq.gz !!!" >>$(@D)/$(*F).bam_check; \
			echo "# BCFTOOLS STATS for $<" >>$(@D)/$(*F).bam_check;  \
			$(SAMTOOLS) stats $< | grep ^SN >>$(@D)/$(*F).bam_check; \
			$(SAMTOOLS) idxstats $< >>$(@D)/$(*F).bam_check; \
			#exit 1; \
		else \
			# OK \
			echo "# Number of reads OK between $< and $$(cat $@.sample_pattern).R1.fastq.gz $$(cat $@.sample_pattern).R2.fastq.gz"  >>$(@D)/$(*F).bam_check; \
		fi; \
	else \
		# Error \
		echo "#[ERROR] No $$(cat $@.sample_pattern).R1.fastq.gz or $$(cat $@.sample_pattern).R2.fastq.gz" >>$(@D)/$(*F).bam_check; \
	fi;
	cat $(@D)/$(*F).bam_check;
	echo "#[$$(date)] BAM Check done" > $@ ;
	#-rm -f $@.sample_pattern



%.bam.metrics/metrics.bam_check_fromUBAM: %.bam %.bam.bai #$*.unaligned.bam $*.unaligned.bam.bai
	# CHECK NUMBER of READS in BAM \
	# Create directory if not created
	mkdir -p $(@D);
	# Pattern
	echo $$(dirname $*)/$$(basename $$(dirname $*)) > $@.sample_pattern;
	# Header
	echo "# Check number of reads between $< and $$(cat $@.sample_pattern).unaligned.bam to ensure no reads removing (only mark reads)" >$(@D)/$(*F).bam_check;
	echo "# (possible new reads marked as secondary aligned)" >>$(@D)/$(*F).bam_check;
	# Check
	if [ -s $$(cat $@.sample_pattern).unaligned.bam ] && [ -s $$(cat $@.sample_pattern).unaligned.bam.bai ]; then \
		echo "# Files $$(cat $@.sample_pattern).unaligned.bam and $$(cat $@.sample_pattern).unaligned.bam.bai exists" >>$(@D)/$(*F).bam_check; \
		# Compare number of reads \
		#grep "Number of reads for /home1/IRC/DATA/DEV/RES/ALL/160923_M01656_0124_000000000-ATJUM/T_Horizon/T_Horizon.unaligned.bam" /home1/IRC/DATA/DEV/RES/ALL/160923_M01656_0124_000000000-ATJUM/T_Horizon/T_Horizon.bwamem.bam.metrics/T_Horizon.bwamem.bam_check | cut -d" " -f2 \
		echo "# $$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $<) reads for $<" >>$(@D)/$(*F).bam_check; \
		echo "# $$($(SAMTOOLS) view -c -F 0x0100 -@ $(THREADS_SAMTOOLS) $$(cat $@.sample_pattern).unaligned.bam) reads for $$(cat $@.sample_pattern).unaligned.bam" >>$(@D)/$(*F).bam_check; \
		grep ".* reads for $<" $(@D)/$(*F).bam_check | cut -d" " -f2; \
		grep ".* reads for $$(cat $@.sample_pattern).unaligned.bam" $(@D)/$(*F).bam_check | cut -d" " -f2; \
		#grep "# Number of reads for $<: " $(@D)/$(*F).bam_check | sed "s/^# Number of reads for $$(cat $@.sample_pattern).unaligned.bam): \(.*\)$/\1/"; \
		#grep "# Number of reads for $$(cat $@.sample_pattern).unaligned.bam): " $(@D)/$(*F).bam_check | sed "s/^# Number of reads for $$(cat $@.sample_pattern).unaligned.bam): \(.*\)$/\1/"; \
		#echo $$(grep "Number of reads for $$(cat $@.sample_pattern).unaligned.bam)" $(@D)/$(*F).bam_check) | sed "s/^.*Number of reads for $$(cat $@.sample_pattern).unaligned.bam): \(.*\)/\1/"); \
		#echo $$(grep "Number of reads for $<: " $(@D)/$(*F).bam_check | sed "s/^.*Number of reads for $<\(.*\)/\1/)"; \
		if [ "$$(grep ".* reads for $<" $(@D)/$(*F).bam_check | cut -d" " -f2)" != "$$(grep ".* reads for $$(cat $@.sample_pattern).unaligned.bam" $(@D)/$(*F).bam_check | cut -d" " -f2)" ]; then \
			# Error \
			echo "#[ERROR] KO Number of reads between $< and $$(cat $@.sample_pattern).unaligned.bam !!!" >>$(@D)/$(*F).bam_check; \
			echo "# BCFTOOLS STATS for $<" >>$(@D)/$(*F).bam_check;  \
			$(SAMTOOLS) stats $< | grep ^SN >>$(@D)/$(*F).bam_check; \
			$(SAMTOOLS) idxstats $< >>$(@D)/$(*F).bam_check; \
			echo "# BCFTOOLS STATS for $$(cat $@.sample_pattern).unaligned.bam" >>$(@D)/$(*F).bam_check;  \
			$(SAMTOOLS) stats $$(cat $@.sample_pattern).unaligned.bam | grep ^SN >>$(@D)/$(*F).bam_check; \
			$(SAMTOOLS) idxstats $$(cat $@.sample_pattern).unaligned.bam >>$(@D)/$(*F).bam_check; \
			#exit 1; \
		else \
			# OK \
			echo "# Number of reads OK between $< and $$(cat $@.sample_pattern).unaligned.bam"  >>$(@D)/$(*F).bam_check; \
		fi; \
	else \
		# Error \
		echo "#[ERROR] No $$(cat $@.sample_pattern).unaligned.bam or $$(cat $@.sample_pattern).unaligned.bam.bai" >>$(@D)/$(*F).bam_check; \
	fi;
	cat $(@D)/$(*F).bam_check;
	echo "#[$$(date)] BAM Check done" > $@ ;
	-rm $@.sample_pattern


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# MAIN RULES '$(MK_RELEASE)' : basicaly to manage VCF, BAM, FASTQ, METRICS, MANIFEST, INTERVAL, BED... using SAMTOOLS, FASTQC GATK, PICARD, FATBAM, BWA, TABIX, IGV..."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )
