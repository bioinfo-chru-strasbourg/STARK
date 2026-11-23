############################
# Main Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE=1.0.0
MK_DATE=17/11/2026

# Release note
# 17/11/2026-1.0.0: Initial release of BAM rules for STARK version 1.0.0 (extracted from main.rules.mk)


# Parameters
MPILEUP_OPTIONS?= -E -d 10000000 -C 50 -Q 10 -q 1
CAP_SOFTCLIPTOQ0?=$(CAP_FOLDER)/CAP.SoftClipToQ0.pl

## BAM/FASTQ Files
####################

# BAM Indexing - BAI-format
%.bam.bai: %.bam
	$(SAMTOOLS) index $< -@ $(THREADS_SAMTOOLS)


# BAM from SAM
# sorting sam file by coordinate and output a bam file
%.bam : %.sam
	if ((! $$($(SAMTOOLS) view -H $< | grep "^@HD.*VN:.*SO:" | wc -l))) || (($$($(SAMTOOLS) view -H $< | grep "^@HD.*VN:.*SO:unsorted" | wc -l))) ; then \
		$(SAMTOOLS) sort $< -o $@ -O BAM -l 1 -T $<.SAMTOOLS_PREFIX -@ $(THREADS_SAMTOOLS); \
	else \
		$(SAMTOOLS) view -o $@ -b -1 -S -T $(GENOME) $< -@ $(THREADS_SAMTOOLS); \
	fi;
	-if [ $(REMOVE_INTERMEDIATE_SAM) -eq 1 ]; then rm -f $<; fi;


# CRAM from BAM
%.cram: %.bam
	echo "test BAM to CRAM: $^"
	$(SAMTOOLS) view -o $@ -O CRAM -S -T $(GENOME) $*.bam -@ $(THREADS_SAMTOOLS);

# BAM Indexing
%.cram.crai: %.cram
	$(SAMTOOLS) index $<


# BAM Compress
%.bam: %.compress.bam
	# Check if sorting is needed
	if ((! $$($(SAMTOOLS) view -H $< | grep "^@HD.*VN:.*SO:" | wc -l))) || (($$($(SAMTOOLS) view -H $< | grep "^@HD.*VN:.*SO:unsorted" | wc -l))) ; then \
		$(SAMTOOLS) sort $< -o $@ -l $(BAM_COMPRESSION) -T $<.SAMTOOLS_PREFIX -@ $(THREADS_SAMTOOLS); \
	else \
		$(SAMTOOLS) view $< -o $@ -O BAM,level=$(BAM_COMPRESSION) -@ $(THREADS_SAMTOOLS); \
	fi;
	rm -rf $*.compress.bam $*.compress.bai;


# BAM Sorting
%.bam : %.sorting.bam
	if ((! $$($(SAMTOOLS) view -H $< | grep "^@HD.*VN:.*SO:" | wc -l))) || (($$($(SAMTOOLS) view -H $< | grep "^@HD.*VN:.*SO:unsorted" | wc -l))) ; then \
		$(SAMTOOLS) sort $< -o $@ -l 1 -T $<.SAMTOOLS_PREFIX -@ $(THREADS_SAMTOOLS); \
	else \
		mv $< $@; \
	fi;
	rm -f $<;

# BAM Remove secondary reads
%.bam : %.removesecondary.bam
	$(SAMTOOLS) view -Fx100 $< -o $@ -O BAM,level=$(BAM_COMPRESSION) -@ $(THREADS_SAMTOOLS); \
	rm -f $< $<.bai;


%.bam.chroms: %.bam %.bam.bai %.design.bed
	# Generate the list of chromosomes to be used for variant calling based on the BAM file and design BED file
	# This ensures that only chromosomes present in both the BAM file and the design BED file are considered for variant calling
	# If no region are specified in the design BED file, all chromosomes with reads in the BAM file will be considered.
	# The resulting file will list chromosomes ordered by the number of reads for each chromosome (descending order).
	# Generate idxstats for the BAM file
	$(SAMTOOLS) idxstats $< > $@.tmp.idxstats
	# Extract chromosomes with reads from idxstats
	awk -F'\t' '$$1!="*" && $$3+$$4>0 {print $$1"\t"$$3+$$4}' $@.tmp.idxstats > $@.tmp.idxstats_chroms
	# Chromosomes of the design (empty file = empty design)
	awk -F'\t' '!/^#/ && !/^(track|browser)/ && NF{print $$1}' $*.design.bed | sort -u > $@.tmp.design_chroms
	# Chromosomes to call (BAM header order): design chromosomes with reads, or all chromosomes with reads if the design is empty
	if [ -s $@.tmp.design_chroms ]; then \
		cut -f1 $@.tmp.idxstats_chroms | grep -Fx -f $@.tmp.design_chroms > $@.tmp.chroms || { rc=$$?; [ $$rc -eq 1 ] || exit $$rc; }; \
	else \
		cp $@.tmp.idxstats_chroms $@.tmp.chroms; \
	fi;
	# Execution order: most reads first (merge order stays chroms.txt)
	awk -F'\t' 'NR==FNR{k[$$1];next} $$1 in k{print $$1"\t"$$3+$$4}' $@.tmp.chroms $@.tmp.idxstats \
		| sort -k2,2nr > $@
	# Clean
	rm -f $@.tmp*


# FASTQ Compression with GZIP
%.fastq.gz: %.fastq
	$(GZ) --fast $< -c > $@;


# FASTQ(s) from BAM
%.R1.fastq %.R2.fastq: %.bam
	$(JAVA) -jar $(PICARD) SamToFastq -INPUT $< -FASTQ $*.R1.fastq -SECOND_END_FASTQ $*.R2.fastq;


# BAM reduction
GATKRR_FLAGS=
%.reduced.bam: %.bam %.bam.bai
	$(JAVA8) $(JAVA_FLAGS) -jar $(GATK3) $(GATKRR_FLAGS) -T ReduceReads -R $(GENOME) -I $< -o $@;



# SAMTOOLS mpileup
####################

# Rule: Generate mpileup from BAM file and design BED file
%.bam.mpileup: %.bam %.bam.bai %.design.bed
	# Clean and create split directory
	rm -rf $@.split && mkdir -p $@.split
	# Generate idxstats for the BAM file (contigs with reads only)
	$(SAMTOOLS) idxstats $< | awk '$$3>0' > $@.split/idxstats.txt
	# List of chromosomes (BAM header order), restricted to the BED if it has data lines
	if (( $$(grep -v -e "^#" -e "^$$" -c $*.design.bed) )); then \
		grep -v -e "^#" -e "^$$" $*.design.bed | cut -f1 | sort -u > $@.split/bedchroms.txt; \
		cut -f1 $@.split/idxstats.txt | grep -Fx -f $@.split/bedchroms.txt > $@.split/chroms.txt || echo "No chromosomes"; \
	else \
		cut -f1 $@.split/idxstats.txt > $@.split/chroms.txt; \
	fi;
	# Execution order: most reads first (concatenation order stays chroms.txt)
	awk 'NR==FNR{k[$$1];next} $$1 in k{print $$3"\t"$$1}' $@.split/chroms.txt $@.split/idxstats.txt \
		| sort -k1,1nr | cut -f2 > $@.split/order.txt;
	# Generate the sub Makefile (flags are defined in the same shell that writes it)
	{ \
	VIEW_L=""; MP_L=""; \
	if [ -s $@.split/bedchroms.txt ]; then VIEW_L="-L $*.design.bed"; MP_L="-l $*.design.bed"; fi; \
	printf '.DELETE_ON_ERROR:\n'; \
	printf 'all: %s\n' "$$(sed 's|.*|$@.split/&.part|' $@.split/order.txt | tr '\n' ' ')"; \
	printf '$@.split/%%.part:\n'; \
	printf '\t$(SAMTOOLS) view $< "$$*" %s -h | perl $(CAP_SOFTCLIPTOQ0) -v1 | $(SAMTOOLS) mpileup - -f $(GENOME) %s $(MPILEUP_OPTIONS) > $$@\n' "$$VIEW_L" "$$MP_L"; \
	} > $@.split/sub.mk
	# Process all chromosome parts in parallel
	+make -f $@.split/sub.mk all
	# Concatenate the parts in BAM order
	sed 's|.*|$@.split/&.part|' $@.split/chroms.txt | xargs -r cat > $@
	rm -rf $@.split


# SoftClip to Quality 0
#########################

# Change SoftClip to Quality 0
%.softclippedtoq0.sam: %.bam %.bam.bai %.manifest
	# Clipping
	$(SAMTOOLS) view $< -h | perl $(CAP_SOFTCLIPTOQ0) -v0 > $@
