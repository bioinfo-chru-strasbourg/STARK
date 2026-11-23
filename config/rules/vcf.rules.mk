############################
# Main Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE=1.0.0
MK_DATE=17/11/2026

# Release note
# 17/11/2026-1.0.0: Initial release of VCF rules for STARK (extract from main.rules.mk)


## VCF files
##############

# VCF Compression
%.vcf.gz: %.vcf %.empty.vcf
	# If no VCF or empty file, create an empty VCF
	if [ ! -s $< ]; then cp $*.empty.vcf $<; fi;
	# VCF Sorting
	mkdir -p $@.SAMTOOLS_PREFIX
	$(BCFTOOLS) sort -T $@.SAMTOOLS_PREFIX $< > $@.sorted
	rm -rf $@.SAMTOOLS_PREFIX
	# VCF compression with BGZIP
	$(BGZIP) -f $@.sorted -c > $@
	rm -f $@.sorted
	# remove files
	-rm -f $*.empty.vcf*


# VCF Indexing with TABIX
%.vcf.gz.tbi: %.vcf.gz
	# Remove index if exists
	-rm -f $@
	# Indexing with TABIX
	$(TABIX) -f -p vcf $<


# VCF Indexing with IGVTOOLS
%.vcf.idx: %.vcf %.empty.vcf
	# If no VCF or empty file, create an empty VCF
	-if [ ! -s $< ]; then cp $*.empty.vcf $<; fi;
	# Indexing with IGVTOOLS
	$(IGVTOOLS) index $<
	# Empty index if fail
	if [ ! -e $@ ]; then touch $@; fi;
	# remove files
	-rm -f $*.empty.vcf*


# EMPTY VCF
%.empty.vcf:
	-mkdir -p $(@D)
	# Header
	echo "##fileformat=VCFv4.1" > $@
	echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> $@
	# Head first line
	echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	"$$(echo $$(basename $@) | awk -F"." '{print $$1}') >> $@


# VCF SORTING
%.vcf: %.sorting.vcf
	# Force replace dictionary of contig in the VCF file
	$(JAVA) $(JAVA_FLAGS_GATK4) -jar $(GATK4) \
		UpdateVCFSequenceDictionary \
		--sequence-dictionary $(DICT) \
		--replace True \
		--variant $< \
		--output $@.with_contigs.vcf;
	# Sort the VCF file by chromosome and position
	$(BCFTOOLS) sort -T $@.with_contigs.SAMTOOLS_PREFIX $@.with_contigs.vcf -o $@
	rm -rf $@.with_contigs.*


# VCF NORMALIZATION with BCFTOOLS
%.vcf: %.normalization.vcf
	# Add or replace dictionary of contig in the VCF file
	if $(BCFTOOLS) view -h "$<" | grep -q '^##contig='; then \
		ln -s $< $@.with_contigs.vcf; \
	else \
		$(JAVA) $(JAVA_FLAGS_GATK4) -jar $(GATK4) \
			UpdateVCFSequenceDictionary \
			--sequence-dictionary $(DICT) \
			--replace True \
			--variant $< \
			--output $@.with_contigs.vcf; \
	fi;
	# Normalize the VCF file with BCFTOOLS
	$(BCFTOOLS) norm -m- --multi-overlaps 0 -f $(GENOME) $@.with_contigs.vcf -o $@
	rm -f $@.with_contigs.*


# MERGE SNP and InDel VCF
%.vcf: %.SNP.vcf %.InDel.vcf
	$(JAVA) $(JAVA_FLAGS_GATK4) -jar $(GATK4) \
		MergeVcfs \
		-I $*.SNP.vcf \
		-I $*.InDel.vcf \
		--CREATE_INDEX false \
		--SEQUENCE_DICTIONARY $(DICT) \
		-O $@;


# VCF to tab delimiter
%.tsv: %.vcf
	$(HOWARD) convert $(HOWARD_CONFIG_OPTIONS) --input=$< --output=$@;
	touch $@
	if [ ! -e $@ ]; then touch $@; fi;
