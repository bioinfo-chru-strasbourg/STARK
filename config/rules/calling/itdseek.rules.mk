############################
# ITDSeek Calling Rules
# Release: 1.0.0
# Date: 09/11/2026
# Author: Antony Le Bechec
############################

# Release note
# 1.0.0-09/11/2026: Initial release of ITDSeek calling rules.


############
# ITDSEEK  #
############

# ITDSEEK filters Flags
# This variable can be defined directly on Apps
ITDSEEK_FILTERS_FLAGS?= \
	--exclude-genotypes '' \
	--exclude 'FORMAT/VAF <= 0.00' \
	--include 'FORMAT/DP>=20' \
	--filter '' \
	--sort


%.itdseek$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf 
	# Generate VCF with ITDSeek
	-$(ITDSEEK) $< $(GENOME) $(SAMTOOLS) > $@.tmp;
	# Check for errors in ITDSeek output
	if (( $$(grep -c ^ERROR $@.tmp) )); then \
		echo "[ERROR] ITDSeek failed: "; \
		grep ^ERROR $@.tmp; \
		cp $*.empty.vcf $@; \
	else \
		# Add header \
		grep "^##" $@.tmp > $@.tmp2; \
		echo '##INFO=<ID=ITD,Number=.,Type=String,Description="ITD detection">' >> $@.tmp2; \
		echo '##FORMAT=<ID=GT,Number=.,Type=String,Description="Genotype">' >> $@.tmp2; \
		echo '##FORMAT=<ID=DP,Number=.,Type=Integer,Description="Read depth">' >> $@.tmp2; \
		echo '##FORMAT=<ID=VAF,Number=.,Type=Float,Description="VAF Variant Frequency, from ITD allele fraction">' >> $@.tmp2; \
		# Add contig \
		echo '##contig=<ID=chr13,assembly=hg19,length=115169878>' >> $@.tmp2; \
		echo '##reference=file://'$(GENOME) >> $@.tmp2; \
		# Add sample columns FORMAT and SAMPLE \
		grep "^#CHROM" $@.tmp | sed "s/INFO$$/INFO\tFORMAT\t"`echo $$(basename $$(dirname $@))`"/" >> $@.tmp2; \
		mkdir -p $@.tmp"_VCF_sort"; \
		touch $@.tmp3; \
		old_IFS=$$IFS; IFS=$$'\n'; for l in $$(grep -v "^#" $@.tmp); do VAF=$$(echo $$l | cut -d" " -f8 | cut -d";" -f5 | cut -d"=" -f2); AD=$$(echo $$l | cut -d" " -f8 | cut -d";" -f1 | cut -d"=" -f2); DP=$$(echo $$l | cut -d" " -f8 | cut -d";" -f1 | cut -d"=" -f2 | awk -F"," '{print $$1+$$2}'); echo -e $$l";ITD=1 GT:DP:VAF 0/1:$$DP:$$VAF" | tr " " "\t" >> $@.tmp3; done; IFS=$$old_IFS; \
		cat $@.tmp2 > $@.tmp2.vcf; \
		cat $@.tmp3 | sort -k1,1V -k2,2n -T $@.tmp"_VCF_sort" >> $@.tmp2.vcf; \
	fi;
	# Filters
	bash $(STARK_FOLDER_BIN)/filter_vcf.sh --input $@.tmp2.vcf --output $@ --threads $(THREADS_BY_CALLER) --bcftools $(BCFTOOLS) $(ITDSEEK_FILTERS_FLAGS);
	# Cleaning
	rm -rf $@.tmp*


# ITDSEEK Release documentation

RELEASE_COMMENT := "\#\# CALLING ITDSeek: ITDSeek tool identify FLT3 ITD variants from aligned BAM"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# CALLING: ITDSeek FLT3 ITD detection, generate *.itdseek.vcf files with with parameters. VCF file is modified to include GT/Genotype information as 0/1 by default, and DP and VAF. Parameters ITDSEEK_FILTERS_FLAGS='$(ITDSEEK_FILTERS_FLAGS)' "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:itdseek:ITDSeek - FLT3 ITD detection:Filter on DP and VAF, GT info 0/1 added, ITDSEEK_FILTERS_FLAGS='$(ITDSEEK_FILTERS_FLAGS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
