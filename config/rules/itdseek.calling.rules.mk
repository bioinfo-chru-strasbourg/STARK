############################
# ITDSeek Calling Rules
# Release: 0.9 
# Date: 23/10/2017
# Author: Antony Le Bechec
############################



############
# ITDSEEK  #
############

# Minimum coverage for a variant called by ITDSeek
DPMIN_ITDSEEK=20
VAFMIN_ITDSEEK=0.00
	
%.itdseek.vcf: %.bam %.bam.bai %.empty.vcf %.genome
	
	# Generate VCF with ITDSeek
	-$(ITDSEEK) $< `cat $*.genome` $(SAMTOOLS) > $@.tmp;
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
		#echo '##FORMAT=<ID=AD,Number=.,Type=String,Description="Allelic depths for the ref and alt alleles in the order listed">' >> $@.tmp2; \
		echo '##FORMAT=<ID=VAF,Number=.,Type=Float,Description="VAF Variant Frequency, from ITD allele fraction">' >> $@.tmp2; \
		#Add contig \
		echo '##contig=<ID=chr13,assembly=hg19,length=115169878>' >> $@.tmp2; \
		echo '##reference=file://'`cat $*.genome` >> $@.tmp2; \
		# Add sample columns FORMAT and $SAMPLE \
		grep "^#CHROM" $@.tmp | sed "s/INFO$$/INFO\tFORMAT\t"`echo $$(basename $$(dirname $@))`"/" >> $@.tmp2; \
		#old_IFS=$$IFS; IFS=$$'\n'; for l in $$(grep -v "^#" $@.tmp); do VAF=$$(echo $$l | cut -d" " -f8 | cut -d";" -f5 | cut -d"=" -f2); AD=$$(echo $$l | cut -d" " -f8 | cut -d";" -f1 | cut -d"=" -f2); DP=$$(echo $$l | cut -d" " -f8 | cut -d";" -f1 | cut -d"=" -f2 | awk -F"," '{print $$1+$$2}'); echo -e $$l";ITD=1 GT:DP:AD:VAF 0/1:$$DP:$$AD:$$VAF" | tr " " "\t" >> $@.tmp2; done; IFS=$$old_IFS; \
		old_IFS=$$IFS; IFS=$$'\n'; for l in $$(grep -v "^#" $@.tmp); do VAF=$$(echo $$l | cut -d" " -f8 | cut -d";" -f5 | cut -d"=" -f2); AD=$$(echo $$l | cut -d" " -f8 | cut -d";" -f1 | cut -d"=" -f2); DP=$$(echo $$l | cut -d" " -f8 | cut -d";" -f1 | cut -d"=" -f2 | awk -F"," '{print $$1+$$2}'); echo -e $$l";ITD=1 GT:DP:VAF 0/1:$$DP:$$VAF" | tr " " "\t" >> $@.tmp2; done; IFS=$$old_IFS; \
		# merge norm allele \
		#$(BCFTOOLS) norm $@.tmp2 -m +any > $@.tmp3; \
		# Filter on DP \
		#$(BCFTOOLS) view -e " FORMAT/DP[*] < $(DPMIN_ITDSEEK) " -e " FORMAT/VAF[*] <= $(VAFMIN_ITDSEEK) "  $@.tmp2 > $@; \
		$(BCFTOOLS) view -e " FORMAT/DP[*] < $(DPMIN_ITDSEEK) && FORMAT/VAF[*] <= $(VAFMIN_ITDSEEK) "  $@.tmp2 > $@; \
	fi;
	# debug
	#-cp $@.tmp $@.debug.tmp
	#-cp $@.tmp2 $@.debug.tmp2
	#-cp $@.tmp3 $@.debug.tmp3
	#-cp $@ $@.debug
	# Cleaning
	rm -f $@.tmp*
	
	
	
	# filter on DP, cause ITDSeek is too relax
	#$(VCFUTILS) varFilter -d $(DPMIN_ITDSEEK) $@.tmp2 > $@;
	

%.itdseek_ok_old.vcf: %.bam %.bam.bai %.empty.vcf %.genome
	
	# Generate VCF with ITDSeek
	-$(ITDSEEK) $< `cat $*.genome` $(SAMTOOLS) > $@.tmp;
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
		#echo '##FORMAT=<ID=AD,Number=.,Type=String,Description="Allelic depths for the ref and alt alleles in the order listed">' >> $@.tmp2; \
		echo '##FORMAT=<ID=VAF,Number=.,Type=Float,Description="VAF Variant Frequency, from ITD allele fraction">' >> $@.tmp2; \
		#Add contig \
		echo '##contig=<ID=chr13,assembly=hg19,length=115169878>' >> $@.tmp2; \
		echo '##reference=file://'`cat $*.genome` >> $@.tmp2; \
		# Add sample columns FORMAT and $SAMPLE \
		grep "^#CHROM" $@.tmp | sed "s/INFO$$/INFO\tFORMAT\t"`echo $$(basename $$(dirname $@))`"/" >> $@.tmp2; \
		#old_IFS=$$IFS; IFS=$$'\n'; for l in $$(grep -v "^#" $@.tmp); do VAF=$$(echo $$l | cut -d" " -f8 | cut -d";" -f5 | cut -d"=" -f2); AD=$$(echo $$l | cut -d" " -f8 | cut -d";" -f1 | cut -d"=" -f2); DP=$$(echo $$l | cut -d" " -f8 | cut -d";" -f1 | cut -d"=" -f2 | awk -F"," '{print $$1+$$2}'); echo -e $$l";ITD=1 GT:DP:AD:VAF 0/1:$$DP:$$AD:$$VAF" | tr " " "\t" >> $@.tmp2; done; IFS=$$old_IFS; \
		old_IFS=$$IFS; IFS=$$'\n'; for l in $$(grep -v "^#" $@.tmp); do VAF=$$(echo $$l | cut -d" " -f8 | cut -d";" -f5 | cut -d"=" -f2); AD=$$(echo $$l | cut -d" " -f8 | cut -d";" -f1 | cut -d"=" -f2); DP=$$(echo $$l | cut -d" " -f8 | cut -d";" -f1 | cut -d"=" -f2 | awk -F"," '{print $$1+$$2}'); echo -e $$l";ITD=1 GT:DP:VAF 0/1:$$DP:$$VAF" | tr " " "\t" >> $@.tmp2; done; IFS=$$old_IFS; \
		# merge norm allele \
		$(BCFTOOLS) norm $@.tmp2 -m +any > $@.tmp3; \
		# Filter on DP \
		$(BCFTOOLS) view -i " FORMAT/DP[*] > $(DPMIN_ITDSEEK) " $@.tmp3 > $@; \
	fi;
	# debug
	cp $@.tmp $@.debug.tmp
	cp $@.tmp2 $@.debug.tmp2
	cp $@.tmp3 $@.debug.tmp3
	cp $@ $@.debug
	# Cleaning
	rm -f $@.tmp*
	

RELEASE_CMD := $(shell echo "\#\# ITDSeek: FLT3 ITD detection, generate *.itdseek.vcf files with with parameters: DPMIN_ITDSEEK='$(DPMIN_ITDSEEK)', VAFMIN_ITDSEEK='$(VAFMIN_ITDSEEK)'. VCF file is modified to include GT/Genotype information as 0/1 by default. Variants are excluded if DP lower than $(DPMIN_ITDSEEK) AND VAF lower than $(VAFMIN_ITDSEEK)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:itdseek:ITDSeek - FLT3 ITD detection:no parameters, GT info 0/1 added"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

