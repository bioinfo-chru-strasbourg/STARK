############################
# FLT3ITDEXT Calling Rules
# Release: 1
# Date: 01/01/2023
# Author: 
############################

###############
# FLT3ITDEXT  #
###############

# --minreads, -mr  Minimum number of supporting reads to be included in VCF
# --genome, -g  Genome build (defaults to "hg19"; or can be "hg38")
# --ngstype, -n NGS platform type (defaults to "HC" [hybrid capture]; or can be "amplicon", "NEB", or "Archer")

%.flt3itdext$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf 
	$(FLT3ITDEXT) -b $< -mr 5 -g $ASSEMBLY --ngstype HC -o $@;
	if (( $$(grep -c ^ERROR $@.tmp) )); then \
		echo "[ERROR] ITDSeek failed: "; \
		grep ^ERROR $@.tmp; \
		cp $*.empty.vcf $@; \
	fi;
	# Filter on DP \
	$(BCFTOOLS) view -e " FORMAT/DP[*] < $(DPMIN_ITDSEEK) || FORMAT/VAF[*] <= $(VAFMIN_ITDSEEK) "  $@.tmp2 > $@; \
	# Cleaning
	rm -rf $@.tmp*




RELEASE_CMD := $(shell echo "\#\# FLT3ITDEXT: FLT3 ITD detection, generate *.flt3itdext.vcf files with with parameters: DPMIN_flt3itdext='$(DPMIN_flt3itdext)', VAFMIN_flt3itdext='$(VAFMIN_flt3itdext)'. VCF file is modified to include GT/Genotype information as 0/1 by default. Variants are excluded if DP lower than $(DPMIN_FLT3ITDEXT) AND VAF lower than $(VAFMIN_FLT3ITDEXT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:iflt3itdext:FLT3ITDEXT - FLT3 ITD detection:no parameters, GT info 0/1 added"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
