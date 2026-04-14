############################
# FLT3ITDEXT Calling Rules
# Release: 1
# Date: 01/01/2023
# Author: Thomas LAVAUX
############################

###############
# FLT3ITDEXT  #
###############

# --minreads, -mr  Minimum number of supporting reads to be included in VCF
# --genome, -g  Genome build (defaults to "hg19"; or can be "hg38")
# --ngstype, -n NGS platform type (defaults to "HC" [hybrid capture]; or can be "amplicon", "NEB", or "Archer")
FLT3INDEX?=/STARK/tools/FLT3_ITD_ext/current/bin/FLT3_dna_e14e15

%.flt3itdext$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf 
	$(FLT3ITDEXT) -b $< -mr 5 -g $(ASSEMBLY) --ngstype HC -i $(FLT3INDEX) -o $@.tmp1;
	# add genotype
	(grep "^##" $@.tmp1.vcf && echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' && grep "^#CHROM" $@.tmp1.vcf | grep "^#" -v $@.tmp.results.vcf | awk '{print $0"\tGT\t0/1"}' ) > $@.tmp2.vcf;
	# Filter 
	$(BCFTOOLS) view -e " FORMAT/DP[*] < $(DPMIN_FLT3ITDEXT) || FORMAT/RAF[*] <= $(RAFMIN_FLT3ITDEXT) "  $@.tmp2.vcf > $@;
	# Cleaning
	-rm -rf $@.tmp*;

RELEASE_CMD := $(shell echo "\#\# FLT3ITDEXT: FLT3 ITD detection, generate *.flt3itdext.vcf files with with parameters: DPMIN_flt3itdext='$(DPMIN_FLT3ITDEXT)', VAFMIN_flt3itdext='$(RAFMIN_FLT3ITDEXT)'. VCF file is modified to include GT/Genotype information as 0/1 by default. Variants are excluded if DP lower than $(DPMIN_FLT3ITDEXT) AND VAF lower than $(VAFMIN_FLT3ITDEXT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:iflt3itdext:FLT3ITDEXT - FLT3 ITD detection:no parameters, GT info 0/1 added"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
