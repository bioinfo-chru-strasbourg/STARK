#!/bin/bash
#################################
##
## NGS environment
## VCF manipulation for EQA
##
## version: 0.9
## date: 18/12/2014
## author: Antony Le Bechec
##
#################################

# USAGE: VCFminimal.sh <INPUT_VCF> <INPUT_BED> > <OUTPUT_VCF>

VCF_INPUT=$1 	# INPUT VCF
BED_INPUT=$2 	# INPUT BED

# HEADER
echo '##fileformat=VCFv4.0'
echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
echo '##FILTER=<ID=PASS,Description="All filters passed">'

# MINIMALIZATION
if [ "$BED_INPUT" != "" ]; then
	$VCFTOOLS/vcftools --vcf $VCF_INPUT --bed $BED_INPUT --recode -c | $BCFTOOLS view -i 'FILTER="PASS"' - 2>/dev/null | $BCFTOOLS annotate -x QUAL,INFO,^FORMAT/GT 2>/dev/null | $VCFTOOLS/vcf-sort -c | grep -v ^##
else
	$BCFTOOLS view -i 'FILTER="PASS"' $VCF_INPUT 2>/dev/null | $BCFTOOLS annotate -x QUAL,INFO,^FORMAT/GT 2>/dev/null | $VCFTOOLS/vcf-sort -c | grep -v ^##
fi;


