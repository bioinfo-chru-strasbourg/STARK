#!/bin/bash
#################################
##
## NGS environment
## Integration into Database
##
## version: 0.9
## date: 12/09/2014
## author: Antony Le Bechec
##
#################################

# Usage : SELF <VCF_LIST> <VCF_OUTPUT> <VCF_FOLDER> <BGZIP> <TABIX>

# 1.1. Configuration
#####################

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_DIR/env.sh

echo "# "
echo "# NGS Folder:            "$NGS_FOLDER
echo "# NGS BIN Folder:        "$NGS_BIN
echo "# NGS SCRIPTS Folder:    "$NGS_SCRIPTS
echo "# MISEQ Folder:          "$MISEQ_FOLDER
echo "# DEMULTIPLEXING Folder: "$DEMULTIPLEXING_FOLDER
echo "# RESULTS Folder:        "$RESULTS_FOLDER
echo "# ANALYSIS Folder:       "$ANALYSIS_FOLDER
echo "# CONFIG ini:            "$CONFIG
echo "# "

# 1.2. INPUT
#############

VCF_LIST=$1
VCF_OUTPUT=$2
if [ "$VCF_OUTPUT" == "" ]; then
	VCF_OUTPUT="merge.vcf"
fi;
VCF_FOLDER=$3
if [ "$VCF_FOLDER" == "" ] || [ ! -d $VCF_FOLDER ]; then
	VCF_FOLDER=$VCFTOOLS_PATH
fi;
BGZIP=$4
if [ "$BGZIP" == "" ] || [ ! -e $BGZIP ]; then
	BGZIP=$BGZIP_BIN
fi;
TABIX=$5
if [ "$TABIX" == "" ] || [ ! -e $TABIX ]; then
	TABIX=$TABIX_BIN
fi;

echo "# Input"
echo "# VCF List:  "$VCF_LIST
echo "# VCF-TOOLS: "$VCF_FOLDER
echo "# BGZIP:     "$BGZIP
echo "# "

# BGZIP

echo "# SORT & BGZIP & TABIX"

VCF_LIST_GZ=""
for VCF in $VCF_LIST
do
	# SORT
	echo "$VCF_FOLDER/vcf-sort $VCF > $VCF.sorted.vcf"
	$VCF_FOLDER/vcf-sort $VCF > $VCF.sorted.vcf
	# BGZIP
	echo "$BGZIP $VCF.sorted.vcf -c > $VCF.gz"
	$BGZIP $VCF.sorted.vcf -c > $VCF.gz
	echo "$TABIX $VCF.gz"
	$TABIX $VCF.gz
	VCF_LIST_GZ=$VCF_LIST_GZ" "$VCF.gz
	#rm $VCF.unsorted.gz
done

echo "# Merge"
echo "$VCF_FOLDER/vcf-merge $VCF_LIST_GZ > $VCF_OUTPUT"
$VCF_FOLDER/vcf-merge $VCF_LIST_GZ > $VCF_OUTPUT

