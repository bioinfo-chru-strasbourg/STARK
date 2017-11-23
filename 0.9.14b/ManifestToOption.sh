#!/bin/bash
#################################
##
## NGS environment
## Manifest to option
##
RELEASE=0.9
DATE=20140717
AUTHOR="Antony Le Bechec"
##
#################################
# HOW to clip a aligned bam using a manifest
# $1: manifest

NGS=/NGS
NGS_SCRIPTS=$NGS/scripts
NGS_BIN=$NGS/bin
NGS_GENOMES=$NGS/genomes
TOOL_ManifestToBED=$NGS_SCRIPTS/FATBAM.ManifestToBED.pl
TOOL_BEDToOPTION=$NGS_SCRIPTS/FATBAM.BEDToOption.pl
TOOL_fastaFromBed=$NGS_BIN/bedtools/fastaFromBed
#FATBAM=/NGS/scripts/FATBAM.pl
REF=$NGS_GENOMES/hg19/hg19.fa

# INPUT
MANIFEST=$1			# /media/IRCV2/V2/RES/ALL/140703_M01656_0014_000000000-D02HD/P1335/P1335.manifest	
if [ "$MANIFEST" == "" ] || [ ! -e $MANIFEST ]; then echo "# ERROR: Manifest '$MANIFEST' doesn't exist "; exit 1; fi;
BED=$MANIFEST".bed"
FASTA=$MANIFEST".fasta"
OPTION=$MANIFEST".option"
rm $BED $OPTION

echo "################################################"
echo "# Illumina Manifest file to option for FATBAM "
echo "# Release: $RELEASE"
echo "# Date: $DATE"
echo "# Author: $AUTHOR"
echo "################################################"
echo "#"
echo "# Manifest file: $1"

# Create BED file
$TOOL_ManifestToBED --input=$MANIFEST --output=$BED --type=PCR #--debug
# Create FASTA File
#echo "$TOOL_fastaFromBed -fi $REF -bed $BED -fo $FASTA"
#$TOOL_fastaFromBed -fi $REF -bed $BED -fo $FASTA
# Create Option
#echo "$TOOL_BEDToOPTION --input=$BED --output=$OPTION"
$TOOL_BEDToOPTION --input=$BED --output=$OPTION #--debug


# OUTPUT
echo ""
echo "# BED"
cat $BED
echo ""
echo "# OPTION"
cat $OPTION
echo ""

