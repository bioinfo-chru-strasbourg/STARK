#!/bin/bash

DATA_FOLDER=/STARK/data/

TEST_NUM=TEST_21


STARK_CMD="STARK"

OUTLYZER=/STARK/tools/outlyzer/current/bin/outLyzer.py
OUTLYZER_NORM=/STARK/tools/stark/current/bin/outlyzer_norm.awk


SAMPLE=HORIZON2

DIR=/STARK/data/SAMPLE/$SAMPLE
BAM=$DIR/HORIZON_R1.bwamem.bam
BED=$DIR/HORIZON_R1.bed
REF=/STARK/databases/genomes/hg19/hg19.fa
OUTPUT=/STARK/data/test_outlyzer.dir/$SAMPLE
OUTPUT_TMP=/STARK/data/test_outlyzer.dir/$SAMPLE/tmp
VCF=$OUTPUT/$SAMPLE.outlyzer.vcf
VCF_ANN=$OUTPUT/$SAMPLE.outlyzer.howard.vcf

CORE=4

mkdir -p $OUTPUT $OUTPUT_TMP

if ((0)); then
	$OUTLYZER calling -bed $BED -bam $BAM -ref $REF -output $OUTPUT_TMP// -core $CORE -verbose 1 -cut $CORE
fi

VCF_OUTLYZER=$(ls $OUTPUT_TMP/*.vcf)

echo "VCF_OUTLYZER="$VCF_OUTLYZER
tail $VCF_OUTLYZER

#cp $VCF_OUTLYZER $VCF

#grep "^##" $VCF_OUTLYZER > $VCF
#echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> $VCF
#grep "^#CHROM" $VCF_OUTLYZER >> $VCF

cat $OUTPUT_TMP/*.vcf | awk -f $OUTLYZER_NORM $VCF_OUTLYZER > $VCF


echo "VCF="$VCF
tail $VCF


if ((1)); then
	source /tool/config/config.app
	HOWARD --input=$VCF --output=$VCF_ANN --annotation=hgvs,outcome,location --verbose --env=/tool/config/tools.app --annovar_folder=/STARK/tools/annovar/current/bin/ --threads=$CORE
fi

echo "VCF_ANN="$VCF_ANN
tail $VCF_ANN


#
