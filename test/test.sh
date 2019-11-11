#!/bin/bash

DATA_FOLDER=/STARK/data/

TEST_NUM=TEST_23


STARK_CMD="STARK"
#STARK_CMD="docker exec --rm STARK_DEV"


# OUTLYZER
$STARK_CMD --analysis_name=$TEST_NUM"_TEST_OUTLYZER_BED_GENES" --application=OUTLYZER+NO_ANNOTATION --reads=$DATA_FOLDER/SAMPLE/TEST/TEST.bam --design=$DATA_FOLDER/SAMPLE/TEST/TEST.bed --genes=$DATA_FOLDER/SAMPLE/TEST/TEST.genes --verbose


# MUTECT
$STARK_CMD --analysis_name=$TEST_NUM"_TEST_MUTECT_BED_GENES" --application=MUTECT+NO_ANNOTATION --reads=$DATA_FOLDER/SAMPLE/TEST/TEST.bam --design=$DATA_FOLDER/SAMPLE/TEST/TEST.bed --genes=$DATA_FOLDER/SAMPLE/TEST/TEST.genes --verbose


# TEST
$STARK_CMD --analysis_name=$TEST_NUM"_TEST_SOMATIC_BED_GENES" --application=SOMATIC+NO_ANNOTATION --reads=$DATA_FOLDER/SAMPLE/TEST/TEST.bam --design=$DATA_FOLDER/SAMPLE/TEST/TEST.bed --genes=$DATA_FOLDER/SAMPLE/TEST/TEST.genes --verbose

$STARK_CMD --analysis_name=$TEST_NUM"_TEST_simple" --reads=$DATA_FOLDER/SAMPLE/TEST/TEST.bam --verbose
$STARK_CMD --analysis_name=$TEST_NUM"_TEST_BED" --reads=$DATA_FOLDER/SAMPLE/TEST/TEST.bam --design=$DATA_FOLDER/SAMPLE/TEST/TEST.bed --verbose
$STARK_CMD --analysis_name=$TEST_NUM"_TEST_MANIFEST" --reads=$DATA_FOLDER/SAMPLE/TEST/TEST.bam --design=$DATA_FOLDER/SAMPLE/TEST/TEST.BRCA.manifest --verbose
$STARK_CMD --analysis_name=$TEST_NUM"_TEST_GENES" --reads=$DATA_FOLDER/SAMPLE/TEST/TEST.bam --genes=$DATA_FOLDER/SAMPLE/TEST/TEST.genes --verbose
$STARK_CMD --analysis_name=$TEST_NUM"_TEST_BED_GENES" --reads=$DATA_FOLDER/SAMPLE/TEST/TEST.bam --design=$DATA_FOLDER/SAMPLE/TEST/TEST.bed --genes=$DATA_FOLDER/SAMPLE/TEST/TEST.genes --verbose
$STARK_CMD --analysis_name=$TEST_NUM"_TEST_BED_GENES2" --reads=$DATA_FOLDER/SAMPLE/TEST/TEST.bam --design=$DATA_FOLDER/SAMPLE/TEST/TEST.bed --genes=$DATA_FOLDER/SAMPLE/TEST/TEST.genes+$DATA_FOLDER/SAMPLE/TEST/TEST.2.genes --verbose
$STARK_CMD --analysis_name=$TEST_NUM"_TEST_BED_GENES_SELECTED" --reads=$DATA_FOLDER/SAMPLE/TEST/TEST.bam --design=$DATA_FOLDER/SAMPLE/TEST/TEST.bed --genes=$DATA_FOLDER/SAMPLE/TEST/TEST.selected.genes --verbose
$STARK_CMD --analysis_name=$TEST_NUM"_TEST_BED_GENES3" --reads=$DATA_FOLDER/SAMPLE/TEST/TEST.bam --design=$DATA_FOLDER/SAMPLE/TEST/TEST.bed --genes=$DATA_FOLDER/SAMPLE/TEST/TEST.genes+$DATA_FOLDER/SAMPLE/TEST/TEST.2.genes+$DATA_FOLDER/SAMPLE/TEST/TEST.selected.genes --verbose
$STARK_CMD --analysis_name=$TEST_NUM"_TEST_BED_GENES3_TRANSCRIPTS" --reads=$DATA_FOLDER/SAMPLE/TEST/TEST.bam --design=$DATA_FOLDER/SAMPLE/TEST/TEST.bed --genes=$DATA_FOLDER/SAMPLE/TEST/TEST.genes+$DATA_FOLDER/SAMPLE/TEST/TEST.2.genes+$DATA_FOLDER/SAMPLE/TEST/TEST.selected.genes --transcripts=$DATA_FOLDER/SAMPLE/TEST/TEST.transcripts --verbose

$STARK_CMD --analysis_name=$TEST_NUM"_TEST_BED_REFGENES" --reads=$DATA_FOLDER/SAMPLE/TEST/TEST.bam --design=$DATA_FOLDER/SAMPLE/TEST/TEST.bed --genes=/STARK/databases/refGene/refGene.hg19.genes --verbose
$STARK_CMD --analysis_name=$TEST_NUM"_TEST_REFGENES" --reads=$DATA_FOLDER/SAMPLE/TEST/TEST.bam --design=$DATA_FOLDER/SAMPLE/TEST/TEST.bed --genes=/STARK/databases/refGene/refGene.hg19.genes --verbose

# HORIZON
$STARK_CMD --analysis_name=$TEST_NUM"_HORIZON_FQ_MANIFEST_NO_ANNOTATION" --reads=$DATA_FOLDER/SAMPLE/HORIZON/HORIZON.R1.fastq.gz --reads2=$DATA_FOLDER/SAMPLE/HORIZON/HORIZON.R2.fastq.gz --application=SOLIDTUMOR+NO_ANNOTATION --design=$DATA_FOLDER/SAMPLE/HORIZON/HORIZON.manifest --verbose
#$STARK_CMD --analysis_name=$TEST_NUM"_HORIZON_FQ_MANIFEST" --reads=$DATA_FOLDER/SAMPLE/HORIZON/HORIZON.R1.fastq.gz --reads2=$DATA_FOLDER/SAMPLE/HORIZON/HORIZON.R2.fastq.gz --application=SOLIDTUMOR --design=$DATA_FOLDER/SAMPLE/HORIZON/HORIZON.manifest --verbose

# HORIZON2
#$STARK_CMD --analysis_name=$TEST_NUM"_HORIZON2_simple" --reads=$DATA_FOLDER/SAMPLE/HORIZON2/HORIZON_R1.archive.cram --application=SOLIDTUMOR+NO_ANNOTATION --verbose
#$STARK_CMD --analysis_name=$TEST_NUM"_HORIZON2_MANIFEST" --reads=$DATA_FOLDER/SAMPLE/HORIZON2/HORIZON_R1.archive.cram --application=SOLIDTUMOR+NO_ANNOTATION --design=$DATA_FOLDER/SAMPLE/HORIZON2/HORIZON.manifest --verbose
$STARK_CMD --analysis_name=$TEST_NUM"_HORIZON2_MANIFEST_GENES" --reads=$DATA_FOLDER/SAMPLE/HORIZON2/HORIZON_R1.archive.cram --application=SOLIDTUMOR+NO_ANNOTATION --design=$DATA_FOLDER/SAMPLE/HORIZON2/HORIZON2.manifest --genes=$DATA_FOLDER/SAMPLE/HORIZON2/HORIZON.bed.genes --verbose

# RUN TEST
$STARK_CMD --run="RUN_TEST_TAG:"$TEST_NUM"_RUN_TEST_TAG" --verbose

# RUN TEST
#$STARK_CMD --run="RUN_TEST_TAG:"$TEST_NUM"_RUN_TEST_TAG" --sample_filter=P1439 --verbose

# RUN TINY
#$STARK_CMD --run="TINY:"$TEST_NUM"_TINY" --verbose

# SAMPLES
$STARK_CMD --analysis_name=$TEST_NUM"_P1439_FQ_MANIFEST" --reads=$DATA_FOLDER/SAMPLE/P1439/P1439.R1.fastq.gz --reads2=$DATA_FOLDER/SAMPLE/P1439/P1439.R2.fastq.gz --application=SOLIDTUMOR+NO_ANNOTATION --design=$DATA_FOLDER/SAMPLE/P1439/P1439.manifest --verbose


# RUNS SOLIDTUMOR SAMPLE P6153
$STARK_CMD --analysis_name=$TEST_NUM"_191001_M01656_0336_000000000-CLH59" --run="191001_M01656_0336_000000000-CLH59" --sample_filter=P6153 --verbose
$STARK_CMD --analysis_name=$TEST_NUM"_191001_M01656_0336_000000000-CLH59" --run="191001_M01656_0336_000000000-CLH59" --sample_filter=M8302,M8303,Z_THS0226_HORIZON --verbose

# M8302
# Z_THS0226_HORIZON

# RUNS HEMATOLOGY
#$STARK_CMD --analysis_name=$TEST_NUM"_190927_M01656_0335_000000000-CLGGY" --run="190927_M01656_0335_000000000-CLGGY" --verbose

# RUNS UMI
#$STARK_CMD --analysis_name=$TEST_NUM"_190830_M01656_0328_000000000-CLK2J" --run="190830_M01656_0328_000000000-CLK2J" --verbose






#
