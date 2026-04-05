# Application Configuration

This document describes the configuration variables for the STARK application, as defined in `config/default.app`.

Example of application configuration file:

```bash

#!/bin/bash
## STARK application EXOME

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP

# APPLICATION INFOS
#####################
APP_NAME="EXOME"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect germline mutations in exome sequencing data"
APP_GROUP="GENETIC"
APP_PROJECT="EXOME"

# ANALYSIS PARAMETERS
#######################

# PIPELINES
PIPELINES="bwamem.gatkHC_EXOME.howard bwamem.gatkUG_EXOME.howard"

# INTERVAL_PADDING / add some padding to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
INTERVAL_PADDING=100

# HOWARD ANNOTATION/PRIOTITIZATION/TRANSLATION CONFIGURATION

# ANNOTATION
# Default annotation with HOWARD for intermediate VCF (for each caller) used by default with annotation rule "howard"
HOWARD_ANNOTATION=""
# Default annotation with HOWARD for minimal VCF annotation (rule howard_minimal)
HOWARD_ANNOTATION_MINIMAL=""
# Default annotation with HOWARD for report
HOWARD_ANNOTATION_REPORT="core,frequency,score,annotation,prediction,snpeff_split"


```

_Application for Exome analysis: From default application, information of application are defind, as well as specific pipelines, interval padding and variant annotation_

## Application Infos

### APP_NAME

Default application name.

- Default: `"DEFAULT"`

### APP_RELEASE

Application release version.

- Default: `"1.2"`

### APP_DESCRIPTION

Description of the application.

- Default: `"Default application"`

### APP_GROUP

Application group.

- Default: `""`

### APP_PROJECT

Application project.

- Default: `""`

## Folders

### STARK_FOLDER_MAIN

Main STARK folder.

- Default: `"/STARK"`

### FOLDER_TOOLS

Folder for all tools needed for STARK, and more, including STARK.

- Default: `"$STARK_FOLDER_MAIN/tools"`

### FOLDER_DATABASES

Folder for all reference genomes and databases. The format for genomes is `$FOLDER_GENOMES/$ASSEMBLY/$ASSEMBLY.fa`. It also contains mandatory databases for variant calling and recalibration.

- Default: `"$STARK_FOLDER_MAIN/databases"`

### FOLDER_INPUT

Input folder.

- Default: `"$STARK_FOLDER_MAIN/input"`

### `FOLDER_RUN`

Illumina Sequencer repository Folder. Subfolder as runs.

- Default: `"$FOLDER_INPUT/runs"` (Commented out by default)

### `FOLDER_MANIFEST`

Illumina Manifests repository. Files to provide in the SampleSheet of each run.

- Default: `"$FOLDER_INPUT/manifests"` (Commented out by default)

### `FOLDER_PEDIGREE`

Pedigree repository.

- Default: `"$FOLDER_INPUT/pedigree"` (Commented out by default)

### FOLDER_OUTPUT

Folder where all results will be generated.

- Default: `"$STARK_FOLDER_MAIN/output"`

### `FOLDER_RESULTS`

RUN files such as BAM, VCF, metrics.

- Default: `"$FOLDER_OUTPUT/results"` (Commented out by default)

### `FOLDER_DEMULTIPLEXING`

Demultiplexing folder.

- Default: `"$FOLDER_OUTPUT/demulitplexing"` (Commented out by default)

### `FOLDER_LOG`

Log files folder.

- Default: `"$FOLDER_OUTPUT/log"` (Commented out by default)

### `FOLDER_TMP`

Temporary files folder.

- Default: `"$FOLDER_OUTPUT/tmp"` (Commented out by default)

### FOLDER_REPOSITORY

Repository folder for results data. Leave blank for no copy.

- Default: `"$FOLDER_OUTPUT/repository"`

### FOLDER_ARCHIVES

Archives folder for results data. Leave blank for no copy.

- Default: `"$FOLDER_OUTPUT/archives"`

### FOLDER_FAVORITES

Favorites folder for results data.

- Default: `""`

- **Configurations**:
  - To NOT use favorites folder: `FOLDER_FAVORITES=`
  - To copy favorites within repository folder: `FOLDER_FAVORITES=$FOLDER_REPOSITORY`
  - To copy favorites within default favorites folder: `FOLDER_FAVORITES=$FOLDER_OUTPUT/favorites`

## Parameters

### `RULES_APP`

Add specific rules to load. These files will be added to the list of rules files from APPS folder. For inheritance, use `RULES_APP="$RULES_APP MYGROUP/*.rules.mk"`.

- Default: `""`

- Example: `RULES_APP="MYGROUP/*.rules.mk" "$APP_FOLDER/*.rules.mk"`

### ASSEMBLY

Assembly to use (e.g., hg19, hg38). Automatically detected from manifest file if available.

- Default: `"hg19"`

### PIPELINES

Pipelines to use for the analysis. If variables `ALIGNERS`, `CALLERS`, and `ANNOTATORS` are defined, this variable will be automatically generated. This variable can be an additional pipeline to those defined by the combination of `ALIGNERS`, `CALLERS`, and `ANNOTATORS`. If no pipeline is finally defined, the default pipeline will be applied.

- Default: `"bwamem.gatkHC.howard"`

- **Format**: `"ALIGNER1.CALLER1.ANNOTATOR ALIGNER1.CALLER2.ANNOTATOR1 ALIGNER2.CALLER1.ANNOTATOR1"`

### `ALIGNERS`

Aligners to use for the analysis.

- Default: `""` (Commented out by default)

- **Example of available aligners**: `bwamem`, `bwasw`, `bwaaln`

- Example: `ALIGNERS="bwamem"`

### `CALLERS`

Callers to use for the analysis.

- Default: `""` (Commented out by default)

- **Example of available callers**: `gatkHC`, `gatkUG`, `VarScan`, `samtools`

- Example: `CALLERS="gatkHC"`

### `ANNOTATORS`

Annotators to use for the analysis.

- Default: `""` (Commented out by default)

- **Example of available annotators**: `howard`, `snpeff`

- Example: `ANNOTATORS="howard"`

### BLANK

Blank samples used to reject for CNV analysis.

- Default: `"BlcADN,blanc,BlcPCR,blcPCR,T_NTC,Z_NTC"`

### BARCODE_MISMATCHES

Number of mismatches allowed for demultiplexing.

- Default: `1`

### BAM_METRICS

Performs BAM metrics (1/TRUE/YES/Y or 0/FALSE/NO/N). Time and space consuming. Switch off for exome/genome for better performances.

- Default: `0`

### BAM GENE COVERAGE METRICS (default 1/TRUE/YES/Y)

Performs BAM GENE COVERAGE METRICS (1/TRUE/YES/Y or 0/FALSE/NO/N) using RNASeQC. Time and space consuming.

- Default: `1`

### METRICS_MINIMUM_MAPPING_QUALITY

Minimum mapping quality to consider in the BAM metrics.

- Default: `10`

### METRICS_MINIMUM_BASE_QUALITY

Minimum base quality to consider in the BAM metrics.

- Default: `10`

### CLIP_OVERLAPPING_READS

From PICARD: For paired reads, soft clip the 3' end of each read if necessary so that it does not extend past the 5' end of its mate.

- Default: `1`

### METRICS_FLAGS

Flagged reads in the metrics BAM (mpileup format).

- Default: `"UNMAP,SECONDARY,QCFAIL,DUP"`

### `SAMTOOLS_METRICS_FLAG_PARAM`

Flagged reads in the metrics BAM (samtools format). Generated from mpileup format if empty.

- Default: `""` (Commented out by default)

- Example: `SAMTOOLS_METRICS_FLAG_PARAM=" -F 0x4 -F 0x100 -F 0x200 -F 0x400"`

### INTERVAL_PADDING

Add some “padding” to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp).

- Default: `0`

### COVERAGE_CRITERIA

For gene coverage metrics, the criteria to calculate the percent of bases over X coverage (e.g., 30 for 30X).

- Default: `"1,5,10,20,30,50,100,200,300"`

- Example: `"1,30"`

### SEQUENCING_DEPTH

Sequencing depth threshold for gene coverage metrics.

- Default: `"1"`

### SEQUENCING_COVERAGE_THRESHOLD

Sequencing coverage threshold for gene coverage metrics.

- Default: `"1"`

### MINIMUM_DEPTH

Fail DP threshold for gene coverage metrics.

- Default: `"30"`

### EXPECTED_DEPTH

Warn DP threshold for gene coverage metrics.

- Default: `"100"`

### DEPTH_COVERAGE_THRESHOLD

Threshold percentage of bases over the DP threshold for gene coverage metrics.

- Default: `"0.95"`

### NB_BASES_AROUND

For gene coverage metrics, the number of bases to look around the exons from the given bed file.

- Default: `0`

### GENESCOVERAGE_PRECISION

Genes Coverage calculation precision.

- Default: `2`

### BEDFILE_GENES

For gene coverage metrics, the bed file containing the 5'UTR, 3'UTR and genomic coding coordinates.

- Default: `""`

### VARANK_ANALYSIS

Performs VARANK ANALYSIS with Alamut (1 for true or 0 for false).

- Default: `0`

### `VARANK_FOLDER`

If `VARANK_ANALYSIS` is enabled, this defines the folder to store the results.

- Default: `"$FOLDER_RESULTS/VARANK"` (Commented out by default)

### BAM_CHECK_STEPS

Check BAM for each manipulation step (clipping, realignment...). Time consuming, but will stop the analysis in case of missing reads. A Metrics on each BAM check the BAM discrepancy in any case.

- Default: `0`

### METRICS_SNPEFF

Generate snpEff variant metrics from VCF. Only for Report final VCF.

- Default: `0`

### PRIORITIZE_PIPELINES_LIST

List of pipelines to prioritize for the report (final.vcf).

- Default: `""`

## FASTQ Processing

### STARK_DEMULTIPLEXING_BASES_MASK

Set mask for demultiplexing.

- Default: `""` (auto from SampleSheet)

- **Examples**: `"Y150,I10,Y10,Y150"`, `"Y150,I8,Y10,Y150"` (UMI index2)

### STARK_DEMULTIPLEXING_MASK_SHORT_ADAPTATER_READ

Set short read size for demultiplexing. If demultiplexing UMI within a read must be set to 0 for Agilent XTHS kits.

- Default: `""` (auto)

### STARK_DEMULTIPLEXING_READS_MAPPING

Redefine FASTQ files in order to identify R1, R2, I1 and I2. Order: R1 I1 I2 R2.

- Default: `"R1 I1 I2 R2"`

- Example: `"R1 I1 R2 R3"` (UMI index2)

### ADAPTER_STRINGENCY

Demultiplexing adapter stringency for BCL2FASTQ.

- Default: `0.9`

### STARK_DEMULTIPLEXING_BCL2FASTQ_OPTIONS

Demultiplexing options for BCL2FASTQ.

- Default: `"--no-lane-splitting --create-fastq-for-index-reads"`

### FASTQ_DEMULTIPLEXING_COMPRESSION_LEVEL

Zlib compression level (1-9) used for FASTQ files during demultiplexing (used by BCL2FASTQ). If `FASTQ_DEMULTIPLEXING_KEEP=1`, a high level of compression (at least 5) is suggested.

- Default: `1`

### FASTQ_COMPRESSION_LEVEL

Zlib compression level (1-9) used for main FASTQ files (used by FASTP).

- Default: `1`

### ENABLE_ADAPTER_TRIMMING

Trim adapter and autodetect adapter for paired end (0 or 1).

- Default: `0` (adapter trimming is disabled)

### FASTQ_QUALITY_FILTERING

Read quality threshold. Reads with quality below this will be removed.

- Default: `0` (disabled)

### POLY_G_MIN_LEN

Force polyG tail trimming. By default, trimming is automatically enabled for Illumina NextSeq/NovaSeq data. This is the minimum length to detect polyG in the read tail.

- Default: `0` (disabled)

### READ_LENGTH_REQUIRED

Reads shorter than this length will be discarded.

- Default: `0` (disabled)

### UMI_LOC

Set the UMI location. If not null, UMI extraction and analysis will be performed. See FASTP/UMI TOOLS documentation for more information.

- Default: `""`

- **Available locations**:
  - `index1`: the first index is used as UMI. If the data is PE, this UMI will be used for both read1/read2.
  - `index2`: the second index is used as UMI. PE data only, this UMI will be used for both read1/read2.
  - `read1`: the head of read1 is used as UMI. If the data is PE, this UMI will be used for both read1/read2.
  - `read2`: the head of read2 is used as UMI. PE data only, this UMI will be used for both read1/read2.
  - `per_index`: read1 will use UMI extracted from index1, read2 will use UMI extracted from index2.
  - `per_read`: read1 will use UMI extracted from the head of read1, read2 will use UMI extracted from the head of read2.

- Example: `UMI_LOC="index2"`

### UMI_BARCODE_PATTERN

Set the UMI Barcode pattern. If not null, STARK will prepare fastq containing UMIs +/- cell barcodes for alignment. If `UMI_LOC` is `per_index` or `per_read`, and `UMI_BAR_CODE_PATTERN` is defined as "NNNNN", it will be redefined as "NNNNN-NNNNN". Only the length of the first part of the duplex barcode will be used with FASTP. See FASTP/UMI TOOLS documentation for more information.

- Default: `""`

- **Examples**: `UMI_BAR_CODE_PATTERN="NNNNNNNNNN"` (simplex), `UMI_BAR_CODE_PATTERN="NNNNN-NNNNN"` (duplex)

### BARCODE_TAG

Barcode to use for Mark Duplicates. If not null, Mark Duplicates will consider this tag. See PICARD documentation for more information.

- Default: `""`

- **Examples**: `BARCODE_TAG="BC"` (10X Genomics), `BARCODE_TAG="BX"` (UMI)

### PICARD_MARKDUP_OPTICAL_DEDUP

Set to `"READ_NAME_REGEX=null"` to disable optical deduplication in PICARD MarkDuplicates.

- Default: `""`

### FASTQ_DEMULTIPLEXING_KEEP

Keep demultiplexed FASTQ files or files from input reads/reads2.

- Default: `0`

### SEQUENCING_DEMULTIPLEXING_FOLDER

Folder for demultiplexing FASTQ within `$SAMPLE.sequencing` folder.

- Default: `"demultiplexing"`

### FASTP_ADDITIONAL_OPTIONS

Additional options for FASTP. See FASTP documentation.

- Default: `""`

### FASTQ_PROCESSING_STEPS

All steps to process input FASTQ files, after sequencing and demultiplexing (if any).

- Default: `"fastq_reheader sort fastp fastq_clean_header compress"`

- **Available steps**:
  - `fastq_reheader`: FASTQ reheader to integrate index within FASTQ comment Illumina tag (e.g. 1:N:0:xxx). Nothing done if already integrated.
  - `fastq_clean_header`: FASTQ read head formatting, especially SAMTOOLS tags. Nothing done if not needed.
  - `compress`: FASTQ files compression (see `FASTQ_COMPRESSION_LEVEL`).
  - `sort`: Sort FASTQ using read name.
  - `fastp`: Process FASTP algorithm and report, UMI extraction, quality filtration...
  - `umi_tools`: Process UMITools algorithm for UMI extraction.

- **Examples**:
  - `"fastq_reheader sort fastp fastq_clean_header compress"` (for UMI technology)
  - `"sort compress"` (for sorting and compression only)

## Pipeline Steps

### POST_SEQUENCING_STEPS

All steps and before alignment. This sequence corresponds to the FASTQ file processing before the alignment (trimming...). The steps are defined as makefiles rules. Check available steps by using the command: `STARK --pipelines_infos`.

- Default: `""` (nothing to do)

### POST_ALIGNMENT_STEPS

All steps after alignment and before calling. This sequence corresponds to the BAM file generated just after the alignment. The steps are defined as makefiles rules. Check available steps by using the command: `STARK --pipelines_infos`.

- Default: `"sorting markduplicates realignment recalibration compress"`

- **Available steps (not up-to-date)**:
  - `sorting`: BAM sorting
  - `compress`: BAM compression (see `$BAM_COMPRESSION` variable)
  - `realignment`: local realignment
  - `recalibration`: reads recalibration
  - `gencore`: gencore is a tool for fast and powerful deduplication for paired-end next-generation sequencing
  - `markduplicates`: Mark duplicated reads in BAM with PICARD MarkDuplicates. Use `BARCODE_TAG` to specify tag.
  - `clipping`: BAM Clipping according to primer definition in manifest file, if any.

- **Examples**:
  - `"sorting realignment clipping compress"` (for Amplicon technology)
  - `"sorting markduplicates realignment compress"` (for Capture technology)
  - `"sorting gencore markduplicates realignment compress"` (for UMI technology)

### POST_CALLING_STEPS

All steps after calling. This sequence corresponds to the VCF file generated just after the calling. The steps are defined as makefiles rules. Check available steps by using the command: `STARK --pipelines_infos`.

- Default: `" "`

- **Available steps (not up-to-date)**:
  - `sorting`: VCF sort
  - `normalization`: VCF normalization
  - `variantrecalibration`: VCF recalibration (using GATK4). Includes variantfiltration if no recalibration possible.
  - `variantfiltration`: VCF filtration (using GATK4).

- **Examples**:
  - `" "` (to avoid at this pipeline step, see `POST_CALLING_MERGING_STEPS`)
  - `"normalization variantfiltration"` (for gene panel)
  - `"normalization variantrecalibration"` (for exome or genome)

### POST_CALLING_MERGING_STEPS

All steps after merging calling VCFs. This sequence corresponds to the VCF file generated after the merge of VCF calling. The steps are defined as makefiles rules. Check available steps by using the command: `STARK --pipelines_infos`.

- Default: `"sorting normalization variantrecalibration variantfiltration"`

- **Available steps (not up-to-date)**:
  - `sorting`: VCF sort
  - `normalization`: VCF normalization
  - `variantrecalibration`: VCF recalibration (using GATK4). Includes variantfiltration if no recalibration possible.
  - `variantfiltration`: VCF filtration (using GATK4).

- **Examples**:
  - `"sorting normalization variantrecalibration"` (for exome or genome)
  - `"sorting normalization variantfiltration"` (for gene panel)

### POST_ANNOTATION_STEPS

All steps after annotation. This sequence corresponds to the VCF file generated just after the annotation. The steps are defined as makefiles rules. Check available steps by using the command: `STARK --pipelines_infos`.

- Default: `" "`

- **Available steps (not up-to-date)**:
  - `sorting`: VCF sorting

- **Examples**:
  - `" "` (to avoid at this pipeline step)
  - `"sorting"`

## Compression

### BAM_COMPRESSION

Final BAM compression level (`ALIGNER.bam`).

- Default: `9`

### BAM_VALIDATION_COMPRESSION

Validation BAM compression level (`ALIGNER.validation.bam`).

- Default: `5`

## Gencore

### GENCORE_SUP_READS

Number of supporting reads to keep clusters. Set to 2 for ultrasensitive filter, 1 to replace Picard Markduplicates.

- Default: `""` (corresponds to `--supporting_reads 1`)

### GENCORE_SCORE_THREESHOLD

Score threshold for deduplication. Set to 8 recommended for dup-rate < 50% if you want to keep all the DNA fragments.

- Default: `""` (corresponds to `--score_threshold 6`)

### GENCORE_RATIO_THREESHOLD

If the ratio of the major base in a cluster is less than this threshold, it will be further compared to the reference. Value should be 0.5~1.0.

- Default: `""` (corresponds to `--ratio_threshold 0.8`)

### GENCORE_DIFF_THREESHOLD

If two reads with identical mapping position have UMI difference <= this threshold, they will be merged to generate a consensus read.

- Default: `""` (corresponds to `--umi_diff_threshold 2`)

### GENCORE_QUAL_THREESHOLD

Quality thresholds for base quality (`--high_qual`, `--moderate_qual`, `--low_qual`).

- Default: `""` (corresponds to default qualities: Q30, Q20, Q15)

- Example: `"--moderate_qual 20"`, `"--high_qual 20 --moderate_qual 15 --low_qual 10"`

### GENCORE_COVERAGE_SAMPLING

The sampling rate for genome scale coverage statistics. For statistics purpose, can be reduced to gain performance.

- Default: `""` (corresponds to `--coverage_sampling=10000`)

## GATK Calling

### USE_VCFDBSNP_WITH_GATK

Use VCF DBSNP for GATK calling and annotation. This option allows GATK to use the VCF DBSNP if the file is available and not empty.

- Default: `0`

### VCFDBSNP

Path to the dbSNP database as VCF used for GATK calling and annotation. This variable is automatically defined in `databases.app`.

- Default: `""`

- Example: `VCFDBSNP=$FOLDER_DATABASES/dbsnp_138.hg19.vcf.gz`

## CRAM

### CRAM_OPTIONS

Final CRAM options for compression (`archive.cram`).

- Default: `"version=3.0,level=9,no_ref"`

- Example: `"version=3.0,level=9,no_ref,use_lzma,use_bzip2,use_fqz,seqs_per_slice=100000"`

### CRAM_REMOVE_TAGS

Final CRAM options for tags to remove (`archive.cram`).

- Default: `"BD,BI"`

- Example: `"BD,BI,OQ"`

## VCF Report

### VCF_MISSING_GENOTYPE

Replace `0/0` or `0|0` genotypes by `./.` and `.|.`.

- Default: `"missing_clean"`

- **Options**:
  - `"missing"`: Replace genotypes.
  - `"missing_clean"`: Replace genotypes and clean by removing other FORMAT fields.

## Resource Management

### THREADS

Number of threads to use for the analysis. `AUTO` will consider CORE-1 threads. The number is auto-adjusted if the value is incorrect.

- Default: `"AUTO"`

### `THREADS_LOADING`

Number of threads used for loading demultiplexed data.

- Default: `THREADS` (Commented out by default)

### `THREADS_WRITING`

Number of threads used for writing demultiplexed data.

- Default: `THREADS` (Commented out by default)

### `THREADS_COPY`

Number of threads used for copying files in repositories.

- Default: `1` (Commented out by default)

### `MEMORY`

Total memory to use by thread. Default is `MemTotal` from `/proc/meminfo` divided by number of threads.

- Default: `AUTO` (Commented out by default)

### `JAVA_MEMORY`

Maximum memory to use for Java.

- Default: `MEMORY` (Commented out by default)

### MAX_VALIDATION_BAM_SIZE

If validation bam is smaller than this size (in Kb), launch `CollectHsMetrics` the classic way. Otherwise, increase RAM and limit concurrent launches.

- Default: `1000000000`

- Example: `MAX_VALIDATION_BAM_SIZE=2097152` (2Go)

### MAX_CONCURRENT_HSMETRICS

Limit concurrent launches of `CollectHsMetrics`.

- Default: `1`

- Example: `MAX_CONCURRENT_HSMETRICS=1` (No concurrence)

### MAX_CONCURRENT_HSMETRICS_RAM

RAM to allocate for `CollectHsMetrics` when BAM size is large.

- Default: `"16g"`

## HOWARD Configuration

### `HOWARD_CONFIG_ANNOTATION`

HOWARD Configuration file for Annotation.

- Default: `"$HOWARD_FOLDER_CONFIG/config.annotation.ini"` (Commented out by default)

### `HOWARD_CONFIG_PRIORITIZATION`

HOWARD Configuration file for Prioritization.

- Default: `"$HOWARD_FOLDER_CONFIG/config.prioritization.ini"` (Commented out by default)

### `HOWARD_CONFIG_DEJAVU_ANNOTATION`

HOWARD DEJAVU Configuration file for Annotation.

- Default: `"$HOWARD_FOLDER_CONFIG/config.annotation.ini"` (Commented out by default)

### HOWARD_CONFIG

HOWARD configuration file.

- Default: `"$HOWARD_FOLDER_CONFIG/config.json"`

### `HOWARD_PARAM`

Default HOWARD parameters.

- Default: `"$HOWARD_FOLDER_CONFIG/param.json"` (Commented out by default)

### `HOWARD_PARAM_MINIMAL`

Default HOWARD parameters for minimal VCF annotation.

- Default: `"$HOWARD_FOLDER_CONFIG/param.json"` (Commented out by default)

### HOWARD_PARAM_REPORT

Default HOWARD parameters for report.

- Default: `"$HOWARD_FOLDER_CONFIG/param.json"`

### `HOWARD_PARAM_ANALYSIS`

Default HOWARD parameters for whole analysis.

- Default: `"$HOWARD_FOLDER_CONFIG/param.json"` (Commented out by default)

### HOWARD_PRIORITIZATION_CONFIG

HOWARD prioritization parameters file.

- Default: `"$HOWARD_FOLDER_CONFIG/prioritization_profiles.json"`

### INFO_TO_FORMAT_ANNOTATIONS

Transfers INFO annotation to FORMAT annotation. Useful for annotations on full VCF to final VCF on each sample.

- Default: `""`

## GATK4 Recalibrator and Filtration

### `VARIANTFILTRATION_OPTIONS`

Variant Filtration main option. See documentation guide for more info.

- Default: `""` (Commented out by default)

### VARIANTFILTRATION_SNP_FILTER_OPTION

One or more expression used with INFO fields to filter SNP. See documentation guide for more info.

- Default: `'--filter-name "SNP_filter_QD" --filter-expression "QD < 2.0" ...'`

### VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION

One or more expressions used with FORMAT (sample/genotype-level) fields to filter SNP. See documentation guide for more info.

- Default: `'--genotype-filter-expression "GQ == 0" --genotype-filter-name "genotype_GQ_filter_VeryVeryLow" ...'`

### VARIANTFILTRATION_INDEL_FILTER_OPTION

One or more expressions used with INFO fields to filter INDEL. See documentation guide for more info.

- Default: `'--filter-name "INDEL_filter_QD" --filter-expression "QD < 2.0" ...'`

### VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION

One or more expressions used with FORMAT (sample/genotype-level) fields to filter INDEL. See documentation guide for more info.

- Default: `'--genotype-filter-expression "GQ == 0" --genotype-filter-name "genotype_GQ_filter_VeryVeryLow" ...'`

### VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS

Remove previous filters applied to the VCF.

- Default: `1`

### `VARIANTRECALIBRATOR_OPTIONS`

Variant Recalibrator main option. See documentation guide for more info.

- Default: `""` (Commented out by default)

### VARIANTRECALIBRATION_SNP_RESOURCES

Variant Recalibrator SNP resources option. These resources need to be available on STARK Databases folder for GATK.

- Default: Depends on `ASSEMBLY`. For hg19: `"-resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg19.sites.vcf.gz ..."`

### VARIANTRECALIBRATION_INDEL_RESOURCES

Variant Recalibrator INDEL resources option. These resources need to be available on STARK Databases folder for GATK.

- Default: Depends on `ASSEMBLY`. For hg19: `"-resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz ..."`

### VARIANTRECALIBRATION_SNP_ANNOTATIONS

Variant Recalibrator SNP annotations option.

- Default: `"-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP"`

### VARIANTRECALIBRATION_INDEL_ANNOTATIONS

Variant Recalibrator INDEL annotations option.

- Default: `"-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum"`

### VARIANTRECALIBRATION_SNP_TRANCHES

Variant Recalibrator SNP tranches option.

- Default: `"-tranche 100.0 -tranche 99.95 -tranche 99.9 ... -tranche 90.0"`

### VARIANTRECALIBRATION_INDEL_TRANCHES

Variant Recalibrator INDEL tranches option.

- Default: `"-tranche 100.0 -tranche 99.95 -tranche 99.9 ... -tranche 90.0"`

### VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_OPTION

Optional Variant Filtration for SNP if Variant Recalibrator failed. Use empty value to switch off.

- Default: `$VARIANTFILTRATION_SNP_FILTER_OPTION`

### VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION

Optional Variant Filtration for SNP (genotype-level) if Variant Recalibrator failed.

- Default: `$VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION`

### VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_OPTION

Optional Variant Filtration for INDEL if Variant Recalibrator failed.

- Default: `$VARIANTFILTRATION_INDEL_FILTER_OPTION`

### VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION

Optional Variant Filtration for INDEL (genotype-level) if Variant Recalibrator failed.

- Default: `$VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION`

## GATK BAM Realignment and Recalibration

### GATK_REALIGNMENT_KNOWN_OPTIONS

Known sites for GATK3 BAM realignment. Realignment is a process of correcting misalignments around indels. It is not performed with GATK4 Haplotype Caller and MuTect2 as they perform local realignment during calling.

- Default: Dynamically generated from `VARIANTRECALIBRATION_INDEL_RESOURCES`.

### GATK_RECALIBRATION_KNOWN_OPTIONS

Known sites for GATK4 BAM recalibration. Recalibration is a process of correcting base quality scores.

- Default: Dynamically generated from `VARIANTRECALIBRATION_SNP_RESOURCES`.

## Report

### REPORT_SECTIONS

List of sections to show in the report.

- Default: `"ALL"`

- **Sections**: `results_summary`, `sequencing_mapping`, `depth`, `coverage`, `variant_calling`, `variant_stats`

- **Annex Sections**: `annex_coverage`, `annex_depth`, `annex_genes_coverage`, `annex_variants`, `annex_annotations`

### REPORT_VARIANTS_FULL

Generate variants files from run with full VCF (include all calling information).

- Default: `0`

## Repository and Archives

### REPOSITORY_FILE_PATTERNS

Repository files patterns to add.

- Default: `' $SAMPLE.*.validation.bam $SAMPLE.*.validation.bam.bai $SAMPLE.*.bam.metrics/$SAMPLE.*.validation.flags.Design.bed $SAMPLE.reports/$SAMPLE.full.Design.vcf.gz $SAMPLE.reports/$SAMPLE.full.Design.vcf.gz.tbi $SAMPLE.reports/$SAMPLE.full.Design.tsv '`

### REPOSITORY_FILE_SUBFOLDER_PATTERNS

Repository files patterns to exclude on results subfolder (usually 'STARK'). Useful to reduce storage and exclude repeated files.

- Default: `' $SAMPLE.*fastq.gz $SAMPLE.*.validation.bam $SAMPLE.*.validation.bam.bai '`

### ARCHIVES_FILE_PATTERNS

Archives files patterns to add.

- Default: `' $SAMPLE.reports/$SAMPLE.full.vcf.gz $SAMPLE.reports/$SAMPLE.full.vcf.gz.tbi $SAMPLE.reports/$SAMPLE.final.tsv $SAMPLE.*.bam.metrics/$SAMPLE.*.validation.flags.*.bed '`

### FAVORITES_FILE_PATTERNS

Favorites files patterns to add.

- Default: `''`

## IGV Session

The IGV session XML file (`*igv_session.xml`) is generated with Samples and Runs files (BAM, VCF, BED). Parameters select files through patterns and folder depth.

### IGV Session Patterns (Sample)

- **`REPOSITORY_SAMPLE_IGV_SESSION_PATTERNS`**: `'$IGV_SESSION_PATTERNS_SAMPLE_DEFAULT_APP'`

- **`REPOSITORY_SAMPLE_IGV_SESSION_MINDEPTH`**: `$IGV_SESSION_MINDEPTH_SAMPLE_DEFAULT_APP`

- **`REPOSITORY_SAMPLE_IGV_SESSION_MAXDEPTH`**: `$IGV_SESSION_MAXDEPTH_SAMPLE_DEFAULT_APP`

- ... (similar variables for `ARCHIVES` and `FAVORITES`)

### IGV Session Patterns (Run)

- **`REPOSITORY_RUN_IGV_SESSION_PATTERNS`**: `'$IGV_SESSION_PATTERNS_RUN_DEFAULT_APP'`

- **`REPOSITORY_RUN_IGV_SESSION_MINDEPTH`**: `$IGV_SESSION_MINDEPTH_RUN_DEFAULT_APP`

- **`REPOSITORY_RUN_IGV_SESSION_MAXDEPTH`**: `$IGV_SESSION_MAXDEPTH_RUN_DEFAULT_APP`

- ... (similar variables for `ARCHIVES` and `FAVORITES`)

### DISPLAYMODE_BAM

IGV Display mode for BAM files.

- Default: `"SQUISHED"`

- **Options**: `SQUISHED`, `COLLAPSED`, `EXPANDED`

### DISPLAYMODE_VCF

IGV display mode for VCF files.

- Default: `"COLLAPSED"`

### DISPLAYMODE_BED

IGV display mode for BED files.

- Default: `"COLLAPSED"`

### IGV_SESSION_RESSOURCES

Additional databases for IGV session (see IGV doc).

- Default: `'<Resource index="https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz.tbi" name="Refseq Genes" path="https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz" type="refgene"/>'`

### IGV_SESSION_DATAPANEL

Additional data panel for IGV session.

- Default: `''`

### IGV_SESSION_FEATUREPANEL

Additional feature panel for IGV session.

- Default: `'<Track attributeKey="Refseq Genes" .../>'`

### IGV_SESSION_DAS

Data as a service URL for IGV session.

- Default: `"http://localhost:4201/static/data/public/repositories"`

### IGV_SESSION_DAS_REPOSITORIES

DAS repositories for IGV session.

- Default: `"$IGV_SESSION_DAS/Repository $IGV_SESSION_DAS/Archives $IGV_SESSION_DAS/Favorites"`

### `IGV_SESSION_JSON_TRACKS_ADDITIONAL`

Additional tracks for IGV session in JSON format.

- Default: `""`

- Example: `'{ "type": "bed", "url": "...", "name": "Gencode V18" }'`
