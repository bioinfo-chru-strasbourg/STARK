#!/bin/bash
#################################
## STARK environment
#################################
# DATABASES
#############
DATABASES_LIST=""
DATABASES_CONFIG_LIST=""

#####################
# ASSEMBLY & GENOME #
#####################

#GENOME_REGEX="--download-genomes-contig-regex='chr[0-9XYM]+$'"
GENOME_REGEX=""
export GENOME_REGEX

if [ -z $ASSEMBLY ] || [ "$ASSEMBLY" == "" ]; then
	ASSEMBLY=hg19
fi;
export ASSEMBLY

if [ -z $GENOME ] || [ "$GENOME" == "" ]; then
	GENOME=$DBFOLDER/genomes/current/$ASSEMBLY/$ASSEMBLY.fa
fi;
export GENOME

if [ -z $DICT ] || [ "$DICT" == "" ]; then
	DICT=$DBFOLDER/genomes/current/$ASSEMBLY/$ASSEMBLY.dict
fi;
export DICT

# REF_CACHE_FOLDER and REF_CACHE
#export REF_CACHE_FOLDER=$GENOME.hts-ref;
#export REF_CACHE="$REF_CACHE_FOLDER/%2s/%2s/%s";

# Indexing genome with BWA2
BWA2_INDEX=0 # 1 = true
export BWA2_INDEX

###########################
# Databases Configuration #
###########################

# Arriba
ARRIBA_CURRENT="https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz";
export ARRIBA_CURRENT

if [ ! -z $FOLDER_DATABASES_ARRIBA ] && [ "$FOLDER_DATABASES_ARRIBA" != "" ]; then
	ARRIBA_DATABASES=$FOLDER_DATABASES_ARRIBA
else
	ARRIBA_DATABASES=$DBFOLDER/arriba/current
fi;
export ARRIBA_DATABASES
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" ARRIBA_DATABASES"

# STAR Fusion (CTAT Lib)
if [ $ASSEMBLY == "hg19" ] ; then CTAT_CURRENT="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz"; fi;
if [ $ASSEMBLY == "hg38" ] ; then CTAT_CURRENT="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz"; fi;
CTAT_PM="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/AnnotFilterRule.pm"
export CTAT_CURRENT
export CTAT_PM

if [ ! -z $FOLDER_DATABASES_CTAT ] && [ "$FOLDER_DATABASES_CTAT" != "" ]; then
	CTAT_DATABASES=$FOLDER_DATABASES_CTAT
else
	CTAT_DATABASES=$DBFOLDER/CTAT_LIB/current
fi;
export CTAT_DATABASES
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" CTAT_DATABASES"

# Gencode
# for hg19 the last gencode version is v19
if [ $ASSEMBLY == "hg19" ] ; then
	GENCODE_VERSION="19"
	GENCODE_CURRENT="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$GENCODE_VERSION/gencode.v$GENCODE_VERSION.annotation.gtf.gz";
fi;
# for hg38 the first gencode version is v20 ; current version (10/2023) is v44
if [ $ASSEMBLY == "hg38" ] ; then 
	GENCODE_VERSION="44"
	GENCODE_CURRENT="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$GENCODE_VERSION/gencode.v$GENCODE_VERSION.primary_assembly.annotation.gtf.gz";
fi;
export GENCODE_VERSION
export GENCODE_CURRENT

if [ ! -z $FOLDER_DATABASES_GENCODE ] && [ "$FOLDER_DATABASES_GENCODE" != "" ]; then
	GENCODE_DATABASES=$FOLDER_DATABASES_GENCODE
else
	GENCODE_DATABASES=$DBFOLDER/gencode/current
fi;
export GENCODE_DATABASES
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" GENCODE_DATABASES"

# Needed for genome indexing with STAR
DBFOLDER_GENCODE=$(dirname $GENCODE_DATABASES)
export DBFOLDER_GENCODE

# ANNOVAR
if [ ! -z $FOLDER_DATABASES_ANNOVAR ] && [ "$FOLDER_DATABASES_ANNOVAR" != "" ]; then
	ANNOVAR_DATABASES=$FOLDER_DATABASES_ANNOVAR
else
	ANNOVAR_DATABASES=$DBFOLDER/annovar/current
fi;
export ANNOVAR_DATABASES
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" ANNOVAR_DATABASES"
# ANNOVAR Files to download with HOWARD
ANNOVAR_FILES="refGene,gnomad_exome,cosmic70,dbnsfp42a,clinvar_202*,nci60"
export ANNOVAR_FILES

# DEJAVU ANNOVAR
if [ ! -z $FOLDER_DATABASES_DEJAVU_ANNOVAR ] && [ "$FOLDER_DATABASES_DEJAVU_ANNOVAR" != "" ]; then
	DEJAVU_ANNOVAR_DATABASES=$FOLDER_DATABASES_DEJAVU_ANNOVAR
else
	DEJAVU_ANNOVAR_DATABASES=$DBFOLDER/annovar/current
fi;
export DEJAVU_ANNOVAR_DATABASES
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" DEJAVU_ANNOVAR_DATABASES"

# SNPEFF
if [ ! -z $FOLDER_DATABASES_SNPEFF ] && [ "$FOLDER_DATABASES_SNPEFF" != "" ] && [ -d "$FOLDER_DATABASES_SNPEFF" ]; then
	SNPEFF_DATABASES=$FOLDER_DATABASES_SNPEFF
elif  [ ! -z $DBFOLDER/snpeff/$SNPEFF_VERSION ] && [ "$DBFOLDER/snpeff/$SNPEFF_VERSION" != "" ] && [ -d "$DBFOLDER/snpeff/$SNPEFF_VERSION" ]; then
	SNPEFF_DATABASES=$DBFOLDER/snpeff/$SNPEFF_VERSION
else
	SNPEFF_DATABASES=$DBFOLDER/snpeff/current
fi;
export SNPEFF_DATABASES
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" SNPEFF_DATABASES"

# dbSNP ; mandatory DB needed for GATK variant recalibration
if [ ! -z $FOLDER_DATABASES_DBSNP ] && [ "$FOLDER_DATABASES_DBSNP" != "" ]; then
	DBSNP_DATABASES=$FOLDER_DATABASES_DBSNP
else
	DBSNP_DATABASES=$DBFOLDER/dbsnp/current
fi;
export DBSNP_DATABASES
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" DBSNP_DATABASES"
DBFOLDER_DBSNP=$DBFOLDER/dbsnp
export VCFDBSNP=$DBFOLDER_DBSNP/current/$ASSEMBLY/dbsnp.vcf.gz

if [ ! -e $VCFDBSNP ]; then
	echo "#[WARNING] No VCFDBSNP '$VCFDBSNP' in the database. Calling step impossible. Please check '$DBFOLDER' folder or configuration file" >>/dev/stderr
fi;
DATABASES_LIST=$DATABASES_LIST" VCFDBSNP"

# dbNSFP
if [ ! -z $FOLDER_DATABASES_DBNSFP ] && [ "$FOLDER_DATABASES_DBNSFP" != "" ]; then
	DBNSFP_DATABASES=$FOLDER_DATABASES_DBNSFP
else
	DBNSFP_DATABASES=$DBFOLDER/dbnsfp/current
fi;
export DBNSFP_DATABASES
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" DBNSNFP_DATABASES"

# refGene/refSeq
DBFOLDER_REFGENE=$DBFOLDER/refGene
export DBFOLDER_REFGENE
export REFSEQ_GENES=$DBFOLDER_REFGENE/current/$ASSEMBLY/refGene.$ASSEMBLY.bed

##################
# GATK DATABASES #
##################

# GATK VARIANT RECALIBRATION URLs
if [ $ASSEMBLY == "hg38" ]; then
	DBFOLDER_GATK_URL_DEFAULT="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0"
	DBFOLDER_GATK_URL_DBSNP="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"
elif [ $ASSEMBLY == "hg19" ]; then
	DBFOLDER_GATK_URL="https://data.broadinstitute.org/snowman/hg19/variant_calling/vqsr_resources/Exome/v2"
fi;
export DBFOLDER_GATK_URL_DEFAULT
export DBFOLDER_GATK_URL_DBSNP

# Gatk resources folder
DBFOLDER_GATK=$DBFOLDER/gatk

# Check resources
GATK_DATABASES_SNP_LIST=$(echo $VARIANTRECALIBRATION_SNP_RESOURCES | tr '\t' ' ' | tr "-" "\n" | sed 's/resource:\([^,]*\),[^ ]* \(.*\)/\1:\2/')
GATK_DATABASES_INDEL_LIST=$(echo $VARIANTRECALIBRATION_INDEL_RESOURCES | tr '\t' ' ' | tr "-" "\n" | sed 's/resource:\([^,]*\),[^ ]* \(.*\)/\1:\2/')
GATK_DATABASES_LIST=$(echo "$GATK_DATABASES_SNP_LIST$GATK_DATABASES_INDEL_LIST" | sort -u)
VARIANTRECALIBRATION_CHECK=1
for GATK_RESOURCE in $GATK_DATABASES_LIST; do
    if [ ! -e $DBFOLDER_GATK/current/$ASSEMBLY/$(echo $GATK_RESOURCE | cut -d: -f2) ]; then
        echo "#[WARNING] No GATK DATABASES '$GATK_RESOURCE' in the database. Recalibration step impossible. Please check '$DBFOLDER_GATK' folder or configuration file" >>/dev/stderr
	    VARIANTRECALIBRATION_CHECK=0
    fi;
done;
export VARIANTRECALIBRATION_CHECK

# Create options
VARIANTRECALIBRATION_SNP_RESOURCES_OPTION=""
if ! (($VARIANTRECALIBRATION_CHECK)); then
    echo "#[WARNING] Missing GATK DATABASES for recalibration in the database. Recalibration step impossible. Please check '$DBFOLDER_GATK' folder or configuration file" >>/dev/stderr
else
    VARIANTRECALIBRATION_SNP_RESOURCES_OPTION=$VARIANTRECALIBRATION_SNP_RESOURCES
    for GATK_RESOURCE in $GATK_DATABASES_SNP_LIST; do
        VARIANTRECALIBRATION_SNP_RESOURCES_OPTION=$(echo $VARIANTRECALIBRATION_SNP_RESOURCES_OPTION | sed 's#'$(echo $GATK_RESOURCE | cut -d: -f2)'#'$DBFOLDER_GATK/current/$ASSEMBLY/$(echo $GATK_RESOURCE | cut -d: -f2)'#')
    done;
    VARIANTRECALIBRATION_INDEL_RESOURCES_OPTION=$VARIANTRECALIBRATION_INDEL_RESOURCES
    for GATK_RESOURCE in $GATK_DATABASES_INDEL_LIST; do
        VARIANTRECALIBRATION_INDEL_RESOURCES_OPTION=$(echo $VARIANTRECALIBRATION_INDEL_RESOURCES_OPTION | sed 's#'$(echo $GATK_RESOURCE | cut -d: -f2)'#'$DBFOLDER_GATK/current/$ASSEMBLY/$(echo $GATK_RESOURCE | cut -d: -f2)'#')
    done;
fi;
export VARIANTRECALIBRATION_SNP_RESOURCES_OPTION
export VARIANTRECALIBRATION_INDEL_RESOURCES_OPTION


##########
# HOWARD #
##########

# Main Folder for HOWARD configuration
export HOWARD_FOLDER_CONFIG=$STARK_FOLDER_CONFIG/howard

if [ -z $HOWARD_CONFIG ] || [ ! -e $HOWARD_CONFIG ]; then
	HOWARD_CONFIG=$HOWARD_FOLDER_CONFIG/config.ini			# INI
fi;
export HOWARD_CONFIG
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" HOWARD_CONFIG"

if [ -z $HOWARD_CONFIG_ANNOTATION ] || [ ! -e $HOWARD_CONFIG_ANNOTATION ]; then
	HOWARD_CONFIG_ANNOTATION=$HOWARD_FOLDER_CONFIG/config.annotation.ini			# INI
fi;
export HOWARD_CONFIG_ANNOTATION
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" HOWARD_CONFIG_ANNOTATION"

if [ -z $HOWARD_CONFIG_PRIORITIZATION ] || [ ! -e $HOWARD_CONFIG_PRIORITIZATION ]; then
	HOWARD_CONFIG_PRIORITIZATION=$HOWARD_FOLDER_CONFIG/config.prioritization.ini			# INI
fi;
export HOWARD_CONFIG_PRIORITIZATION
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" HOWARD_CONFIG_PRIORITIZATION"

# DEJAVU ANNOVAR annotation
if [ "$HOWARD_CONFIG_DEJAVU_ANNOTATION" == "" ] || [ -z $HOWARD_CONFIG_DEJAVU_ANNOTATION ] || [ ! -e $HOWARD_CONFIG_DEJAVU_ANNOTATION ]; then
	HOWARD_CONFIG_DEJAVU_ANNOTATION=$HOWARD_FOLDER_CONFIG/config.annotation.ini			# INI
fi;
export HOWARD_CONFIG_DEJAVU_ANNOTATION
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" HOWARD_CONFIG_DEJAVU_ANNOTATION"

# SNPEFF config
if [ -z $SNPEFF_CONFIG ] || [ ! -e $SNPEFF_CONFIG ]; then
	SNPEFF_CONFIG=$SNPEFF_FOLDER/snpeff.config
fi;
export SNPEFF_CONFIG
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" SNPEFF_CONFIG"