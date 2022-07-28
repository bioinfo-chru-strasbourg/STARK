#!/bin/bash
#################################
## STARK environment
#################################

# DATABASES
#############

DATABASES_LIST=""
DATABASES_CONFIG_LIST=""

# refGene
############

DBFOLDER_REFGENE=$DBFOLDER/refGene
#export REFSEQ_GENES=$DBFOLDER_REFGENE/refGene.$ASSEMBLY.bed
export REFSEQ_GENES=$DBFOLDER_REFGENE/current/refGene.$ASSEMBLY.bed



# dbSNP and other variant sets
#################################

# Mandatory DB (for calling with GATK variant recalibration)

DBFOLDER_DBSNP=$DBFOLDER/dbsnp

# BDSNP DB (used for calling, espacially with GATK)
#export VCFDBSNP=$DBFOLDER_DBSNP/dbsnp.$ASSEMBLY.vcf.gz	#snp138.vcf.gz
export VCFDBSNP=$DBFOLDER_DBSNP/current/dbsnp.$ASSEMBLY.vcf.gz	#snp138.vcf.gz
if [ ! -e $VCFDBSNP ]; then
	echo "#[WARNING] No VCFDBSNP '$VCFDBSNP' in the database. Calling step impossible. Please check '$DBFOLDER' folder or configuration file" >>/dev/stderr
	#exit 1;
fi;
DATABASES_LIST=$DATABASES_LIST" VCFDBSNP"



# GATK DATABASES
##################

# Gatk resources folder

DBFOLDER_GATK=$DBFOLDER/gatk

# SNP

if [ -z "$VARIANTRECALIBRATION_SNP_RESOURCES" ]; then
    VARIANTRECALIBRATION_SNP_RESOURCES="
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.vcf.gz
        -resource:omni,known=false,training=true,truth=true,prior=12.0 1000G_omni2.5.b37.vcf.gz
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.b37.vcf.gz 
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.b37.vcf.gz
        "
    echo "#[WARNING] No GATK DATABASES defined for Variant Recalibration for SNP. Default configuration: $VARIANTRECALIBRATION_SNP_RESOURCES" >>/dev/stderr
fi;
export VARIANTRECALIBRATION_SNP_RESOURCES

GATK_DATABASES_SNP_LIST=$(echo $VARIANTRECALIBRATION_SNP_RESOURCES | tr '\t' ' ' | tr "-" "\n" | sed 's/resource:\([^,]*\),[^ ]* \(.*\)/\1:\2/')

if [ -z "$VARIANTRECALIBRATION_SNP_ANNOTATIONS" ]; then
    #VARIANTRECALIBRATION_SNP_ANNOTATIONS="-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP"
    #VARIANTRECALIBRATION_SNP_ANNOTATIONS="-an ReadPosRankSum -an FS -an SOR -an DP"
    VARIANTRECALIBRATION_SNP_ANNOTATIONS="-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP"
fi;
export VARIANTRECALIBRATION_SNP_ANNOTATIONS

if [ -z "$VARIANTRECALIBRATION_SNP_TRANCHES" ]; then
    VARIANTRECALIBRATION_SNP_TRANCHES="-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.1 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0"
fi;
export VARIANTRECALIBRATION_SNP_TRANCHES

# INDEL

if [ -z "$VARIANTRECALIBRATION_INDEL_RESOURCES" ]; then
    VARIANTRECALIBRATION_INDEL_RESOURCES="
        -resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.b37.vcf.gz
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.b37.vcf.gz
    "
    echo "#[WARNING] No GATK DATABASES defined for Variant Recalibration for INDEL. Default configuration: $VARIANTRECALIBRATION_INDEL_RESOURCES" >>/dev/stderr
fi;
export VARIANTRECALIBRATION_INDEL_RESOURCES

GATK_DATABASES_INDEL_LIST=$(echo $VARIANTRECALIBRATION_INDEL_RESOURCES | tr '\t' ' ' | tr "-" "\n" | sed 's/resource:\([^,]*\),[^ ]* \(.*\)/\1:\2/')

if [ -z "$VARIANTRECALIBRATION_INDEL_ANNOTATIONS" ]; then
    #VARIANTRECALIBRATION_INDEL_ANNOTATIONS="-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP"
    #VARIANTRECALIBRATION_INDEL_ANNOTATIONS="-an ReadPosRankSum -an FS -an SOR -an DP"
    VARIANTRECALIBRATION_INDEL_ANNOTATIONS="-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum"
fi;
export VARIANTRECALIBRATION_INDEL_ANNOTATIONS

if [ -z "$VARIANTRECALIBRATION_INDEL_TRANCHES" ]; then
    VARIANTRECALIBRATION_INDEL_TRANCHES="-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.1 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0"
fi;
export VARIANTRECALIBRATION_INDEL_TRANCHES

# Check resources

GATK_DATABASES_LIST=$(echo "$GATK_DATABASES_SNP_LIST$GATK_DATABASES_INDEL_LIST" | sort -u)
VARIANTRECALIBRATION_CHECK=1

for GATK_RESOURCE in $GATK_DATABASES_LIST; do
    if [ ! -e $DBFOLDER_GATK/current/$(echo $GATK_RESOURCE | cut -d: -f2) ]; then
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
        VARIANTRECALIBRATION_SNP_RESOURCES_OPTION=$(echo $VARIANTRECALIBRATION_SNP_RESOURCES_OPTION | sed 's#'$(echo $GATK_RESOURCE | cut -d: -f2)'#'$DBFOLDER_GATK/current/$(echo $GATK_RESOURCE | cut -d: -f2)'#')
    done;
    VARIANTRECALIBRATION_INDEL_RESOURCES_OPTION=$VARIANTRECALIBRATION_INDEL_RESOURCES
    for GATK_RESOURCE in $GATK_DATABASES_INDEL_LIST; do
        VARIANTRECALIBRATION_INDEL_RESOURCES_OPTION=$(echo $VARIANTRECALIBRATION_INDEL_RESOURCES_OPTION | sed 's#'$(echo $GATK_RESOURCE | cut -d: -f2)'#'$DBFOLDER_GATK/current/$(echo $GATK_RESOURCE | cut -d: -f2)'#')
    done;
fi;
export VARIANTRECALIBRATION_SNP_RESOURCES_OPTION
export VARIANTRECALIBRATION_INDEL_RESOURCES_OPTION



# HOWARD ANNOTATION/PRIOTITIZATION/TRANSLATION CONFIGURATION
##############################################################

# Configuration
#################

# Main Folder for HOWARD configuration
export HOWARD_FOLDER_CONFIG=$STARK_FOLDER_CONFIG/howard

# HOWARD
if [ -z $HOWARD_CONFIG ] || [ ! -e $HOWARD_CONFIG ]; then
	HOWARD_CONFIG=$HOWARD_FOLDER_CONFIG/config.ini			# INI
fi;
export HOWARD_CONFIG
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" HOWARD_CONFIG"

# ANNOVAR
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

# DEJAVU ANNOVAR
if [ $HOWARD_CONFIG_DEJAVU_ANNOTATION == "" ] || [ -z $HOWARD_CONFIG_DEJAVU_ANNOTATION ] || [ ! -e $HOWARD_CONFIG_DEJAVU_ANNOTATION ]; then
	HOWARD_CONFIG_DEJAVU_ANNOTATION=$HOWARD_FOLDER_CONFIG/config.annotation.ini			# INI
fi;
export HOWARD_CONFIG_DEJAVU_ANNOTATION
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" HOWARD_CONFIG_DEJAVU_ANNOTATION"

# SNPEFF
if [ -z $SNPEFF_CONFIG ] || [ ! -e $SNPEFF_CONFIG ]; then
	SNPEFF_CONFIG=$SNPEFF_FOLDER/snpeff.config
fi;
export SNPEFF_CONFIG
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" SNPEFF_CONFIG"


# ANNOVAR / SNPEFF Databases  Configuration
#############################################

# ANNOVAR
if [ ! -z $FOLDER_DATABASES_ANNOVAR ] && [ "$FOLDER_DATABASES_ANNOVAR" != "" ]; then
	ANNOVAR_DATABASES=$FOLDER_DATABASES_ANNOVAR
else
	ANNOVAR_DATABASES=$DBFOLDER/annovar/current
fi;
export ANNOVAR_DATABASES
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" ANNOVAR_DATABASES"


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

