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
DBFOLDER_GATK=$DBFOLDER/gatk
GATK_VARIANT_RECALIBRATION_CHECK=1

if [ -z "$GATK_VARIANT_RECALIBRATION_SNP_RESOURCES" ]; then
    GATK_VARIANT_RECALIBRATION_SNP_RESOURCES="
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.vcf.gz
        -resource:omni,known=false,training=true,truth=true,prior=12.0 1000G_omni2.5.b37.vcf.gz
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.b37.vcf.gz 
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.b37.vcf.gz
        "
    echo "#[WARNING] No GATK DATABASES defined for Variant Recalibration for SNP. Default configuration: $GATK_VARIANT_RECALIBRATION_SNP_RESOURCES" >>/dev/stderr
fi;
export GATK_VARIANT_RECALIBRATION_SNP_RESOURCES

if [ -z "$GATK_VARIANT_RECALIBRATION_INDEL_RESOURCES" ]; then
    GATK_VARIANT_RECALIBRATION_INDEL_RESOURCES="
        -resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.b37.vcf.gz
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.b37.vcf.gz
    "
    echo "#[WARNING] No GATK DATABASES defined for Variant Recalibration for INDEL. Default configuration: $GATK_VARIANT_RECALIBRATION_INDEL_RESOURCES" >>/dev/stderr
fi;
export GATK_VARIANT_RECALIBRATION_INDEL_RESOURCES

GATK_DATABASES_SNP_LIST=$(echo $GATK_VARIANT_RECALIBRATION_SNP_RESOURCES | tr '\t' ' ' | tr "-" "\n" | sed 's/resource:\([^,]*\),[^ ]* \(.*\)/\1:\2/')
GATK_DATABASES_INDEL_LIST=$(echo $GATK_VARIANT_RECALIBRATION_INDEL_RESOURCES | tr '\t' ' ' | tr "-" "\n" | sed 's/resource:\([^,]*\),[^ ]* \(.*\)/\1:\2/')
GATK_DATABASES_LIST=$(echo "$GATK_DATABASES_SNP_LIST$GATK_DATABASES_INDEL_LIST" | sort -u)

for GATK_RESOURCE in $GATK_DATABASES_LIST; do
    #echo $DBFOLDER_GATK/current/$ASSEMBLY/$(echo $GATK_RESOURCE | cut -d: -f2)
    if [ ! -e $DBFOLDER_GATK/current/$ASSEMBLY/$(echo $GATK_RESOURCE | cut -d: -f2) ]; then
        echo "#[WARNING] No GATK DATABASES '$GATK_RESOURCE' in the database. Recalibration step impossible. Please check '$DBFOLDER_GATK' folder or configuration file" >>/dev/stderr
	    GATK_VARIANT_RECALIBRATION_CHECK=0
    fi;
done;
export GATK_VARIANT_RECALIBRATION_CHECK

if ! (($GATK_VARIANT_RECALIBRATION_CHECK)); then
    echo "#[WARNING] Missing GATK DATABASES for recalibration in the database. Recalibration step impossible. Please check '$DBFOLDER_GATK' folder or configuration file" >>/dev/stderr
fi;

# # VCF DB for recalibration
# export VCFDBSNP_RECALIBRATION=$VCFDBSNP
# if (($VARIANT_RECALIBRATION)) && [ ! -e $VCFDBSNP_RECALIBRATION ]; then
# 	echo "#[ERROR] No VCFDBSNP_RECALIBRATION '$VCFDBSNP_RECALIBRATION' in the database. Calling step impossible. Please check '$DBFOLDER' folder or configuration file" >>/dev/stderr
# 	exit 1;
# fi;
# DATABASES_LIST=$DATABASES_LIST" VCFDBSNP_RECALIBRATION"

# # Other databases
# DBFOLDER_OTHERS=$DBFOLDER/others
# export VCFMILLS1000G=$DBFOLDER_OTHERS/Mills_and_1000G_gold_standard.indels.$ASSEMBLY.sites.vcf
# export COSMIC=$DBFOLDER_OTHERS/cosmic.$ASSEMBLY.vcf
# export KNOWN_ALLELES=$VCFMILLS1000G
# export HAPMAP=$DBFOLDER_OTHERS/hapmap_3.3.$ASSEMBLY.sites.vcf
# export OMNI=$DBFOLDER_OTHERS/1000G_omni2.5.$ASSEMBLY.vcf
# export PHASE1_1000G=$DBFOLDER_OTHERS/1000G_phase1.snps.high_confidence.$ASSEMBLY.sites.vcf
# DATABASES_LIST=$DATABASES_LIST" VCFMILLS1000G"
# DATABASES_LIST=$DATABASES_LIST" COSMIC"
# DATABASES_LIST=$DATABASES_LIST" KNOWN_ALLELES"
# DATABASES_LIST=$DATABASES_LIST" HAPMAP"
# DATABASES_LIST=$DATABASES_LIST" OMNI"
# DATABASES_LIST=$DATABASES_LIST" PHASE1_1000G"



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

