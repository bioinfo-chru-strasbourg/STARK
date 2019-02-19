#!/bin/bash
#################################
## STARK environment
#################################

# DATABASES
#############


# refGene
############

DBFOLDER_REFGENE=$DBFOLDER/refGene
export REFSEQ_GENES=$DBFOLDER_REFGENE/refGene.$ASSEMBLY.bed



# dbSNP and other variant sets
#################################

# Mandatory DB (for calling with GATK variant recalibration)

DBFOLDER_DBSNP=$DBFOLDER/dbsnp

# BDSNP DB (used for calling, espacially with GATK)
export VCFDBSNP=$DBFOLDER_DBSNP/dbsnp.$ASSEMBLY.vcf.gz	#snp138.vcf.gz
if [ ! -e $VCFDBSNP ]; then
	echo "#[WARNING] No VCFDBSNP '$VCFDBSNP' in the database. Calling step impossible. Please check '$DBFOLDER' folder or configuration file"
	#exit 1;
fi;

# VCF DB for recalibration
export VCFDBSNP_RECALIBRATION=$VCFDBSNP
if (($VARIANT_RECALIBRATION)) && [ ! -e $VCFDBSNP_RECALIBRATION ]; then
	echo "#[ERROR] No VCFDBSNP_RECALIBRATION '$VCFDBSNP_RECALIBRATION' in the database. Calling step impossible. Please check '$DBFOLDER' folder or configuration file"
	exit 1;
fi;

# Other databases
DBFOLDER_OTHERS=$DBFOLDER/others
export VCFMILLS1000G=$DBFOLDER_OTHERS/Mills_and_1000G_gold_standard.indels.$ASSEMBLY.sites.vcf
export COSMIC=$DBFOLDER_OTHERS/cosmic.$ASSEMBLY.vcf
export KNOWN_ALLELES=$VCFMILLS1000G
export HAPMAP=$DBFOLDER_OTHERS/hapmap_3.3.$ASSEMBLY.sites.vcf
export OMNI=$DBFOLDER_OTHERS/1000G_omni2.5.$ASSEMBLY.vcf
export PHASE1_1000G=$DBFOLDER_OTHERS/1000G_phase1.snps.high_confidence.$ASSEMBLY.sites.vcf



# HOWARD ANNOTATION/PRIOTITIZATION/TRANSLATION CONFIGURATION
##############################################################

# Configuration
#################

if [ -z $HOWARD_CONFIG ] || [ ! -e $HOWARD_CONFIG ]; then
	HOWARD_CONFIG=$HOWARDDIR/config.ini			# INI
fi;
export HOWARD_CONFIG

if [ -z $HOWARD_CONFIG_PRIORITIZATION ] || [ ! -e $HOWARD_CONFIG_PRIORITIZATION ]; then
	HOWARD_CONFIG_PRIORITIZATION=$HOWARDDIR/config.prioritization.ini			# INI
fi;
export HOWARD_CONFIG_PRIORITIZATION

if [ -z $HOWARD_CONFIG_ANNOTATION ] || [ ! -e $HOWARD_CONFIG_ANNOTATION ]; then
	HOWARD_CONFIG_ANNOTATION=$HOWARDDIR/config.annotation.ini			# INI
fi;
export HOWARD_CONFIG_ANNOTATION

# SNPEFF
if [ -z $SNPEFF_CONFIG ] || [ ! -e $SNPEFF_CONFIG ]; then
	SNPEFF_CONFIG=$SNPEFF_FOLDER/snpeff.config
fi;
export SNPEFF_CONFIG


# ANNOVAR / SNPEFF Databases  Configuration
#############################################

# ANNOVAR
if [ -z $ANNOVAR_DATABASES ] || [ ! -e $ANNOVAR_DATABASES ]; then
	ANNOVAR_DATABASES=$DBFOLDER/annovar
fi;
export ANNOVAR_DATABASES

#export SNPEFF_CONFIG=$SNPEFF_FOLDER/snpeff.config	# CONFIG # NOT USED !!! # CHANGE CONFIG FILE in SNPEFF TOOL if necessary
if [ -z $SNPEFF_DATABASES ] || [ ! -e $SNPEFF_DATABASES ]; then
	SNPEFF_DATABASES=$DBFOLDER/snpeff/$SNPEFF_VERSION
fi;
export SNPEFF_DATABASES
#export SNPEFF_DATABASES=$DBFOLDER/snpeff_sources/4.3t	# DATA # NOT USED !!! # CHANGE DATABASE location in CONFIG FILE in SNPEFF TOOL if necessary




