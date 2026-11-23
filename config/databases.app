#!/bin/bash
#################################
## STARK environment
#################################
# DATABASES
#############
DATABASES_LIST=""
DATABASES_CONFIG_LIST=""

#ASSEMBLY=hg38

#####################
# ASSEMBLY & GENOME #
#####################

# Genome regex for chromosome names.
# Default: 'chr[0-9XYM]+\$$' (for canonical contig names like chr1, chr2, chrX, chrY, chrM)
GENOME_REGEX="'chr[0-9XYM]+\$$'"
export GENOME_REGEX

# if [ -z $ASSEMBLY ] || [ "$ASSEMBLY" == "" ]; then
# 	ASSEMBLY=hg19
# fi;
# export ASSEMBLY

# if [ -z $GENOME ] || [ "$GENOME" == "" ]; then
# 	GENOME=$DBFOLDER/genomes/current/$ASSEMBLY/$ASSEMBLY.fa
# fi;
# export GENOME

# if [ -z $DICT ] || [ "$DICT" == "" ]; then
# 	DICT=$DBFOLDER/genomes/current/$ASSEMBLY/$ASSEMBLY.dict
# fi;
# export DICT

# Assembly
# Default: hg19 ; can be changed by application inheritance (e.g. GENOME.app inherit from default.app)
if [ -z "${ASSEMBLY:-}" ]; then
	ASSEMBLY=hg19
fi
export ASSEMBLY

# Genome
# Default: depend on ASSEMBLY ; can be changed by application inheritance (e.g. GENOME.app inherit from default.app)
if [ -z "${GENOME+x}" ] || [ "$ASSEMBLY" != "$_GENOME_ASSEMBLY" ]; then
	GENOME="$DBFOLDER/genomes/current/$ASSEMBLY/$ASSEMBLY.fa"
 	_GENOME_ASSEMBLY="$ASSEMBLY"
fi
export GENOME

# Dict
# Default: depend on ASSEMBLY ; can be changed by application inheritance (e.g. GENOME.app inherit from default.app)
if [ -z "${DICT+x}" ] || [ "$ASSEMBLY" != "$_DICT_ASSEMBLY" ]; then
	DICT="$DBFOLDER/genomes/current/$ASSEMBLY/$ASSEMBLY.dict"
	_DICT_ASSEMBLY="$ASSEMBLY"
fi
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
if [ -z $ARRIBA_URL ] || [ "$ARRIBA_URL" == "" ]; then
	ARRIBA_URL="https://github.com/suhrig/arriba/releases/download"
fi;
export ARRIBA_URL
if [ -z $ARRIBA_RELEASE ] || [ "$ARRIBA_RELEASE" == "" ]; then
	ARRIBA_RELEASE="v2.4.0"
fi;
export ARRIBA_RELEASE
if [ ! -z $FOLDER_DATABASES_ARRIBA ] && [ "$FOLDER_DATABASES_ARRIBA" != "" ]; then
	ARRIBA_DATABASES=$FOLDER_DATABASES_ARRIBA
else
	ARRIBA_DATABASES=$DBFOLDER/arriba/current
fi;
export ARRIBA_DATABASES
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" ARRIBA_DATABASES"

# PFAM
if [ -z $PFAM_URL ] || [ "$PFAM_URL" == "" ]; then
	PFAM_URL="ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases";
fi;
export PFAM_URL
if [ -z $PFAM_RELEASE ] || [ "$PFAM_RELEASE" == "" ]; then
	PFAM_RELEASE="38.2"
fi;
export PFAM_RELEASE
if [ ! -z $FOLDER_DATABASES_PFAM ] && [ "$FOLDER_DATABASES_PFAM" != "" ]; then
	PFAM_DATABASES=$FOLDER_DATABASES_PFAM
else
	PFAM_DATABASES=$DBFOLDER/PFAM/current
fi;
export PFAM_DATABASES
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" PFAM_DATABASES"

# DFAM
if [ -z $DFAM_URL ] || [ "$DFAM_URL" == "" ]; then
	DFAM_URL="https://dfam.org/releases";
fi;
export DFAM_URL
if [ -z $DFAM_RELEASE ] || [ "$DFAM_RELEASE" == "" ]; then
	DFAM_RELEASE="3.1"
fi;
export DFAM_RELEASE
if [ ! -z $FOLDER_DATABASES_DFAM ] && [ "$FOLDER_DATABASES_DFAM" != "" ]; then
	DFAM_DATABASES=$FOLDER_DATABASES_DFAM
else
	DFAM_DATABASES=$DBFOLDER/DFAM/current
fi;
export DFAM_DATABASES
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" DFAM_DATABASES"


# DEEPVARIANT
# Databases for DeepVariant

if [ ! -z $FOLDER_DATABASES_DEEPVARIANT ] && [ "$FOLDER_DATABASES_DEEPVARIANT" != "" ]; then
	DEEPVARIANT_DATABASES=$FOLDER_DATABASES_DEEPVARIANT
else
	DEEPVARIANT_DATABASES=$DBFOLDER/deepvariant/current
fi;
export DEEPVARIANT_DATABASES
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" DEEPVARIANT_DATABASES"


# Default Haploid contigs
# Heterozygous variants in these contigs will be re-genotyped as the most likely (e.g., hom ref or hom alt) because of the haploid nature of these contigs.
# This option is useful to avoid false positive heterozygous calls in haploid contigs, especially in PAR regions.
# For human genome, but can be different for other organisms.
# Default: "chrX,chrY" 
if [ -z $DEEPVARIANT_HAPLOID_CONTIGS ] || [ "$DEEPVARIANT_HAPLOID_CONTIGS" == "" ]; then
	DEEPVARIANT_HAPLOID_CONTIGS="chrX,chrY"
fi;
export DEEPVARIANT_HAPLOID_CONTIGS

# Default PAR BED file for DeepVariant
# If --haploid_contigs is set, then this can be used to provide PAR regions to be excluded from genotype adjustment.
# It will look for a unique BED file in the DeepVariant databases folder for the given assembly.

# Find Default PAR BED
_NEW_DEEPVARIANT_PAR_BED_DEFAULT=
shopt -s nullglob
PAR_BED_FILES=("$DEEPVARIANT_DATABASES/$ASSEMBLY"/*.bed)
PAR_BED_FILES_PAR=("$DEEPVARIANT_DATABASES/$ASSEMBLY"/*_PAR.bed)
shopt -u nullglob
# Determine the default PAR BED file if there is exactly one match in the database folder
if [ ${#PAR_BED_FILES[@]} -eq 1 ]; then
	_NEW_DEEPVARIANT_PAR_BED_DEFAULT=${PAR_BED_FILES[0]}
elif [ ${#PAR_BED_FILES_PAR[@]} -eq 1 ]; then
	_NEW_DEEPVARIANT_PAR_BED_DEFAULT=${PAR_BED_FILES_PAR[0]}
fi;

# Check if DEEPVARIANT_PAR_BED is set and if it matches the default
if [ -z "${DEEPVARIANT_PAR_BED+x}" ] || [[ -z "$DEEPVARIANT_PAR_BED" ]] || [ "${DEEPVARIANT_PAR_BED}" = "${_DEEPVARIANT_PAR_BED_DEFAULT:-}" ]; then
    DEEPVARIANT_PAR_BED="$_NEW_DEEPVARIANT_PAR_BED_DEFAULT"
fi
_DEEPVARIANT_PAR_BED_DEFAULT="$_NEW_DEEPVARIANT_PAR_BED_DEFAULT"
export DEEPVARIANT_PAR_BED

#echo "DEEPVARIANT_PAR_BED is set to: $DEEPVARIANT_PAR_BED"


# _NEW_GENOME_RNA_DEFAULT="$GENOME.ctat/$CTAT_DATABASES_GENE_SOURCE/ref_genome.fa"
# if [ -z "${GENOME_RNA+x}" ] || [ "${GENOME_RNA}" = "${_GENOME_RNA_DEFAULT:-}" ]; then
#     GENOME_RNA="$_NEW_GENOME_RNA_DEFAULT"
# fi
# _GENOME_RNA_DEFAULT="$_NEW_GENOME_RNA_DEFAULT"
# export GENOME_RNA


# CTAT
#######

# STAR Fusion (CTAT lib source). This is the source of the CTAT lib used for STAR fusion. It can be used to build a custom CTAT lib with different gene annotation (refGene or gencode) or different genome version (hg19 or hg38). By default, the CTAT lib used for STAR fusion will be based on refGene annotation and will be located in $CTAT_URL_SOURCE/GRCh37_gencode_v19_CTAT_lib_Mar012021.source.tar.gz for hg19 and $CTAT_URL_SOURCE/GRCh38_gencode_v44_CTAT_lib_Oct292023.source.tar.gz for hg38. If you want to use a different CTAT lib for STAR fusion, you can set the CTAT_LIB_SOURCE variable to the URL of the CTAT lib source to use for STAR fusion. Note that if you change the CTAT lib source for STAR fusion, you should also change the gene source for CTAT lib (CTAT_DATABASES_GENE_SOURCE) to the same value to ensure consistency between fusion detection and annotation.
CTAT_URL_SOURCE="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB"
export CTAT_URL_SOURCE
if [ $ASSEMBLY == "hg19" ] ; then CTAT_LIB_SOURCE="$CTAT_URL_SOURCE/GRCh37_gencode_v19_CTAT_lib_Mar012021.source.tar.gz"; fi;
if [ $ASSEMBLY == "hg38" ] ; then CTAT_LIB_SOURCE="$CTAT_URL_SOURCE/GRCh38_gencode_v44_CTAT_lib_Oct292023.source.tar.gz"; fi;
export CTAT_LIB_SOURCE

# CTAT_DATABASES_GENE_SOURCE can be "refgene" or "gencode". If "refgene" the gene annotation used for CTAT lib will be based on refGene (NCBI RefSeq) annotation. If "gencode" the gene annotation used for CTAT lib will be based on gencode annotation. The choice of gene annotation can impact the fusion detection results, especially for fusions involving genes that are not well annotated in refGene but are annotated in gencode. By default, CTAT_DATABASES_GENE_SOURCE is set to "refgene" because refGene annotation is more conservative and less likely to include dubious gene models that could lead to false positive fusion calls. However, if you want to use gencode annotation for CTAT lib, you can set CTAT_DATABASES_GENE_SOURCE to "gencode". Note that if you change the gene source for CTAT lib, you should also change the gene source for STAR fusion annotation (STAR_FUSION_GENE_SOURCE) to the same value to ensure consistency between fusion detection and annotation.
# Either "refgene" or "gencode"
if [ -z $CTAT_DATABASES_GENE_SOURCE ] || [ "$CTAT_DATABASES_GENE_SOURCE" == "" ]; then
	CTAT_DATABASES_GENE_SOURCE="refgene"
fi;
export CTAT_DATABASES_GENE_SOURCE

# Genome for STAR alignement. By default, the genome used for STAR alignement in STAR fusion will be the same as the genome used for DNA alignement (GENOME variable). However, if you want to use a different genome for STAR alignement in STAR fusion, you can set the GENOME_RNA variable to the path of the genome fasta file to use for STAR alignement in STAR fusion. Note that if you change the genome for STAR alignement in STAR fusion, you should also change the genome for STAR fusion annotation (STAR_FUSION_GENOME) to the same value to ensure consistency between fusion detection and annotation. This genome also depend on the gene source for CTAT lib (CTAT_DATABASES_GENE_SOURCE) because the genome used for STAR alignement in STAR fusion should be the same as the genome used for the gene annotation in CTAT lib to ensure consistency between fusion detection and annotation. By default, GENOME_RNA is set to the same genome as GENOME variable, but with a different path that depends on the gene source for CTAT lib (CTAT_DATABASES_GENE_SOURCE). If CTAT_DATABASES_GENE_SOURCE is set to "refgene", the genome used for STAR alignement in STAR fusion will be based on refGene annotation and will be located in $GENOME.ctat/ref_genome.fa. If CTAT_DATABASES_GENE_SOURCE is set to "gencode", the genome used for STAR alignement in STAR fusion will be based on gencode annotation and will be located in $GENOME.ctat/gencode/ref_genome.fa.
# Default genome for STAR alignement in STAR fusion is based on refGene annotation. Use to switch path if gene source change by application inheritance (e.g. RNASEQ.app inherit from default.app)

_NEW_GENOME_RNA_DEFAULT="$GENOME.ctat/$CTAT_DATABASES_GENE_SOURCE/ref_genome.fa"
if [ -z "${GENOME_RNA+x}" ] || [ "${GENOME_RNA}" = "${_GENOME_RNA_DEFAULT:-}" ]; then
    GENOME_RNA="$_NEW_GENOME_RNA_DEFAULT"
fi
_GENOME_RNA_DEFAULT="$_NEW_GENOME_RNA_DEFAULT"
export GENOME_RNA


# Gencode
##########

# for hg19 the last gencode version is v19
if [ $ASSEMBLY == "hg19" ] ; then
	# GENCODE_VERSION="19"
	# GENCODE_CURRENT="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$GENCODE_VERSION/gencode.v$GENCODE_VERSION.annotation.gtf.gz";
	GENCODE_VERSION="49"
	GENCODE_CURRENT="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/GRCh37_mapping/gencode.v${GENCODE_VERSION}lift37.annotation.gtf.gz"; #gencode.v$GENCODE_VERSION.annotation.gtf.gz";
fi;
# for hg38 the first gencode version is v20 ; current version (10/2023) is v44
if [ $ASSEMBLY == "hg38" ] ; then 
	GENCODE_VERSION="49"
	GENCODE_CURRENT="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.primary_assembly.annotation.gtf.gz";
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
##########

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
##################

if [ ! -z $FOLDER_DATABASES_DEJAVU_ANNOVAR ] && [ "$FOLDER_DATABASES_DEJAVU_ANNOVAR" != "" ]; then
	DEJAVU_ANNOVAR_DATABASES=$FOLDER_DATABASES_DEJAVU_ANNOVAR
else
	DEJAVU_ANNOVAR_DATABASES=$DBFOLDER/annovar/current
fi;
export DEJAVU_ANNOVAR_DATABASES
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" DEJAVU_ANNOVAR_DATABASES"

# SNPEFF
##########

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
#########

if [ ! -z $FOLDER_DATABASES_DBSNP ] && [ "$FOLDER_DATABASES_DBSNP" != "" ]; then
	DBSNP_DATABASES=$FOLDER_DATABASES_DBSNP
else
	DBSNP_DATABASES=$DBFOLDER/dbsnp/current
fi;
export DBSNP_DATABASES

# Version of DBSNP to download ex "b156"
DBSNP_VERSION_DOWNLOAD="b157"
export DBSNP_VERSION_DOWNLOAD


# Version of DBSNP to use for GATK tools (realignement, recalibration, calling)
DBSNP_BUILDID="151"

# export DBSNP_VERSION
# export dbSNPBuildID
export DBSNP_BUILDID

DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" DBSNP_DATABASES"
DBFOLDER_DBSNP=$DBFOLDER/dbsnp

# VCF dbSNP for GATK tools (realignement, recalibration, calling)
if [ -z "${VCFDBSNP+x}" ] || [ "$ASSEMBLY" != "$_VCFDBSNP_ASSEMBLY" ]; then
	VCFDBSNP="$DBFOLDER_DBSNP/current/$ASSEMBLY/$DBSNP_VERSION_DOWNLOAD/dbsnp.b$DBSNP_BUILDID.vcf.gz"
 	_VCFDBSNP_ASSEMBLY="$ASSEMBLY"
fi
export VCFDBSNP

if [ ! -e $VCFDBSNP ]; then
	echo "#[WARNING] No VCFDBSNP '$VCFDBSNP' in the database. Calling step impossible. Please check '$DBFOLDER' folder or configuration file" >>/dev/stderr
fi;
DATABASES_LIST=$DATABASES_LIST" VCFDBSNP"



# dbNSFP
##########

if [ ! -z $FOLDER_DATABASES_DBNSFP ] && [ "$FOLDER_DATABASES_DBNSFP" != "" ]; then
	DBNSFP_DATABASES=$FOLDER_DATABASES_DBNSFP
else
	DBNSFP_DATABASES=$DBFOLDER/dbnsfp/current
fi;
export DBNSFP_DATABASES
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" DBNSFP_DATABASES"

# refGene/refSeq
#################

DBFOLDER_REFGENE=$DBFOLDER/refGene
export DBFOLDER_REFGENE
if [ ! -z $FOLDER_DATABASES_REFGENE ] && [ "$FOLDER_DATABASES_REFGENE" != "" ]; then
	REFGENE_DATABASES=$FOLDER_DATABASES_REFGENE
else
	REFGENE_DATABASES=$DBFOLDER_REFGENE/current
fi;
export REFGENE_DATABASES
#export REFSEQ_GENES=$DBFOLDER_REFGENE/current/$ASSEMBLY/refGene.$ASSEMBLY.bed
export REFSEQ_GENES=$REFGENE_DATABASES/$ASSEMBLY/ncbiRefSeq.bed
export REFSEQ_GENES_GTF=$REFGENE_DATABASES/$ASSEMBLY/ncbiRefSeq.gtf

##################
# GATK DATABASES #
##################

# # GATK VARIANT RECALIBRATION URLs
DBFOLDER_GATK_URL_DEFAULT=ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle
export DBFOLDER_GATK_URL_DEFAULT


# Gatk resources folder
DBFOLDER_GATK=$DBFOLDER/gatk

# Main Folder for GATK configuration
GATK_FOLDER_CONFIG=$STARK_FOLDER_CONFIG/gatk
export GATK_FOLDER_CONFIG


# VARIANTRECALIBRATION_SNP_RESOURCES
VARIANTRECALIBRATION_SNP_RESOURCES=""
# Read file containing GATK variant recalibration resources for the current assembly
if [ -f $GATK_FOLDER_CONFIG/resources.$ASSEMBLY.snp.tsv ]; then
    # Read each line
	VARIANTRECALIBRATION_SNP_RESOURCES=""
	while IFS= read -r line; do
		VARIANTRECALIBRATION_SNP_RESOURCES="$VARIANTRECALIBRATION_SNP_RESOURCES -resource:$line"
	done < $GATK_FOLDER_CONFIG/resources.$ASSEMBLY.snp.tsv
fi;
export VARIANTRECALIBRATION_SNP_RESOURCES

# VARIANTRECALIBRATION_INDEL_RESOURCES
VARIANTRECALIBRATION_INDEL_RESOURCES=""
# Read file containing GATK variant recalibration resources for the current assembly
if [ -f $GATK_FOLDER_CONFIG/resources.$ASSEMBLY.indel.tsv ]; then
    # Read each line
	VARIANTRECALIBRATION_INDEL_RESOURCES=""
	while IFS= read -r line; do
		VARIANTRECALIBRATION_INDEL_RESOURCES="$VARIANTRECALIBRATION_INDEL_RESOURCES -resource:$line"
	done < $GATK_FOLDER_CONFIG/resources.$ASSEMBLY.indel.tsv
fi;
export VARIANTRECALIBRATION_INDEL_RESOURCES



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


# GATK BAM Realignment
# Realignment is a process of correcting misalignments around indels. It is recommended to perform realignment before variant calling, especially for indel calling, to improve the accuracy of variant calling (not for GATK4 Haplotype Caller and MuTect2).
# This option allow GATK to perform realignment, if the file is available and not empty. If not, GATK will be used without realignment.
# Realignemnt is performed with GATK3 RealignerTargetCreator and IndelRealigner, and not performed with GATK4 (GATK4 Haplotype Caller and MuTect2 perform local realignment during calling, so GATK4 do not include Realignemnt tools).
# This option is used with the GATK3 RealignerTargetCreator command to identify regions to realign, and with the GATK3 IndelRealigner command to perform realignment.
# The VCF files include Indels and can be the same as for VCF recalibration (see VARIANTRECALIBRATION_INDEL_RESOURCES).
# Exemple:
# - GATK_REALIGNMENT_KNOWN_OPTIONS="--known $VCFDBSNP"
# default:
GATK_REALIGNMENT_KNOWN_OPTIONS=$(echo -e "$VARIANTRECALIBRATION_INDEL_RESOURCES" | sed "s#-resource.* \(.*\)#-known "$DBFOLDER"/gatk/current/"$ASSEMBLY"/\1#gi" | xargs echo)
export GATK_REALIGNMENT_KNOWN_OPTIONS

# GATK BAM Recalibration
# Recalibration is a process of correcting base quality scores. It is recommended to perform recalibration before variant calling to improve the accuracy of variant calling.
# This option allow GATK to perform recalibration, if the file is available and not empty. If not, GATK will be used without recalibration.
# Recalibration is performed with GATK4 BaseRecalibrator and PrintReads
# This option is used with the GATK4 BaseRecalibrator command to identify covariates to recalibrate, and with the GATK4 PrintReads command to perform recalibration.
# The VCF files include SNPs and can be the same as for VCF recalibration (see VARIANTRECALIBRATION_SNP_RESOURCES).
# Exemple:
# - GATK_RECALIBRATION_KNOWN_OPTIONS="--known-sites $VCFDBSNP"
# default:
GATK_RECALIBRATION_KNOWN_OPTIONS=$(echo -e "$VARIANTRECALIBRATION_SNP_RESOURCES" | sed "s#-resource.* \(.*\)#--known-sites "$DBFOLDER"/gatk/current/"$ASSEMBLY"/\1#gi" | xargs echo)
export GATK_RECALIBRATION_KNOWN_OPTIONS



##########
# SNPEFF #
##########

# Main Folder for HOWARD configuration
export SNPEFF_FOLDER_CONFIG=$STARK_FOLDER_CONFIG/snpeff

# SNPEFF config
if [ -z $SNPEFF_CONFIG ] || [ ! -e $SNPEFF_CONFIG ]; then
	if [ -e $SNPEFF_FOLDER_CONFIG/snpeff.config ]; then
		SNPEFF_CONFIG=$SNPEFF_FOLDER_CONFIG/snpEff.config
	elif [ ! -e $SNPEFF_FOLDER/snpeff.config ]; then
		SNPEFF_CONFIG=$SNPEFF_FOLDER/snpEff.config
	elif [ ! -e $(dirname $SNPEFF)/snpEff.config ]; then
		SNPEFF_CONFIG=$(dirname $SNPEFF)/snpEff.config
	fi;
fi;
export SNPEFF_CONFIG
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" SNPEFF_CONFIG"


##################
# VARIANTCONVERT #
##################

# Main Folder for VARIANTCONVERT configuration
export VARIANTCONVERT_FOLDER_CONFIG=$STARK_FOLDER_CONFIG/variantconvert



##########
# HOWARD #
##########

# Main Folder for HOWARD configuration
export HOWARD_FOLDER_CONFIG=$STARK_FOLDER_CONFIG/howard


############### Main Folder for HOWARD configuration
export HOWARD_FOLDER_CONFIG=$STARK_FOLDER_CONFIG/howard

if [ -z $HOWARD_CONFIG ] || [ ! -e $HOWARD_CONFIG ]; then
	HOWARD_CONFIG=$HOWARD_FOLDER_CONFIG/config.json	
fi;
export HOWARD_CONFIG
DATABASES_CONFIG_LIST=$DATABASES_CONFIG_LIST" HOWARD_CONFIG"
