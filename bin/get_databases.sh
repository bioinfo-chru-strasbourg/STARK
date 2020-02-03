#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARKDatabases"
SCRIPT_DESCRIPTION="STARK download and build databases"
SCRIPT_RELEASE="0.9.3b"
SCRIPT_DATE="31/05/2019"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-11/12/2018: Script creation\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.1b-12/12/2018: Change to Makefile\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.2b-21/12/2018: Add update, build, rebuild and threads options. Change dbsnp source\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.3b-31/05/2019: Add APP configuration\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Configuration
ENV_CONFIG=$(find -L $SCRIPT_DIR/.. -name config.app)
source $ENV_CONFIG 1>/dev/null 2>/dev/null

# Header
function header () {
	echo "#######################################";
	echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
	echo "# $SCRIPT_DESCRIPTION ";
	echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © $SCRIPT_LICENCE";
	echo "#######################################";
}

# Release
function release () {
	echo -e $RELEASE_NOTES
}

# Usage
function usage {
	echo "# USAGE: $(basename $0) [options...]";
	echo "# --application=<STRING|FILE>              APP name or APP file configuration of the APPLICATION.";
	echo "#                                          Use 'default' for default application parameters ('APP/default.app').";
	echo "#                                          Default: Default STARK parameters.";
	echo "# --databases=<FOLDER>                     Databases folder (replace APP parameter)";
	echo "#                                          Will generate STARK databases folder structure";
	echo "# --additional_annotations=<STRING>        Additional HOWARD annotations";
	echo "#                                          Will download ANNOVAR and SNPEFF additional databases (beyond originally defined application)";
	echo "#                                          Format: 'annotation1,annotation2,...'";
	echo "# --build                                  Build all databases.";
	echo "# --rebuild                                Force Rebuild all databases.";
	echo "# --update                                 Update databases (latest dbSNP and HOWARD-ANNOVAR-snpEff databases) and build if needed.";
	echo "# --threads                                Number of threads (depend on system/proxy...).";

	echo "# --verbose                                VERBOSE option";
	echo "# --debug                                  DEBUG option";
	echo "# --release                                RELEASE option";
	echo "# --help                                   HELP option";
	echo "#";

}

# header
header;

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "e:bfut:vdnh" --long "env:,app:,application:,databases:,additional_annotations:,build,rebuild,update,threads:,verbose,debug,release,help" -- "$@" 2> /dev/null)
# || [ -z $@ ]

PARAM=$@
#PARAM=$(echo $@ | tr "\n" " ")
#echo $PARAM;
#PARAM=$(echo $ARGS | sed s/--//gi);
#exit 0;

eval set -- "$ARGS"
while true
do
	#echo "$1=$2"
	#echo "Eval opts";
	case "$1" in
		-e|--env|--app|--application)
			APP="$2"
			shift 2
			;;
		--databases)
			DATABASES="$2"
			shift 2
			;;
		--additional_annotations)
			ADDITIONAL_ANNOTATIONS="$2"
			shift 2
			;;
		-v|--verbose)
			VERBOSE=1
			shift 1
			;;
		-b|--build)
			BUILD=1
			shift 1
			;;
		-r|--rebuild)
			REBUILD=1
			shift 1
			;;
		-u|--update)
			UPDATE=1
			shift 1
			;;
		-t|--threads)
			THREADS_INPUT="$2"
			shift 2
			;;
		-d|--debug)
			VERBOSE=1
			DEBUG=1
			shift 1
			;;
		-n|--release)
			release;
			exit 0
			;;
		-h|--help)
			usage
			exit 0
			;;
		--) shift
			break
			;;
		*) 	echo "# Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done


# ACTION
##########

ACTION=0
[ $BUILD ] || [ $REBUILD ] || [ $UPDATE ] && ACTION=1;


# DATABASES FOLDER
####################

[ ! -z $DATABASES ] && [ ! -d $DATABASES ] && mkdir -p $DATABASES && echo "[INFO] Create databases folder '$DATABASES' "


# ENV
#########

#echo "APP=$APP"; exit;
(($VERBOSE)) && [ ! -z "$APP" ] && echo "#[INFO] Search Application '$APP'"

ENV=$(find_app "$APP" "$STARK_FOLDER_APPS")
source_app "$APP" "$STARK_FOLDER_APPS" 1
APP_NAME=$(name_app "$APP" "$STARK_FOLDER_APPS");

export ENV
export APP

(($VERBOSE)) && [ -z "$APP" ] && [ -z "$ENV" ] && echo "#[INFO] No Application provided. STARK default parameters will be used."
(($VERBOSE)) && [ ! -z "$APP" ] && [ ! -z "$ENV" ] && echo "#[INFO] Application '$APP' found ('$ENV')"
(($VERBOSE)) && [ ! -z "$APP" ] && [ -z "$ENV" ] && echo "#[INFO] Application '$APP' NOT found"

if ((0)); then
	echo "FOLDER_DATABASES=$FOLDER_DATABASES"
	echo "ANNOVAR_DATABASES=$ANNOVAR_DATABASES"
	echo "SNPEFF_DATABASES=$SNPEFF_DATABASES"
	echo "FOLDER_DATABASES_ANNOVAR=$FOLDER_DATABASES_ANNOVAR"
	echo "HOWARD_ANNOTATION=$HOWARD_ANNOTATION"
fi;

# CORES
re='^[0-9]+$'
#CORES=$(ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w)
CORES=$(nproc)
if ! [[ $THREADS =~ $re ]] || [ -z "$THREADS" ] || [ "$THREADS" == "" ] || [ $THREADS -gt $CORES ] ; then
	CORES_FREE=0
	THREADS=$(($CORES-$CORES_FREE))
fi;
# THREADS
if [[ $THREADS_INPUT =~ $re ]] && [ "$THREADS_INPUT" != "" ]; then
	THREADS=$THREADS_INPUT;
fi;

# TMP FOLDER_RUN
TMP_DATABASES_DOWNLOAD_FOLDER=$TMP_FOLDER_TMP/$RANDOM$RANDOM
mkdir -p $TMP_DATABASES_DOWNLOAD_FOLDER
TMP_DATABASES_DOWNLOAD_RAM="$(mktemp -d -p /dev/shm/)"
if [ "$TMP_DATABASES_DOWNLOAD_RAM" == "" ]; then
	TMP_DATABASES_DOWNLOAD_RAM=$TMP_DATABASES_DOWNLOAD_FOLDER;
fi;

# DATE
DATE=$(date '+%Y%m%d-%H%M%S')

# COPY MODE
COPY_MODE_DEFAULT="rsync -rtvu"

# variables
#############

#ASSEMBLY=hg38
#REF=/home1/TOOLS/db/genomes/hg38/hg38.fa
#VCFDBSNP=/home1/TOOLS/db/dbsnp/dbsnp.hg38.vcf.gz
#THREADS=3


# output
##########

echo "#"
echo "# ASSEMBY=$ASSEMBLY"
echo "# REF=$REF"
echo "# REFSEQ_GENES=$REFSEQ_GENES"
echo "# VCFDBSNP=$VCFDBSNP"
echo "# ANNOVAR_DATABASES=$ANNOVAR_DATABASES"
echo "# SNPEFF_DATABASES=$SNPEFF_DATABASES"
echo "# THREADS=$THREADS"
echo "#"

if (($REBUILD)); then
	echo "#"
	echo "# REBUILD Databases"
	echo "#"
fi;

# MK parameters
MK=$TMP_DATABASES_DOWNLOAD_FOLDER/mk
MK_LOG=$TMP_DATABASES_DOWNLOAD_FOLDER/mk.log
MK_ERR=$TMP_DATABASES_DOWNLOAD_FOLDER/mk.err
MK_ALL=""


# DBFOLDER
############

if (($REBUILD)); then
	mv -f $REBUILD.V$DATE;
fi;

# Create DBFOLDER
# MK
echo "$DBFOLDER:
	mkdir -p $DBFOLDER
	chmod 0775 $DBFOLDER
" >> $MK


# GENOME
##########

if [ ! -e $REF ]; then

	echo ""
	echo "### DOWNLOAD GENOME $ASSEMBLY"
	echo "#"

	# REF folder
	REF_FOLDER=$(dirname $REF)

	# Retrieve chromosomes
	UCSC_GENOME_URL="ftp://hgdownload.soe.ucsc.edu/goldenPath/$ASSEMBLY/chromosomes/"
	echo "# Check chromosomes on '$UCSC_GENOME_URL'..."
	CHRS=$(curl -s $UCSC_GENOME_URL | awk '{print $9}' | grep fa.gz* | grep -v "_" | sort -u -k1,1 -V)
	if [ "$CHRS" == "" ]; then
		CHRS=$(curl -s $UCSC_GENOME_URL | awk -F".fa.gz\">" '{print $2}' | awk -F"</A>" '{print $1}' | grep fa.gz* | grep -v "_" | sort -u -k1,1 -V)
	fi;
	CHROM_FILES=""

	# Download chromosomes files
	for CHROM_FILE in $CHRS; do

	    echo "# Download GENOME $ASSEMBLY $CHROM_FILE"
	    mkdir -p $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/chromosomes/

		# MK
	    echo "$TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/chromosomes/$CHROM_FILE:
			wget -S -c -O $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/chromosomes/$CHROM_FILE $UCSC_GENOME_URL/$CHROM_FILE
	    " >> $MK
	    CHROM_FILES=$CHROM_FILES" $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/chromosomes/$CHROM_FILE"

	done

	# Construct genome fasta file
	# MK
	echo "$REF: $DBFOLDER $CHROM_FILES
		mkdir -p $REF_FOLDER
		chmod 0775 $REF_FOLDER
		$GZ -dc $CHROM_FILES > $REF
	" >> $MK

	MK_ALL="$MK_ALL $REF"

fi;


# GENOME INDEX
################


## BOWTIE index

if [ ! -e $(dirname $REF)/$ASSEMBLY.rev.1.bt2 ]; then

	if [ "$BOWTIE" != "" ]; then

	    echo ""
	    echo "### BOWTIE GENOME $ASSEMBLY INDEX"
	    echo "#"

	    # MK
	    echo "$(dirname $REF)/$ASSEMBLY.rev.1.bt2: $REF
			$(dirname $BOWTIE)/bowtie2-build --threads $THREADS --packed $REF $(dirname $REF)/$ASSEMBLY ;
	    " >> $MK

		MK_ALL="$MK_ALL $(dirname $REF)/$ASSEMBLY.rev.1.bt2"

	fi;

fi;


## BWA index

if [ ! -e $REF.bwt ]; then

	if [ "$BWA" != "" ]; then

	    echo ""
	    echo "### BWA GENOME $ASSEMBLY INDEX"
	    echo "#"

		# MK
	    echo "$REF.bwt: $REF
			$BWA index -a bwtsw $REF;
	    " >> $MK

		MK_ALL="$MK_ALL $REF.bwt"

	fi;

fi;


## SAMTOOLS index

if [ ! -e $REF.fai ]; then

	if [ "$SAMTOOLS" != "" ]; then

	    echo ""
	    echo "### SAMTOOLS GENOME $ASSEMBLY INDEX"
	    echo "#"

		# MK
	    echo "$REF.fai: $REF
			$SAMTOOLS faidx $REF;
	    " >> $MK

		MK_ALL="$MK_ALL $REF.fai"

	fi;

fi;

## PICARD index

if [ ! -e $(dirname $REF)/$ASSEMBLY.dict ]; then

	if [ "$PICARD" != "" ]; then

	    echo ""
	    echo "### PICARD GENOME $ASSEMBLY INDEX"
	    echo "#"

		# MK
	    echo "$(dirname $REF)/$ASSEMBLY.dict: $REF
			$JAVA -jar $PICARD CreateSequenceDictionary \
				REFERENCE=$REF \
				OUTPUT=$(dirname $REF)/$ASSEMBLY.dict ;
	    " >> $MK

		MK_ALL="$MK_ALL $(dirname $REF)/$ASSEMBLY.dict"

	fi;

fi;


# REFSEQ_GENES
################

if [ ! -e $REFSEQ_GENES ]; then

	echo ""
	echo "### REFSEQ GENES TXT DOWNLOAD & BED TRANSLATION"
	echo "#"

	# RefGenes TXT file
	REFSEQ_GENES_URL="ftp://hgdownload.soe.ucsc.edu/goldenPath/$ASSEMBLY/database"
	REFSEQ_GENES_TXT=$(echo $REFSEQ_GENES | sed 's/bed$/txt/g')
	REFSEQ_GENES_GENES=$(echo $REFSEQ_GENES | sed 's/bed$/genes/g')

	# REFSEQ_GENE folder
	REFSEQ_GENES_FOLDER=$(dirname $REFSEQ_GENES)
	REFSEQ_GENES_TXT_FILE=$(basename $REFSEQ_GENES_TXT)
	REFSEQ_GENES_GENES_FILE=$(basename $REFSEQ_GENES_GENES)
	REFSEQ_GENES_FILE=$(basename $REFSEQ_GENES)

	# Update
	if (($UPDATE)); then
		if [ -e $REFSEQ_GENES_TXT ]; then mv -f $REFSEQ_GENES_TXT $REFSEQ_GENES_TXT.V$DATE; fi;
		if [ -e $REFSEQ_GENES_GENES ]; then mv -f $REFSEQ_GENES_GENES $REFSEQ_GENES_GENES.V$DATE; fi;
		if [ -e $REFSEQ_GENES ]; then mv -f $REFSEQ_GENES $REFSEQ_GENES.V$DATE; fi;
	fi;
	COPY_MODE_REFSEQ_GENES=$COPY_MODE_DEFAULT

	# OUTPUT
	echo "# REFSEQ_GENES_URL=$REFSEQ_GENES_URL"
	echo "# REFSEQ_GENES_TXT=$REFSEQ_GENES_TXT"
	echo "# REFSEQ_GENES_GENES=$REFSEQ_GENES_GENES"

	# Get refGene txt file
	echo "$TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_TXT_FILE: $DBFOLDER
		mkdir -p $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/
		wget -S -c -O $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/refGene.txt.gz $REFSEQ_GENES_URL/refGene.txt.gz ;
		mkdir -p $REFSEQ_GENES_FOLDER
		chmod 0775 $REFSEQ_GENES_FOLDER
		$GZ -dc  $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/refGene.txt.gz > $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_TXT_FILE
	" >> $MK
	# Copy in database folder
	echo "$REFSEQ_GENES_TXT: $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_TXT_FILE
		$COPY_MODE_REFSEQ_GENES $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_TXT_FILE $REFSEQ_GENES_TXT
	" >> $MK

	# Translate to refGene bed file
	echo "$TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_FILE: $DBFOLDER $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_TXT_FILE
		mkdir -p $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY
		cat $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_TXT_FILE | while IFS='' read -r line; do \
			CHROM=\$\$(echo \$\$line | awk '{print \$\$3}' | cut -d\"_\" -f1); \
			NM=\$\$(echo \$\$line | awk '{print \$\$2}'); \
			GENE=\$\$(echo \$\$line | awk '{print \$\$13}'); \
			STRAND=\$\$(echo \$\$line | awk '{print \$\$4}'); \
			echo \$\$line | awk '{print \$\$10}' | tr \",\" \"\\n\" | grep -v '^\$\$' > $TMP_DATABASES_DOWNLOAD_RAM/refGene.unsorted.bed.NM1 ; \
			echo \$\$line | awk '{print \$\$11}' | tr \",\" \"\\n\" | grep -v '^\$\$' > $TMP_DATABASES_DOWNLOAD_RAM/refGene.unsorted.bed.NM2; \
			paste $TMP_DATABASES_DOWNLOAD_RAM/refGene.unsorted.bed.NM1 $TMP_DATABASES_DOWNLOAD_RAM/refGene.unsorted.bed.NM2 | while read SS ; do \
				echo -e \"\$\$CHROM\t\$\$SS\t\$\$STRAND\t\$\$GENE\t\$\$NM\" >> $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/refGene.unsorted.bed; \
			done; \
		done
		cat $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/refGene.unsorted.bed | sort -k1,1V -k2,2n -k3,3n > $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_FILE
		rm -f $REFSEQ_GENES_TXT.NM1 $REFSEQ_GENES_TXT.NM2;
	" >> $MK
	# Copy in database folder
	echo "$REFSEQ_GENES: $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_FILE
		$COPY_MODE_REFSEQ_GENES $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_FILE $REFSEQ_GENES
	" >> $MK


	# Generates refGene genes file from bed file
	echo "$TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_GENES_FILE: $DBFOLDER $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_FILE
		mkdir -p $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY
		awk -F'\t' 'substr(\$\$6,1,2)=='NM' {print \$\$0}' $REFSEQ_GENES | grep "$(grep "@SQ" $(dirname $REF)/$ASSEMBLY.dict | cut -f2 | cut -d: -f2 | tr '\n' ' ' | sed 's/chr/ -e ^chr/gi')" > $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_GENES_FILE
	" >> $MK

	# Generates refGene genes file from bed file
	echo "$REFSEQ_GENES_GENES: $DBFOLDER $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_GENES_FILE
		$COPY_MODE_REFSEQ_GENES $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_GENES_FILE $REFSEQ_GENES_GENES
	" >> $MK

	MK_ALL="$MK_ALL $REFSEQ_GENES"

fi;

# DBSNP
#########

if [ ! -e $VCFDBSNP ] || (($UPDATE)); then

	echo ""
	echo "### DBSNP DOWNLOAD"
	echo "#"

	# DBSNP folder
	VCFDBSNP_FOLDER=$(dirname $VCFDBSNP)
	VCFDBSNP_FILE=$(basename $VCFDBSNP)

	# Update
	if (($UPDATE)); then
		if [ -e $VCFDBSNP ]; then mv -f $VCFDBSNP $VCFDBSNP.V$DATE; fi;
		if [ -e $VCFDBSNP.tbi ]; then mv -f $VCFDBSNP.tbi $VCFDBSNP.tbi.V$DATE; fi;
	fi;
	COPY_MODE_VCFDBSNP=$COPY_MODE_DEFAULT

	# NCBI Assembly translation
	if [ "$ASSEMBLY" == "hg38" ]; then
		ASSEMBLY_NCBI="GRCh38p7";
		SPECIES_NCBI="human_9606";
		NUM_NCBI="38";
	elif [ "$ASSEMBLY" == "hg19" ]; then
		ASSEMBLY_NCBI="GRCh37p13";
	  	SPECIES_NCBI="human_9606";
	  	NUM_NCBI="25";
	else
	  	ASSEMBLY_NCBI="GRCh37p13";
	  	SPECIES_NCBI="human_9606";
	  	NUM_NCBI="xx";
	fi;
	echo "# NCBI SPECIES=$SPECIES_NCBI"
	echo "# NCBI ASSEMBLY=$ASSEMBLY_NCBI"

	# DB release
	DBSNP_RELEASE="pre_build152" # "pre_build152" or "2.0"

	if [ "$DBSNP_RELEASE" == "2.0" ]; then

		# Find latest dbsnp file
		DBSNP_URL="ftp://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF"
		echo "# Check last DBSNP on '$DBSNP_URL'..."
		DBSNP_LATEST=$(curl -s $DBSNP_URL | awk '{print $9}' | grep "$NUM_NCBI.bgz$" | sort -k1,1 | tail -n 1)
		if [ "$DBSNP_LATEST" == "" ]; then
			DBSNP_LATEST=$(curl -s $DBSNP_URL | awk -F"\">" '{print $2}' | awk -F"</A>" '{print $1}' | grep "$NUM_NCBI.bgz$" | sort -k1,1 | tail -n 1)
		fi;

		DBSNP_LATEST_FILE="$DBSNP_URL/$DBSNP_LATEST"
		#echo "# DBSNP URL=$DBSNP_URL"
		#echo "# DBSNP LATEST=$DBSNP_LATEST"
		#echo "# DBSNP LATEST FILE=$DBSNP_LATEST_FILE"

	else

		# Find latest dbsnp file
		DBSNP_URL="ftp://ftp.ncbi.nlm.nih.gov/snp/pre_build152/organisms/"
		DBSNP_SPECIES_ASSEMBLY="$SPECIES_NCBI/$ASSEMBLY_NCBI"
		echo "# Check last DBSNP on '$DBSNP_URL'..."
		DBSNP_LATEST=$(curl -s $DBSNP_URL | awk '{print $9}' | grep $SPECIES_NCBI.*$ASSEMBLY_NCBI | sort -k1,1 | tail -n 1)
		if [ "$DBSNP_LATEST" == "" ]; then
			DBSNP_LATEST=$(curl -s $DBSNP_URL | awk -F"\">" '{print $2}' | awk -F"</A>" '{print $1}' | grep $SPECIES_NCBI.*$ASSEMBLY_NCBI | sort -k1,1 | tail -n 1)
		fi;
		#curl -s $DBSNP_URL -l  #| awk '{print $9}'
		DBSNP_LATEST_FILE=$DBSNP_URL$DBSNP_LATEST"/VCF/00-common_all.vcf.gz" # 00-All.vcf.gz


	fi;

	DBSNP_LATEST_FILE_BASENAME=$(basename $DBSNP_LATEST_FILE)

	echo "# DBSNP URL=$DBSNP_URL"
	echo "# DBSNP RELEASE=$DBSNP_RELEASE"
	echo "# DBSNP SPECIES/ASSEMBLY=$DBSNP_SPECIES_ASSEMBLY"
	echo "# DBSNP NUM=$NUM_NCBI"
	echo "# DBSNP LATEST=$DBSNP_LATEST"
	echo "# DBSNP LATEST FILE=$DBSNP_LATEST_FILE"
	echo "# DBSNP LATEST FILE BASENAME=$DBSNP_LATEST_FILE_BASENAME"

	# MK
	echo "$VCFDBSNP: $DBFOLDER
		mkdir -p $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER;
		chmod 0775 $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER;
		wget -S -c -O $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER/$DBSNP_LATEST.tbi $DBSNP_LATEST_FILE.tbi;
		wget -S -c -O $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER/$DBSNP_LATEST $DBSNP_LATEST_FILE;
		mkdir -p $VCFDBSNP_FOLDER/$ASSEMBLY
		chmod 0775 $VCFDBSNP_FOLDER/$ASSEMBLY
		$COPY_MODE_VCFDBSNP $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER/$DBSNP_LATEST $VCFDBSNP_FOLDER/$ASSEMBLY/$DBSNP_LATEST_FILE_BASENAME;
		$COPY_MODE_VCFDBSNP $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER/$DBSNP_LATEST.tbi $VCFDBSNP_FOLDER/$ASSEMBLY/$DBSNP_LATEST_FILE_BASENAME.tbi;
		$BCFTOOLS view $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER/$DBSNP_LATEST | $SCRIPT_DIR/contig_NC_to_chr.sh > $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER/$VCFDBSNP_FILE.vcf
		$BGZIP -c $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER/$VCFDBSNP_FILE.vcf > $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER/$VCFDBSNP_FILE
		$TABIX $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER/$VCFDBSNP_FILE
		mkdir -p $VCFDBSNP_FOLDER
		chmod 0775 $VCFDBSNP_FOLDER
		$COPY_MODE_VCFDBSNP $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER/$VCFDBSNP_FILE.tbi $VCFDBSNP.tbi;
		$COPY_MODE_VCFDBSNP $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER/$VCFDBSNP_FILE $VCFDBSNP;
		rm -rf $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER;
	" >> $MK

	MK_ALL="$MK_ALL $VCFDBSNP"

fi;


# ANNOVAR and SNPEFF DATABASES
################################

if [ ! -e $DBFOLDER/HOWARD.download.complete ] || (($UPDATE)); then

	echo ""
	echo "### ANNOVAR and SNPEFF DATABASES"
	echo "#"

	# COPY mode
	COPY_MODE_HOWARD=$COPY_MODE_DEFAULT #"cp -rf"

	# Rebuild or Update
	if (($UPDATE)); then
		if [ -e $DBFOLDER/HOWARD.download.complete ]; then mv -f $DBFOLDER/HOWARD.download.complete $DBFOLDER/HOWARD.download.complete.V$DATE; fi;
	fi;

	# DB to DOWNLOAD
	#HOWARD_DB="ALL,snpeff" # ALL, CORE, snpeff $HOWARD_ANNOTATION
	HOWARD_DB="$HOWARD_ANNOTATION,$HOWARD_ANNOTATION_REPORT,$HOWARD_ANNOTATION_MINIMAL,$ADDITIONAL_ANNOTATIONS" # ALL, CORE, snpeff $HOWARD_ANNOTATION

	if [ ! -e "$HOWARD_CONFIG_ANNOTATION" ]; then
		HOWARD_CONFIG_ANNOTATION=$HOWARDDIR/config.annotation.ini
	fi;

	# PROXY
	#Automaticly import system proxy settings
	if [ -n "$http_proxy" ] ; then
	    #secho $http_proxy | grep "@"
	    if [ $(echo $http_proxy | grep "@" -c) -eq 0 ]; then # If variable has username and password, its parse method different
	    PROXY_HOST=$(echo $http_proxy | sed 's/http:\/\/.*@\(.*\):.*/\1/')
	    PROXY_PORT=$(echo $http_proxy | sed 's/http:\/\/.*@.*:\(.*\)/\1/' | tr -d "/")
	    USERNAME=$(echo $http_proxy | sed 's/http:\/\/\(.*\)@.*/\1/'|awk -F: '{print $1}')
	    PASSWORD=$(echo $http_proxy | sed 's/http:\/\/\(.*\)@.*/\1/'|awk -F: '{print $2}')
	    else # If it doesn't have username and password, its parse method this
	    PROXY_HOST=$(echo $http_proxy | sed 's/http:\/\/\(.*\):.*/\1/')
	    PROXY_PORT=$(echo $http_proxy | sed 's/http:\/\/.*:\(.*\)/\1/' | tr -d "/")
	    fi
	fi

	# If you want to use this feature , enables and disables proxy support for proxy service up or down status
	# uncomment these line, if you installed nmap
	# at ubuntu system you can type this command for this future
	# sudo apt-get install nmap
	#STATUS=$(nmap -sT $PROXY_HOST -p $PROXY_PORT 2>/dev/null| grep open |awk '{print $2}')
	#if [ "$STATUS" != "open" ]; then # If service isn't running, disable proxy support
	#    PROXY_HOST=""
	#    PROXY_PORT=""
	#fi

	JAVA_FLAGS=" -Djava.io.tmpdir=$TMP_DATABASES_DOWNLOAD_FOLDER/HOWARD "
	if [ -n "$PROXY_HOST"   -a  -n "$PROXY_PORT" ] ; then
	    JAVA_FLAGS=" -Dhttp.proxyHost=$PROXY_HOST -Dhttp.proxyPort=$PROXY_PORT"
	    if [ -n "$USERNAME" -a -n "$PASSWORD" ]; then
	    JAVA_FLAGS="$JAVA_FLAGS -Dhttp.proxyUser=$USERNAME -Dhttp.proxyPassword=$PASSWORD"
	    fi
	fi
	#echo $JAVA_FLAGS

	if [ "$SNPEFF" != "" ]; then
		SNPEFF_CMD="$JAVA $JAVA_FLAGS -jar $SNPEFF download $ASSEMBLY -dataDir $TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES 2>$TMP_DATABASES_DOWNLOAD_FOLDER/snpeff.err";
		SNPEFF_CMD="$SNPEFF_CMD; if ((\$\$(grep 'ERROR while connecting to' $TMP_DATABASES_DOWNLOAD_FOLDER/snpeff.err -c))); then wget -O $TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES/snpEff.$SNPEFF_VERSION.$ASSEMBLY.zip \$\$(grep 'ERROR while connecting to' err | cut -f2 | cut -d' ' -f5); unzip $TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES/snpEff.$SNPEFF_VERSION.$ASSEMBLY.zip -d $TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES/; mv $TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES/data/* $TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES/; fi"
		echo "# CMD= $SNPEFF_CMD"
		echo "#"
	fi;

	if [ "$HOWARD" != "" ]; then
		if [ -e $HOWARD_FOLDER_DOCS/example.vcf ]; then
			INPUT_VCF=$HOWARD_FOLDER_DOCS/example.vcf;
		elif [ -e $HOWARD_FOLDER/docs/example.vcf ]; then
			INPUT_VCF=$HOWARD_FOLDER/docs/example.vcf
		elif [ -e $(dirname $HOWARD_FOLDER_BIN)/docs/example.vcf ]; then
			INPUT_VCF=$(dirname $HOWARD_FOLDER_BIN)/docs/example.vcf
		elif [ -e $(dirname $HOWARDDIR)/docs/example.vcf ]; then
			INPUT_VCF=$(dirname $HOWARDDIR)/docs/example.vcf
		elif [ -e $HOWARDDIR/docs/example.vcf ]; then
			INPUT_VCF=$HOWARDDIR/docs/example.vcf
		fi;
		#INPUT_VCF=$HOWARDDIR/docs/example.vcf
		HOWARD_CMD="$HOWARD --input=$INPUT_VCF --output=$DBFOLDER/HOWARD.download.vcf --env=$CONFIG_TOOLS --annotation=$HOWARD_DB --annovar_folder=$ANNOVAR --annovar_databases=$TMP_DATABASES_DOWNLOAD_FOLDER/$ANNOVAR_DATABASES --config_annotation=$HOWARD_CONFIG_ANNOTATION --snpeff_jar=$SNPEFF --snpeff_databases=$TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES --snpeff_threads=$THREADS --tmp=$TMP_DATABASES_DOWNLOAD_FOLDER/HOWARD --assembly=$ASSEMBLY --java=$JAVA --java_flags='\"$JAVA_FLAGS\"' --verbose --threads=1 >$DBFOLDER/HOWARD.download.log;"
		#HOWARD_CMD="$HOWARD --input=$INPUT_VCF --output=$DBFOLDER/HOWARD.download.vcf  --annotation=location --annovar_folder=$ANNOVAR --annovar_databases=$TMP_DATABASES_DOWNLOAD_FOLDER/$ANNOVAR_DATABASES --config_annotation=$HOWARD_CONFIG_ANNOTATION --snpeff_jar=$SNPEFF --snpeff_databases=$TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES --snpeff_threads=$THREADS --tmp=$TMP_DATABASES_DOWNLOAD_FOLDER/HOWARD --assembly=$ASSEMBLY --java=$JAVA --java_flags='\"$JAVA_FLAGS\"' --verbose --threads=1 >$DBFOLDER/HOWARD.download.log;"
		echo "# CMD= $HOWARD_CMD"
		echo "#"
	fi;

	# MK
	echo "$DBFOLDER/HOWARD.download.complete: $DBFOLDER
		mkdir -p $TMP_DATABASES_DOWNLOAD_FOLDER/$ANNOVAR_DATABASES;
		chmod 0775 $TMP_DATABASES_DOWNLOAD_FOLDER/$ANNOVAR_DATABASES
		mkdir -p $TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES;
		chmod 0775 $TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES
		$SNPEFF_CMD
		$HOWARD_CMD
		$COPY_MODE_HOWARD $TMP_DATABASES_DOWNLOAD_FOLDER/$ANNOVAR_DATABASES/ $ANNOVAR_DATABASES/;
		$COPY_MODE_HOWARD $TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES/ $SNPEFF_DATABASES/;
		if ! ((\$\$(grep -c '\[ERROR\]' $DBFOLDER/HOWARD.download.log))); then \
			echo '#[INFO] HOWARD download complete' > $DBFOLDER/HOWARD.download.complete ; \
		fi;
		rm -rf $TMP_DATABASES_DOWNLOAD_FOLDER/$ANNOVAR_DATABASES;
		rm -rf $TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES;
	" >> $MK

	# $HOWARD --input=$HOWARDDIR/docs/example.vcf --output=$DBFOLDER/HOWARD.download.vcf  --annotation=$HOWARD_DB --annovar_folder=$ANNOVAR --annovar_databases=$TMP_DATABASES_DOWNLOAD_FOLDER/$ANNOVAR_DATABASES --config_annotation=$HOWARD_CONFIG_ANNOTATION --snpeff_jar=$SNPEFF --snpeff_databases=$TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES --snpeff_threads=$THREADS --tmp=$TMP_DATABASES_DOWNLOAD_FOLDER/HOWARD --assembly=$ASSEMBLY --java=$JAVA --java_flags=\"'$JAVA_FLAGS\"' --verbose --threads=1;
	# $HOWARD --input=$HOWARDDIR/docs/example.vcf --output=$DBFOLDER/HOWARD.download.vcf  --annotation=$HOWARD_DB --annovar_folder=$ANNOVAR --annovar_databases=$TMP_DATABASES_DOWNLOAD_FOLDER/$ANNOVAR_DATABASES --config_annotation=$HOWARD_CONFIG_ANNOTATION --snpeff_jar=$SNPEFF --snpeff_databases=$TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES --snpeff_threads=$THREADS --tmp=$TMP_DATABASES_DOWNLOAD_FOLDER/HOWARD --assembly=$ASSEMBLY --java=$JAVA --java_flags='$JAVA_FLAGS' --verbose --threads=1;

	MK_ALL="$MK_ALL $DBFOLDER/HOWARD.download.complete"

fi;

# ALL
######

echo "all: $MK_ALL
" >> $MK

echo "### MAKEFILE"
echo "#"
echo "# MK=$MK"
echo "# LOG=$MK_LOG"
echo ""

if (($DEBUG)); then
	cat -n $MK;
fi;


# DATABASES INIT
echo "## DATABASES INIT..."
if (($BUILD)) || (($REBUILD)) || (($UPDATE)); then
	if (($VERBOSE)) || (($DEBUG)); then
		if (($DEBUG)); then
			make -k -j $THREADS $MK_OPTION -s -f $MK all;
		elif (($VERBOSE)); then
			make -k -j $THREADS $MK_OPTION -s -f $MK all;
		fi;
	else
		make -k -j $THREADS $MK_OPTION -f $MK all 1>$MK_LOG 2>$MK_ERR;
	fi;
else
	echo "## use --build, --rebuild or --update to init databases"
fi;

# CLEAN
rm -Rf $TMP_DATABASES_DOWNLOAD_FOLDER
rm -Rf $TMP_DATABASES_DOWNLOAD_RAM

exit 0;
