#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARKDatabases"
SCRIPT_DESCRIPTION="STARK download and build databases"
SCRIPT_RELEASE="0.9.1b"
SCRIPT_DATE="12/12/2018"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-11/12/2018: Script creation\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.1b-12/12/2018: Change to Makefile\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

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
	echo "# -e/--env/--app=<FILE>                       ENV file configuration of the APPLICATION, and dabases to use.";
	echo "#                                             Must be in the STARK folder if relative path";
	echo "#                                             Default: defined in the RUN SampleSheet, or env.sh if not defined";
	echo "# -b/--build                                  Build all databases.";
	echo "# -f/--rebuild                                Force Rebuild all databases.";
	echo "# -u/--update                                 Update databases (latest dbSNP and HOWARD-ANNOVAR-snpEff databases) and build if needed.";
	echo "# -t/--threads                                Number of threads (depend on system/proxy...).";

	echo "# -v/--verbose                                VERBOSE option";
	echo "# -d/--debug                                  DEBUG option";
	echo "# -n/--release                                RELEASE option";
	echo "# -h/--help                                   HELP option";
	echo "#";

}

# header
header;

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "e:bfut:vdnh" --long "env:,build,rebuild,update,threads:,verbose,debug,release,help" -- "$@" 2> /dev/null)
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
		-e|--env|--app)
			ENV="$2"
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


# ENV
#########

if [ -s $ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$ENV;
elif [ -s $SCRIPT_DIR/$ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$SCRIPT_DIR/$ENV;
elif [ "$ENV" == "" ] || [ ! -s $ENV ]; then
	if [ -s $SCRIPT_DIR/"env.sh" ]; then
		ENV=$SCRIPT_DIR/"env.sh";
	else
		ENV="";
		echo "#[WARNING] NO ENV defined. Default ENV used."
	fi;
fi;
# SOURCE ENV if exists
if [ ! -z $ENV ] && [ -s $ENV ]; then
	source $ENV;
fi;


re='^[0-9]+$'
CORES=$(ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w)
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
	#mkdir -p $(dirname $REF)
	#chmod 0775 $(dirname $REF)

	# Retrieve chromosomes
	UCSC_GENOME_URL="ftp://hgdownload.soe.ucsc.edu/goldenPath/$ASSEMBLY/chromosomes/"
	echo "# Check chromosomes on '$UCSC_GENOME_URL'..."
	#curl $UCSC_GENOME_URL
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
		-$GZIP -dc $CHROM_FILES > $REF
	" >> $MK

	MK_ALL="$MK_ALL $REF"

fi;


# GENOME INDEX
################

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

# REFSEQ_GENES
################

if [ ! -e $REFSEQ_GENES ]; then

	echo ""
	echo "### REFSEQ GENES TXT DOWNLOAD & BED TRANSLATION"
	echo "#"

	# RefGenes TXT file
	REFSEQ_GENES_URL="ftp://hgdownload.soe.ucsc.edu/goldenPath/$ASSEMBLY/database"
	REFSEQ_GENES_TXT=$(echo $REFSEQ_GENES | sed 's/bed$/txt/g')

	# REFSEQ_GENE folder
	REFSEQ_GENES_FOLDER=$(dirname $REFSEQ_GENES)
	REFSEQ_GENES_TXT_FILE=$(basename $REFSEQ_GENES_TXT)
	REFSEQ_GENES_FILE=$(basename $REFSEQ_GENES)
	#mkdir -p $(dirname $REFSEQ_GENES)
	#chmod 0775 $(dirname $REFSEQ_GENES)
	# Update
	if (($UPDATE)); then
		if [ -e $REFSEQ_GENES_TXT ]; then mv -f $REFSEQ_GENES_TXT $REFSEQ_GENES_TXT.V$DATE; fi;
		if [ -e $REFSEQ_GENES ]; then mv -f $REFSEQ_GENES $REFSEQ_GENES.V$DATE; fi;
	fi;
	COPY_MODE_REFSEQ_GENES=$COPY_MODE_DEFAULT

	# OUTPUT
	echo "# REFSEQ_GENES_URL=$REFSEQ_GENES_URL"
	echo "# REFSEQ_GENES_TXT=$REFSEQ_GENES_TXT"

	# Get refGene txt file
	echo "$REFSEQ_GENES_TXT: $DBFOLDER
		mkdir -p $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/
		wget -S -c -O $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/refGene.txt.gz $REFSEQ_GENES_URL/refGene.txt.gz ;
		mkdir -p $REFSEQ_GENES_FOLDER
		chmod 0775 $REFSEQ_GENES_FOLDER
		-$GZIP -dc  $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/refGene.txt.gz > $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_TXT_FILE
		$COPY_MODE_REFSEQ_GENES $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_TXT_FILE $REFSEQ_GENES_TXT
	" >> $MK

	# Translate to refGene bed file
	echo "$REFSEQ_GENES: $DBFOLDER $REFSEQ_GENES_TXT
		mkdir -p $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY
		cat $REFSEQ_GENES_TXT | while IFS='' read -r line; do \
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
		$COPY_MODE_REFSEQ_GENES $TMP_DATABASES_DOWNLOAD_FOLDER/$ASSEMBLY/$REFSEQ_GENES_FILE $REFSEQ_GENES
		rm -f $REFSEQ_GENES_TXT.NM1 $REFSEQ_GENES_TXT.NM2;
	" >> $MK

	MK_ALL="$MK_ALL $REFSEQ_GENES"

fi;

# DBSNP
#########

if [ ! -e $VCFDBSNP ] || (($UPDATE)); then

	echo ""
	echo "### DBSNP DOWNLOAD"
	echo "#"

	# Create folder
	VCFDBSNP_FOLDER=$(dirname $VCFDBSNP)
	VCFDBSNP_FILE=$(basename $VCFDBSNP)
	#mkdir -p $(dirname $VCFDBSNP)
	#chmod 0775 $(dirname $VCFDBSNP)

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


	# Find latest dbsnp file
	DBSNP_URL="ftp://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF" 
	#DBSNP_SPECIES_ASSEMBLY="$SPECIES_NCBI/$ASSEMBLY_NCBI"
	echo "# Check last DBSNP on '$DBSNP_URL'..."
	DBSNP_LATEST=$(curl -s $DBSNP_URL | awk '{print $9}' | grep "$NUM_NCBI.bgz$" | sort -k1,1 | tail -n 1)
	if [ "$DBSNP_LATEST" == "" ]; then
		DBSNP_LATEST=$(curl -s $DBSNP_URL | awk -F"\">" '{print $2}' | awk -F"</A>" '{print $1}' | grep "$NUM_NCBI.bgz$" | sort -k1,1 | tail -n 1)
	fi;
	#curl -s $DBSNP_URL -l  #| awk '{print $9}'
	DBSNP_LATEST_FILE="$DBSNP_URL/$DBSNP_LATEST" # 00-All.vcf.gz
	echo "# DBSNP URL=$DBSNP_URL"
	#echo "# DBSNP SPECIES/ASSEMBLY=$DBSNP_SPECIES_ASSEMBLY"
	echo "# DBSNP LATEST=$DBSNP_LATEST"
	echo "# DBSNP LATEST FILE=$DBSNP_LATEST_FILE"

	if ((0)); then
		# Find latest dbsnp file
		DBSNP_URL="ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/"
		DBSNP_SPECIES_ASSEMBLY="$SPECIES_NCBI/$ASSEMBLY_NCBI"
		echo "# Check last DBSNP on '$DBSNP_URL'..."
		DBSNP_LATEST=$(curl -s $DBSNP_URL | awk '{print $9}' | grep $SPECIES_NCBI.*$ASSEMBLY_NCBI | sort -k1,1 | tail -n 1)
		if [ "$DBSNP_LATEST" == "" ]; then
			DBSNP_LATEST=$(curl -s $DBSNP_URL | awk -F"\">" '{print $2}' | awk -F"</A>" '{print $1}' | grep $SPECIES_NCBI.*$ASSEMBLY_NCBI | sort -k1,1 | tail -n 1)
		fi;
		#curl -s $DBSNP_URL -l  #| awk '{print $9}'
		DBSNP_LATEST_FILE=$DBSNP_URL$DBSNP_LATEST"/VCF/00-common_all.vcf.gz" # 00-All.vcf.gz
		echo "# DBSNP URL=$DBSNP_URL"
		echo "# DBSNP SPECIES/ASSEMBLY=$DBSNP_SPECIES_ASSEMBLY"
		echo "# DBSNP LATEST=$DBSNP_LATEST/VCF/00-common_all.vcf.gz"
		echo "# DBSNP LATEST FILE=$DBSNP_LATEST_FILE"
	fi;

	# MK
	echo "$VCFDBSNP: $DBFOLDER
		mkdir -p $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER;
		chmod 0775 $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER;
		wget -S -c -O $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER/$VCFDBSNP_FILE.tbi $DBSNP_LATEST_FILE.tbi;
		wget -S -c -O $TMP_DATABASES_DOWNLOAD_FOLDER/$VCFDBSNP_FOLDER/$VCFDBSNP_FILE $DBSNP_LATEST_FILE;
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

	# Create folder
	#mkdir -p $DBFOLDER;
	#chmod 0775 $DBFOLDER

	# COPY mode
	COPY_MODE_HOWARD=$COPY_MODE_DEFAULT #"cp -rf"


	# Rebuild or Update
	if (($UPDATE)); then
		if [ -e $DBFOLDER/HOWARD.download.complete ]; then mv -f $DBFOLDER/HOWARD.download.complete $DBFOLDER/HOWARD.download.complete.V$DATE; fi;
		#COPY_MODE="rsync -rtvuc" # checksum:-rtvuc date/size:-rtvu
	fi;

	# DB to DOWNLOAD
	#HOWARD_DB="CORE,snpeff" # ALL, CORE, snpeff
	#HOWARD_DB="ALL,snpeff" # ALL, CORE, snpeff
	HOWARD_DB="ALL,snpeff" # ALL, CORE, snpeff
	#HOWARD_DB="snpeff" # ALL, CORE, snpeff
	#HOWARD_DB="popfreq" # ALL, CORE, snpeff

	if [ ! -e "$HOWARD_CONFIG_ANNOTATION" ]; then
		HOWARD_CONFIG_ANNOTATION=$HOWARDDIR/config.annotation.ini
	fi;

	# MK
	echo "$DBFOLDER/HOWARD.download.complete: $DBFOLDER
		mkdir -p $TMP_DATABASES_DOWNLOAD_FOLDER/$ANNOVAR_DATABASES;
		chmod 0775 $TMP_DATABASES_DOWNLOAD_FOLDER/$ANNOVAR_DATABASES
		mkdir -p $TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES;
		chmod 0775 $TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES
		$HOWARD --input=$HOWARDDIR/docs/example.vcf --output=$DBFOLDER/HOWARD.download.vcf  --annotation=$HOWARD_DB --annovar_folder=$ANNOVAR --annovar_databases=$TMP_DATABASES_DOWNLOAD_FOLDER/$ANNOVAR_DATABASES --config_annotation=$HOWARD_CONFIG_ANNOTATION --snpeff_jar=$SNPEFF --snpeff_databases=$TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES --snpeff_threads=$THREADS --tmp=$TMP_DATABASES_DOWNLOAD_FOLDER/HOWARD --assembly=$ASSEMBLY --java=$JAVA --verbose --threads=1;
		$COPY_MODE_HOWARD $TMP_DATABASES_DOWNLOAD_FOLDER/$ANNOVAR_DATABASES/ $ANNOVAR_DATABASES/;
		$COPY_MODE_HOWARD $TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES/ $SNPEFF_DATABASES/;
		if ! ((\$\$(grep -c '\*\*\*' $DBFOLDER/HOWARD.download.log))); then \
			echo '#[INFO] HOWARD download complete' > $DBFOLDER/HOWARD.download.complete ; \
		fi;
		rm -rf $TMP_DATABASES_DOWNLOAD_FOLDER/$ANNOVAR_DATABASES;
		rm -rf $TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES;
	" >> $MK
# if (($UPDATE)); then cp -rf $ANNOVAR_DATABASES/* $TMP_DATABASES_DOWNLOAD_FOLDER/$ANNOVAR_DATABASES/; cp -rf $SNPEFF_DATABASES/* $TMP_DATABASES_DOWNLOAD_FOLDER/$SNPEFF_DATABASES/; fi;
# 1>>$DBFOLDER/HOWARD.download.log 2>>$DBFOLDER/HOWARD.download.log

	MK_ALL="$MK_ALL $DBFOLDER/HOWARD.download.complete"

fi;

# ALL
######

#   $REF $REF.bwt $REF.fai $(dirname $REF)/$ASSEMBLY.dict $(dirname $REF)/$ASSEMBLY.1.bt2 $REFSEQ_GENES_TXT $REFSEQ_GENES $VCFDBSNP $DBFOLDER/HOWARD.download.log
# echo all: $REF $REF.bwt $REF.fai $(dirname $REF)/$ASSEMBLY.dict $(dirname $REF)/$ASSEMBLY.1.bt2 $REFSEQ_GENES_TXT $REFSEQ_GENES $VCFDBSNP $DBFOLDER/HOWARD.download.log

echo "all: $MK_ALL
" >> $MK

echo "## MAKEFILE"
echo "#"
echo "# MK=$MK"
echo "# LOG=$MK_LOG"
echo ""

if (($DEBUG)); then
	cat -n $MK;
fi;


# MK OPTIONS
#MK_OPTION=" -o "$(echo $MK_ALL | sed 's/ / -o /g')
#echo $MK_OPTION


# Change MK to be very old
#default_timestamp='1980-01-01 00:00:00'
#touch -c -d "$default_timestamp" $MK

# DATABASES INIT
echo "## DATABASES INIT..."
if (($BUILD)) || (($REBUILD)) || (($UPDATE)); then
	if (($VERBOSE)) || (($DEBUG)); then
		if (($DEBUG)); then
			make -k -j $THREADS $MK_OPTION -d -f $MK all;
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
