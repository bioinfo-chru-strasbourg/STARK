#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARKDatabases"
SCRIPT_DESCRIPTION="STARK download and build databases"
SCRIPT_RELEASE="0.9.4b"
SCRIPT_DATE="17/06/2020"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-11/12/2018: Script creation\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.1b-12/12/2018: Change to Makefile\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.2b-21/12/2018: Add update, build, rebuild and threads options. Change dbsnp source\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.3b-31/05/2019: Add APP configuration\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.4b-17/06/2020: Clarify code, organisation DB/RELEASE, add STARK.description\n";

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

[ ! -z $DATABASES ] && [ ! -d $DATABASES ] && mkdir -p $DATABASES && echo "#[INFO] Create databases folder '$DATABASES' "


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
#COPY_MODE_DEFAULT="rsync -rtvu"
COPY_MODE_DEFAULT="rsync -ar"


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

# JAVA FLAGS
mkdir -p $TMP_DATABASES_DOWNLOAD_FOLDER/JAVA_FLAGS
JAVA_FLAGS=" -Djava.io.tmpdir=$TMP_DATABASES_DOWNLOAD_FOLDER/JAVA_FLAGS "
if [ -n "$PROXY_HOST"   -a  -n "$PROXY_PORT" ] ; then
    JAVA_FLAGS=" -Dhttp.proxyHost=$PROXY_HOST -Dhttp.proxyPort=$PROXY_PORT"
    if [ -n "$USERNAME" -a -n "$PASSWORD" ]; then
    JAVA_FLAGS="$JAVA_FLAGS -Dhttp.proxyUser=$USERNAME -Dhttp.proxyPassword=$PASSWORD"
    fi
fi


# DATABASES RELEASE
#####################

DATABASES_RELEASE=$DATE



# variables
#############

#ASSEMBLY=hg38
#REF=/home1/TOOLS/db/genomes/hg38/hg38.fa
#VCFDBSNP=/home1/TOOLS/db/dbsnp/dbsnp.hg38.vcf.gz
#THREADS=3


# output
##########

echo ""
echo "#[INFO] DATABASES_RELEASE=$DATABASES_RELEASE"
echo "#[INFO] ASSEMBY=$ASSEMBLY"
echo "#[INFO] REF=$REF"
echo "#[INFO] REFSEQ_GENES=$REFSEQ_GENES"
echo "#[INFO] VCFDBSNP=$VCFDBSNP"
echo "#[INFO] ANNOVAR_DATABASES=$ANNOVAR_DATABASES"
echo "#[INFO] SNPEFF_DATABASES=$SNPEFF_DATABASES"
echo "#[INFO] THREADS=$THREADS"
#echo ""

if (($REBUILD)); then
	echo ""
	echo "# REBUILD Databases"
	#echo "#"
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

# INFOS
#########

DOWNLOAD_METHOD="STARK Databases download script [$SCRIPT_RELEASE-$SCRIPT_DATE]"


##########
# GENOME #
##########

if ((1)); then

	# DB
	DATABASE="genomes"
	DATABASE_NAME="Genomes"
	DATABASE_FULLNAME="Reference Genome Sequences Assembly"
	DATABASE_WEBSITE="https://genome.ucsc.edu/"
	DATABASE_DESCRIPTION="Reference sequence was produced by the Genome Reference Consortium, and is composed of genomic sequence, primarily finished clones that were sequenced as part of the Human Genome Project"

	# DB TARGET
	DB_TARGET=$REF;												# /STARK/databases/genomes/current/hg19.fa
	DB_ASSEMBLY=$ASSEMBLY; 										# hg19
	DB_TARGET_FILE=$(basename $DB_TARGET);						# hg19.fa
	DB_TARGET_FOLDER=$(dirname $DB_TARGET);						# /STARK/databases/genomes/current
	#DB_TARGET_RELEASE=$(basename $DB_TARGET_FOLDER);			# current
	DB_TARGET_RELEASE=$DB_ASSEMBLY;								# hg19
	DB_TARGET_DB_FOLDER=$(dirname $DB_TARGET_FOLDER);			# /STARK/databases/genomes

	# DB RELEASE
	DB_RELEASE=$DATABASES_RELEASE;								# DATE
	#[ -s $DB_TARGET_FOLDER ] && DB_RELEASE=$DB_TARGET_RELEASE	# if /STARK/databases/genomes/hg19 exists then keep release
	DB_RELEASE_DATE=$DATABASES_RELEASE;							# DATE
	DB_RELEASE_FILE=$DB_TARGET_FILE;							# hg19.fa
	DB_RELEASE_FOLDER=$DB_TARGET_DB_FOLDER/$DB_RELEASE;			# /STARK/databases/genomes/DATE
	DB_RELEASE_FILE_PATH="$DB_RELEASE_FOLDER/$DB_RELEASE_FILE";	# /STARK/databases/genomes/DATE/hg19.fa

	# TMP files and folders
	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DB_RELEASE
	mkdir -p $DB_TMP


	## Reference genome
	###################

	if [ ! -e $REF ]; then

		# VERBOSE
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DB_RELEASE' for '$DB_TARGET_RELEASE' [$DB_ASSEMBLY]"
		
		# DEBUG
		#(($DEBUG)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DB_RELEASE' for '$DB_TARGET_RELEASE' ";

		# Retrieve chromosomes
		UCSC_GENOME_URL="ftp://hgdownload.soe.ucsc.edu/goldenPath/$DB_ASSEMBLY/chromosomes/"
		UCSC_GENOME_URL_FILE="md5sum.txt";

		(($VERBOSE)) && echo "#[INFO] Check chromosomes files on '$UCSC_GENOME_URL'..."
		CHRS=""
		CHRS=$(curl $UCSC_GENOME_URL/chr*.fa.gz --list-only -s | grep fa.gz$ | grep -v "_" | sort -u -k1,1 -V)
		if [ "$CHRS" == "" ]; then
			CHRS=$(curl -s $UCSC_GENOME_URL | awk -F".fa.gz\">" '{print $2}' | awk -F"</A>" '{print $1}' | grep fa.gz* | grep -v "_" | sort -u -k1,1 -V)
		fi;
		CHRS_LIST=$(echo $CHRS | tr "\n" " ")
		(($VERBOSE)) && echo "#[INFO] Chromosomes files: '$CHRS_LIST'"

		UCSC_GENOME_URL_FILE_DATE=$(curl -s -I $UCSC_GENOME_URL/$UCSC_GENOME_URL_FILE | grep "Last-Modified: " | sed "s/Last-Modified: //g" | sed "s/\r$//g")

		# RELEASE
		DB_RELEASE_FROM_DOWNLOAD=$(date -d "$UCSC_GENOME_URL_FILE_DATE")

		# Check chromosomes files
		for CHROM_FILE in $CHRS; do
			CHROM_FILES=$CHROM_FILES" $DB_TMP/$CHROM_FILE"
			CHROM_FILES_CURL=$CHROM_FILES_CURL" && curl -s -R $UCSC_GENOME_URL/$CHROM_FILE -o $DB_TMP/$CHROM_FILE"
		done;

		# DATABASE Infos
		DB_INFOS_JSON='
		{
			"code": "'$DATABASE'",
			"name": "'$DATABASE_NAME'",
			"fullname": "'$DATABASE_FULLNAME'",
			"website": "'$DATABASE_WEBSITE'",
			"description": "'$DATABASE_DESCRIPTION'"
		}
		';
		echo "$DB_INFOS_JSON" > $DB_TMP/STARK.database

		# DATABASE Release Infos
		DB_RELEASE_INFOS_JSON='
		{
			"release": "'$DB_RELEASE_FROM_DOWNLOAD'",
			"date": "'$DB_RELEASE_DATE'",
			"files": [ "'$DB_TARGET_FILE'" ],
			"assembly": [ "'$DB_ASSEMBLY'" ],
			"download": {
				"methode": "'$DOWNLOAD_METHOD'",
				"URL": "'$UCSC_GENOME_URL'",
				"file": "'$UCSC_GENOME_URL_FILE'",
				"date": "'$UCSC_GENOME_URL_FILE_DATE'",
				"chromosomes": "'$CHRS_LIST'"
			}
		}
		';
		echo "$DB_RELEASE_INFOS_JSON" > $DB_TMP/STARK.database.release

		# Chromosomes tmp folder
		echo "$DB_TMP:
				mkdir -p $DB_TMP/
		    " >> $MK

		# Chromosomes
		echo "$DB_TMP/$DB_RELEASE_FILE: $DB_TMP
				curl $UCSC_GENOME_URL/README.txt -s -R -o $DB_TMP/README.txt
				curl $UCSC_GENOME_URL/md5sum.txt -s -R -o $DB_TMP/md5sum.txt
				((1)) $CHROM_FILES_CURL
				$GZ -dc $CHROM_FILES > $DB_TMP/$DB_RELEASE_FILE
		    " >> $MK

		# Final
		echo "$DB_TARGET: $DB_TMP/$DB_RELEASE_FILE
			# Folder creation
			mkdir -p $DB_RELEASE_FOLDER
			chmod 0775 $DB_RELEASE_FOLDER
			# Copy Files
			cp -p $DB_TMP/$DB_RELEASE_FILE $DB_RELEASE_FILE_PATH
			cp -p $DB_TMP/README.txt $DB_RELEASE_FOLDER/README.txt
			cp -p $DB_TMP/md5sum.txt $DB_RELEASE_FOLDER/md5sum.txt
			# date 
			touch $DB_RELEASE_FILE_PATH -r $DB_RELEASE_FOLDER/md5sum.txt
			# database release info
			[ ! -s $DB_TARGET_DB_FOLDER/STARK.database ] && cp $DB_TMP/STARK.database $DB_TARGET_DB_FOLDER/STARK.database
			cp $DB_TMP/STARK.database.release $DB_RELEASE_FOLDER
			# links
			[ $DB_TARGET_RELEASE != $DB_RELEASE ] && ln -snf $DB_RELEASE/ $DB_TARGET_FOLDER
			[ latest != $DB_RELEASE ] && ln -snf $DB_RELEASE/ $DB_TARGET_DB_FOLDER/latest
			[ current != $DB_RELEASE ] && ln -snf $DB_RELEASE/ $DB_TARGET_DB_FOLDER/current
			# Clear
			rm -rf $DB_TMP;
		" >> $MK

		# Add to MK
		MK_ALL="$MK_ALL $DB_TARGET"

	fi;


	# Genome index
	##############

	## BOWTIE index

	if [ ! -e $DB_TARGET_FOLDER/$ASSEMBLY.rev.1.bt2 ]; then

		if [ "$BOWTIE" != "" ]; then

		    (($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DB_RELEASE' for '$DB_TARGET_RELEASE' BOWTIE Index"

		    # MK
		    echo "$DB_TARGET_FOLDER/$ASSEMBLY.rev.1.bt2: $REF
				$(dirname $BOWTIE)/bowtie2-build --threads $THREADS --packed $REF $DB_TARGET_FOLDER/$ASSEMBLY ;
		    " >> $MK

			MK_ALL="$MK_ALL $DB_TARGET_FOLDER/$ASSEMBLY.rev.1.bt2"

		fi;

	fi;


	## BWA index

	if [ ! -e $DB_TARGET.bwt ]; then

		if [ "$BWA" != "" ]; then

		    (($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DB_RELEASE' for '$DB_TARGET_RELEASE' BWA Index"
		    
			# MK
		    echo "$DB_TARGET.bwt: $DB_TARGET
				$BWA index -a bwtsw $DB_TARGET;
		    " >> $MK

			MK_ALL="$MK_ALL $DB_TARGET.bwt"

		fi;

	fi;

	## SAMTOOLS index

	if [ ! -e $DB_TARGET.fai ]; then

		if [ "$SAMTOOLS" != "" ]; then

		    (($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DB_RELEASE' for '$DB_TARGET_RELEASE' SAMTOOLS Index"

			# MK
		    echo "$DB_TARGET.fai: $DB_TARGET
				$SAMTOOLS faidx $DB_TARGET;
		    " >> $MK

			MK_ALL="$MK_ALL $DB_TARGET.fai"

		fi;

	fi;

	## PICARD index

	if [ ! -e $DB_TARGET_FOLDER/$ASSEMBLY.dict ]; then

		if [ "$PICARD" != "" ]; then

		    (($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DB_RELEASE' for '$DB_TARGET_RELEASE' PICARD Index"

			# MK
		    echo "$DB_TARGET_FOLDER/$ASSEMBLY.dict: $REF
				$JAVA -jar $PICARD CreateSequenceDictionary \
					REFERENCE=$REF \
					OUTPUT=$DB_TARGET_FOLDER/$ASSEMBLY.dict ;
		    " >> $MK

			MK_ALL="$MK_ALL $DB_TARGET_FOLDER/$ASSEMBLY.dict"

		fi;

	fi;

fi;


################
# REFSEQ GENES #
################

if ((1)); then

	# DB
	DATABASE="refGene"
	DATABASE_NAME="RefGene"
	DATABASE_FULLNAME="Reference Genes"
	DATABASE_WEBSITE="https://genome.ucsc.edu/"
	#DATABASE_FULLNAME="Human Variation Sets in VCF Format" # dbsnp
	DATABASE_DESCRIPTION="Known human protein-coding and non-protein-coding genes taken from the NCBI RNA reference sequences collection (RefSeq)"


	# DB TARGET
	DB_TARGET=$REFSEQ_GENES;									# /STARK/databases/refGene/current/refGene.hg19.bed
	DB_ASSEMBLY=$ASSEMBLY; 										# hg19
	DB_TARGET_FILE=$(basename $DB_TARGET);						# refGene.hg19.bed
	DB_TARGET_FOLDER=$(dirname $DB_TARGET);						# /STARK/databases/refGene/current
	DB_TARGET_RELEASE=$(basename $DB_TARGET_FOLDER);			# current
	#DB_TARGET_RELEASE=$DATABASES_RELEASE;						# DATE
	DB_TARGET_DB_FOLDER=$(dirname $DB_TARGET_FOLDER);			# /STARK/databases/refGene

	# DB RELEASE
	DB_RELEASE=$DATABASES_RELEASE;								# DATE
	#[ -s $DB_TARGET_FOLDER ] && DB_RELEASE=$DB_TARGET_RELEASE	# if /STARK/databases/refGene/current exists then keep release
	DB_RELEASE_DATE=$DATABASES_RELEASE;							# DATE
	DB_RELEASE_FILE=$DB_TARGET_FILE;							# refGene.hg19.bed
	DB_RELEASE_FOLDER=$DB_TARGET_DB_FOLDER/$DB_RELEASE;			# /STARK/databases/refGene/DATE
	DB_RELEASE_FILE_PATH="$DB_RELEASE_FOLDER/$DB_RELEASE_FILE";	# /STARK/databases/refGene/DATE/refGene.hg19.bed

	# TMP files and folders
	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DB_RELEASE
	mkdir -p $DB_TMP


	if [ ! -e $DB_TARGET ] || (($UPDATE)); then

		# VERBOSE
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DB_RELEASE' for '$DB_TARGET_RELEASE' [$DB_ASSEMBLY]"


		# RefGenes TXT file
		REFSEQ_GENES_URL="ftp://hgdownload.soe.ucsc.edu/goldenPath/$DB_ASSEMBLY/database"
		REFSEQ_GENES_URL_FILE="refGene.txt.gz"
		
		REFSEQ_GENES_URL_FILE_DATE=$(curl -s -I $REFSEQ_GENES_URL/$REFSEQ_GENES_URL_FILE | grep "Last-Modified: " | sed "s/Last-Modified: //g" | sed "s/\r$//g")

		# RELEASE
		DB_RELEASE_FROM_DOWNLOAD=$(date -d "$REFSEQ_GENES_URL_FILE_DATE")

		DB_TARGET_TXT=$(echo $DB_TARGET | sed 's/bed$/txt/g')
		DB_TARGET_GENES=$(echo $DB_TARGET | sed 's/bed$/genes/g')
		DB_RELEASE_FILE_PATH_TXT=$(echo $DB_RELEASE_FILE_PATH | sed 's/bed$/txt/g')
		DB_RELEASE_FILE_PATH_GENES=$(echo $DB_RELEASE_FILE_PATH | sed 's/bed$/genes/g')
		
		# REFSEQ_GENE folder
		#REFSEQ_GENES_FOLDER=$(dirname $REFSEQ_GENES)
		#REFSEQ_GENES_TXT_FILE=$(basename $REFSEQ_GENES_TXT)
		DB_RELEASE_FILE_TXT=$(basename $DB_RELEASE_FILE_PATH_TXT)
		#REFSEQ_GENES_GENES_FILE=$(basename $REFSEQ_GENES_GENES)
		DB_RELEASE_FILE_GENES=$(basename $DB_RELEASE_FILE_PATH_GENES)
		#REFSEQ_GENES_FILE=$(basename $REFSEQ_GENES)


		# Other DB
		# geneReviews.txt.gz
		# geneReviewsDetail.txt.gz
		# gtexGene.txt.gz ???
		# gwasCatalog.txt.gz


		# Update
		if (($UPDATE)); then
			if [ -e $DB_TARGET_TXT ]; then mv -f $DB_TARGET_TXT $DB_TARGET_TXT.V$DATE; fi;
			if [ -e $DB_TARGET_GENES ]; then mv -f $DB_TARGET_GENES $DB_TARGET_GENES.V$DATE; fi;
			if [ -e $DB_TARGET ]; then mv -f $RDB_TARGET $DB_TARGET.V$DATE; fi;
		fi;
		COPY_MODE_REFSEQ_GENES=$COPY_MODE_DEFAULT

		# DATABASE Infos
		DB_INFOS_JSON='
		{
			"code": "'$DATABASE'",
			"name": "'$DATABASE_NAME'",
			"fullname": "'$DATABASE_FULLNAME'",
			"website": "'$DATABASE_WEBSITE'",
			"description": "'$DATABASE_DESCRIPTION'"
		}
		';
		echo "$DB_INFOS_JSON" > $DB_TMP/STARK.database

		# DATABASE Release Infos
		DB_RELEASE_INFOS_JSON='
		{
			"release": "'$DB_RELEASE_FROM_DOWNLOAD'",
			"date": "'$DB_RELEASE_DATE'",
			"files": [ "'$DB_TARGET_FILE'" ],
			"assembly": [ "'$DB_ASSEMBLY'" ],
			"download": {
				"methode": "'$DOWNLOAD_METHOD'",
				"URL": "'$REFSEQ_GENES_URL'",
				"file": "'$REFSEQ_GENES_URL_FILE'",
				"date": "'$REFSEQ_GENES_URL_FILE_DATE'"
			}
		}
		';
		echo "$DB_RELEASE_INFOS_JSON" > $DB_TMP/STARK.database.release

		# OUTPUT
		(($VERBOSE)) && echo "#[INFO] REFSEQ GENES URL   = $REFSEQ_GENES_URL"
		(($VERBOSE)) && echo "#[INFO] REFSEQ GENES TXT   = $DB_RELEASE_FILE_TXT"
		(($VERBOSE)) && echo "#[INFO] REFSEQ GENES BED   = $DB_RELEASE_FILE"
		(($VERBOSE)) && echo "#[INFO] REFSEQ GENES GENES = $DB_RELEASE_FILE_GENES"

		# Get refGene txt file
		echo "$DB_TMP/$DB_RELEASE_FILE_TXT: $DBFOLDER
			mkdir -p $DB_TMP/
			#wget --no-verbose -S -c -O $DB_TMP/$REFSEQ_GENES_URL_FILE $REFSEQ_GENES_URL/$REFSEQ_GENES_URL_FILE ;
			curl $REFSEQ_GENES_URL/$REFSEQ_GENES_URL_FILE -s -R -o $DB_TMP/$REFSEQ_GENES_URL_FILE;
			mkdir -p $DB_RELEASE_FOLDER
			chmod 0775 $DB_RELEASE_FOLDER
			#$GZ -dc  $DB_TMP/$REFSEQ_GENES_URL_FILE > $DB_TMP/$DB_RELEASE_FILE_TXT
			$GZ -N -d $DB_TMP/$REFSEQ_GENES_URL_FILE 
			if [ '$DB_TMP/refGene.txt' != '$DB_TMP/$DB_RELEASE_FILE_TXT' ]; then \
				cp -a $DB_TMP/refGene.txt $DB_TMP/$DB_RELEASE_FILE_TXT; \
				rm $DB_TMP/refGene.txt; \
			fi;
		" >> $MK
		# Copy in database folder
		echo "$DB_RELEASE_FILE_PATH_TXT: $DB_TMP/$DB_RELEASE_FILE_TXT
			$COPY_MODE_REFSEQ_GENES $DB_TMP/$DB_RELEASE_FILE_TXT $DB_RELEASE_FILE_PATH_TXT
		" >> $MK

		# Translate to refGene bed file
		echo "$DB_TMP/$DB_RELEASE_FILE: $DBFOLDER $DB_TMP/$DB_RELEASE_FILE_TXT
			mkdir -p $DB_TMP
			cat $DB_TMP/$DB_RELEASE_FILE_TXT | while IFS='' read -r line; do \
				CHROM=\$\$(echo \$\$line | awk '{print \$\$3}' | cut -d\"_\" -f1); \
				NM=\$\$(echo \$\$line | awk '{print \$\$2}'); \
				GENE=\$\$(echo \$\$line | awk '{print \$\$13}'); \
				STRAND=\$\$(echo \$\$line | awk '{print \$\$4}'); \
				echo \$\$line | awk '{print \$\$10}' | tr \",\" \"\\n\" | grep -v '^\$\$' > $TMP_DATABASES_DOWNLOAD_RAM/refGene.unsorted.bed.NM1 ; \
				echo \$\$line | awk '{print \$\$11}' | tr \",\" \"\\n\" | grep -v '^\$\$' > $TMP_DATABASES_DOWNLOAD_RAM/refGene.unsorted.bed.NM2; \
				paste $TMP_DATABASES_DOWNLOAD_RAM/refGene.unsorted.bed.NM1 $TMP_DATABASES_DOWNLOAD_RAM/refGene.unsorted.bed.NM2 | while read SS ; do \
					echo -e \"\$\$CHROM\t\$\$SS\t\$\$STRAND\t\$\$GENE\t\$\$NM\" >> $DB_TMP/refGene.unsorted.bed; \
				done; \
			done
			cat $DB_TMP/refGene.unsorted.bed | sort -k1,1V -k2,2n -k3,3n > $DB_TMP/$DB_RELEASE_FILE
			#rm -f $DB_RELEASE_FILE_PATH_TXT.NM1 $DB_RELEASE_FILE_PATH_TXT.NM2;
		" >> $MK

		# Generates refGene genes file from bed file
		echo "$DB_RELEASE_FILE_PATH: $DBFOLDER $DB_TMP/$DB_RELEASE_FILE
			$COPY_MODE_REFSEQ_GENES $DB_TMP/$DB_RELEASE_FILE $DB_RELEASE_FILE_PATH
		" >> $MK
		
		# Generates refGene genes file from bed file
		echo "$DB_TMP/$DB_RELEASE_FILE_GENES: $DBFOLDER $DB_TMP/$DB_RELEASE_FILE "$(dirname $REF)/$ASSEMBLY.dict"
			mkdir -p $DB_TMP
			awk -F'\t' 'substr(\$\$6,1,2)==\"NM\" {print \$\$0}' $DB_TMP/$DB_RELEASE_FILE | grep \$\$(grep \"@SQ\" "$(dirname $REF)/$ASSEMBLY.dict" | cut -f2 | cut -d: -f2 | tr '\n' ' ' | sed 's/chr/ -e ^chr/gi') > $DB_TMP/$DB_RELEASE_FILE_GENES
		" >> $MK

		# Generates refGene genes file from bed file
		echo "$DB_RELEASE_FILE_PATH_GENES: $DBFOLDER $DB_TMP/$DB_RELEASE_FILE_GENES
			$COPY_MODE_REFSEQ_GENES $DB_TMP/$DB_RELEASE_FILE_GENES $DB_RELEASE_FILE_PATH_GENES
		" >> $MK

		# Copy in database folder
		echo "$DB_TARGET: $DB_RELEASE_FILE_PATH $DB_RELEASE_FILE_PATH_TXT $DB_RELEASE_FILE_PATH_GENES
			# folder
			mkdir -p $DB_RELEASE_FOLDER
			chmod 0775 $DB_RELEASE_FOLDER
			# database release info
			[ ! -s $DB_TARGET_DB_FOLDER/STARK.database ] && cp $DB_TMP/STARK.database $DB_TARGET_DB_FOLDER/STARK.database
			cp $DB_TMP/STARK.database.release $DB_RELEASE_FOLDER/
			# links
			[ $DB_TARGET_RELEASE != $DB_RELEASE ] && ln -snf $DB_RELEASE/ $DB_TARGET_FOLDER
			[ latest != $DB_RELEASE ] && ln -snf $DB_RELEASE/ $DB_TARGET_DB_FOLDER/latest
			[ current != $DB_RELEASE ] && ln -snf $DB_RELEASE/ $DB_TARGET_DB_FOLDER/current
			# Clear
			rm -rf $DB_TMP;
		" >> $MK

		MK_ALL="$MK_ALL $DB_TARGET $DB_RELEASE_FILE_PATH $DB_RELEASE_FILE_PATH_TXT $DB_RELEASE_FILE_PATH_GENES"

	fi;

fi;


#########
# DBSNP #
#########

if ((1)); then

	# DB
	DATABASE="dbsnp"
	DATABASE_NAME="dbSNP"
	DATABASE_FULLNAME="Single-nucleotide polymorphism Database"
	DATABASE_WEBSITE="https://www.ncbi.nlm.nih.gov/snp/"
	DATABASE_DESCRIPTION="Human single nucleotide variations, microsatellites, and small-scale insertions and deletions along with publication, population frequency, molecular consequence, and genomic and RefSeq mapping information for both common variations and clinical mutations"

	# DB TARGET
	DB_TARGET=$VCFDBSNP;										# /STARK/databases/dbsnp/current/dbsnp.hg19.vcf.gz
	DB_ASSEMBLY=$ASSEMBLY; 										# hg19
	DB_TARGET_FILE=$(basename $DB_TARGET);						# dbsnp.hg19.vcf.gz
	DB_TARGET_FOLDER=$(dirname $DB_TARGET);						# /STARK/databases/dbsnp/current
	DB_TARGET_RELEASE=$(basename $DB_TARGET_FOLDER);			# current
	DB_TARGET_DB_FOLDER=$(dirname $DB_TARGET_FOLDER);			# /STARK/databases/dbsnp

	# DB RELEASE
	DB_RELEASE=$DATABASES_RELEASE;								# DATE
	#[ -s $DB_TARGET_FOLDER ] && DB_RELEASE=$DB_TARGET_RELEASE	# if /STARK/databases/dbsnp/current exists then keep release
	DB_RELEASE_DATE=$DATABASES_RELEASE;							# DATE
	DB_RELEASE_FILE=$DB_TARGET_FILE;							# dbsnp.hg19.vcf.gz
	DB_RELEASE_FOLDER=$DB_TARGET_DB_FOLDER/$DB_RELEASE;			# /STARK/databases/dbsnp/DATE
	DB_RELEASE_FILE_PATH="$DB_RELEASE_FOLDER/$DB_RELEASE_FILE";	# /STARK/databases/dbsnp/DATE/dbsnp.hg19.vcf.gz

	# TMP files and folders
	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DB_RELEASE
	mkdir -p $DB_TMP

	# $VCFDBSNP > $DB_TARGET
	#echo $DB_TARGET 

	if [ ! -e $DB_TARGET ] || (($UPDATE)); then

		# VERBOSE
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DB_RELEASE' for '$DB_TARGET_RELEASE' [$DB_ASSEMBLY]"

		# DBSNP folder
		#VCFDBSNP_FOLDER=$(dirname $VCFDBSNP) # DB_RELEASE_FOLDER
		VCFDBSNP_FILE=$(basename $DB_TARGET) # DB_RELEASE_FILE

		# Update
		if (($UPDATE)); then
			if [ -e $DB_TARGET ]; then mv -f $DB_TARGET $DB_TARGET.V$DATE; fi;
			if [ -e $DB_TARGET.tbi ]; then mv -f $DB_TARGET.tbi $DB_TARGET.tbi.V$DATE; fi;
		fi;
		#COPY_MODE_VCFDBSNP=$COPY_MODE_DEFAULT
		COPY_MODE_VCFDBSNP="cp -a "

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
		(($VERBOSE)) && echo "#[INFO] NCBI SPECIES=$SPECIES_NCBI"
		(($VERBOSE)) && echo "#[INFO] NCBI ASSEMBLY=$ASSEMBLY_NCBI"

		# DB release
		DBSNP_RELEASE="pre_build152" # "pre_build152" or "2.0"

		if [ "$DBSNP_RELEASE" == "2.0" ]; then

			# Find latest dbsnp file
			DBSNP_URL="ftp://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF"
			(($VERBOSE)) && echo "#[INFO] Check last DBSNP on '$DBSNP_URL'..."
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
			(($VERBOSE)) && echo "#[INFO] Check last DBSNP on '$DBSNP_URL'..."
			DBSNP_LATEST=$(curl -s $DBSNP_URL | awk '{print $9}' | grep $SPECIES_NCBI.*$ASSEMBLY_NCBI | sort -k1,1 | tail -n 1)
			if [ "$DBSNP_LATEST" == "" ]; then
				DBSNP_LATEST=$(curl -s $DBSNP_URL | awk -F"\">" '{print $2}' | awk -F"</A>" '{print $1}' | grep $SPECIES_NCBI.*$ASSEMBLY_NCBI | sort -k1,1 | tail -n 1)
			fi;
			#curl -s $DBSNP_URL -l  #| awk '{print $9}'
			DBSNP_LATEST_FILE=$DBSNP_URL$DBSNP_LATEST"/VCF/00-common_all.vcf.gz" # 00-All.vcf.gz


		fi;

		DBSNP_LATEST_FILE_BASENAME=$(basename $DBSNP_LATEST_FILE)

		DBSNP_URL=$(dirname $DBSNP_LATEST_FILE)
		DBSNP_URL_FILE=$(basename $DBSNP_LATEST_FILE)
		DBSNP_URL_FILE_DATE=$(curl -s -I $DBSNP_URL/$DBSNP_URL_FILE | grep "Last-Modified: " | sed "s/Last-Modified: //g" | sed "s/\r$//g")

		# RELEASE
		DB_RELEASE_FROM_DOWNLOAD=$(date -d "$DBSNP_URL_FILE_DATE")


		# DATABASE Infos
		DB_INFOS_JSON='
		{
			"code": "'$DATABASE'",
			"name": "'$DATABASE_NAME'",
			"fullname": "'$DATABASE_FULLNAME'",
			"website": "'$DATABASE_WEBSITE'",
			"description": "'$DATABASE_DESCRIPTION'"
		}
		';
		echo "$DB_INFOS_JSON" > $DB_TMP/STARK.database

		# DATABASE Release Infos
		DB_RELEASE_INFOS_JSON='
		{
			"release": "'$DB_RELEASE_FROM_DOWNLOAD'",
			"date": "'$DB_RELEASE_DATE'",
			"files": [ "'$DB_TARGET_FILE'" ],
			"assembly": [ "'$DB_ASSEMBLY'" ],
			"download": {
				"methode": "'$DOWNLOAD_METHOD'",
				"URL": "'$DBSNP_URL'",
				"file": "'$DBSNP_URL_FILE'",
				"date": "'$DBSNP_URL_FILE_DATE'"
			}
		}
		';
		echo "$DB_RELEASE_INFOS_JSON" > $DB_TMP/STARK.database.release


		(($VERBOSE)) && echo "#[INFO] DBSNP URL=$DBSNP_URL"
		(($VERBOSE)) && echo "#[INFO] DBSNP RELEASE=$DBSNP_RELEASE"
		(($VERBOSE)) && echo "#[INFO] DBSNP SPECIES/ASSEMBLY=$DBSNP_SPECIES_ASSEMBLY"
		(($VERBOSE)) && echo "#[INFO] DBSNP NUM=$NUM_NCBI"
		(($VERBOSE)) && echo "#[INFO] DBSNP LATEST=$DBSNP_LATEST"
		(($VERBOSE)) && echo "#[INFO] DBSNP LATEST FILE=$DBSNP_LATEST_FILE"
		(($VERBOSE)) && echo "#[INFO] DBSNP LATEST FILE BASENAME=$DBSNP_LATEST_FILE_BASENAME"

		# MK
		echo "$DB_TARGET: $DBFOLDER
			# Download
			mkdir -p $DB_TMP;
			chmod 0775 $DB_TMP;
			#wget --no-verbose -S -c -O $DB_TMP/$DBSNP_LATEST.tbi $DBSNP_LATEST_FILE.tbi;
			#wget --no-verbose -S -c -O $DB_TMP/$DBSNP_LATEST $DBSNP_LATEST_FILE;
			curl -s -R -o $DB_TMP/$DBSNP_LATEST.tbi $DBSNP_LATEST_FILE.tbi;
			curl -s -R -o $DB_TMP/$DBSNP_LATEST $DBSNP_LATEST_FILE;
			# Sources
			mkdir -p $DB_RELEASE_FOLDER/sources
			chmod 0775 $DB_RELEASE_FOLDER/sources
			$COPY_MODE_VCFDBSNP $DB_TMP/$DBSNP_LATEST $DB_RELEASE_FOLDER/sources/$DBSNP_LATEST_FILE_BASENAME;
			$COPY_MODE_VCFDBSNP $DB_TMP/$DBSNP_LATEST.tbi $DB_RELEASE_FOLDER/sources/$DBSNP_LATEST_FILE_BASENAME.tbi;
			# Add contig
			$BCFTOOLS view $DB_TMP/$DBSNP_LATEST | $SCRIPT_DIR/contig_NC_to_chr.sh > $DB_TMP/$DB_RELEASE_FILE.vcf
			touch $DB_TMP/$DB_RELEASE_FILE.vcf -r $DB_RELEASE_FOLDER/sources/$DBSNP_LATEST_FILE_BASENAME
			# BGZip & Tabix
			$BGZIP -c $DB_TMP/$DB_RELEASE_FILE.vcf > $DB_TMP/$DB_RELEASE_FILE
			$TABIX $DB_TMP/$DB_RELEASE_FILE
			# Copy
			mkdir -p $DB_RELEASE_FOLDER
			chmod 0775 $DB_RELEASE_FOLDER
			$COPY_MODE_VCFDBSNP $DB_TMP/$DB_RELEASE_FILE.tbi $DB_RELEASE_FILE_PATH.tbi;
			$COPY_MODE_VCFDBSNP $DB_TMP/$DB_RELEASE_FILE $DB_RELEASE_FILE_PATH;
			# database release info
			[ ! -s $DB_TARGET_DB_FOLDER/STARK.database ] && cp $DB_TMP/STARK.database $DB_TARGET_DB_FOLDER/STARK.database
			cp $DB_TMP/STARK.database.release $DB_RELEASE_FOLDER
			# Links
			[ $DB_TARGET_RELEASE != $DB_RELEASE ] && ln -snf $DB_RELEASE/ $DB_TARGET_FOLDER
			[ latest != $DB_RELEASE ] && ln -snf $DB_RELEASE/ $DB_TARGET_DB_FOLDER/latest
			[ current != $DB_RELEASE ] && ln -snf $DB_RELEASE/ $DB_TARGET_DB_FOLDER/current
			# Clear
			rm -rf $DB_TMP;
		" >> $MK

		MK_ALL="$MK_ALL $DB_TARGET"

	fi;

fi;



###########
# SNPEFF  #
###########

if ((1)); then

	# DB
	DATABASE="snpeff"
	DATABASE_NAME="SnpEff"
	DATABASE_FULLNAME="SnpEff Annotations"
	DATABASE_WEBSITE="http://snpeff.sourceforge.net/"
	DATABASE_DESCRIPTION="Genetic variant annotation and functional effect prediction toolbox"

	# DB TARGET
	DB_TARGET=$SNPEFF_DATABASES;								# /STARK/databases/snpeff/4.3t
	DB_ASSEMBLY=$ASSEMBLY; 										# hg19
	DB_TARGET_FILE=$(basename $DB_TARGET);						# 4.3t
	DB_TARGET_FOLDER=$DB_TARGET;								# /STARK/databases/snpeff/4.3t
	DB_TARGET_RELEASE=$(basename $DB_TARGET_FOLDER);			# 4.3t
	#[ ! -z SNPEFF_VERSION ] && DB_TARGET_RELEASE=$SNPEFF_VERSION	# 4.3t
	DB_TARGET_DB_FOLDER=$(dirname $DB_TARGET_FOLDER);			# /STARK/databases/snpeff

	# DB RELEASE
	DB_RELEASE=$DATABASES_RELEASE;								# DATE
	#[ -s $DB_TARGET_FOLDER ] && DB_RELEASE=$DB_TARGET_RELEASE	# if /STARK/databases/snpeff/4.3t exists then keep release
	DB_RELEASE_DATE=$DATABASES_RELEASE;							# DATE
	DB_RELEASE_FILE=$DB_TARGET_FILE;							# 4.3t
	DB_RELEASE_FOLDER=$DB_TARGET_DB_FOLDER/$DB_RELEASE;			# /STARK/databases/snpeff/DATE
	DB_RELEASE_FILE_PATH="$DB_RELEASE_FOLDER/$DB_RELEASE_FILE";	# /STARK/databases/snpeff/DATE

	# TMP files and folders
	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DB_RELEASE
	mkdir -p $DB_TMP

	
	#if [ ! -d $DB_TARGET ] || [ -z "$(ls -A $DB_TARGET)" ] || (($UPDATE)); then
	if [ ! -e $DB_TARGET ] || [ -z "$(ls -A $DB_TARGET)" ] || (($UPDATE)); then

		# DB folder exists
		if [ -d $DB_TARGET ] && [ -z "$(ls -A $DB_TARGET)" ]; then
			rm -rf $DB_TARGET 
			(($DEBUG)) && echo "#[INFO] folder $DB_TARGET removed because empty"
		elif [ -d $DB_TARGET ] && ! [ -z "$(ls -A $DB_TARGET)" ]; then
			mv $DB_TARGET $DB_TARGET.$DB_RELEASE
			(($DEBUG)) && echo "#[INFO] folder $DB_TARGET moved to $DB_TARGET.$DB_RELEASE because not empty"
		fi;

		# VERBOSE
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DB_RELEASE' for '$DB_TARGET_RELEASE' [$DB_ASSEMBLY]"

		# COPY mode
		COPY_MODE_SNPEFF=$COPY_MODE_DEFAULT #"cp -rf"

		if [ "$SNPEFF" != "" ]; then

			SNPEFF_DATABASES_FOLDER=$DB_TMP/databases


			SNPEFF_CMD="$JAVA $JAVA_FLAGS -jar $SNPEFF download $DB_ASSEMBLY -dataDir $SNPEFF_DATABASES_FOLDER 1>$SNPEFF_DATABASES_FOLDER/STARK.SNPEFF.download.log 2>$SNPEFF_DATABASES_FOLDER/STARK.SNPEFF.download.err";
			SNPEFF_CMD_ALT="if ((\$\$(grep 'ERROR while connecting to' $SNPEFF_DATABASES_FOLDER/STARK.SNPEFF.download.err -c))); then wget --no-verbose -O $DB_TMP/snpEff.$SNPEFF_VERSION.$ASSEMBLY.zip \$\$(grep 'ERROR while connecting to' $SNPEFF_DATABASES_FOLDER/STARK.SNPEFF.download.err | cut -f2 | cut -d' ' -f5); unzip $DB_TMP/snpEff.$SNPEFF_VERSION.$ASSEMBLY.zip -d $SNPEFF_DATABASES_FOLDER/; mv $SNPEFF_DATABASES_FOLDER/data/* $SNPEFF_DATABASES_FOLDER/; fi"
			(($DEBUG)) && echo "#[INFO] SNPEFF CMD = $SNPEFF_CMD"
			(($DEBUG)) && echo "#[INFO] SNPEFF CMD ALT = $SNPEFF_CMD_ALT"
			#(($VERBOSE)) && echo "#["

			# RELEASE
			SNPEFF_TOOL_VERSION=$(java -jar $SNPEFF -version | cut -f2)
			DB_RELEASE_FROM_DOWNLOAD=$DB_TARGET_RELEASE
			if [ ! -z SNPEFF_VERSION ]; then
				DB_RELEASE_FROM_DOWNLOAD=$SNPEFF_VERSION;
			elif [ ! -z SNPEFF_VERSION ]; then
				DB_RELEASE_FROM_DOWNLOAD=$SNPEFF_TOOL_VERSION;
			fi;
			

			# DATABASE Infos
			DB_INFOS_JSON='
			{
				"code": "'$DATABASE'",
				"name": "'$DATABASE_NAME'",
				"fullname": "'$DATABASE_FULLNAME'",
				"website": "'$DATABASE_WEBSITE'",
				"description": "'$DATABASE_DESCRIPTION'"
			}
			';
			echo "$DB_INFOS_JSON" > $DB_TMP/STARK.database

			# DATABASE Release Infos
			DB_RELEASE_INFOS_JSON='
			{
				"release": "'$DB_RELEASE_FROM_DOWNLOAD'",
				"date": "'$DB_RELEASE_DATE'",
				"assembly": [ "'$DB_ASSEMBLY'" ],
				"download": {
					"methode": "'$DOWNLOAD_METHOD'",
					"info": "From SnpEff binary"
				}
			}
			';
			echo "$DB_RELEASE_INFOS_JSON" > $DB_TMP/STARK.database.release


			# MK
			echo "$DB_TARGET/STARK.SNPEFF.download.complete: $DBFOLDER
				mkdir -p $SNPEFF_DATABASES_FOLDER;
				chmod 0775 $SNPEFF_DATABASES_FOLDER;
				-$SNPEFF_CMD
				$SNPEFF_CMD_ALT
				mkdir -p $DB_RELEASE_FOLDER
				$COPY_MODE_SNPEFF $SNPEFF_DATABASES_FOLDER/* $DB_RELEASE_FOLDER/;
				echo '#[INFO] STARK DATABASE '$DATABASE_NAME' release '$DB_RELEASE' for '$DB_TARGET_RELEASE' [$DB_ASSEMBLY] download complete'  > $DB_RELEASE_FOLDER/STARK.SNPEFF.download.complete
				# database release info
				[ ! -s $DB_TARGET_DB_FOLDER/STARK.database ] && cp $DB_TMP/STARK.database $DB_TARGET_DB_FOLDER/STARK.database
				cp $DB_TMP/STARK.database.release $DB_RELEASE_FOLDER
				# Links
				[ $DB_TARGET_RELEASE != $DB_RELEASE ] && ln -snf $DB_RELEASE/ $DB_TARGET_FOLDER
				[ latest != $DB_RELEASE ] && ln -snf $DB_RELEASE/ $DB_TARGET_DB_FOLDER/latest
				[ current != $DB_RELEASE ] && ln -snf $DB_RELEASE/ $DB_TARGET_DB_FOLDER/current
				# Clear
				rm -rf $DB_TMP;
			" >> $MK

			MK_ALL="$MK_ALL $DB_TARGET/STARK.SNPEFF.download.complete"


		fi;

	fi;

fi;



###########
# ANNOVAR #
###########


if ((1)); then

	# DB
	DATABASE="annovar"
	DATABASE_NAME="ANNOVAR"
	DATABASE_FULLNAME="ANNOVAR Annotations"
	DATABASE_WEBSITE="https://doc-openbio.readthedocs.io/projects/annovar/"
	DATABASE_DESCRIPTION="ANNOVAR is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes"

	# DB TARGET
	DB_TARGET=$ANNOVAR_DATABASES;								# /STARK/databases/annovar/current
	DB_ASSEMBLY=$ASSEMBLY; 										# hg19
	DB_TARGET_FILE=$(basename $DB_TARGET);						# current
	DB_TARGET_FOLDER=$DB_TARGET;								# /STARK/databases/annovar/current
	DB_TARGET_RELEASE=$(basename $DB_TARGET_FOLDER);			# current
	DB_TARGET_DB_FOLDER=$(dirname $DB_TARGET_FOLDER);			# /STARK/databases/annovar

	# DB RELEASE
	DB_RELEASE=$DATABASES_RELEASE;								# DATE
	DB_RELEASE_DATE=$DATABASES_RELEASE;							# DATE
	DB_RELEASE_FILE=$DB_TARGET_FILE;							# current
	DB_RELEASE_FOLDER=$DB_TARGET_DB_FOLDER/$DB_RELEASE;			# /STARK/databases/current/DATE
	DB_RELEASE_FILE_PATH="$DB_RELEASE_FOLDER/$DB_RELEASE_FILE";	# /STARK/databases/current/DATE

	# TMP files and folders
	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DB_RELEASE
	mkdir -p $DB_TMP

	
	#if [ ! -d $DB_TARGET ] || [ -z "$(ls -A $DB_TARGET)" ] || (($UPDATE)); then
	if [ ! -e $DB_TARGET ] || [ -z "$(ls -A $DB_TARGET)" ] || (($UPDATE)); then
	
		# DB folder exists
		if [ -d $DB_TARGET ] && [ -z "$(ls -A $DB_TARGET)" ]; then
			rm -rf $DB_TARGET 
			(($DEBUG)) && echo "#[INFO] folder $DB_TARGET removed because empty"
		elif [ -d $DB_TARGET ] && ! [ -z "$(ls -A $DB_TARGET)" ]; then
			mv $DB_TARGET $DB_TARGET.$DB_RELEASE
			(($DEBUG)) && echo "#[INFO] folder $DB_TARGET moved to $DB_TARGET.$DB_RELEASE because not empty"
		fi;

		# VERBOSE
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DB_RELEASE' for '$DB_TARGET_RELEASE' [$DB_ASSEMBLY]"

		# COPY mode
		COPY_MODE_ANNOVAR=$COPY_MODE_DEFAULT #"cp -rf"
		#COPY_MODE_ANNOVAR="mv "

		HOWARD_DB="$HOWARD_ANNOTATION,$HOWARD_ANNOTATION_REPORT,$HOWARD_ANNOTATION_MINIMAL,$ADDITIONAL_ANNOTATIONS" # ALL, CORE, snpeff $HOWARD_ANNOTATION

		if [ ! -e "$HOWARD_CONFIG_ANNOTATION" ]; then
			HOWARD_CONFIG_ANNOTATION=$HOWARDDIR/config.annotation.ini
		fi;
		HOWARD_CONFIG_ANNOTATION_FILE=$(basename $HOWARD_CONFIG_ANNOTATION)
		(($VERBOSE)) && echo "#[INFO] Annotations: $HOWARD_DB"


		if [ "$HOWARD" != "" ]; then

			# VCF
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

			# validation folder
			HOWARD_VALIDATION_FOLDER=$DB_TMP/example.annotation.validation
			HOWARD_DATABASES_FOLDER=$DB_TMP/databases


			#INPUT_VCF=$HOWARDDIR/docs/example.vcf
			HOWARD_CMD="$HOWARD --input=$INPUT_VCF --output=$HOWARD_VALIDATION_FOLDER/HOWARD.download.vcf --env=$CONFIG_TOOLS --annotation=$HOWARD_DB --annovar_folder=$ANNOVAR --annovar_databases=$HOWARD_DATABASES_FOLDER --config_annotation=$HOWARD_CONFIG_ANNOTATION --snpeff_jar=$SNPEFF --snpeff_threads=$THREADS --tmp=$HOWARD_VALIDATION_FOLDER --assembly=$ASSEMBLY --java=$JAVA --java_flags='\"$JAVA_FLAGS\"' --threads=1 --verbose >$HOWARD_DATABASES_FOLDER/HOWARD.download.log 2>$HOWARD_DATABASES_FOLDER/HOWARD.download.err;"
			(($DEBUG)) && echo "#[INFO] CMD= $HOWARD_CMD"
			#echo "#"

			# DATABASE Infos
			DB_INFOS_JSON='
			{
				"code": "'$DATABASE'",
				"name": "'$DATABASE_NAME'",
				"fullname": "'$DATABASE_FULLNAME'",
				"website": "'$DATABASE_WEBSITE'",
				"description": "'$DATABASE_DESCRIPTION'"
			}
			';
			echo "$DB_INFOS_JSON" > $DB_TMP/STARK.database

			# DATABASE Release Infos
			DB_RELEASE_INFOS_JSON='
			{
				"release": "'$DB_RELEASE_DATE'",
				"date": "'$DB_RELEASE_DATE'",
				"files": [ "'$HOWARD_CONFIG_ANNOTATION_FILE'" ],
				"assembly": [ "'$DB_ASSEMBLY'" ],
				"download": {
					"methode": "'$DOWNLOAD_METHOD'",
					"info": "From ANNOVAR binary",
					"databases": "'$HOWARD_DB'",
					"configuration": "'$HOWARD_CONFIG_ANNOTATION_FILE'"
				}
			}
			';
			echo "$DB_RELEASE_INFOS_JSON" > $DB_TMP/STARK.database.release


			# MK
			echo "$DB_TARGET/STARK.ANNOVAR.download.complete: $DBFOLDER
				# Folders
				mkdir -p $HOWARD_DATABASES_FOLDER;
				mkdir -p $HOWARD_VALIDATION_FOLDER;
				mkdir -p $DB_RELEASE_FOLDER;
				chmod 0775 $HOWARD_DATABASES_FOLDER;
				chmod 0775 $HOWARD_VALIDATION_FOLDER;
				chmod 0775 $DB_RELEASE_FOLDER;
				# Download command
				$HOWARD_CMD
				# Copy files
				$COPY_MODE_ANNOVAR $HOWARD_DATABASES_FOLDER/* $DB_RELEASE_FOLDER/;
				$COPY_MODE_ANNOVAR $HOWARD_CONFIG_ANNOTATION $DB_RELEASE_FOLDER/;
				# Download Complete file
				echo '#[INFO] STARK DATABASE '$DATABASE_NAME' release '$DB_RELEASE' for '$DB_TARGET_RELEASE' [$DB_ASSEMBLY] download complete'  > $DB_RELEASE_FOLDER/STARK.ANNOVAR.download.complete
				# Database release info
				[ ! -s $DB_TARGET_DB_FOLDER/STARK.database ] && cp $DB_TMP/STARK.database $DB_TARGET_DB_FOLDER/STARK.database
				cp $DB_TMP/STARK.database.release $DB_RELEASE_FOLDER
				# Links
				[ $DB_TARGET_RELEASE != $DB_RELEASE ] && ln -snf $DB_RELEASE/ $DB_TARGET_FOLDER
				[ latest != $DB_RELEASE ] && ln -snf $DB_RELEASE/ $DB_TARGET_DB_FOLDER/latest
				[ current != $DB_RELEASE ] && ln -snf $DB_RELEASE/ $DB_TARGET_DB_FOLDER/current
				# Clear
				rm -rf $DB_TMP;
			" >> $MK

			MK_ALL="$MK_ALL $DB_TARGET/STARK.ANNOVAR.download.complete"


		fi;

	fi;

fi;




# ALL
######

echo "all: $MK_ALL
	echo '#[INFO] Build release: $DATABASES_RELEASE' >> $DATABASES/STARK.download.releases
" >> $MK



# DATABASES INIT
(($VERBOSE)) && echo ""
#(($VERBOSE)) && echo "#[INFO] DATABASES INIT..."

if ((1)); then
	if (($BUILD)) || (($REBUILD)) || (($UPDATE)); then
		echo "#[INFO] DATABASES DOWNLOADING..."
		if (($VERBOSE)) || (($DEBUG)); then
			if (($DEBUG)); then
				make -k -j $THREADS $MK_OPTION -s -f $MK all;
			elif (($VERBOSE)); then
				make -k -j $THREADS $MK_OPTION -s -f $MK all;
			fi;
		else
			make -k -j $THREADS $MK_OPTION -f $MK all 1>$MK_LOG 2>$MK_ERR;
		fi;
		echo "#[INFO] DATABASES DOWNLOADED"
	else
		echo "## use --build, --rebuild or --update to download databases"
	fi;
fi;



if (($DEBUG)); then
	echo ""
	echo "### MAKEFILE"
	echo "#"
	echo "# MK=$MK"
	echo "# LOG=$MK_LOG"
	echo ""
	cat -n $MK;
	#cat -n $MK_LOG;
	#cat -n $MK_ERR;
fi;


## Release

if ((0)); then
	echo ""
	(($VERBOSE)) && echo "#[INFO] RELEASES CHECK"

	(($DEBUG)) && ls -l $DATABASES

	if ((1)); then
		echo ""
		(($VERBOSE)) && echo "#[INFO] DOWNLOAD GENOME $ASSEMBLY"
		(($VERBOSE)) && echo "#"
	fi;

fi;


# CLEAN
if ((0)); then
rm -Rf $TMP_DATABASES_DOWNLOAD_FOLDER
rm -Rf $TMP_DATABASES_DOWNLOAD_RAM
fi;

exit 0;
