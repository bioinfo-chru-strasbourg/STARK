#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARKDatabases"
SCRIPT_DESCRIPTION="STARK download and build databases"
SCRIPT_RELEASE="0.9.6.0"
SCRIPT_DATE="28/07/2022"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-11/12/2018: Script creation\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.1b-12/12/2018: Change to Makefile\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.2b-21/12/2018: Add update, build, rebuild and threads options. Change dbsnp source\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.3b-31/05/2019: Add APP configuration\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.4b-17/06/2020: Clarify code, organisation DB/RELEASE, add STARK.description\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.5.0-11/04/2021: Change snpEff download, some bugs fixed, option --current, remove option --rebuild\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.6.0-28/07/2022: Change snpEff download, add GATK databases, some bugs fixed\n";
RELEASE_NOTES=$RELEASE_NOTES"# 1.0.0-01/11/2023: Rewrite for STARK 19 new database structure: clean code, HOWARD database python, fix makefile rules (point to file, not directory), add CTAT/arriba, use aria2c\n";


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
	echo "# --databases_list=<STRING>                List of Databases to consider";
	echo "#                                          Format: 'database1,database2,...'";
	echo "#                                          Default: 'ALL' for all available databases";
	echo "#                                          Available databases: 'dbsnp'";
	echo "# --current                                Make new databases as current.";
	echo "# --build                                  Build all databases.";
	echo "# --update                                 Update databases (latest dbSNP databases) and build if needed.";
	echo "# --threads                                Number of threads (depend on system/proxy...).";
	echo "# --verbose                                VERBOSE option";
	echo "# --debug                                  DEBUG option";
	echo "# --release                                RELEASE option";
	echo "# --help                                   HELP option";
	echo "#";
}

header;

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "e:cbut:vdnh" --long "env:,app:,application:,databases:,databases_list:,current,build,update,threads:,verbose,debug,release,help" -- "$@" 2> /dev/null)
PARAM=$@

eval set -- "$ARGS"
while true
do
	case "$1" in
		-e|--env|--app|--application)
			APP="$2"
			shift 2
			;;
		--databases)
			DATABASES="$2"
			shift 2
			;;
		--databases_list)
			DATABASES_LIST_INPUT=$(echo "$2" | tr "," " ")
			shift 2
			;;
		-v|--verbose)
			VERBOSE=1
			shift 1
			;;
		-c|--current)
			CURRENT=1
			shift 1
			;;
		-b|--build)
			BUILD=1
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

# Script folder
echo "# SEARCHING SCRIPTS"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "# DONE"

# Configuration
echo "# SEARCHING CONFIG APPS"
ENV_CONFIG=$(find -L $SCRIPT_DIR/.. -name config.app)
echo "# DONE"
echo "# SOURCE CONFIGS"
echo $ENV_CONFIG
source $ENV_CONFIG 
echo "# DONE"

# FUNCTIONS
#############

in_array () 
{ 
    param=$1;
    shift;
    for elem in "$@";
    do
        [[ "$param" = "$elem" ]] && return 0;
    done;
    return 1
}

ACTION=0
[ $BUILD ] || [ $UPDATE ] && ACTION=1;

[ ! -z $DATABASES ] && [ ! -d $DATABASES ] && mkdir -p $DATABASES && echo "#[INFO] Create databases folder '$DATABASES' "
[ "$DATABASES_LIST_INPUT" == "" ] && DATABASES_LIST_INPUT="ALL"
echo "#[INFO] Databases list: '$DATABASES_LIST_INPUT' "

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
CORES=$(nproc)

if ! [[ $THREADS =~ $re ]] || [ -z "$THREADS" ] || [ "$THREADS" == "" ] || [ $THREADS -gt $CORES ] ; then
	CORES_FREE=0
	THREADS=$(($CORES-$CORES_FREE))
fi;

if [[ $THREADS_INPUT =~ $re ]] && [ "$THREADS_INPUT" != "" ]; then
	THREADS=$THREADS_INPUT;
fi;

TMP_DATABASES_DOWNLOAD_FOLDER=$TMP_FOLDER_TMP/$RANDOM$RANDOM
mkdir -p $TMP_DATABASES_DOWNLOAD_FOLDER
TMP_DATABASES_DOWNLOAD_RAM="$(mktemp -d -p /dev/shm/)"

if [ "$TMP_DATABASES_DOWNLOAD_RAM" == "" ]; then
	TMP_DATABASES_DOWNLOAD_RAM=$TMP_DATABASES_DOWNLOAD_FOLDER;
fi;

DATE=$(date '+%Y%m%d-%H%M%S')

echo ""
echo "#[INFO] DB_RELEASE=$DATE"
echo "#[INFO] ASSEMBLY=$ASSEMBLY"
echo "#[INFO] THREADS=$THREADS"

MK=$TMP_DATABASES_DOWNLOAD_FOLDER/mk
MK_LOG=$TMP_DATABASES_DOWNLOAD_FOLDER/mk.log
MK_ERR=$TMP_DATABASES_DOWNLOAD_FOLDER/mk.err
MK_ALL=""
> $MK

DOWNLOAD_METHOD="STARK Databases downloading script [$SCRIPT_RELEASE-$SCRIPT_DATE]"

###########
# GENOMES #
###########
DATABASE="genomes"
DATABASE_NAME="Genomes"
DATABASE_FULLNAME="Reference Genome Sequences Assembly"
DATABASE_WEBSITE="https://genome.ucsc.edu/"
DATABASE_DESCRIPTION="Reference sequence was produced by the Genome Reference Consortium, and is composed of genomic sequence, primarily finished clones that were sequenced as part of the Human Genome Project"

if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then
	DBFOLDER_GENOME=$DATABASES/genomes
	if [ ! -e $DBFOLDER_GENOME/current ]; then
		mkdir -p $DBFOLDER_GENOME/current;
	fi;
	
	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_GENOME/current/$ASSEMBLY ] || (($UPDATE)); then
		
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for ' [$ASSEMBLY]"
		
		if (($UPDATE)); then
			if [ -e $DBFOLDER_GENOME/current/$ASSEMBLY ]; then mv -f $DBFOLDER_GENOME/current/$ASSEMBLY $DBFOLDER_GENOME.V$DATE; fi;
		fi;
		
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

		echo "$GENOME: $DBFOLDER_GENOME
			howard databases --assembly='$ASSEMBLY' --download-genomes=$DBFOLDER_GENOME/$DATE --download-genomes-contig-regex=$GENOME_REGEX;
			[ ! -e $DBFOLDER_GENOME/current/$ASSEMBLY ] || unlink $DBFOLDER_GENOME/current/$ASSEMBLY;
			ln -snf $DBFOLDER_GENOME/$DATE/$ASSEMBLY $DBFOLDER_GENOME/current/$ASSEMBLY;
			-[ ! -s $DBFOLDER_GENOME/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_GENOME/STARK.database && chmod o+r $DBFOLDER_GENOME/STARK.database;
			mv $GENOME $GENOME.tmp;
			mv $GENOME.fai $GENOME.fai.tmp;
			cut -f1 $GENOME.fai.tmp|sort -k1V|parallel -k '$SAMTOOLS faidx $GENOME.tmp {}' > $GENOME;
			rm -rf $DB_TMP;
			rm -rf  $DBFOLDER_GENOME/$DATE/$ASSEMBLY/*.tmp;
			rm -rf  $DBFOLDER_GENOME/$DATE/$ASSEMBLY/*.tmp.fai;
		" >> $MK
		MK_ALL="$MK_ALL $GENOME"
	fi;

	# Samtools hts-ref
	if [ ! -d $GENOME.hts-ref ]; then
		if [ "$SAMTOOLS" != "" ]; then
			echo "$GENOME.hts-ref/done: $GENOME
				mkdir -p $GENOME.hts-ref;
				perl $(dirname $SAMTOOLS)/seq_cache_populate.pl -root $GENOME.hts-ref $GENOME;
			" >> $MK
			MK_ALL="$MK_ALL $GENOME.hts-ref/done"
		fi;
	fi;

	# Samtools index
	if [ ! -e $GENOME.fai ]; then
		if [ "$SAMTOOLS" != "" ]; then
		    echo "$GENOME.fai: $GENOME
				$SAMTOOLS faidx $GENOME;
		    " >> $MK
			MK_ALL="$MK_ALL $GENOME.fai"
		fi;
	fi;

	## BOWTIE index
	if [ ! -e $(dirname $GENOME)/$ASSEMBLY.rev.1.bt2 ]; then
		if [ "$BOWTIE" != "" ]; then
			echo "$(dirname $GENOME)/$ASSEMBLY.rev.1.bt2: $GENOME
				$(dirname $BOWTIE)/bowtie2-build --threads $THREADS --packed $GENOME $(dirname $GENOME)/$ASSEMBLY;
			" >> $MK
			MK_ALL="$MK_ALL $(dirname $GENOME)/$ASSEMBLY.rev.1.bt2"
		fi;
	fi;

	## BWA index
	if [ ! -e $GENOME.bwt ]; then
		if [ "$BWA" != "" ]; then
			echo "$GENOME.bwt: $GENOME
				$BWA index -a bwtsw $GENOME;
			" >> $MK
			MK_ALL="$MK_ALL $GENOME.bwt"
		fi;
	fi;

	## BWA2 index
	if [ ! -e $GENOME.bwt.2bit.64 ]; then 
		if [ "$BWA2" != "" ]; then
			if [ BWA2_INDEX == "1" ]; then
				echo "$GENOME.bwt.2bit.64: $GENOME
					$BWA2 index $GENOME;
				" >> $MK
				MK_ALL="$MK_ALL $GENOME.bwt.2bit.64"
			fi;
		fi;
	fi;

	## PICARD index
	if [ ! -e $(dirname $GENOME)/$ASSEMBLY.dict ]; then
		if [ "$PICARD" != "" ]; then
			echo "$(dirname $GENOME)/$ASSEMBLY.dict: $GENOME
				$JAVA -jar $PICARD CreateSequenceDictionary \
					-REFERENCE $GENOME \
					-OUTPUT $(dirname $GENOME)/$ASSEMBLY.dict;
			" >> $MK
			MK_ALL="$MK_ALL $(dirname $GENOME)/$ASSEMBLY.dict"
		fi;
	fi;

	## GATK IMG
	if [ ! -e $GENOME.img ]; then
		if [ "$PICARD" != "" ]; then
			echo "$GENOME.img: $GENOME
				$JAVA -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS -jar $GATK4 BwaMemIndexImageCreator \
					--input $GENOME \
					--output $GENOME.img;
			" >> $MK
			MK_ALL="$MK_ALL $GENOME.img"
		fi;	
	fi;

	## STAR index 
	if [ ! -e $(dirname $GENOME)/$(basename $GENOME).star.idx ]; then
		if [ "$STAR" != "" ]; then
			echo "$(dirname $GENOME)/$(basename $GENOME).star.idx/done: $GENOME $DBFOLDER_GENCODE/current/$ASSEMBLY/gencode.v$GENCODE_VERSION.annotation.gtf
				mkdir -p $(dirname $GENOME)/$(basename $GENOME).star.idx;
				STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $(dirname $GENOME)/$(basename $GENOME).star.idx --genomeFastaFiles $GENOME --sjdbGTFfile $DBFOLDER_GENCODE/current/$ASSEMBLY/gencode.v$GENCODE_VERSION.annotation.gtf;
			" >> $MK
			MK_ALL="$MK_ALL $(dirname $GENOME)/$(basename $GENOME).star.idx/done"
		fi;
	fi;

fi;

###############################
# GATK4 VARIANT RECALIBRATION #
###############################
DATABASE="gatk4"
DATABASE_NAME="GATK4"
DATABASE_FULLNAME="GATK4 Databases"
DATABASE_WEBSITE="https://www.broadinstitute.org/"
DATABASE_DESCRIPTION="Databases for GATK4 Variant Recalibration"
# Ressource documentation
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
# What is variant recalibration
# https://gatk.broadinstitute.org/hc/en-us/articles/13832765070875-VariantRecalibrator
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then

	DBFOLDER_GATK_URL_FULL_COPY=0
	GATK_RESOURCE_NB=0
	MK_DBFOLDER_GATK_ALL=""
	> $MK.existing_gatk_db

	if [ ! -e $DBFOLDER_GATK/current ]; then
		mkdir -p $DBFOLDER_GATK/current;
	fi;

	for GATK_RESOURCE in $GATK_DATABASES_LIST; do
		DB_TARGET_GATK=$DBFOLDER_GATK/current/$ASSEMBLY/$(echo $GATK_RESOURCE | cut -d: -f2);		# /STARK/databases/gatk/current/hg19/1000G_omni2.5.b37.vcf.gz
		DB_RELEASE_FILE=$(basename $DB_TARGET_GATK);												# 1000G_omni2.5.b37.vcf.gz
		DB_RELEASE_FILE_PATH="$DBFOLDER_GATK/$DATE/$ASSEMBLY/$DB_RELEASE_FILE";						# /STARK/databases/gatk/DATE/hg19/1000G_omni2.5.b37.vcf.gz

		DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
		mkdir -p $DB_TMP

		if (($UPDATE)); then
			if [ -e $DBFOLDER_GATK/current/$ASSEMBLY ]; then mv -f $DBFOLDER_GATK/current/$ASSEMBLY $DBFOLDER_GATK.V$DATE; fi;
		fi;

		if [ -e $DB_TARGET_GATK ]; then
			echo "$DB_RELEASE_FILE_PATH: $DBFOLDER $GENOME
				mkdir -p $DBFOLDER_GATK/$DATE/$ASSEMBLY
				rsync -ar $DB_TARGET_GATK $DB_RELEASE_FILE_PATH
				rsync -ar $DB_TARGET_GATK.tbi $DB_RELEASE_FILE_PATH.tbi
			" >> $MK.existing_gatk_db
			MK_DBFOLDER_GATK_ALL_existing_gatk_db="$MK_DBFOLDER_GATK_ALL_existing_gatk_db $DB_RELEASE_FILE_PATH"
			MK_ALL_existing_gatk_db="$MK_ALL_existing_gatk_db $DB_RELEASE_FILE_PATH" 

		else
			((GATK_RESOURCE_NB++))
			(($VERBOSE)) && echo ""
			(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME/$GATK_RESOURCE' release '$DATE' for '$ASSEMBLY'"

			DBFOLDER_GATK_URL_FILE="$(basename $DB_TARGET_GATK)"
			DB_TARGET_FILE_LIST="$DB_TARGET_FILE_LIST $(basename $DB_TARGET_GATK)"

			if [[ "$DBFOLDER_GATK_URL_FILE" =~ .*"dbsnp138"*. ]]; then
				DBFOLDER_GATK_URL=$DBFOLDER_GATK_URL_DBSNP
			else
				DBFOLDER_GATK_URL=$DBFOLDER_GATK_URL_DEFAULT
			fi;

			DBFOLDER_GATK_URL_FILE_DATE=$(curl -s -I $DBFOLDER_GATK_URL/$DBFOLDER_GATK_URL_FILE | grep "Last-Modified: " | sed "s/Last-Modified: //g" | sed "s/\r$//g")
	
			DB_RELEASE_FROM_DOWNLOAD=$(date -d "$DBFOLDER_GATK_URL_FILE_DATE")
			DB_RELEASE_FILE=$(basename $DB_RELEASE_FILE_PATH)

			(($VERBOSE)) && echo "#[INFO] GATK Resource URL   = $DBFOLDER_GATK_URL"
			(($VERBOSE)) && echo "#[INFO] GATK Resource File  = $DBFOLDER_GATK_URL_FILE"
			DBFOLDER_GATK_DOWNLOAD_MULTITHREAD=1

			if ! (($DBFOLDER_GATK_URL_FULL_COPY)); then
				DBFOLDER_GATK_URL_PATH=$(echo $DBFOLDER_GATK_URL | sed 's#^https://##gi'  | sed 's#^http://##gi')
				if (($DBFOLDER_GATK_DOWNLOAD_MULTITHREAD)); then
					echo "$DB_TMP/original/chr_name_conv.txt: $GENOME
						mkdir -p $DB_TMP/original
						cat $GENOME.fai | awk '{a=\$\$1; gsub(\"^chr\",\"\",a); print a\" \"\$\$1}' > $DB_TMP/original/chr_name_conv.txt
				
					" >> $MK
				else
					echo "$DB_TMP/original/chr_name_conv.txt: $GENOME
						mkdir -p $DB_TMP
						wget -q -r --no-parent $DBFOLDER_GATK_URL --directory-prefix=$DB_TMP/original
						mv $DB_TMP/original/$DBFOLDER_GATK_URL_PATH/*vcf.gz* $DB_TMP/original/
						cat $GENOME.fai | awk '{a=\$\$1; gsub(\"^chr\",\"\",a); print a\" \"\$\$1}' > $DB_TMP/original/chr_name_conv.txt
					" >> $MK
				fi;
				DBFOLDER_GATK_URL_FULL_COPY=1
			fi;

			echo "$DB_TMP/$DB_RELEASE_FILE: $GENOME $DB_TMP/original/chr_name_conv.txt
				mkdir -p $DB_TMP/original
				if (($DBFOLDER_GATK_DOWNLOAD_MULTITHREAD)); then \
					curl $DBFOLDER_GATK_URL/$DBFOLDER_GATK_URL_FILE -s -R -o $DB_TMP/original/$DBFOLDER_GATK_URL_FILE; \
				else \
					rsync -ar $DB_TMP/original/$DBFOLDER_GATK_URL_FILE $DB_TMP/$DBFOLDER_GATK_URL_FILE; \
				fi;
				$BGZIP -dc -@$THREADS $DB_TMP/original/$DBFOLDER_GATK_URL_FILE | sed 's/\\t\$\$//gi' | $BGZIP -l1 -@$THREADS -c > $DB_TMP/$DBFOLDER_GATK_URL_FILE.tmp.vcf.gz
				$TABIX $DB_TMP/$DBFOLDER_GATK_URL_FILE.tmp.vcf.gz
				$BCFTOOLS reheader --fai $GENOME.fai --threads $THREADS $DB_TMP/$DBFOLDER_GATK_URL_FILE.tmp.vcf.gz > $DB_TMP/$DBFOLDER_GATK_URL_FILE.tmp2.vcf.gz
				$TABIX $DB_TMP/$DBFOLDER_GATK_URL_FILE.tmp2.vcf.gz
				$BCFTOOLS annotate --rename-chrs $DB_TMP/original/chr_name_conv.txt --threads $THREADS $DB_TMP/$DBFOLDER_GATK_URL_FILE.tmp2.vcf.gz | grep -v '^##contig=<ID=[^,]*>' | $BGZIP -l1 -@$THREADS > $DB_TMP/$DB_RELEASE_FILE.tmp3.vcf.gz
				$TABIX $DB_TMP/$DBFOLDER_GATK_URL_FILE.tmp3.vcf.gz
				$BCFTOOLS view --threads $THREADS -r \$\$(cat $DB_TMP/original/chr_name_conv.txt | cut -d' ' -f2 | tr '\\n' ',') $DB_TMP/$DB_RELEASE_FILE.tmp3.vcf.gz | $BGZIP -@$THREADS > $DB_TMP/$DB_RELEASE_FILE 
				$TABIX $DB_TMP/$DBFOLDER_GATK_URL_FILE
				rm -f $DB_TMP/$DBFOLDER_GATK_URL_FILE.tmp*
				mkdir -p $DBFOLDER_GATK/$DATE/$ASSEMBLY
				chmod 0775 $DBFOLDER_GATK/$DATE/$ASSEMBLY
			" >> $MK

			echo "$DB_RELEASE_FILE_PATH: $DB_TMP/$DB_RELEASE_FILE
				mkdir -p $DBFOLDER_GATK/$DATE/$ASSEMBLY/original
				rsync -ar $DB_TMP/original/$DBFOLDER_GATK_URL_FILE $DBFOLDER_GATK/$DATE/$ASSEMBLY/original/
				rsync -ar $DB_TMP/$DB_RELEASE_FILE $DB_RELEASE_FILE_PATH
				rsync -ar $DB_TMP/$DB_RELEASE_FILE.tbi $DB_RELEASE_FILE_PATH.tbi
			" >> $MK
			MK_DBFOLDER_GATK_ALL="$MK_DBFOLDER_GATK_ALL $DB_RELEASE_FILE_PATH"
			MK_ALL="$MK_ALL $DB_RELEASE_FILE_PATH" 
		fi;
	done;

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

	DB_RELEASE_INFOS_JSON='
	{
		"release": "'$DB_RELEASE_FROM_DOWNLOAD'",
		"date": "'$DATE'",
		"files": [ "'$(echo $DB_TARGET_FILE_LIST | sed 's/ /", "/gi')'" ],
		"assembly": [ "'$DB_ASSEMBLY'" ],
		"download": {
			"methode": "'$DOWNLOAD_METHOD'",
			"URL": "'$DBFOLDER_GATK_URL'",
			"file": "'$(echo $DB_TARGET_FILE_LIST | sed 's/ /,/gi')'",
			"date": "'$DBFOLDER_GATK_URL_FILE_DATE'"
		}
	}
	';
	echo "$DB_RELEASE_INFOS_JSON" > $DB_TMP/STARK.database.release

	if (($GATK_RESOURCE_NB)); then
		cat $MK.existing_gatk_db >> $MK
		MK_DBFOLDER_GATK_ALL="$MK_DBFOLDER_GATK_ALL $MK_DBFOLDER_GATK_ALL_existing_gatk_db"
		MK_ALL="$MK_ALL $MK_ALL_existing_gatk_db" 
		echo "$DBFOLDER_GATK/current/$ASSEMBLY/done: $MK_DBFOLDER_GATK_ALL
			mkdir -p $DBFOLDER_GATK/$DATE/$ASSEMBLY/original
			chmod 0775 $DBFOLDER_GATK/$DATE/$ASSEMBLY -R
			-[ ! -s $DBFOLDER_GATK/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_GATK/STARK.database
			cp $DB_TMP/STARK.database.release $DBFOLDER_GATK/$DATE/$ASSEMBLY/
			chmod o+r $DBFOLDER_GATK/STARK.database $DBFOLDER_GATK/$DATE/$ASSEMBLY/STARK.database.release
			[ ! -e $DBFOLDER_GATK/current/$ASSEMBLY ] || unlink $DBFOLDER_GATK/current/$ASSEMBLY
			ln -snf $DBFOLDER_GATK/$DATE/$ASSEMBLY $DBFOLDER_GATK/current/$ASSEMBLY
			rm -rf $DB_TMP;
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_GATK/current/$ASSEMBLY/done" 
	fi;
fi;


##########
# SNPEFF #
##########
DATABASE="snpeff"
DATABASE_NAME="SnpEff"
DATABASE_FULLNAME="SnpEff Annotations"
DATABASE_WEBSITE="http://snpeff.sourceforge.net/"
DATABASE_DESCRIPTION="Genetic variant annotation and functional effect prediction toolbox"

if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then

	DBFOLDER_SNPEFF=$(dirname $SNPEFF_DATABASES) 
	if [ ! -e $DBFOLDER_SNPEFF/current ]; then
		mkdir -p $DBFOLDER_SNPEFF/current;
	fi;

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_SNPEFF/current/$ASSEMBLY ] || (($UPDATE)); then
		
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for ' [$ASSEMBLY]"
		
		if (($UPDATE)); then
			if [ -e $DBFOLDER_SNPEFF/current/$ASSEMBLY ]; then mv -f $DBFOLDER_SNPEFF/current/$ASSEMBLY $DBFOLDER_SNPEFF.V$DATE; fi;
		fi;

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

		echo "$DBFOLDER_SNPEFF/done: $DBFOLDER
			howard databases --assembly='$ASSEMBLY' --download-snpeff=$DBFOLDER_SNPEFF/$DATE
			-[ ! -s $DBFOLDER_SNPEFF/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_SNPEFF/STARK.database && chmod o+r $DBFOLDER_SNPEFF/STARK.database 
			[ ! -e $DBFOLDER_SNPEFF/current/$ASSEMBLY ] || unlink $DBFOLDER_SNPEFF/current/$ASSEMBLY
			ln -snf $DBFOLDER_SNPEFF/$DATE/$ASSEMBLY $DBFOLDER_SNPEFF/current/$ASSEMBLY
			rm -rf $DB_TMP;
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_SNPEFF/done"
	fi;
fi;

###########
# ANNOVAR #
###########
DATABASE="annovar"
DATABASE_NAME="ANNOVAR"
DATABASE_FULLNAME="ANNOVAR Annotations"
DATABASE_WEBSITE="https://doc-openbio.readthedocs.io/projects/annovar/"
DATABASE_DESCRIPTION="ANNOVAR is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes"

if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then
	
	DBFOLDER_ANNOVAR=$(dirname $ANNOVAR_DATABASES)
	if [ ! -e $DBFOLDER_ANNOVAR/current ]; then
		mkdir -p $DBFOLDER_ANNOVAR/current;
	fi;

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_ANNOVAR/current/$ASSEMBLY ] || (($UPDATE)); then
		
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for ' [$ASSEMBLY]"

		if (($UPDATE)); then
			if [ -e $DBFOLDER_ANNOVAR/current/$ASSEMBLY ]; then mv -f $DBFOLDER_ANNOVAR/current/$ASSEMBLY $DBFOLDER_ANNOVAR.V$DATE; fi;
		fi;

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
	
		echo "$DBFOLDER_ANNOVAR/done: $DBFOLDER
			howard databases --assembly='$ASSEMBLY' --download-annovar=$DBFOLDER_ANNOVAR/$DATE --download-annovar-files='$ANNOVAR_FILES'
			-[ ! -s $DBFOLDER_ANNOVAR/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_ANNOVAR/STARK.database && chmod o+r $DBFOLDER_ANNOVAR/STARK.database 
			[ ! -e $DBFOLDER_ANNOVAR/current/$ASSEMBLY ] || unlink $DBFOLDER_ANNOVAR/current/$ASSEMBLY
			ln -snf $DBFOLDER_ANNOVAR/$DATE/$ASSEMBLY $DBFOLDER_ANNOVAR/current/$ASSEMBLY
			rm -rf $DB_TMP;
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_ANNOVAR/done"
	fi;
fi;

###########
# REFGENE #
###########
DATABASE="refGene"
DATABASE_NAME="RefGene"
DATABASE_FULLNAME="Reference Genes"
DATABASE_WEBSITE="https://genome.ucsc.edu/"
DATABASE_DESCRIPTION="Known human protein-coding and non-protein-coding genes taken from the NCBI RNA reference sequences collection (RefSeq)"

if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then
	
	DBFOLDER_REFGENE=$DBFOLDER/refGene
	if [ ! -e $DBFOLDER_REFGENE/current ]; then
		mkdir -p $DBFOLDER_REFGENE/current;
	fi;

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_REFGENE/current/$ASSEMBLY ] || (($UPDATE)); then
		
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for ' [$ASSEMBLY]"
		
		if (($UPDATE)); then
			if [ -e $DBFOLDER_REFGENE/current/$ASSEMBLY ]; then mv -f $DBFOLDER_REFGENE/current/$ASSEMBLY $DBFOLDER_REFGENE.V$DATE; fi;
		fi;

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

		echo "$DBFOLDER_REFGENE/done: $DBFOLDER
			howard databases --assembly='$ASSEMBLY' --download-refseq=$DBFOLDER_REFGENE/$DATE --download-refseq-format-file='ncbiRefSeq.txt' ;
			mv DBFOLDER_REFGENE/$DATE/$ASSEMBLY/ncbiRefSeq.bed $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/refGene.$ASSEMBLY.bed;
			-[ ! -s $DBFOLDER_REFGENE/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_REFGENE/STARK.database && chmod o+r $DBFOLDER_REFGENE/STARK.database;
			[ ! -e $DBFOLDER_REFGENE/current/$ASSEMBLY ] || unlink $DBFOLDER_REFGENE/current/$ASSEMBLY;
			ln -snf $DBFOLDER_REFGENE/$DATE/$ASSEMBLY $DBFOLDER_REFGENE/current/$ASSEMBLY;
			" >> $MK

		MK_ALL="$MK_ALL $DBFOLDER_REFGENE/done"
	fi;
fi;

#########
# dbSNP #
#########
DATABASE="dbsnp"
DATABASE_NAME="dbSNP"
DATABASE_FULLNAME="Single-nucleotide polymorphism Database"
DATABASE_WEBSITE="https://www.ncbi.nlm.nih.gov/snp/"
DATABASE_DESCRIPTION="Human single nucleotide variations, microsatellites, and small-scale insertions and deletions along with publication, population frequency, molecular consequence, and genomic and RefSeq mapping information for both common variations and clinical mutations"

if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then
	
	DBFOLDER_DBSNP=$(dirname $DBSNP_DATABASES)
	if [ ! -e $DBFOLDER_DBSNP/current ]; then
		mkdir -p $DBFOLDER_DBSNP/current;
	fi;

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_DBSNP/current/$ASSEMBLY ] || (($UPDATE)); then
		
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for ' [$ASSEMBLY]"

		if (($UPDATE)); then
			if [ -e $DBFOLDER_DBSNP/current/$ASSEMBLY ]; then mv -f $DBFOLDER_DBSNP/current/$ASSEMBLY $DBFOLDER_DBSNP.V$DATE; fi;
		fi;

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
		
		echo "$DBFOLDER_DBSNP/done: $DBFOLDER
			howard databases --assembly='$ASSEMBLY' --genomes-folder=$DBFOLDER_GENOME/current/ --download-dbsnp=$DBFOLDER_DBSNP/$DATE --download-dbsnp-releases='$DBSNP_VERSION_DOWNLOAD' --download-dbsnp-vcf --download-dbsnp-parquet —memory=$MEMORY --threads=$THREADS;
			$BCFTOOLS sort $DBFOLDER_DBSNP/$DATE/$ASSEMBLY/$DBSNP_VERSION/dbsnp.vcf.gz -o $DBFOLDER_DBSNP/$DATE/$ASSEMBLY/dbsnp.$DBSNP_VERSION.vcf.gz;
			$TABIX $DBFOLDER_DBSNP/$DATE/$ASSEMBLY/dbsnp.$DBSNP_VERSION.vcf.gz;
			rm -rf $DBFOLDER_DBSNP/$DATE/$ASSEMBLY/$DBSNP_VERSION/;
			mv $DBFOLDER_DBSNP/$DATE/$ASSEMBLY/$DBSNP_VERSION/dbsnp.parquet $DBFOLDER_DBSNP/$DATE/$ASSEMBLY/$DBSNP_VERSION/dbsnp.$DBSNP_VERSION.parquet;
			howard query --input=$DBFOLDER_DBSNP/$DATE/$ASSEMBLY/dbsnp.$DBSNP_VERSION.parquet --query='SELECT "#CHROM", POS, ID, REF, ALT, QUAL, FILTER, INFO FROM variants WHERE dbSNPBuildID<=$DBSNP_VERSION' --include_header --output=$DBFOLDER_DBSNP/$DATE/$ASSEMBLY/dbsnp.$DBSNP_VERSION.vcf.gz —memory=$MEMORY --threads=$THREADS;
			-[ ! -s $DBFOLDER_DBSNP/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_DBSNP/STARK.database && chmod o+r $DBFOLDER_DBSNP/STARK.database;
			[ ! -e $DBFOLDER_DBSNP/current/$ASSEMBLY ] || unlink $DBFOLDER_DBSNP/current/$ASSEMBLY;
			ln -snf $DBFOLDER_DBSNP/$DATE/$ASSEMBLY $DBFOLDER_DBSNP/current/$ASSEMBLY;
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_DBSNP/done"
	fi;
fi;

#########
# dbNSFP #
#########
DATABASE="dbNSFP"
DATABASE_NAME="dbNSFP"
DATABASE_FULLNAME="Non-synonymous single-nucleotide variants database"
DATABASE_WEBSITE="https://dbnsfp.s3.amazonaws.com"
DATABASE_DESCRIPTION="dbNSFP is a database developed for functional prediction and annotation of all potential non-synonymous single-nucleotide variants (nsSNVs) in the human genome."

if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then
	
	DBFOLDER_DBNSFP=$(dirname $DBNSFP_DATABASES)
	if [ ! -e $DBFOLDER_DBNSFP/current ]; then
		mkdir -p $DBFOLDER_DBNSFP/current;
	fi;

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_DBNSFP/current/$ASSEMBLY ] || (($UPDATE)); then
		
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for ' [$ASSEMBLY]"
		
		if (($UPDATE)); then
			if [ -e $DBFOLDER_DBNSFP/current/$ASSEMBLY ]; then mv -f $DBFOLDER_DBNSFP/current/$ASSEMBLY $DBFOLDER_DBNSFP.V$DATE; fi;
		fi;

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

		echo "$DBFOLDER_DBNSFP/done: $DBFOLDER
			howard databases --assembly='$ASSEMBLY' --genomes-folder=$DBFOLDER_GENOME/current/ --download-dbnsfp=$DBFOLDER_DBNSFP/$DATE --download-dbnsfp-vcf --download-dbnsfp-parquet —memory=$MEMORY --threads=$THREADS;
			-[ ! -s $DBFOLDER_DBNSFP/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_DBNSFP/STARK.database && chmod o+r $DBFOLDER_DBNSFP/STARK.database;
			[ ! -e $DBFOLDER_DBNSFP/current/$ASSEMBLY ] || unlink $DBFOLDER_DBNSFP/current/$ASSEMBLY;
			ln -snf $DBFOLDER_DBNSFP/$DATE/$ASSEMBLY $DBFOLDER_DBNSFP/current/$ASSEMBLY;
			rm -rf $DB_TMP;
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_DBNSFP/done"
	fi;
fi;

##########
# ARRIBA #
##########
DATABASE="arriba"
DATABASE_NAME="arriba"
DATABASE_FULLNAME="arriba"
DATABASE_WEBSITE="https://github.com/suhrig/arriba/"
DATABASE_DESCRIPTION="Arriba is a command-line tool for the detection of gene fusions from RNA-Seq data"

if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then
	
	DBFOLDER_ARRIBA=$(dirname $ARRIBA_DATABASES)
	if [ ! -e $DBFOLDER_ARRIBA/current ]; then
		mkdir -p $DBFOLDER_ARRIBA/current;
	fi;

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP


	if [ ! -e $DBFOLDER_ARRIBA/current/$ASSEMBLY ] || (($UPDATE)); then
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for ' [$ASSEMBLY]"

		if (($UPDATE)); then
			if [ -e $DBFOLDER_ARRIBA/current/$ASSEMBLY ]; then mv -f $DBFOLDER_ARRIBA/current/$ASSEMBLY $DBFOLDER_ARRIBA.V$DATE; fi;
		fi;
		
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

		DB_RELEASE_INFOS_JSON='
		{
			"release": "'$DATE'",
			"date": "'$DATE'",
			"files": [ "'$DBFOLDER_ARRIBA'" ],
			"assembly": [ "'$ASSEMBLY'" ],
			"download": {
				"methode": "'$DOWNLOAD_METHOD'",
				"URL": "'$(dirname $ARRIBA_CURRENT)'",
				"file": "'$(basename $ARRIBA_CURRENT)'",
				"date": "'2023-02-08'"
			}
		}
		';
		echo "$DB_RELEASE_INFOS_JSON" > $DB_TMP/STARK.database.release

		(($VERBOSE)) && echo "#[INFO] ARRIBA URL=$ARRIBA_CURRENT"
		(($VERBOSE)) && echo "#[INFO] ARRIBA RELEASE=$DBFOLDER_ARRIBA/$DATE/$ASSEMBLY"

		echo "$DBFOLDER_ARRIBA/done: $DBFOLDER
			aria2c -c -s 16 -x 16 -k 1M -j 1 $ARRIBA_CURRENT -d $DB_TMP;
			mkdir -p $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY;
			chmod 0775 $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY;
			tar -xzf  $DB_TMP/$(basename $ARRIBA_CURRENT) -C $DB_TMP --strip-components=1;
			mv $DB_TMP/database/blacklist_$ASSEMBLY* $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY;
			mv $DB_TMP/database/cytobands_$ASSEMBLY* $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY;
			mv $DB_TMP/database/known_fusions_$ASSEMBLY* $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY;
			mv $DB_TMP/database/protein_domains_$ASSEMBLY* $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY;
			-[ ! -s $DBFOLDER_ARRIBA/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_ARRIBA/STARK.database && chmod o+r $DBFOLDER_ARRIBA/STARK.database; 
			-[ ! -s $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY/STARK.database.release ] && cp $DB_TMP/STARK.database.release $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY/STARK.database.release && chmod o+r $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY/STARK.database.release;
			[ ! -e $DBFOLDER_ARRIBA/current/$ASSEMBLY ] || unlink $DBFOLDER_ARRIBA/current/$ASSEMBLY;
			ln -snf $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY $DBFOLDER_ARRIBA/current/$ASSEMBLY;
			rm -rf $DB_TMP;
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_ARRIBA/done"
	fi;
fi;

########
# CTAT #
########
DATABASE="ctat"
DATABASE_NAME="ctat"
DATABASE_FULLNAME=" CTAT Genome Lib"
DATABASE_WEBSITE="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/"
DATABASE_DESCRIPTION=" CTAT Genome Lib is a resource collection used by the Trinity Cancer Transcriptome Analysis Toolkit (CTAT). This CTAT-genome-lib-builder system is leveraged for preparing a target genome and annotation set for use with Trinity CTAT tools, including fusion transcript detection and cancer mutation discovery"

if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then

	DBFOLDER_CTAT=$(dirname $CTAT_DATABASES)
	if [ ! -e $DBFOLDER_CTAT/current ]; then
		mkdir -p $DBFOLDER_CTAT/current;
	fi;

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_CTAT/current/$ASSEMBLY ] || (($UPDATE)); then
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for [$ASSEMBLY]"

		CTAT_DATE=$(curl -s -I $CTAT_CURRENT | grep "Last-Modified: " | sed "s/Last-Modified: //g" | sed "s/\r$//g");
		CTAT_DATE_RELEASE=$(date -d "$CTAT_DATE");

		if (($UPDATE)); then
			if [ -e $DBFOLDER_CTAT/current/$ASSEMBLY ]; then mv -f $DBFOLDER_CTAT/current/$ASSEMBLY $DBFOLDER_CTAT.V$DATE; fi;
		fi;
		
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

		DB_RELEASE_INFOS_JSON='
		{
			"release": "'$CTAT_DATE'",
			"date": "'$CTAT_DATE_RELEASE'",
			"files": [ "'$CTAT_CURRENT'" ],
			"assembly": [ "'$ASSEMBLY'" ],
			"download": {
				"methode": "'$DOWNLOAD_METHOD'",
				"URL": "'$(dirname $CTAT_CURRENT)'",
				"file": "'$(basename $CTAT_CURRENT)'",
				"date": "'$CTAT_DATE_RELEASE'"
			}
		}
		';
		echo "$DB_RELEASE_INFOS_JSON" > $DB_TMP/STARK.database.release

		(($VERBOSE)) && echo "#[INFO] CTAT URL=$CTAT_CURRENT"
		(($VERBOSE)) && echo "#[INFO] CTAT RELEASE=$DATE"

		echo "$DBFOLDER_CTAT/done: $DBFOLDER
			aria2c -c -s 16 -x 16 -k 1M -j 1 $CTAT_CURRENT -d $DB_TMP;
			wget --progress=bar:force:noscroll $CTAT_PM -P $DB_TMP;
			mkdir -p $DBFOLDER_CTAT/$DATE/$ASSEMBLY;
			chmod 0775 $DBFOLDER_CTAT/$DATE/$ASSEMBLY;
			tar -xzf  $DB_TMP/$(basename $CTAT_CURRENT) -C  $DB_TMP --strip-components=1;
			$JAVA -jar $PICARD CreateSequenceDictionary -REFERENCE $DB_TMP/ctat_genome_lib_build_dir/ref_genome.fa -OUTPUT $DB_TMP/ctat_genome_lib_build_dir/ref_genome.dict;
			cp -R $DB_TMP/ctat_genome_lib_build_dir/* $DBFOLDER_CTAT/$DATE/$ASSEMBLY;
			\cp $DB_TMP/AnnotFilterRule.pm $DBFOLDER_CTAT/$DATE/$ASSEMBLY;
			-[ ! -s $DBFOLDER_CTAT/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_CTAT/STARK.database && chmod o+r $DBFOLDER_CTAT/STARK.database;
			-[ ! -s $DBFOLDER_CTAT/$DATE/$ASSEMBLY/STARK.database.release ] && cp $DB_TMP/STARK.database.release $DBFOLDER_CTAT/$DATE/$ASSEMBLY/STARK.database.release && chmod o+r $DBFOLDER_CTAT/$DATE/$ASSEMBLY/STARK.database.release;
			[ ! -e $DBFOLDER_CTAT/current/$ASSEMBLY ] || unlink $DBFOLDER_CTAT/current/$ASSEMBLY;
			ln -snf $DBFOLDER_CTAT/$DATE/$ASSEMBLY $DBFOLDER_CTAT/current/$ASSEMBLY;
			rm -rf $DB_TMP;
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_CTAT/done"
	fi;
fi;

###########
# GENCODE #
###########
DATABASE="gencode"
DATABASE_NAME="gencode"
DATABASE_FULLNAME="GENCODE"
DATABASE_WEBSITE="https://www.gencodegenes.org/"
DATABASE_DESCRIPTION=" The goal of the GENCODE project is to identify and classify all gene features in the human and mouse genomes with high accuracy based on biological evidence, and to release these annotations for the benefit of biomedical research and genome interpretation"

if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then

	DBFOLDER_GENCODE=$(dirname $GENCODE_DATABASES)
	if [ ! -e $DBFOLDER_GENCODE/current ]; then
		mkdir -p $DBFOLDER_GENCODE/current;
	fi;

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_GENCODE/current/$ASSEMBLY ] || (($UPDATE)); then
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for [$ASSEMBLY]"

		GENCODE_DATE=$(curl -s -I $GENCODE_CURRENT | grep "Last-Modified: " | sed "s/Last-Modified: //g" | sed "s/\r$//g");
		GENCODE_DATE_RELEASE=$(date -d "$GENCODE_DATE");

		if (($UPDATE)); then
			if [ -e $DBFOLDER_GENCODE/current/$ASSEMBLY ]; then mv -f $DBFOLDER_GENCODE/current/$ASSEMBLY $DBFOLDER_GENCODE.V$DATE; fi;
		fi;
		
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

		DB_RELEASE_INFOS_JSON='
		{
			"release": "'$GENCODE_DATE'",
			"date": "'$GENCODE_DATE_RELEASE'",
			"files": [ "'$(basename $GENCODE_CURRENT)'" ],
			"assembly": [ "'$ASSEMBLY'" ],
			"download": {
				"methode": "'$DOWNLOAD_METHOD'",
				"URL": "'$(dirname $GENCODE_CURRENT)'",
				"file": "'$(basename $GENCODE_CURRENT)'",
				"date": "'$GENCODE_DATE_RELEASE'"
			}
		}
		';
		echo "$DB_RELEASE_INFOS_JSON" > $DB_TMP/STARK.database.release

		(($VERBOSE)) && echo "#[INFO] GENCODE URL=$GENCODE_CURRENT"
		(($VERBOSE)) && echo "#[INFO] GENCODE RELEASE=$DATE"

		echo "$DBFOLDER_GENCODE/current/$ASSEMBLY/gencode.v$GENCODE_VERSION.annotation.gtf.gz: $DBFOLDER
			mkdir -p $DBFOLDER_GENCODE/$DATE/$ASSEMBLY;
			chmod 0775 $DBFOLDER_GENCODE/$DATE/$ASSEMBLY;
			aria2c -c -s 16 -x 16 -k 1M -j 1 $GENCODE_CURRENT -d $DB_TMP;
			cp $DB_TMP/$(basename $GENCODE_CURRENT) $DBFOLDER_GENCODE/$DATE/$ASSEMBLY/gencode.v$GENCODE_VERSION.annotation.gtf.gz;
			gzip -d $DBFOLDER_GENCODE/$DATE/$ASSEMBLY/gencode.v$GENCODE_VERSION.annotation.gtf.gz;
			-[ ! -s $DBFOLDER_GENCODE/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_GENCODE/STARK.database && chmod o+r $DBFOLDER_GENCODE/STARK.database;
			-[ ! -s $DBFOLDER_GENCODE/$DATE/$ASSEMBLY/STARK.database.release ] && cp $DB_TMP/STARK.database.release $DBFOLDER_GENCODE/$DATE/$ASSEMBLY/STARK.database.release && chmod o+r $DBFOLDER_GENCODE/$DATE/$ASSEMBLY/STARK.database.release;
			[ ! -e $DBFOLDER_GENCODE/current/$ASSEMBLY ] || unlink $DBFOLDER_GENCODE/current/$ASSEMBLY;
			ln -snf $DBFOLDER_GENCODE/$DATE/$ASSEMBLY $DBFOLDER_GENCODE/current/$ASSEMBLY;
			rm -rf $DB_TMP;
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_GENCODE/current/$ASSEMBLY/gencode.v$GENCODE_VERSION.annotation.gtf.gz"
	fi;
fi;

if [ ! -z "$MK_ALL" ]; then
	echo "$DBFOLDER:
		mkdir -p $DBFOLDER;
		chmod 0775 $DBFOLDER;
	" >> $MK

	echo "all: $MK_ALL
		echo '#[INFO] Build release: $DATE' >> $DATABASES/STARK.download.releases
	" >> $MK

	(($VERBOSE)) && echo "DATABASE INIT"

	if ((1)); then
		if (($BUILD)) || (($UPDATE)); then
			echo "#[INFO] DATABASES DOWNLOADING..."
			echo "#[INFO] THAT CAN TAKE SOME TIME..."
			if (($VERBOSE)) || (($DEBUG)); then
				if (($DEBUG)); then
					make -k -j $THREADS $MK_OPTION -f $MK all;
				elif (($VERBOSE)); then
					make -k -j $THREADS $MK_OPTION -f $MK all;
				fi;
			else
				make -k -j $THREADS $MK_OPTION -f $MK all 1>$MK_LOG 2>$MK_ERR;
				if (($(cat $MK_LOG $MK_ERR | grep "\*\*\*" -c))); then
					echo "#[ERROR] Databases download failed"
					exit 1
				fi;
			fi;
			echo "#[INFO] DATABASES DOWNLOADED"
		else
			echo "## use --build or --update to download or update databases"
		fi;
	fi;
else
	(($VERBOSE)) && echo ""
	(($VERBOSE)) && echo "#[INFO] Nothing to download"
fi;

if (($DEBUG)); then
	echo ""
	echo "### MAKEFILE"
	echo "#"
	echo "# MK=$MK"
	echo "# LOG=$MK_LOG"
	echo ""
	cat -n $MK
fi;

if ((0)); then
	echo ""
	(($VERBOSE)) && echo "#[INFO] RELEASES CHECK"
	(($DEBUG)) && ls -l $DATABASES

	if ((1)); then
		echo ""
		(($VERBOSE)) && echo "#[INFO] DOWNLOAD GENOME [$ASSEMBLY]"
		(($VERBOSE)) && echo "#"
	fi;
fi;

if ((0)); then
	rm -Rf $TMP_DATABASES_DOWNLOAD_FOLDER
	rm -Rf $TMP_DATABASES_DOWNLOAD_RAM
fi;

exit 0;