#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARKDatabases"
SCRIPT_DESCRIPTION="STARK download and build databases"
SCRIPT_RELEASE="1.2.0"
SCRIPT_DATE="27/03/2026"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

# Release note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-11/12/2018: Script creation\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.1b-12/12/2018: Change to Makefile\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.2b-21/12/2018: Add update, build, rebuild and threads options. Change dbsnp source\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.3b-31/05/2019: Add APP configuration\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.4b-17/06/2020: Clarify code, organisation DB/RELEASE, add STARK.description\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.5.0-11/04/2021: Change snpEff download, some bugs fixed, option --current, remove option --rebuild\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.6.0-28/07/2022: Change snpEff download, add GATK databases, some bugs fixed\n";
RELEASE_NOTES=$RELEASE_NOTES"# 1.0.0-01/11/2023: Rewrite for STARK 19 new database structure: clean code, HOWARD database python, fix makefile rules (point to file, not directory), add CTAT/arriba, use aria2c\n";
RELEASE_NOTES=$RELEASE_NOTES"# 1.1.0-08/04/2025: Update for STARK 19: clear code, fixes and improves\n";
RELEASE_NOTES=$RELEASE_NOTES"# 1.1.1-17/02/2026: Update GATK4 resources, add assembly option\n";
RELEASE_NOTES=$RELEASE_NOTES"# 1.2.0-27/03/2026: Fix bugs, genome is mandatory for any download\n";


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
	echo "# --assembly=<ASSEMBLY>                    Assembly to use (e.g. 'hg19', 'hg38') (default APP configuration).";
	echo "# --databases=<FOLDER>                     Databases folder (replace APP parameter)";
	echo "#                                          Will generate STARK databases folder structure";
	echo "# --databases_list=<STRING>                List of Databases to consider";
	echo "#                                          Format: 'database1,database2,...'";
	echo "#                                          Default: 'ALL' for all available databases";
	echo "#                                          Available databases: 'dbsnp'";
	echo "# --release                                Make new databases as a specific release (default 'current').";
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
ARGS=$(getopt -o "e:cbut:vdnh" --long "env:,app:,application:,assembly:,databases:,databases_list:,current,release:,assembly:,build,update,threads:,verbose,debug,release,help" -- "$@" 2> /dev/null)
PARAM=$@

eval set -- "$ARGS"
while true
do
	case "$1" in
		-e|--env|--app|--application)
			APP="$2"
			shift 2
			;;
		--assembly)
			ASSEMBLY="$2"
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
		--release)
			RELEASE="$2"
			shift 2
			;;
		# --assembly)
		# 	ASSEMBLY_INPUT="$2"
		# 	shift 2
		# 	;;
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
(($DEBUG)) && echo "#[DEBUG] SEARCHING SCRIPTS"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
(($DEBUG)) && echo "#[DEBUG] DONE"

# Configuration
(($DEBUG)) && echo "#[DEBUG] SEARCHING CONFIG APPS"
ENV_CONFIG=$(find -L $SCRIPT_DIR/.. -name config.app)
(($DEBUG)) && echo "#[DEBUG] DONE"
(($DEBUG)) && echo "#[DEBUG] SOURCE CONFIGS"
#echo $ENV_CONFIG
source $ENV_CONFIG 
(($DEBUG)) && echo "#[DEBUG] DONE"

# Memory
if [ "$MEMORY" == "" ]; then
	MEMORY=1
else
	MEMORY=$(echo $MEMORY | sed 's/[^0-9]//g')
fi;
MEMORYG=$MEMORY"G" 


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

# ARIA
ARIA_CMD="aria2c -c -s $THREADS -x 16 -k 1M --async-dns=false -j $THREADS"

TMP_DATABASES_DOWNLOAD_FOLDER=$TMP_FOLDER_TMP/$RANDOM$RANDOM
mkdir -p $TMP_DATABASES_DOWNLOAD_FOLDER
TMP_DATABASES_DOWNLOAD_RAM="$(mktemp -d -p /dev/shm/)"

if [ "$TMP_DATABASES_DOWNLOAD_RAM" == "" ]; then
	TMP_DATABASES_DOWNLOAD_RAM=$TMP_DATABASES_DOWNLOAD_FOLDER;
fi;

# DATABASES build release
DATE=$(date '+%Y%m%d-%H%M%S')

# DATABASES release
if [ "$RELEASE"	== "" ]; then
	RELEASE=current
fi;

# if [ "$ASSEMBLY_INPUT" == "" ]; then
# 	ASSEMBLY=$ASSEMBLY # From APP
# else
# 	ASSEMBLY=$ASSEMBLY_INPUT
# fi;

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

# GENOME for all other databases
GENOME=$DATABASES/genomes/$RELEASE/$ASSEMBLY/$ASSEMBLY.fa

if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT || [ ! -e $GENOME ]; then

	DBFOLDER_GENOME=$DATABASES/genomes
	if [ ! -e $DBFOLDER_GENOME/$RELEASE ]; then
		mkdir -p $DBFOLDER_GENOME/$RELEASE;
	fi;
	
	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_GENOME/$RELEASE/$ASSEMBLY ] || (($UPDATE)); then
		
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for ' [$ASSEMBLY]"
		
		if (($UPDATE)); then
			if [ -e $DBFOLDER_GENOME/$RELEASE/$ASSEMBLY ]; then mv -f $DBFOLDER_GENOME/$RELEASE/$ASSEMBLY $DBFOLDER_GENOME.V$DATE; fi;
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
			$HOWARD databases --assembly='$ASSEMBLY' --download-genomes=$DBFOLDER_GENOME/$DATE --download-genomes-contig-regex=$GENOME_REGEX;
			[ ! -e $DBFOLDER_GENOME/$RELEASE/$ASSEMBLY ] || unlink $DBFOLDER_GENOME/$RELEASE/$ASSEMBLY;
			ln -snf ../$DATE/$ASSEMBLY $DBFOLDER_GENOME/$RELEASE/$ASSEMBLY;
			-[ ! -s $DBFOLDER_GENOME/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_GENOME/STARK.database && chmod o+r $DBFOLDER_GENOME/STARK.database;
			#mv $GENOME $GENOME.tmp;
			#mv $GENOME.fai $GENOME.fai.tmp;
			#cut -f1 $GENOME.fai.tmp|sort -k1V|parallel -k '$SAMTOOLS faidx $GENOME.tmp {}' > $GENOME;
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
# Source: ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle
if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then

	DBFOLDER_GATK_URL_FULL_COPY=0
	GATK_RESOURCE_NB=0
	MK_DBFOLDER_GATK_ALL=""
	> $MK.existing_gatk_db

	if [ ! -e $DBFOLDER_GATK/$RELEASE ]; then
		mkdir -p $DBFOLDER_GATK/$RELEASE;
	fi;

	for GATK_RESOURCE in $GATK_DATABASES_LIST; do
		DB_TARGET_GATK=$DBFOLDER_GATK/$RELEASE/$ASSEMBLY/$(echo $GATK_RESOURCE | cut -d: -f2);		# /STARK/databases/gatk/current/hg19/1000G_omni2.5.b37.vcf.gz
		DB_RELEASE_FILE=$(basename $DB_TARGET_GATK);												# 1000G_omni2.5.b37.vcf.gz
		DB_RELEASE_FILE_PATH="$DBFOLDER_GATK/$DATE/$ASSEMBLY/$DB_RELEASE_FILE";						# /STARK/databases/gatk/DATE/hg19/1000G_omni2.5.b37.vcf.gz

		DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
		mkdir -p $DB_TMP

		if (($UPDATE)); then
			if [ -e $DBFOLDER_GATK/$RELEASE/$ASSEMBLY ]; then mv -f $DBFOLDER_GATK/$RELEASE/$ASSEMBLY $DBFOLDER_GATK.V$DATE; fi;
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

			# URL for download
			DBFOLDER_GATK_URL=$DBFOLDER_GATK_URL_DEFAULT/$ASSEMBLY

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
		echo "$DBFOLDER_GATK/$RELEASE/$ASSEMBLY/done: $MK_DBFOLDER_GATK_ALL $GENOME
			mkdir -p $DBFOLDER_GATK/$DATE/$ASSEMBLY/original
			chmod 0775 $DBFOLDER_GATK/$DATE/$ASSEMBLY -R
			-[ ! -s $DBFOLDER_GATK/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_GATK/STARK.database
			cp $DB_TMP/STARK.database.release $DBFOLDER_GATK/$DATE/$ASSEMBLY/
			chmod o+r $DBFOLDER_GATK/STARK.database $DBFOLDER_GATK/$DATE/$ASSEMBLY/STARK.database.release
			[ ! -e $DBFOLDER_GATK/$RELEASE/$ASSEMBLY ] || unlink $DBFOLDER_GATK/$RELEASE/$ASSEMBLY
			#ln -snf $DBFOLDER_GATK/$DATE/$ASSEMBLY $DBFOLDER_GATK/$RELEASE/$ASSEMBLY
			ln -snf ../$DATE/$ASSEMBLY $DBFOLDER_GATK/$RELEASE/$ASSEMBLY
			rm -rf $DB_TMP;
			touch $DBFOLDER_GATK/$RELEASE/$ASSEMBLY/done
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_GATK/$RELEASE/$ASSEMBLY/done" 
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
	if [ ! -e $DBFOLDER_SNPEFF/$RELEASE ]; then
		mkdir -p $DBFOLDER_SNPEFF/$RELEASE;
	fi;

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_SNPEFF/$RELEASE/$ASSEMBLY ] || (($UPDATE)); then
		
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for ' [$ASSEMBLY]"
		
		if (($UPDATE)); then
			if [ -e $DBFOLDER_SNPEFF/$RELEASE/$ASSEMBLY ]; then mv -f $DBFOLDER_SNPEFF/$RELEASE/$ASSEMBLY $DBFOLDER_SNPEFF.V$DATE; fi;
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

		echo "$DBFOLDER_SNPEFF/done: $DBFOLDER $GENOME
			$HOWARD databases --assembly='$ASSEMBLY' --download-snpeff=$DBFOLDER_SNPEFF/$DATE --config=$HOWARD_CONFIG
			-[ ! -s $DBFOLDER_SNPEFF/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_SNPEFF/STARK.database && chmod o+r $DBFOLDER_SNPEFF/STARK.database 
			[ ! -e $DBFOLDER_SNPEFF/$RELEASE/$ASSEMBLY ] || unlink $DBFOLDER_SNPEFF/$RELEASE/$ASSEMBLY
			#ln -snf $DBFOLDER_SNPEFF/$DATE/$ASSEMBLY $DBFOLDER_SNPEFF/$RELEASE/$ASSEMBLY
			ln -snf ../$DATE/$ASSEMBLY $DBFOLDER_SNPEFF/$RELEASE/$ASSEMBLY
			rm -rf $DB_TMP;
			touch $DBFOLDER_SNPEFF/done
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_SNPEFF/done"
		(($VERBOSE)) && cat $MK
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
	if [ ! -e $DBFOLDER_ANNOVAR/$RELEASE ]; then
		mkdir -p $DBFOLDER_ANNOVAR/$RELEASE;
	fi;

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_ANNOVAR/$RELEASE/$ASSEMBLY ] || (($UPDATE)); then
		
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for ' [$ASSEMBLY]"

		if (($UPDATE)); then
			if [ -e $DBFOLDER_ANNOVAR/$RELEASE/$ASSEMBLY ]; then mv -f $DBFOLDER_ANNOVAR/$RELEASE/$ASSEMBLY $DBFOLDER_ANNOVAR.V$DATE; fi;
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
	
		echo "$DBFOLDER_ANNOVAR/done: $DBFOLDER $GENOME
			$HOWARD databases --assembly='$ASSEMBLY' --download-annovar=$DBFOLDER_ANNOVAR/$DATE --download-annovar-files='$ANNOVAR_FILES'
			-[ ! -s $DBFOLDER_ANNOVAR/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_ANNOVAR/STARK.database && chmod o+r $DBFOLDER_ANNOVAR/STARK.database 
			[ ! -e $DBFOLDER_ANNOVAR/$RELEASE/$ASSEMBLY ] || unlink $DBFOLDER_ANNOVAR/$RELEASE/$ASSEMBLY
			ln -snf ../$DATE/$ASSEMBLY $DBFOLDER_ANNOVAR/$RELEASE/$ASSEMBLY
			rm -rf $DB_TMP;
			touch $DBFOLDER_ANNOVAR/done
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
	
	#DBFOLDER_REFGENE=$DBFOLDER/refGene
	DBFOLDER_REFGENE=$(dirname $REFGENE_DATABASES)
	if [ ! -e $DBFOLDER_REFGENE/$RELEASE ]; then
		mkdir -p $DBFOLDER_REFGENE/$RELEASE;
	fi;

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_REFGENE/$RELEASE/$ASSEMBLY ] || (($UPDATE)); then
		
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for ' [$ASSEMBLY]"
		
		if (($UPDATE)); then
			if [ -e $DBFOLDER_REFGENE/$RELEASE/$ASSEMBLY ]; then mv -f $DBFOLDER_REFGENE/$RELEASE/$ASSEMBLY $DBFOLDER_REFGENE.V$DATE; fi;
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
			"files": [ "ncbiRefSeq.*" ],
			"assembly": [ "'$ASSEMBLY'" ],
			"download": {
				"methode": "HOWARD and STARK scripts",
				"URL": "http://hgdownload.soe.ucsc.edu/goldenPath",
				"file": "ncbiRefSeq.*",
				"date": "'$DATE'"
			}
		}
		';
		echo "$DB_RELEASE_INFOS_JSON" > $DB_TMP/STARK.database.release

		echo "$DBFOLDER_REFGENE/$RELEASE/$ASSEMBLY: $DBFOLDER $GENOME
			# Download ncbiRefSeq TXT and BED
			$HOWARD databases --assembly='$ASSEMBLY' --download-refseq=$DBFOLDER_REFGENE/$DATE --download-refseq-format-file='ncbiRefSeq.txt' ;
			# Convert into GTF
			awk -F '\t' -v OFS='\t' -f $STARK_FOLDER_BIN/refSeq_to_gtf.awk $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/ncbiRefSeq.txt > $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/ncbiRefSeq.gtf
			# Download refSeq in GTF
			$PYTHON $STARK_FOLDER_BIN/get_refGene.py --assembly $ASSEMBLY --patch latest --output $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/refSeq.gtf
			# Collapse GTF
			$PYTHON $STARK_FOLDER_BIN/collapse_annotation.py $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/refSeq.gtf $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/refSeq.collapsed.gtf;
			# README.md
			echo '# refSeq files' > $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '## ncbiRefSeq.txt' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '- Source: NCBI RefSeq database downloaded by HOWARD' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '- URL: https://www.ncbi.nlm.nih.gov/refseq/' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '- Downloaded file: ncbiRefSeq.txt' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '- Download date: '$DATE >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '## ncbiRefSeq.bed' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '- Source: NCBI RefSeq database converted by HOWARD' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '## ncbiRefSeq.gtf' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '- Source: NCBI RefSeq database converted by STARK refSeq_to_gtf.awk script' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '## refSeq.gtf' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '- Source: NCBI RefSeq database in GTF downloaded by get_refGene.py' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '- URL: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '- Release: latest' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '- Download date: '$DATE >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '## refSeq.collapsed.gtf' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			echo '- Source: NCBI RefSeq database in GTF converted by collapse_annotation.py' >> $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/README.md;
			# Links
			-[ ! -s $DBFOLDER_REFGENE/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_REFGENE/STARK.database && chmod o+r $DBFOLDER_REFGENE/STARK.database;
			-[ ! -s $DBFOLDER_REFGENE/$RELEASE/$ASSEMBLY/STARK.database.release ] && cp $DB_TMP/STARK.database.release $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/STARK.database.release && chmod o+r $DBFOLDER_REFGENE/$DATE/$ASSEMBLY/STARK.database.release;
			[ ! -e $DBFOLDER_REFGENE/$RELEASE/$ASSEMBLY ] || unlink $DBFOLDER_REFGENE/$RELEASE/$ASSEMBLY;
			ln -snf ../$DATE/$ASSEMBLY $DBFOLDER_REFGENE/$RELEASE/$ASSEMBLY;
			touch $DBFOLDER_REFGENE/$RELEASE/$ASSEMBLY
			" >> $MK

		MK_ALL="$MK_ALL $DBFOLDER_REFGENE/$RELEASE/$ASSEMBLY"
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
	if [ ! -e $DBFOLDER_GENCODE/$RELEASE ]; then
		mkdir -p $DBFOLDER_GENCODE/$RELEASE;
	fi;

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$GENCODE_VERSION
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_GENCODE/$RELEASE/$ASSEMBLY ] || (($UPDATE)); then
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$GENCODE_VERSION' for [$ASSEMBLY]"

		GENCODE_DATE=$(curl -s -I $GENCODE_CURRENT | grep "Last-Modified: " | sed "s/Last-Modified: //g" | sed "s/\r$//g");
		GENCODE_DATE_RELEASE=$(date -d "$GENCODE_DATE");

		if (($UPDATE)); then
			if [ -e $DBFOLDER_GENCODE/$RELEASE/$ASSEMBLY ]; then mv -f $DBFOLDER_GENCODE/$RELEASE/$ASSEMBLY $DBFOLDER_GENCODE.V$DATE; fi;
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
			"release": "'$GENCODE_VERSION'",
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
		(($VERBOSE)) && echo "#[INFO] GENCODE RELEASE=$GENCODE_VERSION"

		echo "$GENCODE_DATABASES/$ASSEMBLY/done: $DBFOLDER $GENOME
			mkdir -p $DBFOLDER_GENCODE/$GENCODE_VERSION/$ASSEMBLY;
			chmod 0775 $DBFOLDER_GENCODE/$GENCODE_VERSION/$ASSEMBLY;
			$ARIA_CMD $GENCODE_CURRENT -d $DB_TMP;
			cp $DB_TMP/$(basename $GENCODE_CURRENT) $DBFOLDER_GENCODE/$GENCODE_VERSION/$ASSEMBLY/gencode.gtf.gz;
			gzip -d $DBFOLDER_GENCODE/$GENCODE_VERSION/$ASSEMBLY/gencode.gtf.gz;
			$PYTHON $STARK_FOLDER_BIN/collapse_annotation.py $DBFOLDER_GENCODE/$GENCODE_VERSION/$ASSEMBLY/gencode.gtf $DBFOLDER_GENCODE/$GENCODE_VERSION/$ASSEMBLY/gencode.collapsed.gtf;
			-[ ! -s $DBFOLDER_GENCODE/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_GENCODE/STARK.database && chmod o+r $DBFOLDER_GENCODE/STARK.database;
			-[ ! -s $DBFOLDER_GENCODE/$GENCODE_VERSION/$ASSEMBLY/STARK.database.release ] && cp $DB_TMP/STARK.database.release $DBFOLDER_GENCODE/$GENCODE_VERSION/$ASSEMBLY/STARK.database.release && chmod o+r $DBFOLDER_GENCODE/$GENCODE_VERSION/$ASSEMBLY/STARK.database.release;
			[ ! -e $DBFOLDER_GENCODE/$RELEASE/$ASSEMBLY ] || unlink $DBFOLDER_GENCODE/$RELEASE/$ASSEMBLY;
			ln -snf ../$GENCODE_VERSION/$ASSEMBLY $DBFOLDER_GENCODE/$RELEASE/$ASSEMBLY;
			rm -rf $DB_TMP;
			touch $GENCODE_DATABASES/$ASSEMBLY/done;
		" >> $MK
		MK_ALL="$MK_ALL $GENCODE_DATABASES/$ASSEMBLY/done"
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
	if [ ! -e $DBFOLDER_DBSNP/$RELEASE ]; then
		mkdir -p $DBFOLDER_DBSNP/$RELEASE;
	fi;

	if [ "$DBFOLDER_GENOME" == "" ]; then
		DBFOLDER_GENOME=$DATABASES/genomes
	fi;

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_DBSNP/$RELEASE/$ASSEMBLY ] || (($UPDATE)); then
		
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for ' [$ASSEMBLY]"

		if (($UPDATE)); then
			if [ -e $DBFOLDER_DBSNP/$RELEASE/$ASSEMBLY ]; then mv -f $DBFOLDER_DBSNP/$RELEASE/$ASSEMBLY $DBFOLDER_DBSNP.V$DATE; fi;
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
		
		echo "$DBFOLDER_DBSNP/done: $DBFOLDER $GENOME
			if [ ! -e $DBFOLDER_DBSNP/$DATE/$ASSEMBLY/$DBSNP_VERSION_DOWNLOAD/dbsnp.parquet ]; then \
				$HOWARD databases --assembly=$ASSEMBLY --genomes-folder=$DBFOLDER_GENOME/$RELEASE --download-dbsnp=$DBFOLDER_DBSNP/$DATE --download-dbsnp-releases=$DBSNP_VERSION_DOWNLOAD --download-dbsnp-parquet --threads=$THREADS; \
			fi;
			if [ ! -e $DBFOLDER_DBSNP/$DATE/$ASSEMBLY/$DBSNP_VERSION_DOWNLOAD/dbsnp.vcf.gz ]; then \
				howard convert --input=$DBFOLDER_DBSNP/$DATE/$ASSEMBLY/$DBSNP_VERSION_DOWNLOAD/dbsnp.parquet --output=$DBFOLDER_DBSNP/$DATE/$ASSEMBLY/$DBSNP_VERSION_DOWNLOAD/dbsnp.vcf.gz --threads=$THREADS --access=RO; \
			fi;
			howard query --input=$DBFOLDER_DBSNP/$DATE/$ASSEMBLY/$DBSNP_VERSION_DOWNLOAD/dbsnp.parquet --query=\"SELECT \\\"#CHROM\\\", POS, ID, REF, ALT, '.' AS QUAL, '.' AS FILTER, INFO FROM variants WHERE COMMON\" --output=$DBFOLDER_DBSNP/$DATE/$ASSEMBLY/$DBSNP_VERSION_DOWNLOAD/dbsnp.COMMON.vcf.gz --memory=$MEMORYG --threads=$THREADS --access=RO
			$TABIX $DBFOLDER_DBSNP/$DATE/$ASSEMBLY/$DBSNP_VERSION_DOWNLOAD/dbsnp.COMMON.vcf.gz
			howard query --input=$DBFOLDER_DBSNP/$DATE/$ASSEMBLY/$DBSNP_VERSION_DOWNLOAD/dbsnp.parquet --query=\"SELECT \\\"#CHROM\\\", POS, ID, REF, ALT, QUAL, FILTER, INFO FROM variants WHERE dbSNPBuildID<=$DBSNP_BUILDID\" --output=$DBFOLDER_DBSNP/$DATE/$ASSEMBLY/$DBSNP_VERSION_DOWNLOAD/dbsnp.b$DBSNP_BUILDID.vcf.gz --memory=$MEMORYG --threads=$THREADS --access=RO
			$TABIX -f $DBFOLDER_DBSNP/$DATE/$ASSEMBLY/$DBSNP_VERSION_DOWNLOAD/dbsnp.COMMON.vcf.gz
			$TABIX -f $DBFOLDER_DBSNP/$DATE/$ASSEMBLY/$DBSNP_VERSION_DOWNLOAD/dbsnp.b$DBSNP_BUILDID.vcf.gz
			-[ ! -s $DBFOLDER_DBSNP/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_DBSNP/STARK.database && chmod o+r $DBFOLDER_DBSNP/STARK.database;
			[ ! -e $DBFOLDER_DBSNP/$RELEASE/$ASSEMBLY ] || unlink $DBFOLDER_DBSNP/$RELEASE/$ASSEMBLY;
			ln -snf ../$DATE/$ASSEMBLY $DBFOLDER_DBSNP/$RELEASE/$ASSEMBLY;
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_DBSNP/done"
		(($VERBOSE)) && cat $MK
	fi;
fi;

# #########
# # dbNSFP #
# #########
# DATABASE="dbNSFP"
# DATABASE_NAME="dbNSFP"
# DATABASE_FULLNAME="Non-synonymous single-nucleotide variants database"
# DATABASE_WEBSITE="https://dbnsfp.s3.amazonaws.com"
# DATABASE_DESCRIPTION="dbNSFP is a database developed for functional prediction and annotation of all potential non-synonymous single-nucleotide variants (nsSNVs) in the human genome."

# if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then
	
# 	DBFOLDER_DBNSFP=$(dirname $DBNSFP_DATABASES)
# 	if [ ! -e $DBFOLDER_DBNSFP/current ]; then
# 		mkdir -p $DBFOLDER_DBNSFP/current;
# 	fi;

# 	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
# 	mkdir -p $DB_TMP
# 	chmod 0775 $DB_TMP;

# 	if [ ! -e $DBFOLDER_DBNSFP/current/$ASSEMBLY ] || (($UPDATE)); then
		
# 		(($VERBOSE)) && echo ""
# 		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for ' [$ASSEMBLY]"
		
# 		if (($UPDATE)); then
# 			if [ -e $DBFOLDER_DBNSFP/current/$ASSEMBLY ]; then mv -f $DBFOLDER_DBNSFP/current/$ASSEMBLY $DBFOLDER_DBNSFP.V$DATE; fi;
# 		fi;

# 		DB_INFOS_JSON='
# 		{
# 			"code": "'$DATABASE'",
# 			"name": "'$DATABASE_NAME'",
# 			"fullname": "'$DATABASE_FULLNAME'",
# 			"website": "'$DATABASE_WEBSITE'",
# 			"description": "'$DATABASE_DESCRIPTION'"
# 		}
# 		';
# 		echo "$DB_INFOS_JSON" > $DB_TMP/STARK.database

# 		echo "$DBFOLDER_DBNSFP/done: $DBFOLDER
# 			$HOWARD databases --assembly='$ASSEMBLY' --genomes-folder=$DBFOLDER_GENOME/current/ --download-dbnsfp=$DBFOLDER_DBNSFP/$DATE --download-dbnsfp-vcf --download-dbnsfp-parquet --memory=$MEMORYG --threads=$THREADS;
# 			-[ ! -s $DBFOLDER_DBNSFP/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_DBNSFP/STARK.database && chmod o+r $DBFOLDER_DBNSFP/STARK.database;
# 			[ ! -e $DBFOLDER_DBNSFP/current/$ASSEMBLY ] || unlink $DBFOLDER_DBNSFP/current/$ASSEMBLY;
# 			ln -snf $DBFOLDER_DBNSFP/$DATE/$ASSEMBLY $DBFOLDER_DBNSFP/current/$ASSEMBLY;
# 			rm -rf $DB_TMP;
# 		" >> $MK
# 		MK_ALL="$MK_ALL $DBFOLDER_DBNSFP/done"
# 	fi;
# fi;

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
	if [ ! -e $DBFOLDER_ARRIBA/$RELEASE ]; then
		mkdir -p $DBFOLDER_ARRIBA/$RELEASE;
	fi;

	# Arriba database URL and release
	ARRIBA_CURRENT=$ARRIBA_URL/$ARRIBA_RELEASE"/arriba_$ARRIBA_RELEASE.tar.gz";

	# Arriba database release date
	DATE_RELEASE=$ARRIBA_RELEASE

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE_RELEASE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP


	if [ ! -e $DBFOLDER_ARRIBA/$RELEASE/$ASSEMBLY ] || (($UPDATE)); then
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for ' [$ASSEMBLY]"

		if (($UPDATE)); then
			if [ -e $DBFOLDER_ARRIBA/$RELEASE/$ASSEMBLY ]; then mv -f $DBFOLDER_ARRIBA/$RELEASE/$ASSEMBLY $DBFOLDER_ARRIBA.V$DATE; fi;
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
			"release": "'$DATE_RELEASE'",
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

		echo "$DBFOLDER_ARRIBA/$RELEASE/$ASSEMBLY: $DBFOLDER $GENOME
			$ARIA_CMD $ARRIBA_CURRENT -d $DB_TMP;
			mkdir -p $DBFOLDER_ARRIBA/$DATE_RELEASE/$ASSEMBLY;
			chmod 0775 $DBFOLDER_ARRIBA/$DATE_RELEASE/$ASSEMBLY;
			tar -xzf  $DB_TMP/$(basename $ARRIBA_CURRENT) -C $DB_TMP --strip-components=1;
			mv $DB_TMP/database/blacklist_$ASSEMBLY* $DBFOLDER_ARRIBA/$DATE_RELEASE/$ASSEMBLY;
			mv $DB_TMP/database/cytobands_$ASSEMBLY* $DBFOLDER_ARRIBA/$DATE_RELEASE/$ASSEMBLY;
			mv $DB_TMP/database/known_fusions_$ASSEMBLY* $DBFOLDER_ARRIBA/$DATE_RELEASE/$ASSEMBLY;
			mv $DB_TMP/database/protein_domains_$ASSEMBLY* $DBFOLDER_ARRIBA/$DATE_RELEASE/$ASSEMBLY;
			-[ ! -s $DBFOLDER_ARRIBA/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_ARRIBA/STARK.database && chmod o+r $DBFOLDER_ARRIBA/STARK.database; 
			-[ ! -s $DBFOLDER_ARRIBA/$DATE_RELEASE/$ASSEMBLY/STARK.database.release ] && cp $DB_TMP/STARK.database.release $DBFOLDER_ARRIBA/$DATE_RELEASE/$ASSEMBLY/STARK.database.release && chmod o+r $DBFOLDER_ARRIBA/$DATE_RELEASE/$ASSEMBLY/STARK.database.release;
			[ ! -e $DBFOLDER_ARRIBA/$RELEASE/$ASSEMBLY ] || unlink $DBFOLDER_ARRIBA/$RELEASE/$ASSEMBLY;
			ln -snf ../$DATE_RELEASE/$ASSEMBLY $DBFOLDER_ARRIBA/$RELEASE/$ASSEMBLY;
			rm -rf $DB_TMP;
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_ARRIBA/$RELEASE/$ASSEMBLY"
	fi;
fi;


########
# PFAM #
########
DATABASE="pfam"
DATABASE_NAME="pfam"
DATABASE_FULLNAME="PFAM Database of protein families"
DATABASE_WEBSITE="https://www.pfam.org/"
DATABASE_DESCRIPTION="PFAM is a database of protein families, providing a comprehensive collection of protein domains and families."

if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then

	#[ -z "$PFAM_DATABASES" ] && PFAM_DATABASES=$DATABASES/pfam/current

	DBFOLDER_PFAM=$(dirname $PFAM_DATABASES)
	if [ ! -e $DBFOLDER_PFAM/$RELEASE ]; then
		mkdir -p $DBFOLDER_PFAM/$RELEASE;
	fi;
	(($DEBUG)) && echo "[DEBUG] DBFOLDER_PFAM=$DBFOLDER_PFAM"
	(($DEBUG)) && echo "[DEBUG] RELEASE=$RELEASE"

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$PFAM_RELEASE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP



	if [ ! -e $DBFOLDER_PFAM/$RELEASE/$ASSEMBLY ] || (($UPDATE)); then
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$PFAM_RELEASE' for [$ASSEMBLY]"

		PFAM_DATE=$(curl -s -I $PFAM_URL/Pfam$PFAM_RELEASE/Pfam-A.hmm.gz | grep "Last-Modified: " | sed "s/Last-Modified: //g" | sed "s/\r$//g");

		if (($UPDATE)); then
			if [ -e $DBFOLDER_PFAM/$RELEASE/$ASSEMBLY ]; then mv -f $DBFOLDER_PFAM/$RELEASE/$ASSEMBLY $DBFOLDER_PFAM.V$DATE; fi;
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
			"release": "'$PFAM_RELEASE'",
			"date": "'$DATE'",
			"files": [ "'$DBFOLDER_PFAM'" ],
			"assembly": [ "'$ASSEMBLY'" ],
			"download": {
				"methode": "'$DOWNLOAD_METHOD'",
				"URL": "'$PFAM_URL/Pfam$PFAM_RELEASE'",
				"file": "'Pfam-A.hmm.gz'",
				"date": "'$PFAM_DATE'"
			}
		}
		';
		echo "$DB_RELEASE_INFOS_JSON" > $DB_TMP/STARK.database.release

		echo "$PFAM_DATABASES/$ASSEMBLY: $DBFOLDER $GENOME
			mkdir -p $DBFOLDER_PFAM
			$ARIA_CMD $PFAM_URL/Pfam${PFAM_RELEASE}/Pfam-A.hmm.gz -d $DB_TMP;
			mkdir -p $DBFOLDER_PFAM/$PFAM_RELEASE/$ASSEMBLY;
			chmod 0775 $DBFOLDER_PFAM/$PFAM_RELEASE/$ASSEMBLY;
			$GZ -d $DB_TMP/Pfam-A.hmm.gz
			mv $DB_TMP/Pfam-A.hmm* $DBFOLDER_PFAM/$PFAM_RELEASE/$ASSEMBLY;
			# Pre index the Pfam-A.hmm file with hmmpress to speed up the ctat-genome-lib-builder step
			$DOCKER_RUN --rm --name ctat_prep_genome_lib-$ASSEMBLY-PFAM-"$(date +%Y%m%d-%H%M%S)" trinityctat/starfusion hmmpress $DBFOLDER_PFAM/$PFAM_RELEASE/$ASSEMBLY/Pfam-A.hmm;
			-[ ! -s $DBFOLDER_PFAM/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_PFAM/STARK.database && chmod o+r $DBFOLDER_PFAM/STARK.database; 
			-[ ! -s $DBFOLDER_PFAM/$PFAM_RELEASE/$ASSEMBLY/STARK.database.release ] && cp $DB_TMP/STARK.database.release $DBFOLDER_PFAM/$PFAM_RELEASE/$ASSEMBLY/STARK.database.release && chmod o+r $DBFOLDER_PFAM/$PFAM_RELEASE/$ASSEMBLY/STARK.database.release;
			[ ! -e $DBFOLDER_PFAM/$RELEASE/$ASSEMBLY ] || unlink $DBFOLDER_PFAM/$RELEASE/$ASSEMBLY;
			ln -snf ../$PFAM_RELEASE/$ASSEMBLY $DBFOLDER_PFAM/$RELEASE/$ASSEMBLY;
			rm -rf $DB_TMP;
			touch $PFAM_DATABASES/$ASSEMBLY
		" >> $MK

		MK_ALL="$MK_ALL $PFAM_DATABASES/$ASSEMBLY"
	fi;
fi;


########
# DFAM #
########
DATABASE="dfam"
DATABASE_NAME="dfam"
DATABASE_FULLNAME="DFAM Database of transposable element families"
DATABASE_WEBSITE="https://www.dfam.org/"
DATABASE_DESCRIPTION="DFAM is a database of transposable element families, providing a comprehensive collection of transposable element sequences and annotations."

if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then

	#[ -z "$DFAM_DATABASES" ] && DFAM_DATABASES=$DATABASES/dfam/current

	DBFOLDER_DFAM=$(dirname $DFAM_DATABASES)
	if [ ! -e $DBFOLDER_DFAM/$RELEASE ]; then
		mkdir -p $DBFOLDER_DFAM/$RELEASE;
	fi;
	(($DEBUG)) && echo "[DEBUG] DBFOLDER_DFAM=$DBFOLDER_DFAM"
	(($DEBUG)) && echo "[DEBUG] RELEASE=$RELEASE"

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DFAM_RELEASE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP



	if [ ! -e $DBFOLDER_DFAM/$RELEASE/$ASSEMBLY ] || (($UPDATE)); then
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DFAM_RELEASE' for [$ASSEMBLY]"

		DFAM_DATE=$(curl -s -I $DFAM_URL/Dfam$DFAM_RELEASE | grep "Last-Modified: " | sed "s/Last-Modified: //g" | sed "s/\r$//g");

		if (($UPDATE)); then
			if [ -e $DBFOLDER_DFAM/$RELEASE/$ASSEMBLY ]; then mv -f $DBFOLDER_DFAM/$RELEASE/$ASSEMBLY $DBFOLDER_DFAM.V$DFAM_RELEASE; fi;
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
			"release": "'$DFAM_RELEASE'",
			"date": "'$DATE'",
			"files": [ "'$DBFOLDER_DFAM'" ],
			"assembly": [ "'$ASSEMBLY'" ],
			"download": {
				"methode": "'$DOWNLOAD_METHOD'",
				"URL": "'$DFAM_URL/Dfam$DFAM_RELEASE'",
				"file": "'Dfam-A.hmm.gz'",
				"date": "'$DFAM_DATE'"
			}
		}
		';
		echo "$DB_RELEASE_INFOS_JSON" > $DB_TMP/STARK.database.release

		echo "$DFAM_DATABASES/$ASSEMBLY: $DBFOLDER $GENOME
			mkdir -p $DBFOLDER_DFAM
			$ARIA_CMD $DFAM_URL/Dfam_${DFAM_RELEASE}/infrastructure/dfamscan/homo_sapiens_dfam.hmm -d $DB_TMP;
			$ARIA_CMD $DFAM_URL/Dfam_${DFAM_RELEASE}/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3f -d $DB_TMP;
			$ARIA_CMD $DFAM_URL/Dfam_${DFAM_RELEASE}/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3i -d $DB_TMP;
			$ARIA_CMD $DFAM_URL/Dfam_${DFAM_RELEASE}/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3m -d $DB_TMP;
			$ARIA_CMD $DFAM_URL/Dfam_${DFAM_RELEASE}/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3p -d $DB_TMP;
			mkdir -p $DBFOLDER_DFAM/$DFAM_RELEASE/$ASSEMBLY;
			chmod 0775 $DBFOLDER_DFAM/$DFAM_RELEASE/$ASSEMBLY;
			mv $DB_TMP/*.hmm* $DBFOLDER_DFAM/$DFAM_RELEASE/$ASSEMBLY;
			-[ ! -s $DBFOLDER_DFAM/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_DFAM/STARK.database && chmod o+r $DBFOLDER_DFAM/STARK.database; 
			-[ ! -s $DBFOLDER_DFAM/$DFAM_RELEASE/$ASSEMBLY/STARK.database.release ] && cp $DB_TMP/STARK.database.release $DBFOLDER_DFAM/$DFAM_RELEASE/$ASSEMBLY/STARK.database.release && chmod o+r $DBFOLDER_DFAM/$DFAM_RELEASE/$ASSEMBLY/STARK.database.release;
			[ ! -e $DBFOLDER_DFAM/$RELEASE/$ASSEMBLY ] || unlink $DBFOLDER_DFAM/$RELEASE/$ASSEMBLY;
			ln -snf ../$DFAM_RELEASE/$ASSEMBLY $DBFOLDER_DFAM/$RELEASE/$ASSEMBLY;
			rm -rf $DB_TMP;
			touch $DFAM_DATABASES/$ASSEMBLY
		" >> $MK

		MK_ALL="$MK_ALL $DFAM_DATABASES/$ASSEMBLY"
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

	#DBFOLDER_CTAT=$(dirname $CTAT_DATABASES)
	DBFOLDER_CTAT_ROOT=$GENOME.ctat
	DBFOLDER_CTAT=$DBFOLDER_CTAT_ROOT/$CTAT_DATABASES_GENE_SOURCE
	# if [ ! -e $DBFOLDER_CTAT ]; then
	# 	mkdir -p $DBFOLDER_CTAT;
	# fi;
	(($DEBUG)) && echo "#[DEBUG] DBFOLDER_CTAT=$DBFOLDER_CTAT"

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_CTAT ] || (($UPDATE)); then
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for [$ASSEMBLY]"

		# CTAT_DATE=$(curl -s -I $CTAT_CURRENT | grep "Last-Modified: " | sed "s/Last-Modified: //g" | sed "s/\r$//g");
		# CTAT_DATE_RELEASE=$(date -d "$CTAT_DATE");

		if (($UPDATE)); then
			if [ -e $DBFOLDER_CTAT ]; then mv -f $DBFOLDER_CTAT $DBFOLDER_CTAT.V$DATE; fi;
			#mkdir -p $DBFOLDER_CTAT;
		fi;
		
		# DB_INFOS_JSON='
		# {
		# 	"code": "'$DATABASE'",
		# 	"name": "'$DATABASE_NAME'",
		# 	"fullname": "'$DATABASE_FULLNAME'",
		# 	"website": "'$DATABASE_WEBSITE'",
		# 	"description": "'$DATABASE_DESCRIPTION'"
		# }
		# ';
		# echo "$DB_INFOS_JSON" > $DB_TMP/STARK.database

		# DB_RELEASE_INFOS_JSON='
		# {
		# 	"release": "'$RELEASE'",
		# 	"date": "'$RELEASE'",
		# 	"files": [ "'$DBFOLDER_CTAT'" ],
		# 	"assembly": [ "'$ASSEMBLY'" ],
		# 	"download": {
		# 		"methode": "'$DOWNLOAD_METHOD'",
		# 		"URL": "'$(dirname $CTAT_CURRENT)'",
		# 		"file": "'$(basename $CTAT_CURRENT)'",
		# 		"date": "'$CTAT_DATE_RELEASE'"
		# 	}
		# }
		# ';
		# echo "$DB_RELEASE_INFOS_JSON" > $DB_TMP/STARK.database.release

		if [ "$CTAT_DATABASES_GENE_SOURCE" == "refgene" ]; then
			CTAT_REF_GENE_SOURCE=$REFGENE_DATABASES/$ASSEMBLY/ncbiRefSeq.gtf
		elif [ "$CTAT_DATABASES_GENE_SOURCE" == "gencode" ]; then
			CTAT_REF_GENE_SOURCE=$GENCODE_DATABASES/$ASSEMBLY/gencode.gtf
		else
			echo "ERROR: CTAT_DATABASES_GENE_SOURCE '$CTAT_DATABASES_GENE_SOURCE' not supported for CTAT database preparation"
			exit 1;
		fi
		CTAT_REF_GENE_SOURCE_FOLDER=$(dirname $CTAT_REF_GENE_SOURCE)
		

		# Fusion anotation lib preparation
		echo "$DBFOLDER_CTAT_ROOT/fusion_lib: $DBFOLDER $GENOME
			mkdir -p $DBFOLDER_CTAT_ROOT/fusion_lib;
			$ARIA_CMD $CTAT_LIB_SOURCE -d $DBFOLDER_CTAT_ROOT;
			tar -xzf $DBFOLDER_CTAT_ROOT/$(basename $CTAT_LIB_SOURCE) -C $DBFOLDER_CTAT_ROOT/fusion_lib --wildcards --no-anchored "fusion_lib.\*gz"  --transform='s:.*/::';
			touch $DBFOLDER_CTAT_ROOT/fusion_lib/done;
		" >> $MK

		# refGene preparation genome lib
		echo "$DBFOLDER_CTAT: $DBFOLDER $GENOME $DFAM_DATABASES/$ASSEMBLY $PFAM_DATABASES/$ASSEMBLY $CTAT_REF_GENE_SOURCE_FOLDER $DBFOLDER_CTAT_ROOT/fusion_lib
			mkdir -p $DBFOLDER_CTAT;
			# Prepare annotation GTF
			awk 'NR==FNR {chroms[\$\$1]; next} \$\$1 in chroms' <(cut -f1 $GENOME.fai) $CTAT_REF_GENE_SOURCE > $DBFOLDER_CTAT/ref_gene.gtf;
			ls $DBFOLDER_CTAT_ROOT/fusion_lib/fusion_lib*gz | head -n1;
			# Prepare genome lib
			$DOCKER_RUN --rm --name ctat_prep_genome_lib-$ASSEMBLY-$CTAT_DATABASES_GENE_SOURCE-"$(date +%Y%m%d-%H%M%S)" trinityctat/starfusion /usr/local/src/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl --genome_fa $GENOME --gtf $DBFOLDER_CTAT/ref_gene.gtf --fusion_annot_lib=\$\$(ls $DBFOLDER_CTAT_ROOT/fusion_lib/fusion_lib*gz | head -n1) --dfam_db $DFAM_DATABASES/$ASSEMBLY/homo_sapiens_dfam.hmm --pfam_db $PFAM_DATABASES/$ASSEMBLY/Pfam-A.hmm --output_dir $DBFOLDER_CTAT --CPU $THREADS 1> $DBFOLDER_CTAT/prep_genome_lib.log 2> $DBFOLDER_CTAT/prep_genome_lib.err;
			$JAVA -jar $PICARD CreateSequenceDictionary -REFERENCE $DBFOLDER_CTAT/ref_genome.fa -OUTPUT $DBFOLDER_CTAT/ref_genome.dict;
			echo "CTAT for $CTAT_DATABASES_GENE_SOURCE generated" > $DBFOLDER_CTAT/done
		" >> $MK

		MK_ALL="$MK_ALL $DBFOLDER_CTAT"

	fi;
fi;


# ncbiRefSeq.gtf

# ########
# # CTAT OLD #
# ########
# DATABASE="ctat_old"
# DATABASE_NAME="ctat"
# DATABASE_FULLNAME=" CTAT Genome Lib"
# DATABASE_WEBSITE="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/"
# DATABASE_DESCRIPTION=" CTAT Genome Lib is a resource collection used by the Trinity Cancer Transcriptome Analysis Toolkit (CTAT). This CTAT-genome-lib-builder system is leveraged for preparing a target genome and annotation set for use with Trinity CTAT tools, including fusion transcript detection and cancer mutation discovery"

# if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then

# 	DBFOLDER_CTAT=$(dirname $CTAT_DATABASES)
# 	if [ ! -e $DBFOLDER_CTAT/$RELEASE ]; then
# 		mkdir -p $DBFOLDER_CTAT/$RELEASE;
# 	fi;

# 	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
# 	mkdir -p $DB_TMP
# 	chmod 0775 $DB_TMP;

# 	if [ ! -e $DBFOLDER_CTAT/$RELEASE/$ASSEMBLY ] || (($UPDATE)); then
# 		(($VERBOSE)) && echo ""
# 		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for [$ASSEMBLY]"

# 		CTAT_DATE=$(curl -s -I $CTAT_CURRENT | grep "Last-Modified: " | sed "s/Last-Modified: //g" | sed "s/\r$//g");
# 		CTAT_DATE_RELEASE=$(date -d "$CTAT_DATE");

# 		if (($UPDATE)); then
# 			if [ -e $DBFOLDER_CTAT/$RELEASE/$ASSEMBLY ]; then mv -f $DBFOLDER_CTAT/$RELEASE/$ASSEMBLY $DBFOLDER_CTAT.V$DATE; fi;
# 		fi;
		
# 		DB_INFOS_JSON='
# 		{
# 			"code": "'$DATABASE'",
# 			"name": "'$DATABASE_NAME'",
# 			"fullname": "'$DATABASE_FULLNAME'",
# 			"website": "'$DATABASE_WEBSITE'",
# 			"description": "'$DATABASE_DESCRIPTION'"
# 		}
# 		';
# 		echo "$DB_INFOS_JSON" > $DB_TMP/STARK.database

# 		DB_RELEASE_INFOS_JSON='
# 		{
# 			"release": "'$CTAT_DATE'",
# 			"date": "'$CTAT_DATE_RELEASE'",
# 			"files": [ "'$CTAT_CURRENT'" ],
# 			"assembly": [ "'$ASSEMBLY'" ],
# 			"download": {
# 				"methode": "'$DOWNLOAD_METHOD'",
# 				"URL": "'$(dirname $CTAT_CURRENT)'",
# 				"file": "'$(basename $CTAT_CURRENT)'",
# 				"date": "'$CTAT_DATE_RELEASE'"
# 			}
# 		}
# 		';
# 		echo "$DB_RELEASE_INFOS_JSON" > $DB_TMP/STARK.database.release

# 		(($VERBOSE)) && echo "#[INFO] CTAT URL=$CTAT_CURRENT"
# 		(($VERBOSE)) && echo "#[INFO] CTAT RELEASE=$DATE"

# 		echo "$DBFOLDER_CTAT/done: $DBFOLDER
# 			$ARIA_CMD $CTAT_CURRENT -d $DB_TMP;
# 			wget --progress=bar:force:noscroll $CTAT_PM -P $DB_TMP;
# 			mkdir -p $DBFOLDER_CTAT/$DATE/$ASSEMBLY;
# 			chmod 0775 $DBFOLDER_CTAT/$DATE/$ASSEMBLY;
# 			tar -xzf  $DB_TMP/$(basename $CTAT_CURRENT) -C  $DB_TMP --strip-components=1;
# 			$JAVA -jar $PICARD CreateSequenceDictionary -REFERENCE $DB_TMP/ctat_genome_lib_build_dir/ref_genome.fa -OUTPUT $DB_TMP/ctat_genome_lib_build_dir/ref_genome.dict;
# 			cp -R $DB_TMP/ctat_genome_lib_build_dir/* $DBFOLDER_CTAT/$DATE/$ASSEMBLY;
# 			\cp $DB_TMP/AnnotFilterRule.pm $DBFOLDER_CTAT/$DATE/$ASSEMBLY;
# 			-[ ! -s $DBFOLDER_CTAT/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_CTAT/STARK.database && chmod o+r $DBFOLDER_CTAT/STARK.database;
# 			-[ ! -s $DBFOLDER_CTAT/$DATE/$ASSEMBLY/STARK.database.release ] && cp $DB_TMP/STARK.database.release $DBFOLDER_CTAT/$DATE/$ASSEMBLY/STARK.database.release && chmod o+r $DBFOLDER_CTAT/$DATE/$ASSEMBLY/STARK.database.release;
# 			[ ! -e $DBFOLDER_CTAT/$RELEASE/$ASSEMBLY ] || unlink $DBFOLDER_CTAT/$RELEASE/$ASSEMBLY;
# 			ln -snf ../$DATE/$ASSEMBLY $DBFOLDER_CTAT/$RELEASE/$ASSEMBLY;
# 			rm -rf $DB_TMP;
# 		" >> $MK
# 		MK_ALL="$MK_ALL $DBFOLDER_CTAT/done"
# 	fi;
# fi;

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