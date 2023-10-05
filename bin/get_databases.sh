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
RELEASE_NOTES=$RELEASE_NOTES"# STARK 19: \n";



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
# function in_array
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
# ACTION
ACTION=0
[ $BUILD ] || [ $UPDATE ] && ACTION=1;

# DATABASES FOLDER
[ ! -z $DATABASES ] && [ ! -d $DATABASES ] && mkdir -p $DATABASES && echo "#[INFO] Create databases folder '$DATABASES' "

# DATABASES LIST
[ "$DATABASES_LIST_INPUT" == "" ] && DATABASES_LIST_INPUT="ALL"
echo "#[INFO] Databases list: '$DATABASES_LIST_INPUT' "

# ENV
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

DATE=$(date '+%Y%m%d-%H%M%S')

# output
echo ""
echo "#[INFO] DB_RELEASE=$DATE"
echo "#[INFO] ASSEMBLY=$ASSEMBLY"
echo "#[INFO] THREADS=$THREADS"

MK=$TMP_DATABASES_DOWNLOAD_FOLDER/mk
MK_LOG=$TMP_DATABASES_DOWNLOAD_FOLDER/mk.log
MK_ERR=$TMP_DATABASES_DOWNLOAD_FOLDER/mk.err
MK_ALL=""
> $MK

# INFOS
DOWNLOAD_METHOD="STARK Databases downloading script [$SCRIPT_RELEASE-$SCRIPT_DATE]"

##############################
# GATK VARIANT RECALIBRATION #
##############################

DATABASE="gatk"
DATABASE_NAME="GATK"
DATABASE_FULLNAME="GATK Databases"
DATABASE_WEBSITE="https://www.broadinstitute.org//"
DATABASE_DESCRIPTION="Databases for GATK Variant Recalibration"

if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then

	DBFOLDER_GATK_URL_FULL_COPY=0
	GATK_RESOURCE_NB=0
	MK_DBFOLDER_GATK_ALL=""
	> $MK.existing_gatk_db

	for GATK_RESOURCE in $GATK_DATABASES_LIST; do
		DB_TARGET=$DBFOLDER_GATK/current/$ASSEMBLY/$(echo $GATK_RESOURCE | cut -d: -f2);	# /STARK/databases/gatk/current/1000G_omni2.5.b37.vcf.gz
		DB_TARGET_FILE=$(basename $DB_TARGET);									# 1000G_omni2.5.b37.vcf.gz
		DB_TARGET_FOLDER=$(dirname $DB_TARGET);									# /STARK/databases/gatk/current
		DB_RELEASE=$DATE;
		DB_RELEASE_FILE=$DB_TARGET_FILE;										# 1000G_omni2.5.b37.vcf.gz
		DB_RELEASE_FOLDER=$DBFOLDER_GATK/$DATE/$ASSEMBLY;						# /STARK/databases/gatk/DATE/hg19
		DB_RELEASE_FILE_PATH="$DB_RELEASE_FOLDER/$DB_RELEASE_FILE";				# /STARK/databases/gatk/DATE/hg19/1000G_omni2.5.b37.vcf.gz

		DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
		mkdir -p $DB_TMP

		if (($UPDATE)); then
			if [ -e $DB_TARGET ]; then mv -f $DB_TARGET $DB_TARGET.V$DATE; fi;
		fi;

		if [ -e $DB_TARGET_FOLDER/$(echo $GATK_RESOURCE | cut -d: -f2) ]; then
			echo "$DB_RELEASE_FILE_PATH: $DBFOLDER $GENOME
				mkdir -p $DB_RELEASE_FOLDER
				rsync -ar $DB_TARGET_FOLDER/$(echo $GATK_RESOURCE | cut -d: -f2) $DB_RELEASE_FILE_PATH
				rsync -ar $DB_TARGET_FOLDER/$(echo $GATK_RESOURCE | cut -d: -f2).tbi $DB_RELEASE_FILE_PATH.tbi
			" >> $MK.existing_gatk_db
			MK_DBFOLDER_GATK_ALL_existing_gatk_db="$MK_DBFOLDER_GATK_ALL_existing_gatk_db $DB_RELEASE_FILE_PATH"
			MK_ALL_existing_gatk_db="$MK_ALL_existing_gatk_db $DB_RELEASE_FILE_PATH" 

		else
			((GATK_RESOURCE_NB++))
			(($VERBOSE)) && echo ""
			(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME/$GATK_RESOURCE' release '$DATE' for '$ASSEMBLY'"

			if [ $ASSEMBLY == "hg38" ]; then
				DBFOLDER_GATK_URL="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0"
			elif [ $ASSEMBLY == "hg19" ]; then
				DBFOLDER_GATK_URL="https://data.broadinstitute.org/snowman/hg19/variant_calling/vqsr_resources/Exome/v2"
			else
				DBFOLDER_GATK_URL="https://data.broadinstitute.org/snowman/$ASSEMBLY/variant_calling/vqsr_resources/Exome/v2"
			fi;

			DBFOLDER_GATK_URL_FILE="$DB_TARGET_FILE"
			DB_TARGET_FILE_LIST="$DB_TARGET_FILE_LIST $DB_TARGET_FILE"
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
				mkdir -p $DB_RELEASE_FOLDER
				chmod 0775 $DB_RELEASE_FOLDER
			" >> $MK

			echo "$DB_RELEASE_FILE_PATH: $DB_TMP/$DB_RELEASE_FILE
				mkdir -p $DB_RELEASE_FOLDER/original
				rsync -ar $DB_TMP/original/$DBFOLDER_GATK_URL_FILE $DB_RELEASE_FOLDER/original/
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
		"date": "'$$DATE'",
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
		echo "$DB_TARGET_FOLDER: $MK_DBFOLDER_GATK_ALL
			mkdir -p $DB_RELEASE_FOLDER/original
			chmod 0775 $DB_RELEASE_FOLDER -R
			-[ ! -s $DBFOLDER_GATK/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_GATK/STARK.database
			cp $DB_TMP/STARK.database.release $DB_RELEASE_FOLDER/
			chmod o+r $DBFOLDER_GATK/STARK.database $DB_RELEASE_FOLDER/STARK.database.release
			[ ! -e $DBFOLDER_GATK/current ] || unlink $DBFOLDER_GATK/current
			ln -snf $DBFOLDER_GATK/$DATE $DBFOLDER_GATK/current
			rm -rf $DB_TMP;
		" >> $MK
		MK_ALL="$MK_ALL $DB_TARGET_FOLDER" 
	fi;
fi;


#############################
# ANNOVAR & SNPEFF & REFSEQ #
#############################
# Install with HOWARD v2

DATABASE="snpeff"
DATABASE_NAME="SnpEff"
DATABASE_FULLNAME="SnpEff Annotations"
DATABASE_WEBSITE="http://snpeff.sourceforge.net/"
DATABASE_DESCRIPTION="Genetic variant annotation and functional effect prediction toolbox"

if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then

	DBFOLDER_SNPEFF=$(dirname $SNPEFF_DATABASES)

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_SNPEFF ] || (($UPDATE)); then
		if (($UPDATE)); then
			if [ -e $DBFOLDER_SNPEFF ]; then mv -f $DBFOLDER_SNPEFF $DBFOLDER_SNPEFF.V$DATE; fi;
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

		echo "$DBFOLDER_SNPEFF: $DBFOLDER
			howard databases --assembly='$ASSEMBLY' --download-snpeff=$DBFOLDER_SNPEFF/$DATE
			
			-[ ! -s $DBFOLDER_SNPEFF/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_SNPEFF/STARK.database && chmod o+r $DBFOLDER_SNPEFF/STARK.database 
			[ ! -e $DBFOLDER_SNPEFF/current ] || unlink $DBFOLDER_SNPEFF/current
			ln -snf $DBFOLDER_SNPEFF/$DATE $DBFOLDER_SNPEFF/current

			rm -rf $DB_TMP;
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_SNPEFF"
	fi;
fi;


if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then
	DATABASE="annovar"
	DATABASE_NAME="ANNOVAR"
	DATABASE_FULLNAME="ANNOVAR Annotations"
	DATABASE_WEBSITE="https://doc-openbio.readthedocs.io/projects/annovar/"
	DATABASE_DESCRIPTION="ANNOVAR is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes"
	
	DBFOLDER_ANNOVAR=$(dirname $ANNOVAR_DATABASES)

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_ANNOVAR ] || (($UPDATE)); then
		if (($UPDATE)); then
			if [ -e $DBFOLDER_ANNOVAR ]; then mv -f $DBFOLDER_ANNOVAR $DBFOLDER_ANNOVAR.V$DATE; fi;
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
	
		# remove dbnsfp42a (memoryleak)
		echo "$DBFOLDER_ANNOVAR: $DBFOLDER
			howard databases --assembly='$ASSEMBLY' --download-annovar=$DBFOLDER_ANNOVAR/$DATE --download-annovar-files='refGene,gnomad_exome,cosmic70,clinvar_202*,nci60' 
			
			-[ ! -s $DBFOLDER_ANNOVAR/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_ANNOVAR/STARK.database && chmod o+r $DBFOLDER_ANNOVAR/STARK.database 
			[ ! -e $DBFOLDER_ANNOVAR/current ] || unlink $DBFOLDER_ANNOVAR/current
			ln -snf $DBFOLDER_ANNOVAR/$DATE $DBFOLDER_ANNOVAR/current

			rm -rf $DB_TMP;
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_ANNOVAR"
	fi;
fi;


DATABASE="refGene"
DATABASE_NAME="RefGene"
DATABASE_FULLNAME="Reference Genes"
DATABASE_WEBSITE="https://genome.ucsc.edu/"
DATABASE_DESCRIPTION="Known human protein-coding and non-protein-coding genes taken from the NCBI RNA reference sequences collection (RefSeq)"

if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then
	
	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_REFGENE ] || (($UPDATE)); then
		if (($UPDATE)); then
			if [ -e $DBFOLDER_REFGENE ]; then mv -f $DBFOLDER_REFGENE $DBFOLDER_REFGENE.V$DATE; fi;
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


		echo "$DB_TMP/ncbiRefSeq.txt: $DBFOLDER
			howard databases --assembly='$ASSEMBLY' --download-refseq=$DB_TMP
			cat $DB_TMP/ncbiRefSeq.txt | while IFS='' read -r line; do \
				CHROM=\$\$(echo \$\$line | awk '{print \$\$3}' | cut -d\"_\" -f1); \
				NM=\$\$(echo \$\$line | awk '{print \$\$2}'); \
				GENE=\$\$(echo \$\$line | awk '{print \$\$13}'); \
				STRAND=\$\$(echo \$\$line | awk '{print \$\$4}'); \
				echo \$\$line | awk '{print \$\$10}' | tr \",\" \"\\n\" | grep -v '^\$\$' > $TMP_DATABASES_DOWNLOAD_RAM/refGene.unsorted.bed.NM1 ; \
				echo \$\$line | awk '{print \$\$11}' | tr \",\" \"\\n\" | grep -v '^\$\$' > $TMP_DATABASES_DOWNLOAD_RAM/refGene.unsorted.bed.NM2; \
				paste $TMP_DATABASES_DOWNLOAD_RAM/refGene.unsorted.bed.NM1 $TMP_DATABASES_DOWNLOAD_RAM/refGene.unsorted.bed.NM2 | while read SS ; do \
					echo -e \"\$\$CHROM\t\$\$SS\t\$\$GENE\t\$\$NM\t\$\$STRAND\" >> $DB_TMP/refGene.unsorted.bed; \
				done; \
			done
			cat $DB_TMP/refGene.unsorted.bed | sort -k1,1V -k2,2n -k3,3n > $DB_TMP/refGene.$ASSEMBLY.bed
			rsync -ar $DB_TMP/refGene.$ASSEMBLY.bed $DBFOLDER_REFGENE/$DATE/refGene.$ASSEMBLY.bed
		" >> $MK

		if false; then
			# Generates refGene genes file (txt file) from bed file with ref genome
			echo "$DB_TMP/refGene.$ASSEMBLY.txt: $DBFOLDER $DB_TMP/refGene.$ASSEMBLY.bed "$(dirname $GENOME)/$ASSEMBLY.dict"
				mkdir -p $DB_TMP
				awk -F'\t' 'substr(\$\$6,1,2)==\"NM\" {print \$\$0}' $DB_TMP/refGene.$ASSEMBLY.bed | grep \$\$(grep \"@SQ\" "$(dirname $GENOME)/$ASSEMBLY.dict" | cut -f2 | cut -d: -f2 | tr '\n' ' ' | sed 's/chr/ -e ^chr/gi') > $DB_TMP/refGene.$ASSEMBLY.txt
				rsync -ar $DB_TMP/refGene.$ASSEMBLY.txt $DBFOLDER_REFGENE/$DATE/refGene.$ASSEMBLY.txt
				" >> $MK
		fi;

		echo "$DBFOLDER_REFGENE: $DBFOLDER
			-[ ! -s $DBFOLDER_REFGENE/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_REFGENE/STARK.database && chmod o+r $DBFOLDER_REFGENE/STARK.database 
			[ ! -e $DBFOLDER_REFGENE/current ] || unlink $DBFOLDER_REFGENE/current
			ln -snf $DBFOLDER_REFGENE/$DATE $DBFOLDER_REFGENE/current
			rm -rf $DB_TMP;
			" >> $MK

		MK_ALL="$MK_ALL $DBFOLDER_REFGENE"
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

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP


	if [ ! -e $DBFOLDER_ARRIBA ] || (($UPDATE)); then
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for ' [$ASSEMBLY]"

		ARRIBA_CURRENT="https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz";

		if (($UPDATE)); then
			if [ -e $DBFOLDER_ARRIBA ]; then mv -f $DBFOLDER_ARRIBA $DBFOLDER_ARRIBA.V$DATE; fi;
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

		echo "$DBFOLDER_ARRIBA: $DBFOLDER
			wget --progress=bar:force:noscroll $ARRIBA_CURRENT -P $DB_TMP
			mkdir -p $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY
			chmod 0775 $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY
			tar -xzf  $DB_TMP/$(basename $ARRIBA_CURRENT) -C $DB_TMP --strip-components=1
			mv $DB_TMP/database/blacklist_$ASSEMBLY* $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY
			mv $DB_TMP/database/cytobands_$ASSEMBLY* $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY
			mv $DB_TMP/database/known_fusions_$ASSEMBLY* $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY
			mv $DB_TMP/database/protein_domains_$ASSEMBLY* $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY

			-[ ! -s $DBFOLDER_ARRIBA/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_ARRIBA/STARK.database && chmod o+r $DBFOLDER_ARRIBA/STARK.database 
			-[ ! -s $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY/STARK.database.release ] && cp $DB_TMP/STARK.database.release $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY/STARK.database.release && chmod o+r $DBFOLDER_ARRIBA/$DATE/$ASSEMBLY/STARK.database.release
			
			[ ! -e $DBFOLDER_ARRIBA/current ] || unlink $DBFOLDER_ARRIBA/current
			ln -snf $DBFOLDER_ARRIBA/$DATE $DBFOLDER_ARRIBA/current
			rm -rf $DB_TMP;
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_ARRIBA"
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

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_CTAT ] || (($UPDATE)); then
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for [$ASSEMBLY]"

	if [ $ASSEMBLY == "hg19" ] ; then CTAT_CURRENT="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz"; fi;
	if [ $ASSEMBLY == "hg38" ] ; then CTAT_CURRENT="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz"; fi;
		CTAT_DATE=$(curl -s -I $CTAT_CURRENT | grep "Last-Modified: " | sed "s/Last-Modified: //g" | sed "s/\r$//g");
		CTAT_DATE_RELEASE=$(date -d "$CTAT_DATE");

		if (($UPDATE)); then
			if [ -e $DBFOLDER_CTAT ]; then mv -f $DBFOLDER_CTAT $DBFOLDER_CTAT.V$DATE; fi;
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

		# DATABASE Release Infos
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

		# MK
		echo "$DBFOLDER_CTAT: $DBFOLDER
			wget --progress=bar:force:noscroll $CTAT_CURRENT -P $DB_TMP;
			mkdir -p $DBFOLDER_CTAT/$DATE/$ASSEMBLY
			chmod 0775 $DBFOLDER_CTAT/$DATE/$ASSEMBLY
			tar -xzf  $DB_TMP/$(basename $CTAT_CURRENT) -C  $DB_TMP --strip-components=1
			cp -R $DB_TMP/ctat_genome_lib_build_dir/* $DBFOLDER_CTAT/$DATE/$ASSEMBLY
			-[ ! -s $DBFOLDER_CTAT/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_CTAT/STARK.database && chmod o+r $DBFOLDER_CTAT/STARK.database 
			-[ ! -s $DBFOLDER_CTAT/$DATE/$ASSEMBLY/STARK.database.release ] && cp $DB_TMP/STARK.database.release $DBFOLDER_CTAT/$DATE/$ASSEMBLY/STARK.database.release && chmod o+r $DBFOLDER_CTAT/$DATE/$ASSEMBLY/STARK.database.release

			[ ! -e $DBFOLDER_CTAT/current ] || unlink $DBFOLDER_CTAT/current
			ln -snf $DBFOLDER_CTAT/$DATE $DBFOLDER_CTAT/current

			rm -rf $DB_TMP;
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_CTAT"
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

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_GENCODE ] || (($UPDATE)); then
		(($VERBOSE)) && echo ""
		(($VERBOSE)) && echo "#[INFO] DATABASE '$DATABASE_NAME' release '$DATE' for [$ASSEMBLY]"
	# for hg19 the last gencode version is v19
	if [ $ASSEMBLY == "hg19" ] ; then GENCODE_CURRENT="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"; fi;
	# for hg38 the first gencode version is v20 ; current version (10/2023) is v44
	if [ $ASSEMBLY == "hg38" ] ; then GENCODE_CURRENT="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz"; fi;
		GENCODE_DATE=$(curl -s -I $GENCODE_CURRENT | grep "Last-Modified: " | sed "s/Last-Modified: //g" | sed "s/\r$//g");
		GENCODE_DATE_RELEASE=$(date -d "$GENCODE_DATE");

		if (($UPDATE)); then
			if [ -e $DBFOLDER_GENCODE ]; then mv -f $DBFOLDER_GENCODE $DBFOLDER_GENCODE.V$DATE; fi;
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

		# DATABASE Release Infos
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

		(($VERBOSE)) && echo "#[INFO] CTAT URL=$GENCODE_CURRENT"
		(($VERBOSE)) && echo "#[INFO] CTAT RELEASE=$DATE"

		# MK
		echo "$DBFOLDER_GENCODE: $DBFOLDER
			mkdir -p $DBFOLDER_GENCODE/$DATE/$ASSEMBLY
			chmod 0775 $DBFOLDER_GENCODE/$DATE/$ASSEMBLY
			wget --progress=bar:force:noscroll $GENCODE_CURRENT -P $DB_TMP
			cp $DB_TMP/$(basename $GENCODE_CURRENT) $DBFOLDER_GENCODE/$DATE/$ASSEMBLY
			gzip -d $DBFOLDER_GENCODE/$DATE/$ASSEMBLY/$(basename $GENCODE_CURRENT)
			-[ ! -s $DBFOLDER_GENCODE/STARK.database ] && cp $DB_TMP/STARK.database $DBFOLDER_GENCODE/STARK.database && chmod o+r $DBFOLDER_GENCODE/STARK.database 
			-[ ! -s $DBFOLDER_GENCODE/$DATE/$ASSEMBLY/STARK.database.release ] && cp $DB_TMP/STARK.database.release $DBFOLDER_GENCODE/$DATE/$ASSEMBLY/STARK.database.release && chmod o+r $DBFOLDER_GENCODE/$DATE/$ASSEMBLY/STARK.database.release

			[ ! -e $DBFOLDER_GENCODE/current ] || unlink $DBFOLDER_GENCODE/current
			ln -snf $DBFOLDER_GENCODE/$DATE $DBFOLDER_GENCODE/current

			rm -rf $DB_TMP;
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_GENCODE"
	fi;
fi;


if in_array $DATABASE $DATABASES_LIST_INPUT || in_array ALL $DATABASES_LIST_INPUT; then
	# GENOME=$DATABASES/genomes/current/$ASSEMBLY/$ASSEMBLY.fa
	DBFOLDER_GENOME=$DATABASES/genomes

	DB_TMP=$TMP_DATABASES_DOWNLOAD_FOLDER/$DATABASE/$DATE
	mkdir -p $DB_TMP
	chmod 0775 $DB_TMP;

	if [ ! -e $DBFOLDER_GENOME ] || (($UPDATE)); then
		if (($UPDATE)); then
			if [ -e $DBFOLDER_GENOME ]; then mv -f $DBFOLDER_GENOME $DBFOLDER_GENOME.V$DATE; fi;
		fi;
		
		echo "$DBFOLDER_GENOME: $DBFOLDER
			howard databases --assembly='$ASSEMBLY' --download-genomes=$DBFOLDER_GENOME/$DATE
			
			[ ! -e $DBFOLDER_GENOME/current ] || unlink $DBFOLDER_GENOME/current
			ln -snf $DBFOLDER_GENOME/$DATE $DBFOLDER_GENOME/current

			rm -rf $DB_TMP;
		" >> $MK
		MK_ALL="$MK_ALL $DBFOLDER_GENOME"
	fi;
fi;

if [ -e $GENOME ] ; then


	# Samtools genome index
	if [ ! -d $GENOME.hts-ref ] ; then
		if [ "$SAMTOOLS" != "" ]; then
		    echo "$GENOME.hts-ref: $GENOME
				mkdir -p $GENOME.hts-ref;
				perl $(dirname $SAMTOOLS)/seq_cache_populate.pl -root $GENOME.hts-ref $GENOME;
				echo 'done.' > $GENOME.hts-ref;
		    " >> $MK
			MK_ALL="$MK_ALL $GENOME.hts-ref"
		fi;
	fi;

	## BOWTIE index
	if [ ! -e $(dirname $GENOME)/$ASSEMBLY.rev.1.bt2 ]; then
		if [ "$BOWTIE" != "" ]; then
		    echo "$(dirname $GENOME)/$ASSEMBLY.rev.1.bt2: $GENOME
				$(dirname $BOWTIE)/bowtie2-build --threads $THREADS --packed $GENOME $(dirname $GENOME)/$ASSEMBLY.rev ;
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

	## BWA2 index TODO make it optional
	if [ ! -e $GENOME.bwt.2bit.64 ]; then
		if [ "$BWA2" != "" ]; then
		    echo "$GENOME.bwt.2bit.64: $GENOME
				$BWA2 index $GENOME;
		    " >> $MK
			MK_ALL="$MK_ALL $GENOME.bwt.2bit.64"
		fi;
	fi;

	## SAMTOOLS index
	if [ ! -e $GENOME.fai ]; then
		if [ "$SAMTOOLS" != "" ]; then
		    echo "$GENOME.fai: $GENOME
				$SAMTOOLS faidx $GENOME;
		    " >> $MK
			MK_ALL="$MK_ALL $GENOME.fai"
		fi;
	fi;

	## PICARD index
	if [ ! -e $(dirname $GENOME)/$ASSEMBLY.dict ]; then
		if [ "$PICARD" != "" ]; then
		    echo "$(dirname $GENOME)/$ASSEMBLY.dict: $GENOME
				$JAVA -jar $PICARD CreateSequenceDictionary \
					-REFERENCE $GENOME \
					-OUTPUT $(dirname $GENOME)/$ASSEMBLY.dict ;
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
		    echo "$(dirname $GENOME)/$(basename $GENOME).star.idx: $GENOME
			STAR --runThreadN $THREADS --runMode genomeGenerate --genomeDir $(dirname $GENOME)/ --genomeFastaFiles $GENOME --sjdbGTFfile $DBFOLDER_GENCODE/current/$ASSEMBLY/gencode.*annotation.gtf ;
		    " >> $MK
			MK_ALL="$MK_ALL $(dirname $GENOME)/$(basename $GENOME).star.idx"
		fi;
	fi;
fi;

if [ ! -z "$MK_ALL" ]; then
	echo "$DBFOLDER:
		mkdir -p $DBFOLDER
		chmod 0775 $DBFOLDER
	" >> $MK

	echo "all: $MK_ALL
		echo '#[INFO] Build release: $DATE' >> $DATABASES/STARK.download.releases
	" >> $MK

	# DATABASES INIT
	(($VERBOSE)) && echo "DATABASE INIT"

	if ((1)); then
		if (($BUILD)) || (($UPDATE)); then
			echo "#[INFO] DATABASES DOWNLOADING..."
			if (($VERBOSE)) || (($DEBUG)); then
				if (($DEBUG)); then
					make -k -j $THREADS $MK_OPTION -f $MK all;
				elif (($VERBOSE)); then
					make -k -j $THREADS $MK_OPTION -f $MK all;
				fi;
			else
				make -k -j $THREADS $MK_OPTION -f $MK all 1>$MK_LOG 2>$MK_ERR;
				if (($(cat $MK_LOG $MK_ERR | grep "\*\*\*" -c))); then
					echo "#[ERROR] Fail download databases"
					exit 1
				fi;
			fi;
			echo "#[INFO] DATABASES DOWNLOADED"
		else
			echo "## use --build or --update to download databases"
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
		(($VERBOSE)) && echo "#[INFO] DOWNLOAD GENOME $ASSEMBLY"
		(($VERBOSE)) && echo "#"
	fi;
fi;

if ((0)); then
	rm -Rf $TMP_DATABASES_DOWNLOAD_FOLDER
	rm -Rf $TMP_DATABASES_DOWNLOAD_RAM
fi;

exit 0;