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
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Configuration
ENV_CONFIG=$(find -L $SCRIPT_DIR/.. -name config.app)
source $ENV_CONFIG 
# 1>/dev/null 2>/dev/null

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

# default Genome
if [ -e $GENOME ] ; then

	# Genome index
	if [ ! -d $GENOME.hts-ref ] || [ ! "$(ls -A $GENOME.hts-ref)" ]; then
		if [ "$SAMTOOLS" != "" ]; then
		    echo "$GENOME.hts-ref/done: $GENOME				mkdir -p $GENOME.hts-ref;
				perl $(dirname $SAMTOOLS)/seq_cache_populate.pl -root $GENOME.hts-ref $GENOME 1>/dev/null 2>/dev/null;
				echo 'done.' > $REF.hts-ref/done;
		    " >> $MK
			MK_ALL="$MK_ALL $GENOME.hts-ref/done"
		fi;
	fi;

	## BOWTIE index
	if [ ! -e $(dirname $GENOME)/$ASSEMBLY.rev.1.bt2 ]; then
		if [ "$BOWTIE" != "" ]; then
		    echo " $(dirname $GENOME)/$ASSEMBLY.rev.1.bt2: $GENOME
				$(dirname $BOWTIE)/bowtie2-build --threads $THREADS --packed $GENOME $ASSEMBLY.rev.1.bt2 ;
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
		    echo "$(dirname $GENOME)/$ASSEMBLY.dict: $REF
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
			STAR --runThreadN $THREADS --runMode genomeGenerate --genomeDir $(dirname $GENOME)/$(basename $GENOME).star.idx --genomeFastaFiles $GENOME --sjdbGTFfile $DBFOLDER_GENCODE/current/$ASSEMBLY/gencode.*annotation.gtf ;
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