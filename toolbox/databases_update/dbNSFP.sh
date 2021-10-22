#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARK_DATABASES_DBNSFP"
SCRIPT_DESCRIPTION="STARK DATABASES dbNSFP databases"
SCRIPT_RELEASE="0.9.0.0"
SCRIPT_DATE="08/04/2021"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="HUS"
SCRIPT_LICENCE="GNU-AGPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.0.0-08/04/2021: Script creation\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Configuration
ENV_CONFIG=$(find -L $SCRIPT_DIR/../.. -name config.app)

source $ENV_CONFIG 1>/dev/null 2>/dev/null

ENV_TOOLS=$(find -L $SCRIPT_DIR/../.. -name tools.app)

source $ENV_TOOLS 1>/dev/null 2>/dev/null


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
	echo "# RELEASE NOTES:";
	echo -e $RELEASE_NOTES
}

# Usage
function usage {
	echo "# USAGE: $(basename $0) [options...]";
	echo "# --dbnsfp_release=<STRING>                     dbNSFP Release";
	echo "#                                               Default: '4.2a'";
	echo "# --dbnsfp_url=<URL>                            dbNSFP URL to download";
	echo "#                                               Will be download if <dbnsfp_zip> does not exists";
	echo "#                                               Default: 'ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP<dbnsfp_release>.zip'";
	echo "# --dbnsfp_url_readme=<URL>                     dbNSFP URL to readme file";
	echo "# --dbnsfp_folder=<FOLDER>                      dbNSFP folder (decompressed)";
	echo "#                                               Default: '/STARK/databases/dbnsfp/<dbnsfp_release>'";
	echo "# --dbnsfp_zip=<FILE>                           dbNSFP zip file (downloaded)";
	echo "#                                               Will be decompressed if <dbnsfp_folder> does not exists";
	echo "#                                               Default: '<dbnsfp_folder>/dbNSFP<dbnsfp_release>.zip'";
	echo "# --dbnsfp_zip_readme=<URL>                     dbNSFP zip readme file (downloaded)";
	echo "#                                               Default: '<dbnsfp_zip>.readme'";
	echo "# --dbnsfp_file=<FILE>                          dbNSFP file (recomposed)";
	echo "#                                               Default: '<dbnsfp_folder>/dbNSFP<dbnsfp_release>_variant.gz'";
	echo "# --dbnsfp_annovar=<FILE>                       dbNSFP annovar (formatted)";
	echo "#                                               Formatted from <dbnsfp_file> with <assembly>";
	echo "#                                               Default: '<dbnsfp_folder>/<assembly>_dbNSFP<dbnsfp_release_annovar_format>_custom.txt'";
	echo "#                                               Note: dbnsfp_relase will be formatted into <dbnsfp_release_annovar_format> to remove '.'";
	echo "#                                               Note: index file <dbnsfp_annovar>.idx will be generated";
	echo "# --dbnsfp_howard_config_annotation=<FILE>      HOWARD config annotation ini file";
	echo "#                                               Default: '<dbnsfp_annovar>.config.annotation.ini'";
	echo "# --assembly=<STRING>                           Assembly";
	echo "#                                               Default: 'hg19'";
	echo "# --verbose                                     VERBOSE";
	echo "# --debug                                       DEBUG";
	echo "# --release                                     RELEASE";
	echo "# --help                                        HELP";
	echo "#";
}

# header
header;


####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "vdnh" --long "dbnsfp_release:,dbnsfp_url:,dbnsfp_url_readme:,dbnsfp_folder:,dbnsfp_zip:,dbnsfp_zip_readme:,dbnsfp_file:,dbnsfp_annovar:,dbnsfp_howard_config_annotation:,assembly:,verbose,debug,release,help" -- "$@" 2> /dev/null)

eval set -- "$ARGS"
while true
do
	case "$1" in
		--dbnsfp_release)
			DBNSFP_RELEASE="$2"
			shift 2
			;;
		--dbnsfp_url)
			DBNSFP_URL="$2";
			shift 2
			;;
		--dbnsfp_url_readme)
			DBNSFP_URL_README="$2";
			shift 2
			;;
		--dbnsfp_folder)
			DBNSFP_FOLDER=$2;
			shift 2
			;;
		--dbnsfp_zip)
			DBNSFP_ZIP="$2";
			shift 2
			;;
		--dbnsfp_zip_readme)
			DBNSFP_ZIP_README="$2";
			shift 2
			;;
		--dbnsfp_file)
			DBNSFP_FILE=$2;
			shift 2
			;;
		--dbnsfp_annovar)
			DBNSFP_ANNOVAR=$2;
			shift 2
			;;
		--dbnsfp_howard_config_annotation)
			DBNSFP_HOWARD_CONFIG_ANNOTATION=$2;
			shift 2
			;;
		--assembly)
			ASSEMBLY=$2;
			shift 2
			;;
		-v|--verbose)
			VERBOSE=1
			shift 1
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




####################################################################################################################################
# Checking the input parameter
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#if [ -z "$DBNSFP_ZIP" ]; then
if ((0)); then
	echo "#[ERROR] Required parameter: --dbnsfp_file. Use --help to display the help." && echo "" && usage && exit 1;
fi
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


## PARAMETERS
##############

#dbnsfp_release:,dbnsfp_ftp:,dbnsfp_folder:,dbnsfp_zip:,dbnsfp_file:,dbnsfp_annovar:,assembly:


# DBNSFP_RELEASE
if [ "$DBNSFP_RELEASE" == "" ]; then
	DBNSFP_RELEASE="4.2a"
fi;
echo "#[INFO] DBNSFP Release '$DBNSFP_RELEASE'"


# ASSEMBLY
if [ "$ASSEMBLY" == "" ]; then
	ASSEMBLY="hg19"
fi;
echo "#[INFO] DBNSFP Assembly '$ASSEMBLY'"


# DBNSFP_RELEASE_FORMATTED
DBNSFP_RELEASE_FORMATTED=$(echo $DBNSFP_RELEASE | sed "s/\.//gi")
echo "#[INFO] DBNSFP Release formatted '$DBNSFP_RELEASE_FORMATTED'"


# DBNSFP_URL
if [ "$DBNSFP_URL" == "" ]; then
	DBNSFP_URL="ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP$DBNSFP_RELEASE.zip"
fi;
echo "#[INFO] DBNSFP URL '$DBNSFP_URL'"


# DBNSFP_URL_README
if [ "$DBNSFP_URL_README" == "" ]; then
	# try to find URL
	#[ "$DBNSFP_RELEASE" == "4.2a" ] && DBNSFP_URL_README="https://drive.google.com/file/d/1Vse3b_qw_E46eDcsuLegF5HxpodAZ6Og/view?usp=sharing"
	[ "$DBNSFP_RELEASE" == "4.2a" ] && DBNSFP_URL_README="https://drive.google.com/uc?export=download&id=1Vse3b_qw_E46eDcsuLegF5HxpodAZ6Og"
	#wget "https://drive.google.com/uc?export=download&id=1Vse3b_qw_E46eDcsuLegF5HxpodAZ6Og"
	#DBNSFP_URL_README=""
fi;
echo "#[INFO] DBNSFP URL README '$DBNSFP_URL_README'"


# DBNSFP_FOLDER
if [ "$DBNSFP_FOLDER" == "" ]; then
	DBNSFP_FOLDER="/STARK/databases/dbnsfp/$DBNSFP_RELEASE"
fi;
if ! mkdir -p $DBNSFP_FOLDER 2>/dev/null; then
	echo "#[WARNING] DBNSFP Folder '$DBNSFP_FOLDER' can NOT be created"
fi;
if [ ! -d "$DBNSFP_FOLDER" ]; then
	echo "#[ERROR] DBNSFP Folder '$DBNSFP_FOLDER' does NOT exists"
	usage
	exit 0
fi;
echo "#[INFO] DBNSFP Folder '$DBNSFP_FOLDER'"


# DBNSFP_ZIP
if [ "$DBNSFP_ZIP" == "" ]; then
	DBNSFP_ZIP="$DBNSFP_FOLDER/dbNSFP$DBNSFP_RELEASE.zip"
fi;
echo "#[INFO] DBNSFP Zip '$DBNSFP_ZIP'"


# DBNSFP_ZIP_README
if [ "$DBNSFP_ZIP_README" == "" ]; then
	DBNSFP_ZIP_README="$DBNSFP_ZIP.readme"
fi;
echo "#[INFO] DBNSFP Readme '$DBNSFP_ZIP_README'"



# DBNSFP_FILE
if [ "$DBNSFP_FILE" == "" ]; then
	DBNSFP_FILE="$DBNSFP_FOLDER/dbNSFP"$DBNSFP_RELEASE"_variant.gz"
fi;
echo "#[INFO] DBNSFP File '$DBNSFP_FILE'"


# DBNSFP_ANNOVAR
if [ "$DBNSFP_ANNOVAR" == "" ]; then
	DBNSFP_ANNOVAR="$DBNSFP_FOLDER/"$ASSEMBLY"_dbNSFP"$DBNSFP_RELEASE_FORMATTED"_custom.txt"
fi;
echo "#[INFO] DBNSFP Annovar file '$DBNSFP_ANNOVAR'"

# DBNSFP_HOWARD_CONFIG_ANNOTATION
if [ "$DBNSFP_HOWARD_CONFIG_ANNOTATION" == "" ]; then
	DBNSFP_HOWARD_CONFIG_ANNOTATION="$DBNSFP_ANNOVAR.config.annotation.ini"
fi;
echo "#[INFO] DBNSFP HOWARD config annotation ini file '$DBNSFP_HOWARD_CONFIG_ANNOTATION'"


###### ZIP ######

if [ ! -e $DBNSFP_ZIP ]; then

	echo "#[INFO] DBNSFP Download '$DBNSFP_ZIP' from '$DBNSFP_URL'"
	COMMAND_URL="curl -s $DBNSFP_URL -o $DBNSFP_ZIP.tmp >/dev/null"
	#(($DEBUG)) && echo $COMMAND_URL
	if eval $COMMAND_URL; then
		mv $DBNSFP_ZIP.tmp $DBNSFP_ZIP
	else
		echo "#[ERROR] DBNSFP Download '$DBNSFP_ZIP' from '$DBNSFP_URL' FAILED!!!"
		exit 0
	fi

else

	echo "#[INFO] DBNSFP Zip '$DBNSFP_ZIP' exists"

fi

if [ ! -e $DBNSFP_ZIP_README ]; then


	if [ "$DBNSFP_URL_README" != "" ]; then

		echo "#[INFO] DBNSFP Download '$DBNSFP_ZIP_README' from '$DBNSFP_URL_README'"
		#COMMAND_URL_README="curl $DBNSFP_URL_README -o $DBNSFP_ZIP_README.tmp"
		COMMAND_URL_README="wget -q '$DBNSFP_URL_README' -O $DBNSFP_ZIP_README.tmp >/dev/null"
		#(($DEBUG)) && echo $COMMAND_URL_README
		if eval $COMMAND_URL_README; then
			cat $DBNSFP_ZIP_README.tmp | sed 's/\r//g' > $DBNSFP_ZIP_README
			rm $DBNSFP_ZIP_README.tmp
		else
			echo "#[ERROR] DBNSFP Download '$DBNSFP_ZIP_README' from '$DBNSFP_URL_README' FAILED!!!"
			exit 0
		fi

	else

		echo "#[WARNING] DBNSFP Download '$DBNSFP_URL_README' failed because no URL"

	fi

else

	echo "#[INFO] DBNSFP Zip '$DBNSFP_ZIP_README' exists"

fi


###### FILE ######

### DEV
DBNSFP_ZIP_LIST_GZ_VARIANT_LIMIT=""
#DBNSFP_ZIP_LIST_GZ_VARIANT_LIMIT=" | head -n100000 "

if [ ! -e $DBNSFP_FILE ]; then

	echo "#[INFO] DBNSFP File creation '$DBNSFP_FILE' from '$DBNSFP_ZIP'"

	# Tmp file
	DBNSFP_FILE_TMP=$DBNSFP_FILE.tmp
	DBNSFP_FILE_TMP_HEAD=$DBNSFP_FILE_TMP.head
	DBNSFP_FILE_TMP_BODY=$DBNSFP_FILE_TMP.body
	> $DBNSFP_FILE_TMP
	> $DBNSFP_FILE_TMP_HEAD
	> $DBNSFP_FILE_TMP_BODY

	### List files
	COMMAND_ZIP_LIST="unzip -l $DBNSFP_ZIP | grep '_variant.chr.*gz$' | awk -F' ' '{print \$4}' > "$DBNSFP_FILE_TMP".files_list"
	#(($DEBUG)) && echo $COMMAND_ZIP_LIST
	if ! eval $COMMAND_ZIP_LIST; then
		echo "#[ERROR] DBNSFP Zip list '$DBNSFP_ZIP'"
		exit 0
	fi
	#(($DEBUG)) && echo "#[INFO] DBNSFP Zip list " && cat $DBNSFP_FILE".files_list" | sort -k1,1V
	#exit 0

	### Extract and recompose file
	
	

	# header
	dbnsfp_file_gz_first=$(cat $DBNSFP_FILE_TMP".files_list" | head -n1)
	#(($DEBUG)) && echo "dbnsfp_file_gz_first=$dbnsfp_file_gz_first"
	COMMAND_ZIP_LIST_GZ_FIRST="unzip -p $DBNSFP_ZIP $dbnsfp_file_gz_first | gzip -dc | head -n1 | sed 's/\r//g' >$DBNSFP_FILE_TMP_HEAD"
	#(($DEBUG)) && echo $COMMAND_ZIP_LIST_GZ_FIRST
	eval $COMMAND_ZIP_LIST_GZ_FIRST
	
	# content
	for dbnsfp_file_gz in $(cat $DBNSFP_FILE_TMP".files_list"); do
		#(($DEBUG)) && echo "dbnsfp_file_gz=$dbnsfp_file_gz" 
		COMMAND_ZIP_LIST_GZ="unzip -p $DBNSFP_ZIP $dbnsfp_file_gz | gzip -dc $DBNSFP_ZIP_LIST_GZ_VARIANT_LIMIT | sed 1d | sed 's/\r//g' >>$DBNSFP_FILE_TMP_BODY"
		#(($DEBUG)) && echo $COMMAND_ZIP_LIST_GZ
		eval $COMMAND_ZIP_LIST_GZ
	done

	cat $DBNSFP_FILE_TMP_HEAD > $DBNSFP_FILE_TMP
	#cat $DBNSFP_FILE_TMP_BODY | sort -k1,1V >> $DBNSFP_FILE_TMP
	cat $DBNSFP_FILE_TMP_BODY >> $DBNSFP_FILE_TMP


	if ! gzip -c $DBNSFP_FILE_TMP > $DBNSFP_FILE; then
		echo "#[ERROR] DBNSFP Annovar file creation '$DBNSFP_ANNOVAR' from '$DBNSFP_FILE' FAILED!!!"
	fi;

	rm -f $DBNSFP_FILE_TMP*

	#(($DEBUG)) && echo "#[INFO] DBNSFP File Tmp " && cat $DBNSFP_FILE


else

	echo "#[INFO] DBNSFP File '$DBNSFP_FILE' exists"

fi


###### ANNOVAR ######

### DEV
DBNSFP_ZIP_LIST_GZ_VARIANT_LIMIT=""
#DBNSFP_ANNOVAR_HEAD_LIMIT=" | cut -f1-20 "

if [ ! -e $DBNSFP_ANNOVAR ]; then

		echo "#[INFO] DBNSFP Annovar file creation '$DBNSFP_ANNOVAR' from '$DBNSFP_FILE'"

		#(($DEBUG)) && cat $DBNSFP_FILE | head -n2 | cut -f1-20 | column -t
		
		DBNSFP_ANNOVAR_TMP=$DBNSFP_ANNOVAR.tmp

		# out: #Chr	Start	End	Ref	Alt	reste...

		# Find columns: Chr	Start	End	Ref	Alt
		COLUMN_INDEX_CHR=1
		COLUMN_INDEX_START=2
		COLUMN_INDEX_END=2
		COLUMN_INDEX_REF=3
		COLUMN_INDEX_ALT=4

		ASSEMBLY_FOR_FILE=$ASSEMBLY

		# chr
		COMMAND_FILE_HEAD_CHR_INDEX="zcat $DBNSFP_FILE | head -n1 $DBNSFP_ANNOVAR_HEAD_LIMIT | tr '\t' '\n' | grep -n '^[#]*"$ASSEMBLY_FOR_FILE"[_]*chr$' | cut -d: -f1"
		COLUMN_INDEX_CHR_FOUND=$(eval $COMMAND_FILE_HEAD_CHR_INDEX)
		#(($DEBUG)) && echo "COLUMN_INDEX_CHR_FOUND=$COLUMN_INDEX_CHR_FOUND"
		[ "$COLUMN_INDEX_CHR_FOUND" != "" ] && COLUMN_INDEX_CHR=$COLUMN_INDEX_CHR_FOUND
		#(($DEBUG)) && echo "COLUMN_INDEX_CHR=$COLUMN_INDEX_CHR"

		# start
		COMMAND_FILE_HEAD_START_INDEX="zcat $DBNSFP_FILE | head -n1 $DBNSFP_ANNOVAR_HEAD_LIMIT | tr '\t' '\n' | grep -n '^[#]*"$ASSEMBLY_FOR_FILE"[_]*pos(1-based)$' | cut -d: -f1"
		COLUMN_INDEX_START_FOUND=$(eval $COMMAND_FILE_HEAD_START_INDEX)
		#(($DEBUG)) && echo "COLUMN_INDEX_START_FOUND=$COLUMN_INDEX_START_FOUND"
		[ "$COLUMN_INDEX_START_FOUND" != "" ] && COLUMN_INDEX_START=$COLUMN_INDEX_START_FOUND
		#(($DEBUG)) && echo "COLUMN_INDEX_START=$COLUMN_INDEX_START"

		# end
		COMMAND_FILE_HEAD_END_INDEX="zcat $DBNSFP_FILE | head -n1 $DBNSFP_ANNOVAR_HEAD_LIMIT | tr '\t' '\n' | grep -n '^[#]*"$ASSEMBLY_FOR_FILE"[_]*pos(1-based)$' | cut -d: -f1"
		COLUMN_INDEX_END_FOUND=$(eval $COMMAND_FILE_HEAD_END_INDEX)
		#(($DEBUG)) && echo "COLUMN_INDEX_END_FOUND=$COLUMN_INDEX_END_FOUND"
		[ "$COLUMN_INDEX_END_FOUND" != "" ] && COLUMN_INDEX_END=$COLUMN_INDEX_END_FOUND
		#(($DEBUG)) && echo "COLUMN_INDEX_END=$COLUMN_INDEX_END"

		# ref
		COMMAND_FILE_HEAD_REF_INDEX="zcat $DBNSFP_FILE | head -n1 $DBNSFP_ANNOVAR_HEAD_LIMIT | tr '\t' '\n' | grep -n '^ref$' | cut -d: -f1"
		COLUMN_INDEX_REF_FOUND=$(eval $COMMAND_FILE_HEAD_REF_INDEX)
		#(($DEBUG)) && echo "COLUMN_INDEX_REF_FOUND=$COLUMN_INDEX_REF_FOUND"
		[ "$COLUMN_INDEX_REF_FOUND" != "" ] && COLUMN_INDEX_REF=$COLUMN_INDEX_REF_FOUND
		#(($DEBUG)) && echo "COLUMN_INDEX_REF=$COLUMN_INDEX_REF"

		# alt
		COMMAND_FILE_HEAD_ALT_INDEX="zcat $DBNSFP_FILE | head -n1 $DBNSFP_ANNOVAR_HEAD_LIMIT | tr '\t' '\n' | grep -n '^alt$' | cut -d: -f1"
		COLUMN_INDEX_ALT_FOUND=$(eval $COMMAND_FILE_HEAD_ALT_INDEX)
		#(($DEBUG)) && echo "COLUMN_INDEX_ALT_FOUND=$COLUMN_INDEX_ALT_FOUND"
		[ "$COLUMN_INDEX_ALT_FOUND" != "" ] && COLUMN_INDEX_ALT=$COLUMN_INDEX_ALT_FOUND
		#(($DEBUG)) && echo "COLUMN_INDEX_ALT=$COLUMN_INDEX_ALT"

		COMMAND_ANNOVAR_HEAD="zcat $DBNSFP_FILE | head -n1 | sed 's/[\(\#\)\.]/_/gi' | tr '-' '_' | tr '+' '_' | awk '{print \"#Chr\tStart\tEnd\tRef\tAlt\t\"\$0}' > $DBNSFP_ANNOVAR_TMP"
		#(($DEBUG)) && echo $COMMAND_ANNOVAR_HEAD
		eval $COMMAND_ANNOVAR_HEAD
		COMMAND_ANNOVAR_BODY="zcat $DBNSFP_FILE | sed 1d | awk '{print \$$COLUMN_INDEX_CHR\"\t\"\$$COLUMN_INDEX_START\"\t\"\$$COLUMN_INDEX_END\"\t\"\$$COLUMN_INDEX_REF\"\t\"\$$COLUMN_INDEX_ALT\"\t\"\$0}' >> $DBNSFP_ANNOVAR_TMP"
		#(($DEBUG)) && echo $COMMAND_ANNOVAR_BODY
		eval $COMMAND_ANNOVAR_BODY


		#mv $DBNSFP_ANNOVAR_TMP $DBNSFP_ANNOVAR

		### annovar sort and idx
		if perl $SCRIPT_DIR/../index_annovar.pl --outfile $DBNSFP_ANNOVAR_TMP.sort_and_index $DBNSFP_ANNOVAR_TMP 2>/dev/null; then
			 if mv $DBNSFP_ANNOVAR_TMP.sort_and_index $DBNSFP_ANNOVAR && mv $DBNSFP_ANNOVAR_TMP.sort_and_index.idx $DBNSFP_ANNOVAR.idx; then
				(($DEBUG)) && echo "#[INFO] DBNSFP Annovar file sor and indexing done."
			 else
				echo "#[ERROR] DBNSFP Annovar file creation '$DBNSFP_ANNOVAR' from '$DBNSFP_FILE' FAILED!!!"
				exit 0
			 fi;
		else
			echo "#[ERROR] DBNSFP Annovar file creation '$DBNSFP_ANNOVAR' from '$DBNSFP_FILE' FAILED!!!"
			exit 0
		fi;

		# remove tmp
		rm -f $DBNSFP_ANNOVAR_TMP*

else

	echo "#[INFO] DBNSFP Annovar file '$DBNSFP_ANNOVAR' exists"

fi
#(($DEBUG)) && head -n100000 $DBNSFP_ANNOVAR


###### HOWARD ######

#DBNSFP_HOWARD_CONFIG_ANNOTATION=$DBNSFP_ANNOVAR.config.annotation.txt

# DEV
#rm -f $DBNSFP_HOWARD_CONFIG_ANNOTATION

if [ ! -e $DBNSFP_HOWARD_CONFIG_ANNOTATION ]; then

	echo "#[INFO] DBNSFP HOWARD Config Annotation ini file creation '$DBNSFP_HOWARD_CONFIG_ANNOTATION' from '$DBNSFP_ANNOVAR'"

	DBNSFP_HOWARD_CONFIG_ANNOTATION_TMP=$DBNSFP_HOWARD_CONFIG_ANNOTATION.tmp

	#$HOWARD_FOLDER/toolbox/dbnsfp_to_config_annotation.sh --input=$DBNSFP_ANNOVAR --output=$DBNSFP_HOWARD_CONFIG_ANNOTATION_TMP --input_readme=$DBNSFP_ZIP_README | head -n 50
	#exit 0

	ANNOVAR_CODE=$(basename $DBNSFP_ANNOVAR | sed 's/^'$ASSEMBLY'_//gi' | sed 's/\.txt$//gi')
	#echo $ANNOVAR_CODE
	#exit 0
	#$SCRIPT_DIR
	if ! $SCRIPT_DIR/dbnsfp_to_config_annotation.sh --input=$DBNSFP_ANNOVAR --output=$DBNSFP_HOWARD_CONFIG_ANNOTATION_TMP --input_name=$ANNOVAR_CODE --input_readme=$DBNSFP_ZIP_README >/dev/null ; then
		echo "#[ERROR] DBNSFP HOWARD Config Annotation ini file creation '$DBNSFP_HOWARD_CONFIG_ANNOTATION' from '$DBNSFP_ANNOVAR' FAILED!!!"
		exit 0	
	fi

	#head -n50 $DBNSFP_HOWARD_CONFIG_ANNOTATION_TMP | sed "s/^\[/[dbNSFP"$DBNSFP_RELEASE_FORMATTED"_custom_/gi"
	if ! cat $DBNSFP_HOWARD_CONFIG_ANNOTATION_TMP | sed "s/^\[/[dbNSFP"$DBNSFP_RELEASE_FORMATTED"_custom_/gi" > $DBNSFP_HOWARD_CONFIG_ANNOTATION; then
		echo "#[ERROR] DBNSFP HOWARD Config Annotation ini file creation '$DBNSFP_HOWARD_CONFIG_ANNOTATION' from '$DBNSFP_ANNOVAR' FAILED!!!"
		exit 0
	fi

	rm -f $DBNSFP_HOWARD_CONFIG_ANNOTATION_TMP*

else

	echo "#[INFO] DBNSFP Annovar file '$DBNSFP_HOWARD_CONFIG_ANNOTATION' exists"

fi


###### HOWARD test #######

if [ -e "$DBNSFP_HOWARD_CONFIG_ANNOTATION" ]; then
	
	TEST_FOLDER=$(dirname $DBNSFP_HOWARD_CONFIG_ANNOTATION)
	ANNOTATION="dbNSFP42a_custom__chr,dbNSFP42a_custom_ref,dbNSFP42a_custom_alt,dbNSFP42a_custom_clinvar_clnsig"
	ANNOTATION="ALL"

	$HOWARD --input=$HOWARD_FOLDER/docs/example.vcf --output=$TEST_FOLDER/HOWARD.validation.vcf --env=/tool/config/tools.app --annotation=$ANNOTATION --annovar_folder=$ANNOVAR --annovar_databases=$TEST_FOLDER --config_annotation=$DBNSFP_HOWARD_CONFIG_ANNOTATION --snpeff_jar=$SNPEFF --snpeff_threads=3 --tmp=/tmp --assembly=hg19 --java=$JAVA --threads=3 --verbose

fi;




###### HOWARD ######

#echo $HOWARD

#echo $ENV_CONFIG



exit 0

















# CONFIG_ANNOTATION
if [ "$CONFIG_ANNOTATION" == "" ] || [ ! -e "$CONFIG_ANNOTATION" ] && ! (($CHECK_UPDATE)); then
    echo "#[ERROR] No Config Annotation file '$CONFIG_ANNOTATION'"
    exit 0
fi;
echo "#[INFO] Config Annotation file '$CONFIG_ANNOTATION'"


# ANNOTATION
if [ "$ANNOTATION" == "" ]; then
    ANNOTATION="ALL"
fi;
echo "#[INFO] Annotation '$ANNOTATION'"


# STARK_MAIN_ANNOTATION
if [ "$STARK_MAIN_FOLDER" == "" ]; then
	STARK_MAIN_FOLDER=$HOME/STARK
fi;
if [ "$STARK_MAIN_FOLDER" == "" ] || [ ! -d "$STARK_MAIN_FOLDER" ]; then
	echo "#[ERROR] No STARK main folder '$STARK_MAIN_FOLDER'"
    exit 0
fi;
echo "#[INFO] STARK main folder '$STARK_MAIN_FOLDER'"


# RELEASE
if [ "$ANNOVAR_NEW_RELEASE" == "" ]; then
	ANNOVAR_NEW_RELEASE=$(date '+%Y%m%d-%H%M%S')
fi;
[ "$ANNOVAR_NEW_RELEASE" == ""  ] && ANNOVAR_NEW_RELEASE=$(date '+%Y%m%d-%H%M%S') 
echo "#[INFO] Release '$ANNOVAR_NEW_RELEASE'"


# RELEASE LATEST
if [ -e $STARK_MAIN_FOLDER/databases/annovar/latest/ ]; then
	ANNOVAR_LATEST_RELEASE=$(basename $(realpath $STARK_MAIN_FOLDER/databases/annovar/latest/)) 
else
	echo "#[INFO] No latest release "
fi;


# FOLDERS
TMP_ID=$RANDOM$ANDOM

# HOST
TMP_FOLDER=$STARK_MAIN_FOLDER/output/tmp/$TMP_ID
HOWARD_UPDATE_RELEASE_FOLDER=$TMP_FOLDER/$ANNOVAR_NEW_RELEASE
HOWARD_UPDATE_RELEASE_FOLDER_DATABASES=$HOWARD_UPDATE_RELEASE_FOLDER/databases
HOWARD_UPDATE_RELEASE_FOLDER_TMP=$HOWARD_UPDATE_RELEASE_FOLDER/tmp

# INNER
STARK_MAIN_FOLDER_INNER=/STARK
TMP_FOLDER_INNER=$STARK_MAIN_FOLDER_INNER/output/tmp/$TMP_ID
HOWARD_UPDATE_RELEASE_FOLDER_INNER=$TMP_FOLDER_INNER/$ANNOVAR_NEW_RELEASE
HOWARD_UPDATE_RELEASE_FOLDER_DATABASES_INNER=$HOWARD_UPDATE_RELEASE_FOLDER_INNER/databases
HOWARD_UPDATE_RELEASE_FOLDER_TMP_INNER=$HOWARD_UPDATE_RELEASE_FOLDER_INNER/tmp


# DEBUG
(($DEBUG)) && echo "#[INFO] ANNOVAR_LATEST_RELEASE=$ANNOVAR_LATEST_RELEASE"
(($DEBUG)) && echo "#[INFO] HOWARD_UPDATE_RELEASE_FOLDER=$HOWARD_UPDATE_RELEASE_FOLDER"
(($DEBUG)) && echo "#[INFO] HOWARD_UPDATE_RELEASE_FOLDER_DATABASES=$HOWARD_UPDATE_RELEASE_FOLDER_DATABASES"
(($DEBUG)) && echo "#[INFO] HOWARD_UPDATE_RELEASE_FOLDER_TMP=$HOWARD_UPDATE_RELEASE_FOLDER_TMP"


# CHECK UPDATE
if (($CHECK_UPDATE)); then

	if [ -e $STARK_MAIN_FOLDER/databases/annovar/latest/hg19_avdblist.txt ]; then

		mkdir -p $TMP_FOLDER

		echo "#[INFO] Check for update..."
		docker exec -ti stark-module-stark-submodule-stark-service-cli bash -c 'source /tool/config/config.app && perl $ANNOVAR/annotate_variation.pl -webfrom annovar -downdb avdblist -buildver hg19 /tmp && grep "^NOTICE" -v /tmp/hg19_avdblist.txt ' > $TMP_FOLDER/hg19_avdblist.txt
		checksum_new=$(grep "^NOTICE" -v  $TMP_FOLDER/hg19_avdblist.txt | sort -u | sha256sum | cut -d" " -f1 )
		checksum_latest=$(grep "^NOTICE" -v  $STARK_MAIN_FOLDER/databases/annovar/latest/hg19_avdblist.txt | sort -u | sha256sum | cut -d" " -f1 )
		if [ "$checksum_new" != "$checksum_latest" ]; then
			echo "#[INFO] Update needed"
			diff -E <(grep "^NOTICE" -v $TMP_FOLDER/hg19_avdblist.txt | sort -u) <(grep "^NOTICE" -v $STARK_MAIN_FOLDER/databases/annovar/latest/hg19_avdblist.txt | sort -u) | grep "^<" | column -t
		else
			echo "#[INFO] NO update needed"
		fi
	
		rm -rf $TMP_FOLDER

	else 

		echo "#[INFO] No latest release. Can not check update "

	fi

	exit 0

fi


if ((1)); then

	mkdir -p $HOWARD_UPDATE_RELEASE_FOLDER $HOWARD_UPDATE_RELEASE_FOLDER_DATABASES $HOWARD_UPDATE_RELEASE_FOLDER_TMP

	# Fichier de configuration d'HOWARD
	# Créer et copier le fichier de configuration config.annotation.ini
	cp $CONFIG_ANNOTATION $HOWARD_UPDATE_RELEASE_FOLDER_DATABASES/config.annotation.ini


	# Téléchargement (par HOWARD dans STARK CLI)
	# Dans STARK CLI (pour l'acces a HOWARD et aux paramétrages de STARK-tools)
	docker exec -ti stark-module-stark-submodule-stark-service-cli bash -c 'cd /STARK && source /tool/config/config.app && $HOWARD --input=$HOWARD_FOLDER/docs/example.vcf --output='$HOWARD_UPDATE_RELEASE_FOLDER_DATABASES_INNER'/HOWARD.download.vcf --env=/tool/config/tools.app --annotation='$ANNOTATION' --annovar_folder=$ANNOVAR --annovar_databases='$HOWARD_UPDATE_RELEASE_FOLDER_DATABASES_INNER' --config_annotation='$HOWARD_UPDATE_RELEASE_FOLDER_DATABASES_INNER'/config.annotation.ini --snpeff_jar=$SNPEFF --snpeff_threads=1 --tmp='$HOWARD_UPDATE_RELEASE_FOLDER_TMP_INNER' --assembly=hg19 --java=$JAVA --threads=1 --verbose' 1>$HOWARD_UPDATE_RELEASE_FOLDER_DATABASES/HOWARD.download.log 2>$HOWARD_UPDATE_RELEASE_FOLDER_DATABASES/HOWARD.download.err

	# Telechargement de la liste des bases de données 'hg19_avdblist.txt'
	#docker exec -ti stark-module-stark-submodule-stark-service-cli bash -c 'cd /STARK && source /tool/config/config.app && perl $ANNOVAR/annotate_variation.pl -webfrom annovar -downdb avdblist -buildver hg19 '$HOWARD_UPDATE_RELEASE_FOLDER_DATABASES_INNER' '
	docker exec -ti stark-module-stark-submodule-stark-service-cli bash -c 'source /tool/config/config.app && perl $ANNOVAR/annotate_variation.pl -webfrom annovar -downdb avdblist -buildver hg19 /tmp && grep "^NOTICE" -v /tmp/hg19_avdblist.txt ' > $HOWARD_UPDATE_RELEASE_FOLDER_DATABASES/hg19_avdblist.txt
	

	# Update (sur le serveur)
	cd $STARK_MAIN_FOLDER
	# créer la version d'annovar
	mkdir -p $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_NEW_RELEASE
	# Copier la nouvelle version d'ANNOVAR dans database
	rsync -rauz $HOWARD_UPDATE_RELEASE_FOLDER_DATABASES/* $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_NEW_RELEASE/
	# Ajouter les fichiers manquants de la précédente version d'ANNOVAR dans database
	if [ "$ANNOVAR_LATEST_RELEASE" != "" ]; then
		rsync -rauz $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_LATEST_RELEASE/* $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_NEW_RELEASE/
	fi;
	# modifier le fichier STARK.database.release au besoin ("download")
	#mv $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_NEW_RELEASE/STARK.database.release $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_NEW_RELEASE/STARK.database.release.previous
	#sed "s/$ANNOVAR_LATEST_RELEASE/$ANNOVAR_NEW_RELEASE/gi" $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_NEW_RELEASE/STARK.database.release.previous > $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_NEW_RELEASE/STARK.database.release
	echo '{
			"release": "'$ANNOVAR_NEW_RELEASE'",
			"date": "'$ANNOVAR_NEW_RELEASE'",
			"files": [ "config.annotation.ini" ],
			"assembly": [ "hg19" ],
			"download": {
					"methode": "'$SCRIPT_DESCRIPTION' ['$SCRIPT_RELEASE'-'$SCRIPT_DATE']",
					"info": "From ANNOVAR binary",
					"databases": "'$ANNOTATION'",
					"configuration": "config.annotation.ini"
			}
	}' > $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_NEW_RELEASE/STARK.database.release


	# Mise a disposition
	# Copier le fichier de configuration dans
	cp $HOWARD_UPDATE_RELEASE_FOLDER_DATABASES/config.annotation.ini $STARK_MAIN_FOLDER/config/howard/config.annotation.$ANNOVAR_NEW_RELEASE.ini
	# Changer latest
	if [ "$ANNOVAR_LATEST_RELEASE" != "" ]; then
		rm -f $STARK_MAIN_FOLDER/databases/annovar/latest
		ln -snf $ANNOVAR_NEW_RELEASE/ $STARK_MAIN_FOLDER/databases/annovar/latest
	fi


fi


### Production
echo "
### Manually execute these command to push ANNOVAR database into production
### !!! possibly destructive !!!
# change current release
rm -f $STARK_MAIN_FOLDER/databases/annovar/current
ln -snf $ANNOVAR_NEW_RELEASE $STARK_MAIN_FOLDER/databases/annovar/current
# change HOWARD configuration
cp $STARK_MAIN_FOLDER/config/howard/config.annotation.ini $STARK_MAIN_FOLDER/config/howard/config.annotation.ini.from_update_$ANNOVAR_NEW_RELEASE
cp $STARK_MAIN_FOLDER/config/howard/config.annotation.$ANNOVAR_NEW_RELEASE.ini $STARK_MAIN_FOLDER/config/howard/config.annotation.ini
"
if [ "$ANNOVAR_LATEST_RELEASE" != "" ]; then
echo "
# Compress previous release
(cd databases/annovar && tar -vczf $ANNOVAR_LATEST_RELEASE.tar.gz $ANNOVAR_LATEST_RELEASE/)
rm -rf databases/annovar/$ANNOVAR_LATEST_RELEASE
"
fi;


rm -rf $TMP_FOLDER

exit 0

