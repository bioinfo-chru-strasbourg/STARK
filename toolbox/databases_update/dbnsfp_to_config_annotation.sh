#!/bin/bash
#################################
## HOWARD 
#################################

SCRIPT_NAME="HOWARD_DBNSFP"
SCRIPT_DESCRIPTION="HOWARD DBNSFPx to config annotation ini file"
SCRIPT_RELEASE="0.9.1"
SCRIPT_DATE="10/01/2019"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="HUS"
SCRIPT_LICENCE="GNU AGPL V3"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-07/11/2017:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tScript creation\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.1-10/01/2019:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tHelp and code cleaning\n";


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
	echo "# RELEASE NOTES:";
	echo -e $RELEASE_NOTES
}

# Usage
function usage {
	echo "# USAGE: $(basename $0) --input=<FILE>  [options...]";
	echo "# Following options are available:";
	echo "# --input=<FILE>             Input DBNSPF file";
	echo "# --output=<FILE>            Output annotation configuration file";
	echo "# --input_name=<STRING>      Annovar code of the input DBNSFP file";
	echo "# --input_release=<STRING>   Release of the input DBNSFP file";
	echo "# --input_date=<STRING>      Date of the input DBNSFP file";
	echo "# --input_readme=<STRING>    Readme of the input DBNSFP file";
	echo "# --input_readme_url=<URL>   Readme URL of the input DBNSFP file";
	echo "# --env=<FILE>               Environment configuration for multithreading (BGZIP, TABIX, BCFTOOLS, VCFTOOLS)";
	#echo "# --force                    Force annotation even if already exists in VCF header";
	echo "# --tmp=<FOLDER>             Temporary folder (default /tmp)";
	echo "# --verbose                  VERBOSE option";
	echo "# --debug                    DEBUG option";
	echo "# --release                  RELEASE option";
	echo "# --help                     HELP option";
	echo "# Example: $(basename $0) --input_readme_url='https://drive.google.com/uc?export=download&id=1Vse3b_qw_E46eDcsuLegF5HxpodAZ6Og'"
	echo "#";
	echo "";

}



# EXAMPLE :
# ./dbnsfp_to_config_annotation.sh --input=dbnsfp33a.txt --output=config.annotation.dbnsfp33a.ini

# header
header;


ARGS=$(getopt -o "i:o:e:a:f:s:r:xt:m:vdnh" --long "input:,output:,input_name:,input_release:,input_date:,input_readme:,input_readme_url:,env:,tmp:,verbose,debug,release,help" -- "$@" 2> /dev/null)
#ARGS=$(getopt --long "input:,output:,annotation:,multithreading,threads:,verbose,debug,release,help" -- "$@" 2> /dev/null)
if [ $? -ne 0 ]; then
	:
	#echo $?
	#usage;
	#exit;
fi;

PARAM=$@

eval set -- "$ARGS"
while true
do
	#echo "$1=$2"
	#echo "Eval opts";
	case "$1" in
		--input)
			if [ ! -e $2 ] || [ "$2" == "" ]; then
				echo "#[ERROR] No file '$2'"
				usage; exit;
			else
				INPUT="$2";
			fi;
			shift 2
			;;
		--output)
			OUTPUT="$2";
			shift 2
			;;
		--input_name)
			INPUT_NAME="$2";
			shift 2
			;;
		--input_release)
			INPUT_RELEASE="$2";
			shift 2
			;;
		--input_date)
			INPUT_DATE="$2";
			shift 2
			;;
		--input_readme)
			INPUT_README="$2";
			shift 2
			;;
		--input_readme_url)
			INPUT_README_URL="$2";
			shift 2
			;;
		--env)
			ENV="$2";
			shift 2
			;;
		--tmp)
			TMP_INPUT="$2"
			shift 2
			;;
		--verbose)
			VERBOSE=1
			shift 1
			;;
		--debug)
			VERBOSE=1
			DEBUG=1
			shift 1
			;;
		--release)
			release;
			exit 0
			;;
		--help)
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
if [ ! -z $ENV ] && [ -s $ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$ENV;
	echo "#[INFO] ENV '$ENV' found."
elif [ -s $SCRIPT_DIR/$ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$SCRIPT_DIR/$ENV;
	echo "#[INFO] ENV '$ENV' found."
elif [ "$ENV" == "" ] || [ ! -s $ENV ]; then
	if [ -s $SCRIPT_DIR/"env.sh" ]; then
		ENV=$SCRIPT_DIR/"env.sh";
		echo "#[INFO] Default ENV '$ENV' used."
	else
		ENV="";
		echo "#[WARNING] NO ENV defined. No ENV used."
	fi;
fi;
if [ -e $ENV ] && [ "$ENV" != "" ]; then
	source $ENV
fi;


if [ "$INPUT_README_URL" == "" ]; then
	[ "$INPUT_RELEASE" == "4.2a" ] && INPUT_README_URL="https://drive.google.com/uc?export=download&id=1Vse3b_qw_E46eDcsuLegF5HxpodAZ6Og";
fi;

if [ "$INPUT_README" == "" ]; then
	if [ "$INPUT_README_URL" != "" ]; then
		INPUT_README=/tmp/dbnsfp.readme
		echo "#[INFO] DBNSFP Download '$INPUT_README' from '$INPUT_README_URL'"
		#COMMAND_URL_README="curl $INPUT_README_URL -o $DBNSFP_ZIP_README.tmp"
		COMMAND_URL_README="wget -q '$INPUT_README_URL' -O $INPUT_README.tmp >/dev/null"
		#(($DEBUG)) && echo $COMMAND_URL_README
		if eval $COMMAND_URL_README; then
			cat $INPUT_README.tmp | sed 's/\r//g' > $INPUT_README
			rm $INPUT_README.tmp
		else
			echo "#[ERROR] DBNSFP Download '$INPUT_README' from '$INPUT_README_URL' FAILED!!!"
			exit 0
		fi
	fi;
fi

if [ -z $INPUT_RELEASE ]; then #  || [ -z $OUTPUT ]; then
	if [ "$INPUT" != "" ] && [ -z $INPUT ]; then
		INPUT_RELEASE=$(basename $INPUT)
	fi;
	if [ -e $INPUT_README ] && [ "$INPUT_README" != "" ]; then
		INPUT_RELEASE=$(cat $INPUT_README | sed 's/\r$//' | grep "^dbNSFP version"  -m1 |cut -d" " -f3)
		echo "#[INFO] RELEASE '$INPUT_RELEASE' found in README file '$INPUT_README'.";
	fi;
	if [ "$INPUT_RELEASE" == "" ] || [ -z $INPUT_RELEASE ]; then
		INPUT_RELEASE="4.2a"
	fi;
	echo "#[WARNING] INPUT RELEASE Missing. default '$INPUT_RELEASE' used";
fi


# Mandatory parameters
if [ -z $INPUT ]; then #  || [ -z $OUTPUT ]; then
	if [ "$INPUT_RELEASE" != "" ] && [ ! -z $INPUT_RELEASE ]; then
		INPUT=/tmp/dbnsfp$(echo $INPUT_RELEASE | tr -d ".")
		curl http://www.openbioinformatics.org/annovar/download/hg19_$(basename $INPUT).txt.gz 2>/dev/null | zcat | head -n1 > $INPUT
	fi;
	if [ -z $INPUT ] || [ "$INPUT" == "" ]; then
		echo "#[ERROR] INPUT Missing";
		usage;
		exit;
	fi;
fi

if [ -z $INPUT_DATE ]; then #  || [ -z $OUTPUT ]; then
	INPUT_DATE=$(basename $INPUT)
	if [ -e $INPUT_README ] && [ "$INPUT_README" != "" ]; then
		INPUT_DATE=$(cat $INPUT_README | sed 's/\r$//' | grep "^Release:" -A1 | tail -n1 | tr -d "\t")
		echo "#[INFO] DATE '$INPUT_DATE' found in README file '$INPUT_README'.";
	fi;
	echo "#[WARNING] INPUT DATE Missing. '$INPUT_DATE' used";
fi

if [ -z $INPUT_NAME ]; then #  || [ -z $OUTPUT ]; then
	INPUT_NAME=$(basename $INPUT)
	if [ -e $INPUT_README ] && [ "$INPUT_README" != "" ]; then
		INPUT_NAME="dbnsfp"$(echo $INPUT_RELEASE | tr -d ".")
	fi;
	echo "#[WARNING] INPUT NAME Missing. default '$INPUT_NAME' used";
fi

echo "#[INFO] INPUT_NAME=$INPUT_NAME";
echo "#[INFO] INPUT_DATE=$INPUT_DATE";
echo "#[INFO] INPUT_RELEASE=$INPUT_RELEASE";


# Mandatory parameters
if [ -z $OUTPUT ] || [ "$OUTPUT" == "" ]; then #  || [ -z $OUTPUT ]; then
	echo "#[WARNING] OUTPUT Missing";
	#OUTPUT=$(echo $INPUT | sed "s/.vcf$/.annotated.vcf/g"); #/\.vcf$/\.output\.vcf/;
	OUTPUT=config.annotation.$(basename $INPUT).ini; #/\.vcf$/\.output\.vcf/;
fi


# TMP
if [ ! -z $TMP_INPUT ] && [ ! -d $TMP_INPUT ]; then
	mkdir -p $TMP_INPUT;
fi;
if [ ! -z $TMP_INPUT ] && [ -d $TMP_INPUT ]; then
	TMP_SYS_FOLDER=$TMP_INPUT;
fi;
if [ -z $TMP_SYS_FOLDER ] || [ ! -d $TMP_SYS_FOLDER ]; then
	TMP_SYS_FOLDER=/tmp
	(($VERBOSE)) && echo "#[WARNING] TMP=$TMP_SYS_FOLDER"
	#usage;
	#exit;
fi



# HEADER
echo "#[INFO] INPUT=$INPUT"
echo "#[INFO] OUTPUT=$OUTPUT"
(($VERBOSE)) && echo "#[INFO] LOG=$LOG";
(($VERBOSE)) && echo "#[INFO] ERR=$ERR";
#echo "##################"

> $OUTPUT

line=6
for ANNOTATION in $(head -n1 $INPUT | cut -f6- | tr "\t" "\n"); do
	annotation_type="annotation"
	((line++))
	otherinfo=$(($line-6))
	if [[ $ANNOTATION =~ ^.*_score$ ]] || [[ $ANNOTATION =~ ^.*_rankscore$ ]] || [[ $ANNOTATION =~ ^.*_phred$ ]]; then
		annotation_type="score"
	elif [[ $ANNOTATION =~ ^.*_pred$ ]]; then
		annotation_type="prediction"
	fi
	if [ -e $INPUT_README ] && [ "$INPUT_README" != "" ]; then
		DESCRIPTION=$(cat $INPUT_README | sed 's/\r$//' | sed ':a;N;$!ba;s/\n\t\t/ /g' | grep "$ANNOTATION: " | cut -d":" -f2- | sed "s/^ //" | tr -d "\"';")
	else
		DESCRIPTION=$ANNOTATION
	fi;
	
	echo "#[INFO] ANNOTATION '$ANNOTATION'";
	
	ANNOTATION_OUTPUT="[$ANNOTATION]
annovar_code=$INPUT_NAME
annovar_annotation_type=filter
release=$INPUT_RELEASE
available=true
core=false
annotation_type=$annotation_type
date=$INPUT_DATE
otherinfo=$otherinfo
description=$DESCRIPTION

"
	
	(($VERBOSE)) && echo "#[INFO] $ANNOTATION_OUTPUT"
	echo "$ANNOTATION_OUTPUT" >> $OUTPUT


	
done;


