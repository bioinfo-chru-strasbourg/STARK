#!/bin/bash
#################################
##
## NGS environment
## Analysis of multiple RUNS
##
## author: Antony Le Bechec
##
#################################

SCRIPT_NAME="Launch synthetic"
SCRIPT_DESCRIPTION="Analysis a synthetic FASTQ"
SCRIPT_RELEASE="0.9d"
SCRIPT_DATE="15/07/2015"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"

echo "#######################################";
echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
echo "# $SCRIPT_DESCRIPTION ";
echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © ";
echo "#######################################";

# Realse note
RELEASE_NOTES="# Release Notes\n"
RELEASE_NOTES=$RELEASE_NOTES"0.9d-15/07/2015: Creation of the script\n"


######################################
# 1. Configuration, Input parameters #
######################################

# 1.1. Configuration
#####################

####################################################################################################################################
# Define the function to print the usage of the script
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    sh launch.validation.synthetic.sh -f1 fastq_R1 -f2 fastq_R2 -b bed -e env [-h]

		Options:
		 	-e, --env
		 	ENV file configuration.
		 	-f1, --fastq_R1
		 	This arg is optional. FASTQ file in .fastq.gz format OR .bam
		 	-f2 --fastq_R2
		 	This arg is optional. FASTQ Read2 file in .fastq.gz format (first FASTQ parameter will be considered as FASTQ Read1 file).
		 	-b --bed
		 	BED (optional): BED file in BED format (no header needed)
		  	-h, --help
		 	Print this message and exit the program.
		__EOF__
}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "e:f1:f2:b:h" --long "env:,fastq_R1:,fastq_R2:,bed:,help" -- "$@" 2> /dev/null)
[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && \
	exit 1
eval set -- "$ARGS"
while true
do
	case "$1" in
		-e|--env)
			ENV="$2"
			shift 2 
			;;
		-f1|--fastq_R1)
			FASTQ_R1="$2"
			shift 2 
			;;
		-f2|--fastq_R2)
			FASTQ_R2="$2"
			shift 2 
			;;
		-b|--bed)
			BED="$2"
			shift 2 
			;;
		-h|--help)
			usage
			exit 0
			;;
		--) shift
			break 
			;;
		*) 	echo "Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# INPUT

SYNTHETIC_DATA_FOLDER=validation/fastq

# FASTQ_R1
#FASTQ_R1_DEFAULT="synthetic.R1.fastq"
FASTQ_R1_DEFAULT="$SYNTHETIC_DATA_FOLDER/synthetic.unaligned.bam"
#FASTQ_R1_DEFAULT="$SYNTHETIC_DATA_FOLDER/synthetic.R1.fastq"
if [ "$1" != "" ] && [ -e $1 ]; then
	FASTQ_R1=$1;
elif [ -s $FASTQ_R1_DEFAULT ]; then
	FASTQ_R1=$FASTQ_R1_DEFAULT;
else
	echo "[ERROR] No input FASTQ_R1 file!";
	exit 0;	
fi;


# FASTQ_R2 
FASTQ_R2=$2
if [ ! -e $FASTQ_R2 ]; then
	FASTQ_R2="";
fi;

# BED
BED=$3
BED_DEFAULT="$SYNTHETIC_DATA_FOLDER/synthetic.bed"
ls $BED_DEFAULT;
if [ "$BED" != "" ] && [ -e $BED ]; then
	BED=$BED
elif [ -s $BED_DEFAULT ]; then
	BED=$BED_DEFAULT;
	echo "[WARNING] No BED file. Default BED file '$BED_DEFAULT' will be used";
else
	BED="";
	echo "[WARNING] No BED file and NO Default BED file! No BED file will be used";
fi;


# ENV
ENV=$4
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
if [ "$ENV" != "" ]; then
	source $ENV;
fi;

echo $ENV
echo $RESULTS_FOLDER


#CORES=$(echo `ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w`" - 1 " | bc)
#THREADS=$CORES
THREADS=1

#PIPELINES="bwamem.gatkHC.howard bwamemUnclipped.gatkHC.howard bwamem.gatkUG.howard bwamemUnclipped.gatkUG.howard"
PIPELINES="bwamem.gatkHC.howard bwamem.gatkUG.howard"

$STARK/launch.sample.sh -e "$ENV" -f "$FASTQ_R1" -q "$FASTQ_R2" -b "$BED" -s "SAMPLE" -r "SYNTHETIC" -p "$PIPELINES"

#echo "make -j $THREADS -e NGSEnv=$NGS_FOLDER ENV=$ENV PARAM=$MAKEFILE_ANALYSIS VALIDATION=1 OUTDIR=$RES_SOURCE -f $NGS_SCRIPTS/NGSWorkflow.mk 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA"
#make -j $THREADS -e NGSEnv=$NGS_FOLDER ENV=$ENV PARAM=$MAKEFILE_ANALYSIS VALIDATION=1 OUTDIR=$RES_SOURCE -f $NGS_SCRIPTS/NGSWorkflow.mk 1>>$LOGFILE_ANA 2>>$LOGFILE_ANA







