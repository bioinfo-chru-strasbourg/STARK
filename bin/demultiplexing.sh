#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="IlluminaRunDemultiplexing"
SCRIPT_DESCRIPTION="Demultiplexing of an Illumina Run"
SCRIPT_RELEASE="0.9.9b"
SCRIPT_DATE="26/09/2016"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"

echo "#######################################";
echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
echo "# $SCRIPT_DESCRIPTION ";
echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © ";
echo "#######################################";

# Realse note
#RELEASE_NOTES="#\n"
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.7b-20/11/2015: Creation/ReFormating script\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.8b-26/11/2015: SampleSheet as input option. Add Manifest default folder to find Manifest file defined in SampleSheet\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.9b-26/09/2016: Cleaning Sample Name in the Sample Sheet by removing spaces.\n";

# NOTES
if [ "${1^^}" == "RELEASE" ] || [ "${1^^}" == "RELEASE_NOTES" ] || [ "${1^^}" == "NOTES" ]; then
	echo "# RELEASE NOTES:";
	echo -e $RELEASE_NOTES
	exit 0;
fi;


# 1.1. Parameters
##################

NB_PARAMETERS=$#
if [ $NB_PARAMETERS -eq 0 ] || [ "${1^^}" == "HELP" ]; then
	echo "# Usage: "`basename $0`" <ENV> <RUN> <DEMULTIPLEXING_DIR> <DEMULTIPLEXING> <SAMPLESHEET> ";
	echo "# ENV                   ENV file configuration";
	echo "# RUN                   RUN to demultiplex";
	echo "# DEMULTIPLEXING_DIR    Demultiplexing result Folder";
	echo "# DEMULTIPLEXING        Launch demultiplexing (not available)";
	echo "# SAMPLESHEET           SampleSheet to use";
	exit 1;
fi;

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# INPUT
ENV=$1
RUN=$2
DEMULTIPLEXING_DIR=$3
DEMULTIPLEXING=$4
SAMPLESHEET=$5

# SET ENV
#SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#source $SCRIPT_DIR/env.sh

# Input

if [ "$DEMULTIPLEXING_DIR" == "" ];
then 
	DEMULTIPLEXING_DIR=$DEMULTIPLEXING_FOLDER
fi;

if [ "$DEMULTIPLEXING" == "" ];
then 
	DEMULTIPLEXING=0
fi;

# ENV
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

# ENV
#source $SCRIPT_DIR/env.sh
if [ "$ENV" != "" ]; then
	source $ENV;
fi;


# Config
NGS_DIR=$NGS_FOLDER
NGS_BIN_DIR=$NGS_DIR/bin
#CASAVA=$NGS_BIN_DIR/casava
#BCL2FASTQ=$NGS_BIN_DIR/bcl2fastq
#BCL2FASTQ=$CASAVA_BCLTOFASTQ #$BCL2FASTQ_BCL2FASTQ
BCL2FASTQ=$BCL2FASTQ_BCLTOFASTQ
#BCL2FASTQ=$CASAVA_BCLTOFASTQ
#MISEQ_DIR=/media/miseq/MSR # /media/IRC/RAW/MSR
#MISEQ_DIR=/media/IRC/RAW/MSR
MISEQ_DIR=$MISEQ_FOLDER
#DEMULTIPLEXING_DIR=/media/miseq/demultiplexing.test # /mediSampleSheeta/IRC/RES/...
#DEMULTIPLEXING_DIR=/media/IRC/RES2/demultiplexing
#DEMULTIPLEXING_DIR=$DEMULTIPLEXING_FOLDER
MSR_SUBFOLDER=Data/Intensities/BaseCalls
LOGFILE=$DEMULTIPLEXING_DIR/$RUN/demultiplexing.configuration.log
CORES=$(ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w ) # AUTO, otherwize #THREAD=12
THREADS=$CORES
CURRENT_DIR=`pwd`

# Input Test
if [ ! -d $MISEQ_DIR/$RUN ]; then
	echo "#[`date`] ERROR: RUN '$RUN' missing! Folder '$MISEQ_DIR/$RUN' doesn't exist"
	exit 1;
fi

SAMPLE_SHEET_MISEQ_ORIGINAL_DEFAULT=$MISEQ_DIR/$RUN/SampleSheet.csv
#SAMPLE_SHEET_MISEQ_ORIGINAL=$MISEQ_DIR/$RUN/SampleSheet.csv
if [ -s $SAMPLESHEET ]; then
	SAMPLE_SHEET_MISEQ_ORIGINAL=$SAMPLESHEET;
elif [ -s $SAMPLE_SHEET_MISEQ_ORIGINAL_DEFAULT ]; then
	SAMPLE_SHEET_MISEQ_ORIGINAL=$SAMPLE_SHEET_MISEQ_ORIGINAL_DEFAULT;
fi;

if [ ! -e $SAMPLE_SHEET_MISEQ_ORIGINAL ]; then
	echo "#[`date`] ERROR: SampleSheet.csv for RUN '$RUN' missing! File '$SAMPLE_SHEET_MISEQ_ORIGINAL' doesn't exist"
	exit 1;
fi


# Main
INPUT_DIR=$MISEQ_DIR/$RUN/$MSR_SUBFOLDER
RUNFOLDER_DIR=$MISEQ_DIR/$RUN
OUTPUT_DIR=$DEMULTIPLEXING_DIR/$RUN

mkdir -p $OUTPUT_DIR

# Create Sample Sheet, and options for BCL2FASTQ
#SAMPLE_SHEET_MISEQ_ORIGINAL=$MISEQ_DIR/$RUN/SampleSheet.csv
SAMPLE_SHEET_MISEQ=$OUTPUT_DIR/SampleSheet.csv
#SAMPLE_SHEET_CASAVA=$OUTPUT_DIR/SampleSheet.casava.csv
SAMPLE_SHEET_BCL2FASTQ=$OUTPUT_DIR/SampleSheet.bcl2fastq.csv
MANIFESTS=$OUTPUT_DIR/manifests.txt
MANIFESTS_LIST=$OUTPUT_DIR/manifests_list.txt
INTERVALS=$OUTPUT_DIR/run.intervals
READS_LENGTH=$OUTPUT_DIR/readsLength.txt
MASK=$OUTPUT_DIR/mask.txt
ADAPTERS=$OUTPUT_DIR/adapters.txt
MAKEFILE=$OUTPUT_DIR/Makefile

echo "#[`date`] Demultiplexing Configuration log" >> $LOGFILE

cp $SAMPLE_SHEET_MISEQ_ORIGINAL $SAMPLE_SHEET_MISEQ

# Cleaning SampleSheet SampleName
DATA_SECTION_LINE=$(grep "^\[Data\]" $SAMPLE_SHEET_MISEQ_ORIGINAL -n | awk -F: '{print $1}')
SAMPLE_SECTION_FIRST_LINE=$(($DATA_SECTION_LINE+1))
# Print before DATA section
head -n $SAMPLE_SECTION_FIRST_LINE $SAMPLE_SHEET_MISEQ_ORIGINAL > $SAMPLE_SHEET_MISEQ
# print DATA section and the rest
tail -n $(($SAMPLE_SECTION_FIRST_LINE-$(grep ^ -c $SAMPLE_SHEET_MISEQ_ORIGINAL)))  $SAMPLE_SHEET_MISEQ_ORIGINAL | while read L; do
  	if [ "$L" == "" ] || [[ $L =~ ^\[.*\] ]]; then
  		echo $L >> $SAMPLE_SHEET_MISEQ;
		#break;
	else
		echo $L | tr " " "_" >> $SAMPLE_SHEET_MISEQ;
	fi;
done;
echo "#[`date`] COMMAND: SAMPLESHEET cleaning (Sample Name bad characters)" >> $LOGFILE
$COMMAND >> $LOGFILE

#cat $SAMPLE_SHEET_MISEQ
#exit 0;


# BCL2FASTQ Sample Sheet
#COMMAND="$NGS_BIN_DIR/scripts/ElectricFancyFox.py $SAMPLE_SHEET_MISEQ $SAMPLE_SHEET_BCL2FASTQ bcl2fastqSampleSheet $RUN"
COMMAND="$SCRIPT_DIR/ElectricFancyFox.py $SAMPLE_SHEET_MISEQ $SAMPLE_SHEET_BCL2FASTQ bcl2fastqSampleSheet $RUN"

echo "#[`date`] COMMAND: $COMMAND" >> $LOGFILE
$COMMAND >> $LOGFILE

# Manifests_list
COMMAND="$SCRIPT_DIR/ElectricFancyFox.py $SAMPLE_SHEET_MISEQ $MANIFESTS_LIST manifests_list $RUN"
echo "#[`date`] COMMAND: $COMMAND" >> $LOGFILE
$COMMAND >> $LOGFILE

# Manifests
COMMAND="$SCRIPT_DIR/ElectricFancyFox.py $SAMPLE_SHEET_MISEQ $MANIFESTS manifests $RUN"
echo "#[`date`] COMMAND: $COMMAND" >> $LOGFILE
$COMMAND >> $LOGFILE

# Copy manifests and create intervals
if [ -e $INTERVALS ]; then
	rm $INTERVALS
fi
if [ -e $OUTPUT_DIR/run.manifest ]; then
	rm $OUTPUT_DIR/run.manifest
fi

if ((1)); then
for MANIFEST in `cat $MANIFESTS | tr " " "|" | tr "," " "`; do
	# Find manifest file
	MANIFEST_FILE=`echo $MANIFEST | tr "|" " "`
	# Copy manifest
	if [ -e $MISEQ_DIR/$RUN/$MANIFEST_FILE ]; then
		COMMAND="cp \"$MISEQ_DIR/$RUN/$MANIFEST_FILE\" $OUTPUT_DIR/  >> $LOGFILE"
	elif [ -e $MANIFEST_FOLDER/$MANIFEST_FILE ]; then
		COMMAND="cp \"$MANIFEST_FOLDER/$MANIFEST_FILE\" $OUTPUT_DIR/  >> $LOGFILE"
	fi;
	echo "#[`date`] COMMAND: $COMMAND" >> $LOGFILE
	eval $COMMAND
	#cp "$MISEQ_DIR/$RUN/$MANIFEST_FILE" $OUTPUT_DIR/ >> $LOGFILE
	
	#if [ -e "$OUTPUT_DIR/$MANIFEST_FILE" ]; then
	#	ls -l "$OUTPUT_DIR/$MANIFEST_FILE" >> $LOGFILE;
	#else
	#	echo "[ERROR] in Manifest '$MANIFEST' copy" 
	#fi;
	if [ ! -e "$OUTPUT_DIR/$MANIFEST_FILE" ]; then
		echo "[ERROR] in Manifest '$MANIFEST' copy" >> $LOGFILE;
	fi;
	# Add manifest to RUN intervals
	COMMAND="$SCRIPT_DIR/manifest.pl --manifest=\"$OUTPUT_DIR/$MANIFEST_FILE\" --intervals=$INTERVALS.tmp"
	echo "#[`date`] COMMAND: $COMMAND" >> $LOGFILE
	$SCRIPT_DIR/manifest.pl --manifest="$OUTPUT_DIR/$MANIFEST_FILE" --intervals=$INTERVALS.tmp >> $LOGFILE
	# Create run.interval
	COMMAND="cat $INTERVALS.tmp >> $INTERVALS"
	echo "#[`date`] COMMAND: $COMMAND" >> $LOGFILE
	cat $INTERVALS.tmp >> $INTERVALS
	# Remove temporary file run.interval.tmp
	COMMAND="rm $INTERVALS.tmp"
	echo "#[`date`] COMMAND: $COMMAND" >> $LOGFILE
	rm $INTERVALS.tmp
	# Create "run.manifest"
	COMMAND="cp \"$OUTPUT_DIR/$MANIFEST_FILE\" $OUTPUT_DIR/run.manifest"
	echo "#[`date`] COMMAND: $COMMAND" >> $LOGFILE
	eval $COMMAND
	
done;
fi;


# Reads Length file
COMMAND="$SCRIPT_DIR/ElectricFancyFox.py $SAMPLE_SHEET_MISEQ $READS_LENGTH readsLength $RUN"
echo "#[`date`] COMMAND: $COMMAND" >> $LOGFILE
$COMMAND >> $LOGFILE

# BCL2FASTQ Mask option file
COMMAND="$SCRIPT_DIR/ElectricFancyFox.py $SAMPLE_SHEET_MISEQ $MASK mask $RUN"
echo "#[`date`] COMMAND: $COMMAND" >> $LOGFILE
$COMMAND >> $LOGFILE

# BCL2FASTQ Adapters
COMMAND="$SCRIPT_DIR/ElectricFancyFox.py $SAMPLE_SHEET_MISEQ $ADAPTERS adapters $RUN"
echo "#[`date`] COMMAND: $COMMAND" >> $LOGFILE
$COMMAND >> $LOGFILE



# Parameters from SampleSheet

RUNS_SAMPLES=""
NB_SAMPLE=0
RUN_Investigator_Name=$(grep "Investigator Name" $SAMPLE_SHEET_MISEQ | tr -d "\r\n" | awk -F, '{print $2}')
DATA_SECTION_LINE=$(grep "^\[Data\]" $SAMPLE_SHEET_MISEQ -n | awk -F: '{print $1}')
SAMPLE_SECTION_FIRST_LINE=$(($DATA_SECTION_LINE+1))
SAMPLE_LINES=$(tail -n $(($SAMPLE_SECTION_FIRST_LINE-$(grep ^ -c $SAMPLE_SHEET_MISEQ)))  $SAMPLE_SHEET_MISEQ)
SAMPLE_PROJECT_COL=$(grep -i ^Sample_ID $SAMPLE_SHEET_MISEQ | tr -d '\r\n' | sed -e 's/,/\n/g' | grep -i -n Sample_Project | cut -d \: -f 1)
for L in $(tail -n $(($SAMPLE_SECTION_FIRST_LINE-$(grep ^ -c $SAMPLE_SHEET_MISEQ)))  $SAMPLE_SHEET_MISEQ);
do
	if [ "$L" == "" ] || [[ $L =~ ^\[.*\] ]]; then
		break;
	else
		((NB_SAMPLE++))
		SAMPLE_ID=$(echo $L | awk -F, '{print $1}')
		SAMPLE_PROJECT=$(echo $L | tr -d '\r\n' | cut -d \, -f $SAMPLE_PROJECT_COL )
		if [ "$SAMPLE_PROJECT" == "" ]; then SAMPLE_PROJECT=$RUN_Investigator_Name; fi;
		RUNS_SAMPLES="$RUNS_SAMPLES$RUN:$SAMPLE_ID:$RUN_Investigator_Name "
	fi;
done;

# --writing-threads
#WRITING_THREADS=1	# default
WRITING_THREADS=$(($NB_SAMPLE>$CORES_TO_USE?$CORES_TO_USE:$NB_SAMPLE));

if [ "$BARCODE_MISMATCHES" == "" ]; then
	BARCODE_MISMATCHES="1"
	echo "No barcode mismatches defined -> default used (1 mismatch)"
fi

# Demultiplexing configuration
ADAPTER_STRINGENCY=0.9 # default=0.9
#COMMAND="$BCL2FASTQ --force --runfolder-dir $RUNFOLDER_DIR  --output-dir $OUTPUT_DIR --sample-sheet $SAMPLE_SHEET_BCL2FASTQ --use-bases-mask `cat $MASK` --adapter-sequence $ADAPTERS --fastq-cluster-count 0 --flowcell-id $RUN"
COMMAND="$BCL2FASTQ --runfolder-dir $RUNFOLDER_DIR  --output-dir $OUTPUT_DIR --sample-sheet $SAMPLE_SHEET_MISEQ  --barcode-mismatches $BARCODE_MISMATCHES --fastq-compression-level 9 --no-lane-splitting -r $CORES_TO_USE -d $CORES_TO_USE  -w $WRITING_THREADS "
#--use-bases-mask `cat $MASK` 
# for bcl2fastq not casava : --adapter-stringency $ADAPTER_STRINGENCY 
#COMMAND="$CASAVA/configureBcl2fastq --force --input-dir $INPUT_DIR  --output-dir $OUTPUT_DIR --sample-sheet $SAMPLE_SHEET_BCL2FASTQ --use-bases-mask `cat $MASK` --adapter-sequence $ADAPTERS --fastq-cluster-count 0 --mismatches 1 --flowcell-id $RUN"
echo "#[`date`] COMMAND: " $COMMAND >> $LOGFILE

$COMMAND 1>>$LOGFILE 2>>$LOGFILE
 
touch $OUTPUT_DIR/Makefile

# Log
#echo "#[`date`] COMMAND END: " $(tail -n 1 $LOGFILE) #>> $LOGFILE


# End file
tail -n 1 $LOGFILE > $OUTPUT_DIR/STARKComplete.txt

#/home1/IRC/TOOLS/tools/bcl2fastq/2.17.1.14/bin/bcl2fastq --runfolder-dir /home1/IRC/DATA/RAW/MSR/151030_M01656_0065_000000000-D0NKK --output-dir /home1/IRC/DATA/DEV/DEM/151030_M01656_0065_000000000-D0NKK --sample-sheet /home1/IRC/DATA/DEV/DEM/151030_M01656_0065_000000000-D0NKK/SampleSheet.csv --use-bases-mask Y251,I8,I8,Y251 --help


# Demultiplexing
#if [ $DEMULTIPLEXING -eq 1 ]; then
#	#COMMAND="make -j 2 --directory=$OUTPUT_DIR -f $OUTPUT_DIR/Makefile"
#	COMMAND="make -j $THREADS --directory=$OUTPUT_DIR -f $OUTPUT_DIR/Makefile"
#	echo "#[`date`] COMMAND: " $COMMAND >> $LOGFILE
#	$COMMAND 1>>$LOGFILE 2>>$LOGFILE
#fi;

exit 0;


