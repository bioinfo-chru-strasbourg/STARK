#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="DBintegration"
SCRIPT_DESCRIPTION="Integration into Database (RUNS, SAMPLES, VCF).
# Depends on the structure RUN/SAMPLE/*.vcf"
SCRIPT_RELEASE="0.9.1beta"
SCRIPT_DATE="11/03/2015"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"

echo "#######################################";
echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
echo "# $SCRIPT_DESCRIPTION ";
echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © ";
echo "#######################################";


# 1.1. Configuration
#####################

#HELP
if [ "${1^^}" == "HELP" ]; then
	echo "# USAGE:  $0 <CONFIG> <RUNS> <SAMPLES>";
	echo "# CONFIG  INI file for DB connexion, folders...";
	echo "# RUNS    List of RUNs to integrate";
	echo "# SAMPLES List of SAMPLEs to integrate";
	#echo "Usage : SELF <CONFIG> <RUNS> <SAMPLES>";
	exit 0;
fi;

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_DIR/env.sh

CONFIG=$1
if [ ! -e $CONFIG ] && [ ! -e $NGS_SCRIPTS/$CONFIG ] || [ "$CONFIG" == "" ]; then
	CONFIG=$NGS_SCRIPTS/config.ini
fi;

echo "# "
echo "# NGS Folder:            "$NGS_FOLDER
echo "# NGS BIN Folder:        "$NGS_BIN
echo "# NGS SCRIPTS Folder:    "$NGS_SCRIPTS
echo "# MISEQ Folder:          "$MISEQ_FOLDER
echo "# DEMULTIPLEXING Folder: "$DEMULTIPLEXING_FOLDER
echo "# RESULTS Folder:        "$RESULTS_FOLDER
echo "# ANALYSIS Folder:       "$ANALYSIS_FOLDER
echo "# CONFIG ini:            "$CONFIG
echo "# "


#DIR=/media/IRCV2/V2/RES/ALL
DIR=$RESULTS_FOLDER
#RUNS=`find $DIR -maxdepth 1 -type d -name [^\.]\* | sed "s:^$DIR::"`
#echo "$RUNS"; 
#RUN="140812_M01656_0017_000000000-AA7WY"
#SAMPLES="14C0002 14C0003"
MANIFEST_REF="BRCA" # To automatize in SampleSheet
PROJECT_REF_DEFAULT="MISC" # To automatize in SampleSheet
#PIPELINES=".bwamem.gatkHC.vap .bwamem.gatkUG.vap .bwamem.gatkHC.trakxs .bwamem.gatkUG.trakxs"
RUNS_LIST=$2 #140812_M01656_0017_000000000-AA7WY
#echo $RUNS_LIST; exit 0;
SAMPLES_LIST=$3 #"14C0002 14C0003"
#SAMPLES_LIST="14C0003"
PIPELINES_LIST= #".bwamem.gatkHC.vap" # TODO 14C0003.bwamem.gatkUG.vap.vcf
VERBOSE=--verbose
DEBUG=--debug
SNAPSHOT=0	# 0 or 1
METRICS_EXTS="flagstat idxstats sample_interval_summary HsMetrics depthbed genomeCoverageBedbed" # fastqc_data.txt  depth TO BIG

## MANIFEST
DIR_MANIFESTS=/media/IRCV2/RAW/Manifests
if [ "$MANIFEST_FOLDER" != "" ]; then
	DIR_MANIFESTS=$MANIFEST_FOLDER
elif [ "$RAW_FOLDER" != "" ]; then
	DIR_MANIFESTS=$RAW_FOLDER/Manifests
fi;

if [ -e $DIR_MANIFESTS ]; then

	MANIFESTS=`find $DIR_MANIFESTS -maxdepth 1 -type f -name \* | sed "s:^$DIR_MANIFESTS::" | sed "s:/::g" | tr " " "|" `

	if [ 1 == 1 ]; then
	for MANIFEST in $MANIFESTS
	do
		MANIFEST_FILE=`echo $MANIFEST | tr "|" " "`
		echo "# MANIFEST: $MANIFEST_FILE"
		if [ -e "$DIR_MANIFESTS/$MANIFEST_FILE" ]; then
			time $NGS_SCRIPTS/DBintegration.pl --config=$CONFIG --type=manifest --data=manifest.ref="$MANIFEST_FILE",manifest="$DIR_MANIFESTS/$MANIFEST_FILE" $VERBOSE $DEBUG
		fi;
	done;
	fi;

fi;

if [ "$RUNS_LIST" == "" ];
then
	RUNS=`find $DIR -maxdepth 1 -type d -name [^\.]\* | sed "s:^$DIR::" | sed "s:/::g"`
else
	RUNS=$RUNS_LIST;
fi;
#echo $RUNS


for RUN in $RUNS
do
	echo "# RUN: $RUN"
	if [ "$SAMPLES_LIST" == "" ];
	then
		SAMPLES=`find $DIR/$RUN -maxdepth 1 -type d -name [^\.]\* | sed "s:$DIR::" | sed "s:$RUN::" | sed "s:/::g"`
	else
		SAMPLES=$SAMPLES_LIST
	fi;
	for SAMPLE in $SAMPLES
	do
		#echo "$RUN/$SAMPLE"
		echo "# SAMPLE: $RUN/$SAMPLE"

		GROUP_REF=""
		INSTITUTE_REF=""
		LAB_REF=""
		PROJECT_REF=""
		USER_REF=""

		# SampleSheet
		#SAMPLESHEET=$DIR/$RUN/$SAMPLE/SampleSheet.csv
		SAMPLESHEET=$DIR/$RUN/$SAMPLE/$SAMPLE.SampleSheet.csv
		# Manifest
		MANIFEST=$DIR/$RUN/$SAMPLE/$SAMPLE.manifest
		MANIFEST_REF=
		#echo "# Manifest found: $MANIFEST (ref: $MANIFEST_REF)"
		# Project
		PROJECT_REF=$PROJECT_REF_DEFAULT	# TODO
		RUN_INVESTIGATOR=`grep -i 'Investigator Name' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2`
		#RUN_INFOS=     `grep -i 'Investigator Name' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2`
		if [[ $RUN_INVESTIGATOR =~ .*-.* ]]
		then
			RUN_GROUP=`grep -i 'Investigator Name' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2 | cut -d \- -f 1` # -s  for without delimiter
			RUN_PROJECT=`grep -i 'Investigator Name' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2 | cut -d \- -f 2`
			RUN_USER=`grep -i 'Investigator Name' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2 | cut -d \- -f 3`
		else
			RUN_GROUP=$RUN_INVESTIGATOR
		fi
		
		#RUN_INSTITUTE=${RUN_GROUP:0:3}
		#RUN_LAB=${RUN_GROUP:3}
		#RUN_DATE=`grep -i '^Date' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2`
		#RUN_APPLICATION=`grep -i '^Application' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2`
		#RUN_WORKFLOW=`grep -i '^Workflow' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2`
		#RUN_ASSAY=`grep -i '^Assay' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2`
		#RUN_DESCRIPTION=`grep -i '^Description' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2`
		#RUN_CHEMISTRY=`grep -i '^Chemistry' $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f 2`
		
		
		
		INDEX_SAMPLEPROJECT=`grep -i ^Sample_ID $SAMPLESHEET | tr -d '\r\n' | sed -e 's/,/\n/g' | grep -i -n Sample_Project | cut -d \: -f 1`
		SAMPLE_INVESTIGATOR=`grep ^$SAMPLE, $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f $INDEX_SAMPLEPROJECT`

		if [[ $SAMPLE_INVESTIGATOR =~ .*-.* ]]
		then
			SAMPLE_GROUP=`grep ^$SAMPLE, $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f $INDEX_SAMPLEPROJECT | cut -d \- -f 1`
			SAMPLE_PROJECT=`grep ^$SAMPLE, $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f $INDEX_SAMPLEPROJECT | cut -d \- -f 2`
			SAMPLE_USER=`grep ^$SAMPLE, $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f $INDEX_SAMPLEPROJECT | cut -d \- -f 3`
		else
			SAMPLE_GROUP=$SAMPLE_INVESTIGATOR
		fi
		#SAMPLE_GROUP=`grep ^$SAMPLE, $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f $INDEX_SAMPLEPROJECT | cut -d \- -s -f 1`
		#SAMPLE_PROJECT=`grep ^$SAMPLE, $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f $INDEX_SAMPLEPROJECT | cut -d \- -s -f 2`
		#SAMPLE_USER=`grep ^$SAMPLE, $SAMPLESHEET | tr -d '\r\n' | cut -d \, -f $INDEX_SAMPLEPROJECT | cut -d \- -s -f 3`
		
		# GROUP
		if [ ! "$SAMPLE_GROUP" == "" ] # && [ -d $RESULTS_BYGROUP_FOLDER/$RUN_GROUP ]
		then
			GROUP_REF=$SAMPLE_GROUP
		elif [ ! "$RUN_GROUP" == "" ] 
		then
			GROUP_REF=$RUN_GROUP
		else
			GROUP_REF="MISC"
			RUN_INSTITUTE="MISC"
			RUN_LAB="MISC"
		fi

		# PROJECT
		if [ ! "$SAMPLE_PROJECT" == "" ] # && [ -d $RESULTS_BYGROUP_FOLDER/$RUN_GROUP ]
		then
			PROJECT_REF=$SAMPLE_PROJECT
		elif [ ! "$RUN_PROJECT" == "" ] # && [ -d $RESULTS_BYGROUP_FOLDER/$RUN_GROUP ]
		then
			PROJECT_REF=$RUN_PROJECT
		elif [ ! "$GROUP_REF" == "" ] # && [ -d $RESULTS_BYGROUP_FOLDER/$RUN_GROUP ]
		then
			PROJECT_REF=$GROUP_REF
		else
			PROJECT_REF="unknown"
		fi

		# USER
		if [ ! "$SAMPLE_USER" == "" ] # && [ -d $RESULTS_BYGROUP_FOLDER/$RUN_GROUP ]
		then
			USER_REF=$SAMPLE_USER
		elif [ ! "$RUN_USER" == "" ] # && [ -d $RESULTS_BYGROUP_FOLDER/$RUN_GROUP ]
		then
			USER_REF=$RUN_USER
		else
			USER_REF="unknown"
		fi


		# Institute and Lab
		if [ "$GROUP_REF" != "MISC" ] && [ "$GROUP_REF" != "" ]; then
			INSTITUTE_REF=${GROUP_REF:0:3}
			LAB_REF=${GROUP_REF:3}
		else
			GROUP_REF="MISC"
			INSTITUTE_REF="MISC"
			LAB_REF="MISC"
		fi;
		#echo "# RUN GROUP: $RUN_GROUP"
		#echo "# RUN PROJECT: $RUN_PROJECT"
		#echo "# RUN USER: $RUN_USER"
		#echo "# SAMPLE GROUP: $SAMPLE_GROUP"
		#echo "# SAMPLE PROJECT: $SAMPLE_PROJECT"
		#echo "# SAMPLE USER: $SAMPLE_USER"
		echo "# GROUP: $GROUP_REF"
		echo "# INSTITUTE: $INSTITUTE_REF"
		echo "# LAB: $LAB_REF"
		echo "# PROJECT: $PROJECT_REF"
		echo "# USER: $USER_REF"
		#GROUP_REF=$RUN_PROJECT
		#continue

		if [ "$VCFS" == "" ];
		then
			#PIPELINES=`find $DIR/$RUN/$SAMPLE -maxdepth 1 -type d -name [^\.]\* | sed "s:$DIR/$RUN/$SAMPLE::" | sed "s:/::g"`
			#PIPELINES=`find $DIR/$RUN/$SAMPLE -maxdepth 1 -type d -name [^\.]\*[vcf] | sed "s:$DIR/$RUN/$SAMPLE::" | sed "s:/::g"`
			#-maxdepth 1 -type f -name 14C0002\*\.vcf
			VCFS=`find $DIR/$RUN/$SAMPLE -maxdepth 1 -type f -name $SAMPLE\*$PIPELINES_LIST\*\.vcf | sed "s:$DIR/$RUN/$SAMPLE::" | sed "s:/::g"`
		fi;


		for VCF in $VCFS
		do
			echo "# VCF: $RUN/$SAMPLE/$VCF";
			# VCF
			VCF_FILE=$DIR/$RUN/$SAMPLE/$VCF
			
			# PIPELINE
			ALIGNER=`echo $VCF | awk -F"." '{print $2}'`
			CALLER=`echo $VCF | awk -F"." '{print $3}'`
			ANNOTATOR=`echo $VCF | awk -F"." '{print $4}'`
			PIPELINE=".$ALIGNER.$CALLER.$ANNOTATOR"
			#echo "PIPELINE: $PIPELINE   ALIGNER: $ALIGNER   CALLER: $CALLER   ANNOTATOR: $ANNOTATOR"

			# BAM
			BAM=""
			if [ -e $DIR/$RUN/$SAMPLE/$SAMPLE.$ALIGNER.bam ]; then
				echo "# Associated BAM '$DIR/$RUN/$SAMPLE/$SAMPLE.$ALIGNER.bam'"
				BAM="$SAMPLE.$ALIGNER.bam"
			fi; 

			if [ -e $SAMPLESHEET ] && [ -e $MANIFEST ] && [ `grep ^ $VCF_FILE -c` -gt 0 ]; then

				# VCF
				DATA="sample.ref=$SAMPLE,run.ref:$RUN,samplesheet:$SAMPLESHEET,manifest.ref=$MANIFEST_REF,manifest=$MANIFEST,project.ref=$PROJECT_REF,group.ref=$GROUP_REF,group.institute=$INSTITUTE_REF,group.lab=$LAB_REF,user.ref=$USER_REF,vcf=$VCF_FILE,vcf.bam=$BAM,vcf.pipeline=$PIPELINE,vcf.aligner=$ALIGNER,vcf.caller=$CALLER,vcf.annotator=$ANNOTATOR"
				time $NGS_SCRIPTS/DBintegration.pl --config=$CONFIG --type=vcf --data=$DATA $VERBOSE $DEBUG
				truc=1
				
				# FASTQC
				FASTQC=`find $DIR/$RUN/$SAMPLE -maxdepth 3 -type f -name fastqc_data.txt`
				if [ -e $FASTQC ] && [ ! "$FASTQC" == "" ]; then
					$NGS_SCRIPTS/DBintegration.pl --config=$CONFIG --type=QC --data=$DATA,QC=$FASTQC,QC.metrics=FASTQC,QC.type=UnalignedBAM $VERBOSE $DEBUG
					truc=1
				fi;
				#echo "FASTQC: "$FASTQC

				# METRICS (usually 1 by Aligner)
				QCS=`find $DIR/$RUN/$SAMPLE -mindepth 1 -maxdepth 1 -type d -name [^\.]\*\.$ALIGNER\.\*.metrics`
				#echo "QCS: "$QCS

				for QC in $QCS
				do
					echo "#   QC: "$QC
					for METRICS_EXT in $METRICS_EXTS
					do
						#METRICS=`find $QC -maxdepth 3 -type f -name [^\.]\*$METRICS_EXT`
						METRICS=`find $QC -maxdepth 3 -type f -name \*$METRICS_EXT`
						#echo "      METRICS: "$METRICS
						if [ -e $METRICS ] && [ "$METRICS" != "" ] && [ `grep ^ $METRICS -c` -gt 0 ]; then
							echo "#      METRICS: "$METRICS
							time $NGS_SCRIPTS/DBintegration.pl --config=$CONFIG --type=QC --data=$DATA,QC=$METRICS,QC.metrics=$METRICS_EXT,QC.type=AlignedBAM $VERBOSE $DEBUG
							truc=1
						fi;
					done;
				done;

				# SNAPSHOT
				SNAPSHOTS_DIR=`find $DIR/$RUN/$SAMPLE -mindepth 1 -maxdepth 1 -type d -name [^\.]\*\.$ALIGNER\.\*.snapshots`
				echo "SNAPSHOTS DIR: "$SNAPSHOTS_DIR

				if [ $SNAPSHOT == 1 ]; then

				for SNAPSHOTS in $SNAPSHOTS_DIR
				do
					echo "#   SNAPSHOTS: "$SNAPSHOTS
					#echo "find $SNAPSHOTS -maxdepth 3 -type f -name *png"
					SNAPSHOTS_FILES=`find $SNAPSHOTS -maxdepth 3 -type f -name \*snapshot\*png`
					#echo "SNAPSHOTS FILES: "$SNAPSHOTS_FILES
					for SNAPSHOTS_FILE in $SNAPSHOTS_FILES
					do
						echo "#      SNAPSHOTS_FILE: "$SNAPSHOTS_FILE
						REF=`basename $SNAPSHOTS_FILE`
						FILENAME_NO_EXT="${REF%.*}"
						SNAPSHOT_INFO="${FILENAME_NO_EXT##*.}"
						#CHROM_POS=`echo $REF | awk -F"." '{print $2}'`
						CHROM=`echo $SNAPSHOT_INFO | awk -F"_" '{print $1}'`
						POS=`echo $SNAPSHOT_INFO | awk -F"_" '{print $2}'`
						START=`echo $SNAPSHOT_INFO | awk -F"_" '{print $3}'`
						STOP=`echo $SNAPSHOT_INFO | awk -F"_" '{print $4}'`
						echo "# REF=$REF CHROM=$CHROM POS=$POS START=$START STOP=$STOP"
						time $NGS_SCRIPTS/DBintegration.pl --config=$CONFIG --type=snapshot --data=$DATA,snapshot=$SNAPSHOTS_FILE,snapshot.ref=$REF,snapshot.chrom=$CHROM,snapshot.pos=$POS,snapshot.start=$START,snapshot.stop=$STOP,snapshot.bam=$BAM $VERBOSE $DEBUG
						truc=1
					done;
				done;

				fi;

			else
				echo "# [ERROR] Missing information: either no SAMPLESHEET '$SAMPLESHEET', no MANIFEST '$MANIFEST', or no variant in VCF_FILE '$VCF_FILE'";

			fi;

		done;
		VCFS="";
	done;
	SAMPLES="";

done;
exit 0;


