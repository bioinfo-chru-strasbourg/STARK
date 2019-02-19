#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="Monitor"
SCRIPT_DESCRIPTION="Monitor Analysis (TMP folder)"
SCRIPT_RELEASE="0.9.6.2b"
SCRIPT_DATE="22/03/2016"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

HEADER="#######################################
# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]
# $SCRIPT_DESCRIPTION 
# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © 
#######################################";
#echo "$HEADER";

# Realse note
#RELEASE_NOTES="#\n"
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.3b-09/10/2015: Add RUN configuration depending on Group and Project defined in RUN SampleSheet\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.4b-13/10/2015: Add Complete and Running flag files\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.5b-16/02/2016: Add Exec time \n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.6b-30/08/2016: Add number of samples to analyse by log. Add colors. Add Tail and errors \n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.6.1b-19/09/2016: Add details options -1\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.6.2b-22/03/2017: Change log view\n";

# NOTES
if [ "${1^^}" == "RELEASE" ] || [ "${1^^}" == "RELEASE_NOTES" ] || [ "${1^^}" == "NOTES" ]; then
	echo "# RELEASE NOTES:";
	echo -e $RELEASE_NOTES
	exit 0;
fi;


# 1.1. Configuration
#####################

#HELP
if [ "${1^^}" == "HELP" ]; then
	echo "# USAGE:  $0 <ENV> <NB_RUN> <NB_LOG> <DETAILS> <TAIL> <FILTER>";
	echo "# ENV             ENVironnement file (define folders, tools...). Use ALL to check all Result folders";
	echo "# NB_RUN		Number of RUN to show (default ALL). Number by ENV if all Result folders checked";
	echo "# NB_LOG		Number of LOG to show by RUN (default 1)";
	echo "# DETAILS		1 for more details, -1 for less details, 0 for standard details (default 0)";
	echo "# TAIL		Number of last log lines to show. default 0";
	echo "# FILTER		Filter RUN on state (RTA_COMPLETE, QUEUED, RUNNING, COMPLETE)";
	echo "# HEADER		1 to show script header, 0 not. Default 1";
	exit 0;
fi;

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# FUNCTIONS
function join { local IFS="$1"; shift; echo "$*"; }

# INPUT
ENV=$1
NB_RUN=$2
NB_LOG=$3
DETAILS=$4
TAIL=$5
FILTER=$6
HEADER=$7


# ENV
if [ "$ENV" == "ALL" ]; then ENV_ALL=1; fi;
if [ -s $ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$ENV;
elif [ -s $APPS/$ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$APPS/$ENV;
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


# NB_RUN
if [ "$NB_RUN" == "" ]; then
	NB_RUN=1000000000000
fi;


# NB_LOG
if [ "$NB_LOG" == "" ]; then
	NB_LOG=1
fi;

# DETAILS
#DETAILS=1
if [ "$DETAILS" == "" ]; then
	DETAILS=0;
#else
#	DETAILS=1;
fi;

#echo $DETAILS; exit 0;
if [ "$HEADER" == "" ]; then HEADER=1; fi;

# HEADER
if (($HEADER)); then
	echo "$HEADER";
fi;
# TAIL
#TAIL=1

#CONFIG=$1
#if [ ! -e $CONFIG ] && [ ! -e $NGS_SCRIPTS/$CONFIG ] || [ "$CONFIG" == "" ]; then
#	CONFIG=$NGS_SCRIPTS/config.ini
#fi;





(($HEADER)) && echo "# ENV '"$(basename $ENV)"'"
if (($ENV_ALL)) && (($HEADER)); then echo "# All Result Folders checked... "; fi;

: '
if (($DETAILS)); then
echo "# "
	echo "# NGS Folder:            "$NGS_FOLDER
	#echo "# NGS BIN Folder:        "$NGS_BIN
	echo "# NGS SCRIPTS Folder:    "$NGS_SCRIPTS
	#echo "# MISEQ Folder:          "$MISEQ_FOLDER
	echo "# DEMULTIPLEXING Folder: "$DEMULTIPLEXING_FOLDER
	echo "# RESULTS Folder:        "$RESULTS_FOLDER
	#echo "# ANALYSIS Folder:       "$ANALYSIS_FOLDER
	echo "# TMP Folder:            "$TMP_FOLDER
	#echo "# CONFIG ini:            "$CONFIG
fi;
'

#echo "# "

#COLORS
#COLOR_END="\e[0m"
COLOR_NULL="\e[0m"
COLOR_GREY="\e[90m"
COLOR_0="\e[93m" #"\e[31m" RED
COLOR_PROCESSING="\e[93m"
COLOR_COMPLETE="\e[92m"	
COLOR_OK="\e[32m"	
COLOR_KO="\e[31m"	


TMP_FOLDER_LIST=$RESULTS_FOLDER
#: '
	if (($ENV_ALL)); then

		TMP_FOLDER_LIST="";
		#echo "TRUC "$(echo $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/env.*sh)
		for ENV_DEF in $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/env.*sh; do
			#echo $ENV_DEF
			#if [ $(basename $ENV_DEF) != env.sh ]; then
			APP=$(echo $(basename $ENV_DEF) | sed "s/^env.//gi" | sed "s/.sh$//gi" | sed "s/sh$//gi")
			if [ "$APP" == "" ]; then
				APP="default";
				ENV=$APP;
			fi
			# SOURCE
			#source $ENV_DEF;
			#RAW_FOLDER=$($SCRIPT_DIR/extract_variable_from_env.sh $ENV_DEF "RAW_FOLDER")
			#MISEQ_FOLDER=$($SCRIPT_DIR/extract_variable_from_env.sh $ENV_DEF "MISEQ_FOLDER");
			#MAIN_FOLDER=$($SCRIPT_DIR/extract_variable_from_env.sh $ENV_DEF "MAIN_FOLDER")
			#PIPELINES=$($SCRIPT_DIR/extract_variable_from_env.sh $ENV_DEF "PIPELINES")
			#TMP_FOLDER_TEST=$($SCRIPT_DIR/extract_variable_from_env.sh $ENV_DEF "TMP_FOLDER") #RESULTS_FOLDER
			TMP_FOLDER_TEST=$($SCRIPT_DIR/extract_variable_from_env.sh $ENV_DEF "RESULTS_FOLDER") #RESULTS_FOLDER
			# OUTPUT
			TMP_FOLDER_LIST=$TMP_FOLDER_LIST" $TMP_FOLDER_TEST"
			#RAW_LIST=$RAW_LIST" $RAW_FOLDER"
			#echo $RAW_LIST
			#fi;
			#echo "$ENV_DEF $TMP_FOLDER_TEST"
		done;
		TMP_FOLDER_LIST=$(echo $TMP_FOLDER_LIST | sed "s/ /\n/gi" | sort | uniq | sed "s/\n/ /gi" );
		

	fi;
#'
#echo "$TMP_FOLDER_LIST"
#exit 0;
for TMP_FOLDER in $TMP_FOLDER_LIST; do 

#echo "## RESULT FOLDER $TMP_FOLDER "
#echo "############################ "

if [ ! -d $TMP_FOLDER ]; then
	continue;
fi;

TMP_FOLDER_AVAIL=$(df -h $TMP_FOLDER | tail -n1 | awk '{print $4}')
TMP_FOLDER_AVAIL_PERCENT=$(df -h $TMP_FOLDER | tail -n1 | awk '{print $5}')

#echo "ls $TMP_FOLDER$RESULTS_SUBFOLDER -c -r | tail -n $NB_RUN"
#echo $NB_RUN; #exit 0;

#RESULTS_SUBFOLDER="/RES";

#echo "$TMP_FOLDER$RESULTS_SUBFOLDER/$RUN"; continue;

for RUN in `ls $TMP_FOLDER$RESULTS_SUBFOLDER -c -r | tail -n $NB_RUN`;
do

	#echo $RUN; continue;
	

	# STATE
	RTA_COMPLETE_FILE=$MISEQ_FOLDER/$RUN/RTAComplete.txt
	STARK_QUEUED_FILE=$TMP_FOLDER$RESULTS_SUBFOLDER/$RUN/$STARK_QUEUED
	STARK_RUNNING_FILE=$TMP_FOLDER$RESULTS_SUBFOLDER/$RUN/$STARK_RUNNING
	STARK_COMPLETE_FILE=$TMP_FOLDER$RESULTS_SUBFOLDER/$RUN/$STARK_COMPLETE

	#echo $MISEQ_FOLDER/$RUN/RTAComplete.txt; continue;

	# FILTER
	if [ "$FILTER" != "" ]; then
		if [ "$FILTER" == "RTA_COMPLETE" ] && [ ! -e $RTA_COMPLETE_FILE ]; then continue; fi;
		if [ "$FILTER" == "QUEUED" ] && [ ! -e $STARK_QUEUED_FILE ]; then continue; fi;
		if [ "$FILTER" == "RUNNING" ] && [ ! -e $STARK_RUNNING_FILE ]; then continue; fi;
		if [ "$FILTER" == "COMPLETE" ] && [ ! -e $STARK_COMPLETE_FILE ]; then continue; fi;
	fi;

	#echo "$TMP_FOLDER$RESULTS_SUBFOLDER/$RUN"; continue;
	# Find the number of Sample
	#NB_SAMPLES=$(find $TMP_FOLDER$RESULTS_SUBFOLDER/$RUN -mindepth 1 -maxdepth 1 -type d  | wc -l)
	
	#TEST
	NB_SAMPLES="?"
	
	#SAMPLE_SHEET_MISEQ_ORIGINAL=$TMP_FOLDER$RESULTS_SUBFOLDER/$RUN/SampleSheet.csv
	SAMPLE_SHEET_MISEQ=$DEMULTIPLEXING_FOLDER/$RUN/SampleSheet.csv
	RUNS_SAMPLES=""
	if [ -s $SAMPLE_SHEET_MISEQ ]; then
		#SAMPLE_SHEET_MISEQ=$TMP_FOLDER$RESULTS_SUBFOLDER/$RUN/SampleSheet.csv
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
				RUNS_SAMPLES="$RUNS_SAMPLES$RUN:$SAMPLE_ID:$SAMPLE_PROJECT "
			fi;
		done;
		NB_SAMPLES=$NB_SAMPLE
	fi;
	if [ $NB_SAMPLES == "" ]; then NB_SAMPLES="?"; fi;
	#echo $NB_SAMPLE
	#NB_SAMPLES="0"


	echo "#";
	#echo "## RUN ########## $RUN [$NB_SAMPLES samples] [$TMP_FOLDER $TMP_FOLDER_AVAIL/$TMP_FOLDER_AVAIL_PERCENT Avail.]";
	echo "## RUN ########## $RUN [$NB_SAMPLES samples] [$TMP_FOLDER $TMP_FOLDER_AVAIL Avail.]";
	#echo $LASTLOG; continue;

	if (($DETAILS>0)); then
		if [ -e $RTA_COMPLETE_FILE ]; then echo "## RTA Complete"; fi
		if [ -e $STARK_QUEUED_FILE ]; then echo "## QUEUED"; fi
		if [ -e $STARK_RUNNING_FILE ]; then echo "## RUNNING..."; fi
		if [ -e $STARK_COMPLETE_FILE ]; then echo "## COMPLETE"; fi
	fi;


	# Find the last analysis
	#LASTLOG=$(ls -r $TMP_FOLDER$RESULTS_SUBFOLDER/$RUN/run_analysis.V*.log 2>/dev/null| grep -E '.*V[0-9]+.*[0-9]+.log' - | head -$NB_LOG)
	LASTLOG=$(ls -r $TMP_FOLDER$RESULTS_SUBFOLDER/$RUN/*analysis.V*.log 2>/dev/null| grep -E '.*V[0-9]+.*[0-9]+.log' - | head -$NB_LOG)

	#echo $LASTLOG; continue;

	#if [ -e $LASTLOG ] && [ "$LASTLOG" != "" ]; then

		for LOG in $LASTLOG; do

			if [[ -r $LOG ]]; then

				#echo $LOG;
				MAKEFILE_ANALYSIS_RUN=$(echo $LOG | sed s/log$/param.mk/)
				#echo $MAKEFILE_ANALYSIS_RUN
				if [ -z "$RUNS_SAMPLES" ]; then
					RUNS_SAMPLES=$(grep "^RUNS_SAMPLES=" $MAKEFILE_ANALYSIS_RUN | cut -d= -f2)
				fi;
				NB_SAMPLE=$(grep "^RUNS_SAMPLES=" $MAKEFILE_ANALYSIS_RUN | wc -w)
				#echo $NB_SAMPLE
				#continue;
				
				

				LASTANA_PATTERN=$(echo $LOG | sed s/.log$//g)
				NB_ERROR=$(grep "make: \*\*\*.*Error.*" $LOG -c);
				STARTED=$(grep "Main Analysis Process for RUN '$RUN' START" $LOG -c);
				#START=$(grep "Main Analysis Process for RUN '.*' START" $LOG | sed 's/.*\[\([^]]*\)\].*/\1/g')
				START=$(grep "Main Analysis Process for .* '.*' START" $LOG | sed 's/.*\[\([^]]*\)\].*/\1/g')
				FINISHED=$(grep "Main Analysis Process for RUN '$RUN' END" $LOG -c);
				#END=$(grep "Main Analysis Process for RUN '.*' END" $LOG | sed 's/.*\[\([^]]*\)\].*/\1/g')
				END=$(grep "Main Analysis Process for .* '.*' END" $LOG | sed 's/.*\[\([^]]*\)\].*/\1/g')

				# EXEC TIME
				TODAY=`date '+%Y%m%d-%H%M%S'`
				if [ "$START" == "" ]; then
						START=$TODAY
				fi;
				if [ "$END" == "" ]; then
						END=$TODAY
				fi;
				let DAYS=(`date +%s -d $(echo $END | cut -d\- -f1)`-`date +%s -d $(echo $START | cut -d\- -f1)`)
				StartDate=$(date -u -d "$(join ":" $(echo $START | cut -d\- -f2 | fold -w2))" +"%s")
				let FinalDate=($(date -u -d "$(join ":" $(echo $END | cut -d\- -f2 | fold -w2))" +"%s")+$DAYS)
				EXEC_TIME=$(date -u -d "0 $FinalDate sec - $StartDate sec" +"%s")
				EXEC_DATETIME=$(date -u -d "0 $FinalDate sec - $StartDate sec" +"%H:%M:%S")

				if [ "$FINISHED" == "0" ] && [ -e $TMP_FOLDER$RESULTS_SUBFOLDER/$RUN/$STARK_RUNNING ]; then
					END="RUNNING [$EXEC_DATETIME]...";
				elif [ "$FINISHED" == "0" ] && [ -e $TMP_FOLDER$RESULTS_SUBFOLDER/$RUN/$STARK_RUNNING ]; then
					END="RUNNING [$EXEC_DATETIME]...";
				else
					END="COMPLETE [$EXEC_DATETIME]";
				fi;
				# TIME
				#grep "Main Analysis Process for RUN '$RUN' END" run_analysis.V20160216-120530.log | sed 's/.*\[\([^]]*\)-\([^]]*\)\].*/\1\2/g'

				# NB Sample for this analysis
				#NB_SAMPLES_TO_ANALYZE=$(echo $RUNS_SAMPLES | wc -w)
			
				NB_SAMPLES_TO_ANALYZE=$(grep "^RUNS_SAMPLES=" $MAKEFILE_ANALYSIS_RUN | wc -w)
				#NB_SAMPLES=$NB_SAMPLES_TO_ANALYZE
				#RUNS_SAMPLES_LIST=$(echo $RUNS_SAMPLES | awk )
				#ls $(echo $RUNS_SAMPLES | cut -d: -f1,2 | sed "s/:/\//gi")


				# Print
				#echo -ne "\e[31m";
				if [ $NB_ERROR == 0 ]; then 
					COLOR_NB_ERROR=$COLOR_OK
				else
					COLOR_NB_ERROR=$COLOR_KO
				fi;
				echo -e "## ANALYSIS ##### "$(basename $LASTANA_PATTERN)" [$NB_SAMPLES_TO_ANALYZE/$NB_SAMPLES Samples] [$COLOR_NB_ERROR$NB_ERROR errors$COLOR_NULL] $END"
			
			
			
				if (($DETAILS<0)); then continue; fi;

				#echo $LOG
				ALIGNERS=""
				CALLERS=""
				ANNOTATORS=""
				PIPELINES=""
				source $LASTANA_PATTERN.param.sh
				
				if [ -z "$ALIGNERS" ]; then
					ALIGNERS=$(echo $PIPELINES | tr " " "\n"  | cut -d. -f1 | sort | uniq | tr "\n" " ")
				fi;
				if [ -z "$CALLERS" ]; then
					CALLERS=$(echo $PIPELINES | tr " " "\n"  | cut -d. -f2 | sort | uniq | tr "\n" " ")
				fi;
				if [ -z "$ANNOTATORS" ]; then
					ANNOTATORS=$(echo $PIPELINES | tr " " "\n"  | cut -d. -f3 | sort | uniq | tr "\n" " ")
				fi;
				
				#echo $ALIGNERS
				NB_ALIGNERS=$(echo $ALIGNERS | wc -w)
				#echo "$NB_ALIGNERS Aligners"
				#echo $CALLERS
				NB_CALLERS=$(echo $CALLERS | wc -w)
				#echo "$NB_CALLERS Callers"
				#echo $ANNOTATORS
				NB_ANNOTATORS=$(echo $ANNOTATORS | wc -w)
				#echo "$NB_ANNOTATORS Annotators"

				# Find the pipelines launched in the last analysis
				# Calculation
				#NB_ALIGNMENTS_DONE=$(ls $TMP_FOLDER$RESULTS_SUBFOLDER/$RUN/*/*unaligned.bam | wc -l )
				#NB_ALIGNMENTS_DONE=0
				#NB_ALIGNMENTS_DONE=$(ls $(for i in $RUNS_SAMPLES; do echo $TMP_FOLDER$RESULTS_SUBFOLDER/$i/*unaligned.bam | cut -d: -f1,2 | sed "s/:/\//gi" | sed "s/$/\/*unaligned.bam/gi" ; done;) 2>/dev/null | wc -l )
				NB_ALIGNMENTS_DONE=$(echo $(for i in $(echo $RUNS_SAMPLES | tr " " "\n" | cut -d: -f1,2 --output-delimiter=/); do echo $TMP_FOLDER$RESULTS_SUBFOLDER/$i/*unaligned.bam ; done;) 2>/dev/null | wc -w )
				#NB_ALIGNMENTS_DONE=$(for i in $(echo $RUNS_SAMPLES | tr " " "\n" | cut -d: -f1,2 --output-delimiter=/); do echo $TMP_FOLDER$RESULTS_SUBFOLDER/$i/*unaligned.bam ; done; 2>/dev/null | wc -l )
				#echo "NB_ALIGNMENTS_DONE=$NB_ALIGNMENTS_DONE"
				#for i in $RUNS_SAMPLES; do echo $TMP_FOLDER$RESULTS_SUBFOLDER/$i/*unaligned.bam | cut -d: -f1,2 | sed "s/:/\//gi" | sed "s/$/\/*unaligned.bam/gi" ; done;
				#echo $(ls $(for i in $RUNS_SAMPLES; do echo $TMP_FOLDER$RESULTS_SUBFOLDER/$i/*unaligned.bam | cut -d: -f1,2 | sed "s/:/\//gi" ; done;))
				#ls $(for i in $TEST; do echo $i"*vcf" | cut -d: -f1,2 | sed "s/:/\//gi" | sed "s/$/*how*/gi"; done;)
				#ls $(for i in $RUNS_SAMPLES; do echo $i | cut -d: -f1,2 | sed "s/:/\//gi"; done;)
				#echo "scale=2 ; ( $NB_ALIGNMENTS_DONE / ( $NB_SAMPLES_TO_ANALYZE ) ) * 100"
				PERCENTAGE_ALIGNMENTS_DONE=$(echo "scale=2 ; ( $NB_ALIGNMENTS_DONE / ( $NB_SAMPLES_TO_ANALYZE ) ) * 100" | bc)
				#echo "# $PERCENTAGE_ALIGNMENTS_DONE% of uBAM"
				COLOR=$COLOR_PROCESSING
				if [ "$PERCENTAGE_ALIGNMENTS_DONE" == "0" ]; then COLOR=$COLOR_0; fi;
				if [ "$PERCENTAGE_ALIGNMENTS_DONE" == "100.00" ]; then COLOR=$COLOR_COMPLETE; fi;
				printf "# $COLOR%6s$COLOR_NULL" $PERCENTAGE_ALIGNMENTS_DONE; echo -n "% BAM unaligned";
				if (($DETAILS)); then echo " ($NB_ALIGNMENTS_DONE out of $NB_SAMPLES_TO_ANALYZE)"; else echo ""; fi;
				#echo "# "$(echo "scale=2 ; ( "$(ls $TMP_FOLDER$RESULTS_SUBFOLDER/$RUN/*/*unfiltered.vcf | wc -l )" / ( $NB_SAMPLES * ( $NB_ALIGNERS * $NB_CALLERS * $NB_ANNOTATORS ) ) ) * 100" | bc)"% of VCF";

				PIPELINES=""
				PIPELINES_CALLERS=""
				for ALIGNER in $ALIGNERS;
				do
					for CALLER in $CALLERS;
					do
						PIPELINES_CALLERS=$PIPELINES_CALLERS" $ALIGNER.$CALLER";
						for ANNOTATOR in $ANNOTATORS;
						do
							PIPELINES=$PIPELINES" $ALIGNER.$CALLER.$ANNOTATOR";
						done;
					done;
				done;

				# ALIGNERS
				NB_ALIGNERS=$(echo $ALIGNERS | wc -w )
				NB_ALIGNERS_DONE=0;
				ALIGNERS_DETAILS=""
				for ALIGNER in $ALIGNERS;
				do
					#echo $PIPELINE
					#NB_ALIGNER_DONE=$(ls $TMP_FOLDER$RESULTS_SUBFOLDER/$RUN/*/*.$ALIGNER.bam 2>/dev/null | wc -l )
					NB_ALIGNER_DONE=$(ls $(for i in $RUNS_SAMPLES; do echo $TMP_FOLDER$RESULTS_SUBFOLDER/$i/*unaligned.bam | cut -d: -f1,2 | sed "s/:/\//gi" | sed "s/$/\/*.$ALIGNER.bam/gi" ; done;) 2>/dev/null | wc -l )
					#NB_ALIGNER_DONE=$(echo $(for i in $(echo $RUNS_SAMPLES | tr " " "\n" | cut -d: -f1,2 --output-delimiter=/); do echo $TMP_FOLDER$RESULTS_SUBFOLDER/$i/*.$ALIGNER.bam ; done;) 2>/dev/null | wc -w )
					NB_ALIGNER_DONE=$(echo $(for i in $(echo $RUNS_SAMPLES | tr " " "\n" | cut -d: -f1,2 --output-delimiter=/); do ls $TMP_FOLDER$RESULTS_SUBFOLDER/$i/*.$ALIGNER.bam 2>/dev/null ; done;) | wc -w)
					#for i in $RUNS_SAMPLES; do echo $TMP_FOLDER$RESULTS_SUBFOLDER/$i/*unaligned.bam | cut -d: -f1,2 | sed "s/:/\//gi" | sed "s/$/\/*.$ALIGNER.bam/gi" ; done;
					NB_ALIGNERS_DONE=$(echo "$NB_ALIGNERS_DONE + $NB_ALIGNER_DONE" | bc )
					PERCENTAGE_ALIGNER_DONE=$(echo "scale=2 ; ( $NB_ALIGNER_DONE / ( $NB_SAMPLES_TO_ANALYZE ) ) * 100" | bc)
					COLOR=$COLOR_PROCESSING
					if [ $PERCENTAGE_ALIGNER_DONE == 0 ]; then COLOR=$COLOR_0; fi;
					if [ $PERCENTAGE_ALIGNER_DONE == 100.00 ]; then COLOR=$COLOR_COMPLETE; fi;
					PERCENTAGE_ALIGNER_PRINT=$(printf "$COLOR%6s$COLOR_NULL" $PERCENTAGE_ALIGNER_DONE)
					ALIGNERS_DETAILS=$ALIGNERS_DETAILS"\n#         $PERCENTAGE_ALIGNER_PRINT% BAM for $ALIGNER";
				done;
				NB_ALIGNERS_FILES=$(echo "$NB_SAMPLES_TO_ANALYZE * $NB_ALIGNERS" | bc)
				PERCENTAGE_ALIGNERS_DONE=$(echo "scale=2 ; ( $NB_ALIGNERS_DONE / $NB_ALIGNERS_FILES ) * 100" | bc)
				#echo "# PERCENTAGE_ALIGNERS_DONE=$PERCENTAGE_ALIGNERS_DONE scale=2 ; ( $NB_ALIGNERS_DONE / $NB_ALIGNERS_FILES ) * 100"
				COLOR=$COLOR_PROCESSING
				if [ $PERCENTAGE_ALIGNERS_DONE == 0 ]; then COLOR=$COLOR_0; fi;
				if [ $PERCENTAGE_ALIGNERS_DONE == 100.00 ]; then COLOR=$COLOR_COMPLETE; fi;
				printf "# $COLOR%6s$COLOR_NULL" $PERCENTAGE_ALIGNERS_DONE; echo -n "% BAM aligned";
				if (($DETAILS)); then echo -e " ($NB_ALIGNERS_DONE out of $NB_ALIGNERS_FILES)$ALIGNERS_DETAILS"; else echo ""; fi;


				if [ "$FINISHED" == "0" ]; then # If running
					# PIPELINES CALLERS
					NB_PIPELINES_CALLERS=$(echo $PIPELINES_CALLERS | wc -w )
					NB_PIPELINES_CALLERS_DONE=0;
					PIPELINES_CALLERS_DETAILS="";
					for PIPELINE in $PIPELINES_CALLERS;
					do
						#echo $PIPELINE
						#NB_PIPELINE_DONE=$(ls $TMP_FOLDER$RESULTS_SUBFOLDER/$RUN/*/*.$PIPELINE.vcf 2>/dev/null | wc -l )
						#NB_PIPELINE_DONE=$(ls $(for i in $RUNS_SAMPLES; do echo $TMP_FOLDER$RESULTS_SUBFOLDER/$i/*unaligned.bam | cut -d: -f1,2 | sed "s/:/\//gi" | sed "s/$/\/*.$PIPELINE.vcf/gi" ; done;) 2>/dev/null | wc -l )
						NB_PIPELINE_DONE=$(echo $(for i in $(echo $RUNS_SAMPLES | tr " " "\n" | cut -d: -f1,2 --output-delimiter=/); do ls $TMP_FOLDER$RESULTS_SUBFOLDER/$i/*.$PIPELINE.vcf 2>/dev/null; done;) 2>/dev/null | wc -w )
						NB_PIPELINES_CALLERS_DONE=$(echo "$NB_PIPELINES_CALLERS_DONE + $NB_PIPELINE_DONE" | bc )
						PERCENTAGE_PIPELINE_DONE=$(echo "scale=2 ; ( $NB_PIPELINE_DONE / ( $NB_SAMPLES_TO_ANALYZE ) ) * 100" | bc)
						COLOR=$COLOR_PROCESSING
						if [ $PERCENTAGE_PIPELINE_DONE == 0 ]; then COLOR=$COLOR_0; fi;
						if [ $PERCENTAGE_PIPELINE_DONE == 100.00 ]; then COLOR=$COLOR_COMPLETE; fi;
						PERCENTAGE_PIPELINE_PRINT=$(printf "$COLOR%6s$COLOR_NULL" $PERCENTAGE_PIPELINE_DONE)
						PIPELINES_CALLERS_DETAILS=$PIPELINES_CALLERS_DETAILS"\n#         $PERCENTAGE_PIPELINE_PRINT% VCF for $PIPELINE";
						#if (($DETAILS)); then printf "# %6s" $PERCENTAGE_PIPELINE_DONE; echo  "% VCF for $PIPELINE"; fi;
					done;
					NB_PIPELINES_CALLERS_FILES=$(echo "( $NB_SAMPLES_TO_ANALYZE * $NB_PIPELINES_CALLERS )" | bc)
					PERCENTAGE_PIPELINES_CALLERS_DONE=$(echo "scale=2 ; ( $NB_PIPELINES_CALLERS_DONE / $NB_PIPELINES_CALLERS_FILES ) * 100" | bc)
					#echo "# $PERCENTAGE_PIPELINES_CALLERS_DONE% of VCF"
					COLOR=$COLOR_PROCESSING
					if [ $PERCENTAGE_PIPELINES_CALLERS_DONE == 0 ]; then COLOR=$COLOR_0; fi;
					if [ $PERCENTAGE_PIPELINES_CALLERS_DONE == 100.00 ]; then COLOR=$COLOR_COMPLETE; fi;
					printf "# $COLOR%6s$COLOR_NULL" $PERCENTAGE_PIPELINES_CALLERS_DONE; echo -n "% VCF non annotated"
					if (($DETAILS)); then echo -e " ($NB_PIPELINES_CALLERS_DONE out of $NB_PIPELINES_CALLERS_FILES)$PIPELINES_CALLERS_DETAILS"; else echo ""; fi;

				fi;

				# PIPELINES
				NB_PIPELINES=$(echo $PIPELINES | wc -w )
				NB_PIPELINES_DONE=0;
				PIPELINES_DETAILS="";
				for PIPELINE in $PIPELINES;
				do
					#echo $PIPELINE
					#NB_PIPELINE_DONE=$(ls $TMP_FOLDER$RESULTS_SUBFOLDER/$RUN/*/*.$PIPELINE.vcf 2>/dev/null | wc -l )
					#NB_PIPELINE_DONE=$(ls $(for i in $RUNS_SAMPLES; do echo $TMP_FOLDER$RESULTS_SUBFOLDER/$i/*unaligned.bam | cut -d: -f1,2 | sed "s/:/\//gi" | sed "s/$/\/*.$PIPELINE.vcf/gi" ; done;) 2>/dev/null | wc -l )
					NB_PIPELINE_DONE=$(echo $(for i in $(echo $RUNS_SAMPLES | tr " " "\n" | cut -d: -f1,2 --output-delimiter=/); do ls $TMP_FOLDER$RESULTS_SUBFOLDER/$i/*.$PIPELINE.vcf 2>/dev/null ; done;) 2>/dev/null | wc -w )
					NB_PIPELINES_DONE=$(echo "$NB_PIPELINES_DONE + $NB_PIPELINE_DONE" | bc )
					PERCENTAGE_PIPELINE_DONE=$(echo "scale=2 ; ( $NB_PIPELINE_DONE / ( $NB_SAMPLES_TO_ANALYZE ) ) * 100" | bc)
					COLOR=$COLOR_PROCESSING
					if [ $PERCENTAGE_PIPELINE_DONE == 0 ]; then COLOR=$COLOR_0; fi;
					if [ $PERCENTAGE_PIPELINE_DONE == 100.00 ]; then COLOR=$COLOR_COMPLETE; fi;
					PERCENTAGE_PIPELINE_PRINT=$(printf "$COLOR%6s$COLOR_NULL" $PERCENTAGE_PIPELINE_DONE)
					PIPELINES_DETAILS=$PIPELINES_DETAILS"\n#         $PERCENTAGE_PIPELINE_PRINT% VCF for $PIPELINE";
					#if (($DETAILS)); then printf "# %6s" $PERCENTAGE_PIPELINE_DONE; echo  "% VCF for $PIPELINE"; fi;
				done;
				NB_PIPELINES_FILES=$(echo "( $NB_SAMPLES_TO_ANALYZE * $NB_PIPELINES )" | bc)
				PERCENTAGE_PIPELINES_DONE=$(echo "scale=2 ; ( $NB_PIPELINES_DONE / $NB_PIPELINES_FILES ) * 100" | bc)
				#echo "# $PERCENTAGE_PIPELINES_DONE% of VCF"
				COLOR=$COLOR_PROCESSING
				if [ $PERCENTAGE_PIPELINES_DONE == 0 ]; then COLOR=$COLOR_0; fi;
				if [ $PERCENTAGE_PIPELINES_DONE == 100.00 ]; then COLOR=$COLOR_COMPLETE; fi;
				printf "# $COLOR%6s$COLOR_NULL" $PERCENTAGE_PIPELINES_DONE; echo -n "% VCF annotated"
				if (($DETAILS)); then echo -e " ($NB_PIPELINES_DONE out of $NB_PIPELINES_FILES)$PIPELINES_DETAILS"; else echo ""; fi;

				if [ "$TAIL" != "" ] && [ "$TAIL" != "0" ]; then 
					echo "#"
					echo "# LOG ($TAIL last lines)"
					echo -n "$COLOR_GREY"
					tail -n$TAIL $LOG | tr "\t" " " | cut -c 1-$(tput cols) ;
					echo -n "$COLOR_NULL"
				
					#NB_ERRORS=$(grep "\*\*\*" $LOG -c)
					ERRORS=$(grep "\*\*\*" $LOG | tail -n$TAIL  | tr "\t" " ")
					if (($NB_ERROR > 0)); then 
						echo "# $COLOR_KO$NB_ERROR ERRORS"
						#echo -ne $ERRORS
						grep "make: \*\*\*.*Error.*" $LOG # | tail -n$TAIL
						echo "$COLOR_NULL"
					fi;

				fi;
			
			else
			
				echo "##[WARNING] LOG file '$( basename $LOG)' not readable. No more information can be found...";
			
			fi;
			
		done;

	#fi;

done;

	
done;


if (($DETAILS)) && ((0)); then
	echo "#";
	echo "## CPU/MEM Use ##";
	# TOP
	#TOP=$(top -b -n 1 | head -n 4)
	#TOP=$(echo -e "$TOP")

	# TOP CPU
	#CPU=$(NUMCPUS=`grep ^proc /proc/cpuinfo | wc -l`; FIRST=`cat /proc/stat | awk '/^cpu / {print $5}'`; sleep 1; SECOND=`cat /proc/stat | awk '/^cpu / {print $5}'`; USED=`echo 2 k 100 $SECOND $FIRST - $NUMCPUS / - p | dc`; echo ${USED})
	CPU=$(HZ=`awk '/\#define HZ/ {print $3}' /usr/include/asm-generic/param.h`; NUMCPUS=`grep ^proc /proc/cpuinfo | wc -l`; FIRST=`cat /proc/stat | awk '/^cpu / {print $5}'`; sleep 1; SECOND=`cat /proc/stat | awk '/^cpu / {print $5}'`; USED=`echo 2 k 100 $SECOND $FIRST - $NUMCPUS $HZ \* / 100 \* - p | dc`; echo ${USED})
	#TOP_CPU=$(top -b -n 1 | grep ^Cpu | awk -F" " '{print $2}' | tr -d "us,")
	echo "# CPU "$CPU"%"

	# TOP MEM

	#TOP_MEM=$(top -b -n 1 | grep ^Mem | awk -F" " '{print "scale=1; "$4" / "$2" * 100 "}' | tr -d "k" | bc)
	MEM=$(free | grep Mem | awk '{print "scale=2; "$4/$2 * 100" /1 "}' | bc)
	echo "# MEM "$MEM"%"

fi;


#echo "";
