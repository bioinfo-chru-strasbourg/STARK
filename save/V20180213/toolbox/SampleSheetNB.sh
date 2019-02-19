#!/bin/bash
#################################
##
## NGS environment
##
#################################

# Number of sample in a SampleSheet
# Assume DATA section is the last one

CMD="ls > log"
eval $CMD

exit 0;

SAMPLESHEET=$1

MAKEFILE_ANALYSIS_RUN=mk.test
RUN=run
ALIGNERS="bwamem"
CALLERS="gatkUG gatkHC"
ANNOTATORS="howard"
PIPELINES="test"

echo "# SAMPLESHEET	$SAMPLESHEET"

if [ -e $SAMPLESHEET ]; then
	echo "# SAMPLESHEET OK"
	RUNS_SAMPLES=""
	NB_SAMPLE=0
	RUN_Investigator_Name=$(grep "Investigator Name" $SAMPLESHEET | tr -d "\r\n" | awk -F, '{print $2}')
	DATA_SECTION_LINE=$(grep "^\[Data\]" $SAMPLESHEET -n | awk -F: '{print $1}')
	SAMPLE_SECTION_FIRST_LINE=$(($DATA_SECTION_LINE+1))
	SAMPLE_LINES=$(tail -n $(($SAMPLE_SECTION_FIRST_LINE-$(grep ^ -c $SAMPLESHEET)))  $SAMPLESHEET)
 	SAMPLE_PROJECT_COL=$(grep -i ^Sample_ID $SAMPLESHEET | tr -d '\r\n' | sed -e 's/,/\n/g' | grep -i -n Sample_Project | cut -d \: -f 1)
	for L in $(tail -n $(($SAMPLE_SECTION_FIRST_LINE-$(grep ^ -c $SAMPLESHEET)))  $SAMPLESHEET);
	do
		if [ "$L" == "" ] || [[ $L =~ ^\[.*\] ]]; then
			break;
		else
			((NB_SAMPLE++))
			SAMPLE_ID=$(echo $L | awk -F, '{print $1}')
			SAMPLE_PROJECT=$(echo $L | tr -d '\r\n' | cut -d \, -f $SAMPLE_PROJECT_COL )
			if [ "$SAMPLE_PROJECT" == "" ]; then SAMPLE_PROJECT=$RUN_Investigator_Name; fi;
			RUNS_SAMPLES="$RUNS_SAMPLES$RUN:$SAMPLE_ID:$RUN_Investigator_Name "
			echo "SAMPLE_ID $SAMPLE_ID"
		fi;
	done;
	echo "NB_SAMPLE=$NB_SAMPLE"
else
	echo "[ERROR] No SampleSheet '$SAMPLESHEET'";
	exit 1;
fi;


SAMPLE_SHEET_MISEQ=$SAMPLESHEET

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

	# Additional parameters to Makefile Header
	echo "" > $MAKEFILE_ANALYSIS_RUN
	echo "#############" >> $MAKEFILE_ANALYSIS_RUN
	echo "# Parameters #" >> $MAKEFILE_ANALYSIS_RUN
	echo "#############" >> $MAKEFILE_ANALYSIS_RUN
	echo "" >> $MAKEFILE_ANALYSIS_RUN
	echo "RUNS=$RUN" >> $MAKEFILE_ANALYSIS_RUN
	echo "RUNS_SAMPLES=$RUNS_SAMPLES" >> $MAKEFILE_ANALYSIS_RUN
	echo "" >> $MAKEFILE_ANALYSIS_RUN
	echo "" >> $MAKEFILE_ANALYSIS_RUN
	echo "#########################" >> $MAKEFILE_ANALYSIS_RUN
	echo "# Additional parameters #" >> $MAKEFILE_ANALYSIS_RUN
	echo "#########################" >> $MAKEFILE_ANALYSIS_RUN
	echo "" >> $MAKEFILE_ANALYSIS_RUN
	echo "RUNS=$RUN" >> $MAKEFILE_ANALYSIS_RUN
	echo "ALIGNERS=$ALIGNERS" >> $MAKEFILE_ANALYSIS_RUN
	echo "CALLERS=$CALLERS" >> $MAKEFILE_ANALYSIS_RUN 
	echo "ANNOTATORS=$ANNOTATORS" >> $MAKEFILE_ANALYSIS_RUN 
	echo "INTERSEC=$INTERSEC" >> $MAKEFILE_ANALYSIS_RUN 
	echo "" >> $MAKEFILE_ANALYSIS_RUN

cat $MAKEFILE_ANALYSIS_RUN


CORES_TO_USE=31
grep "^RUNS_SAMPLES=" $MAKEFILE_ANALYSIS_RUN
NB_SAMPLE=0

# NB SAMPLE
NB_SAMPLE=$(grep "^RUNS_SAMPLES=" $MAKEFILE_ANALYSIS_RUN | wc -w)
echo "NB_SAMPLE=$NB_SAMPLE"

# NB ALIGNERS
NB_ALIGNERS=$(echo $ALIGNERS | wc -w)
echo "NB_ALIGNERS=$NB_ALIGNERS"

# NB CALLERS
NB_CALLERS=$(echo $CALLERS | wc -w)
echo "NB_CALLERS=$NB_CALLERS"

# NB ANNOTATORS
NB_ANNOTATORS=$(echo $ANNOTATORS | wc -w)
echo "NB_ANNOTATORS=$NB_ANNOTATORS"


# NB SAMPLE
# CREATE PIPELINES from ALIGNERS/CALLERS/ANNOTATORS
for ALIGNER in $ALIGNERS; do
	for CALLER in $CALLERS; do
		for ANNOTATOR in $ANNOTATORS; do
			PIPELINES="$PIPELINES $ALIGNER.$CALLER.$ANNOTATOR"
		done;
	done;
done;
PIPELINES=$(echo $PIPELINES | tr " " "\n" | sort | uniq | tr "\n" " ")
NB_PIPELINES=$(echo $PIPELINES | wc -w)
echo "NB_PIPELINES=$NB_PIPELINES"

echo "$THREADS_BY_SAMPLE/$NB_ALIGNERS	"$(($THREADS_BY_SAMPLE/$NB_ALIGNERS))


THREADS_BY_SAMPLE=$(($CORES_TO_USE/$NB_SAMPLE))			# Number of threads by sample
THREADS_BY_PIPELINE=$(($THREADS_BY_SAMPLE/$NB_PIPELINES))	# Number of threads for a pipeline's command
THREADS_BY_ALIGNER=$(($THREADS_BY_SAMPLE/$NB_ALIGNERS))	# Number of threads for a aligner's command
THREADS_BY_CALLER=$(($THREADS_BY_SAMPLE/$NB_CALLERS))	# Number of threads for a caller's command
THREADS_BWA=$THREADS_BY_ALIGNER								# NB of threads for BWA command
THREADS_RTC=$THREADS_BY_ALIGNER								# NB of threads for GATK realigment RealignerTargetCreator command


echo "THREADS_BY_SAMPLE=$THREADS_BY_SAMPLE"
echo "THREADS_BY_PIPELINE=$THREADS_BY_PIPELINE"
echo "THREADS_BY_ALIGNER=$THREADS_BY_ALIGNER"
echo "THREADS_BY_CALLER=$THREADS_BY_CALLER"
echo "THREADS_BWA=$THREADS_BWA"
echo "THREADS_RTC=$THREADS_RTC"



















