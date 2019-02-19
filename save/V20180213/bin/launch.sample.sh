#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="LaunchSample"
SCRIPT_DESCRIPTION="Launch a Sample Analysis"
SCRIPT_RELEASE="0.9.6b"
SCRIPT_DATE="27/10/2016"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

# Realse note
#RELEASE_NOTES="#\n"
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.3b-09/10/2015: Add RUN configuration depending on Group and Project defined in RUN SampleSheet\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.4b-13/10/2015: Add Complete and Running flag files\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.5.1b-11/04/2016: Optimize Threads and modifiy options order\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.5.2b-13/04/2016: Adding reports process\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.6b-27/10/2016: Multi input FASTQ/BAM/CRAM/SAM and associated parameters\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Configuration
ENV_CONFIG=$(find $SCRIPT_DIR/.. -name config.app)
source $ENV_CONFIG


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
	echo "# USAGE: $(basename $0) --fastq=<FILE1,FILE2,...> [options...]";
	echo "# -e|--application|--app|--env=<APP|FILE>     APP name or APP file configuration of the APPLICATION.";
	echo "#                                             Must be in the STARK folder if relative path";
	echo "#                                             Default: env.sh if not defined";
	echo "# -f|--fastq|--fastq_R1=<FILE1,FILE2,...>     List of FASTQ|BAM|SAM|CRAM (mandatory)";
	echo "#                                             Formats: *fastq.gz|*fq.gz|*bam|*ubam|*cram|*ucram|*sam|*usam";
	echo "# -q|--fastq_R2=<FILE1,FILE2,...>             List of corresponding FASTQ Read2 file (beware of correspondance).";
	echo "#                                             Format:  (*fastq.gz|fq.gz)";
	echo "# -b|--bed|--manifest=<FILE1,FILE2,...>       List of corresponding BED files in BED format.";
	echo "#                                             If not *.bed file, considered as Illumina manifest.";
	echo "#                                             Default first BED or empty";
	echo "# -j|--genes=<FILE1,FILE2,...>                List of corresponding GENES files in BED format.";
	echo "#                                             Default first GENES file or empty";
	echo "# -m|--transcripts=<FILE1,FILE2,...>          List of corresponding TRANSCRIPTS files in TSV format (transcript geneID).";
	echo "#                                             Default first GENES file or empty";
	echo "# -s|--sample=<STRING1,STRING2,...>           List of corresponding SAMPLE Name";
	echo "#                                             Automatically detected from input files.";
	echo "# -r|--run=<STRING1,STRING2,...>              List of corresponding RUN Name";
	echo "#                                             Default: date of first RUN of the list).";
	echo "# -o|--output|--results=<FOLDER>              OUTPUT directory to generate OUTPUT|RUN|SAMPLE|* files";
	echo "#                                             Default: Defined in ENV, or first --fastq file folder";
	echo "# -u|--repository=<FOLDER>                    Repository directory to generate REPOSITORY|RUN|SAMPLE|* specific files (TODO)";
	echo "#                                             Default: no copy in a repository";
	echo "# -p|--pipelines=<STRING1,STRING2,...>        PIPELINES to launch, in the format ALIGNER.CALLER.ANNOTATOR, separated by a comma";
	echo "#                                             Default: 'bwamem.gatkHC.howard'";
	echo "# -t|--threads=<INTEGER>                      Number of thread to use";
	echo "#                                             Default: all cores in your system minus one";
	echo "# -g|--by_sample                              Split analysis by SAMPLE, all threads on each sample, one by one";
	echo "# -a|--adapters=<STRING>                      Adapters to trim..";

	echo "# -v|--verbose                                VERBOSE option";
	echo "# -d|--debug                                  DEBUG option";
	echo "# -n|--release                                RELEASE option";
	echo "# -h|--help                                   HELP option";
	echo "#";
}

# header
header;


# USAGE exemple
# ./launch.sample.sh -f "$DATA/SAMPLE/HORIZON/HORIZON.R1.fastq.gz,$DATA/SAMPLE/HORIZON/HORIZON.bam,$DATA/SAMPLE/TEST/TEST.cram" -q "$DATA/SAMPLE/HORIZON/HORIZON.R2.fastq.gz" -s "H1,H2" -r "test" -b "/home1/TOOLS/data/SAMPLE/HORIZON/HORIZON.manifest,/home1/TOOLS/data/SAMPLE/HORIZON/HORIZON.bed,/home1/TOOLS/data/SAMPLE/TEST/TEST.bed" --output="$DATA/RESULTS"
# /home1/IRC/TOOLS/tools/stark/dev/STARK.sh -f "$DATA/SAMPLE/HORIZON/HORIZON.R1.fastq.gz,$DATA/SAMPLE/HORIZON/HORIZON.bam,$DATA/SAMPLE/TEST/TEST.cram" -q "$DATA/SAMPLE/HORIZON/HORIZON.R2.fastq.gz" -s "H1,H2" -r "testH,testH,testT" -b "/home1/TOOLS/data/SAMPLE/HORIZON/HORIZON.manifest,/home1/TOOLS/data/SAMPLE/HORIZON/HORIZON.bed,/home1/TOOLS/data/SAMPLE/TEST/TEST.bed" -e env.HUSTUMSOL.sh --output="$DATA/RESULTS" --repository="$DATA/REPOSITORY"
# /home1/IRC/TOOLS/tools/stark/dev/STARK.sh -f $DATA/SAMPLE/TEST/TEST.cram -b /home1/TOOLS/data/SAMPLE/TEST/TEST.bed --output="$DATA/RESULTS" -r testT

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "e:f:q:b:j:m:s:r:o:u:p:t:ga:vdnh" --long "env:,app:,application:,fastq:,fastq_R1:,fastq_R2:,bed:,manifest:,genes:,transcripts:,sample:,run:,output:,results:,repository:,pipelines:,threads:,by_sample,adapters:,verbose,debug,release,help" -- "$@" 2> /dev/null)
[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && usage && \
	exit 1
eval set -- "$ARGS"

PARAM=$@

while true
do
	case "$1" in
		-e|--env|--app|--application)
			APP="$2"
			shift 2
			;;
		-f|--fastq|--fastq_R1)
			FASTQ="$2"
			# transform RUNS list
			FASTQ=$(echo $FASTQ | tr "," " ")
			shift 2
			;;
		-q|--fastq_R2)
			FASTQ_R2="$2"
			# transform RUNS list
			FASTQ_R2=$(echo $FASTQ_R2 | tr "," " ")
			shift 2
			;;
		-b|--bed|--manifest)
			BED_INPUT="$2"
			BED_INPUT=$(echo $BED_INPUT | tr "," " ")
			shift 2
			;;
		-j|--genes)
			BEDFILE_GENES_INPUT="$2"
			BEDFILE_GENES_INPUT=$(echo $BEDFILE_GENES_INPUT | tr "," " ")
			shift 2
			;;
		-m|--transcripts)
			TRANSCRIPTS_INPUT="$2"
			TRANSCRIPTS_INPUT=$(echo $TRANSCRIPTS_INPUT | tr "," " ")
			shift 2
			;;
		-s|--sample)
			SAMPLE="$2"
			SAMPLE=$(echo $SAMPLE | tr "," " ")
			shift 2
			;;
		-r|--run)
			RUN="$2"
			RUN=$(echo $RUN | tr "," " ")
			shift 2
			;;
		-o|--output|--results)
			OUTPUT="$2"
			shift 2
			;;
		-u|--repository)
			REPOSITORY="$2"
			shift 2
			;;
		-p|--pipelines)
			PIPELINES_INPUT="$2"
			PIPELINES_INPUT=$(echo $PIPELINES_INPUT | tr "," " ")
			shift 2
			;;
		-t|--threads)
			THREADS_INPUT="$2"
			shift 2
			;;
		-g|--by_sample)
			BY_SAMPLE=1
			shift 1
			;;
		-a|--adapters)
			ADAPTERS="$2"
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
		*) 	echo "Option $1 is not recognized. " "Use -h or --help to display the help." && usage && \
			exit 1
			;;
	esac
done
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Checking the input parameter
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if [ "$FASTQ" == "" ]
then
	echo "Option --fastq is required. " "Use -h or --help to display the help." && usage && exit 1;
fi
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


COMMAND_COPY="rsync -aucqpAXoghi --no-links --no-perms --no-owner --no-group" # "cp -auv" or "rsync -auv" # auvpAXog
COMMAND_LINK="ln " # "cp -auv" or "rsync -auv" # auvpAXog
PERMS="a+rwx"


# FASTQ
for F in $FASTQ; do
	if [ "$F" != "" ] && [ -e "$F" ]; then
		F=$F;
		(($DEBUG)) && echo "F=$F";
		if ! (($(echo "$F" | grep ".fastq.gz$\|.fq.gz$\|.bam$\|.ubam$\|.cram$\|.ucram$" -c))); then
			echo "[ERROR]! Format of input file '$F' Unknown! Please check file format (.fastq.gz|.fq.gz|.bam|.ubam|.cram|.ucram)";
			exit 0;
		fi
	else
		echo "[ERROR] No input FASTQ '$F' file!";
		exit 0;
	fi;
done;


for F in $FASTQ_R2; do
	if [ "$F" != "" ] && [ -e "$F" ]; then
		F=$F;
		(($DEBUG)) && echo "F=$F";
		if ! (($(echo "$F" | grep ".fastq.gz$\|.fq.gz$\|.bam$\|.ubam$\|.cram$\|.ucram$" -c))); then
			echo "[ERROR]! Format of input file '$F' Unknown! Please check file format (.fastq.gz|.fq.gz|.bam|.ubam|.cram|.ucram)";
			exit 0;
		fi
	else
		echo "[ERROR] No input FASTQ R2 '$F' file!";
		exit 0;
	fi;
done;



ENV=$(find_app "$APP" "$STARK_FOLDER_APPS")

# SOURCE ENV if exists
if [ ! -z $ENV ] && [ -s $ENV ]; then
	echo "#[INFO] APPLICATION file '"$(echo $ENV | sed s#$STARK_FOLDER_APPS#APPS#)"' found."
	source $ENV;
else
	echo "#[ERROR] NO APPLICATION file '$ENV' found."
	exit 1
fi;


# FASTQ
#if [ "$FASTQ" != "" ] && [ -e "$FASTQ" ]; then
#	FASTQ=$FASTQ;
#	if ! (($(echo "$FASTQ" | grep ".fastq.gz$\|.fq.gz$\|.bam$\|.ubam$\|.cram$\|.ucram$" -c))); then
#		echo "[ERROR] Format of input file '$FASTQ' Unknown! Please check file format (.fastq.gz|.fq.gz|.bam|.ubam|.cram|.ucram)";
#		exit 0;
#	fi
#else
#	echo "[ERROR] No input FASTQ '$FASTQ' file!";
#	exit 0;
#fi;


# FASTQ_R2
#if [ ! -e $FASTQ_R2 ]; then
#	FASTQ_R2="";
#fi;

# BED
BED_DEFAULT="";
if [ "$BED_INPUT" != "" ]; then
	# Test BED files exist
	for B in $BED_INPUT; do
		if [ ! -e $B ]; then
			echo "[ERROR] BED '$B' file DOES NOT exist!";
			exit 0;
		fi;
	done;
	BED_DEFAULT=$(echo $BED_INPUT | cut -d" " -f1) #`date '+%Y%m%d-%H%M%S'`

fi;
# Test BED default file exists if not null
#if [ "$BED_DEFAULT" != "" ] && [ ! -e "$BED_DEFAULT" ]; then
#	echo "[ERROR] BED '$BED_DEFAULT' file DOES NOT exist!";
#	exit 0;
#fi;

# BEDFILE_GENES
if [ "$BEDFILE_GENES_INPUT" != "" ] && [ -s $BEDFILE_GENES_INPUT ]; then
	BEDFILE_GENES=$BEDFILE_GENES_INPUT;
fi;

BEDFILE_GENES_DEFAULT="";
if [ "$BEDFILE_GENES" != "" ]; then
	# Test BED files exist
	for G in $BEDFILE_GENES; do
		if [ ! -e $G ]; then
			echo "[ERROR] BEDFILE_GENES '$G' file DOES NOT exist!";
			exit 0;
		fi;
	done;
	BEDFILE_GENES_DEFAULT=$(echo $BEDFILE_GENES | cut -d" " -f1) #`date '+%Y%m%d-%H%M%S'`
fi;

# TRANSCRIPTS
if [ "$TRANSCRIPTS_INPUT" != "" ] && [ -s $TRANSCRIPTS_INPUT ]; then
	TRANSCRIPTS=$TRANSCRIPTS_INPUT;
fi;

TRANSCRIPTS_DEFAULT="";
if [ "$TRANSCRIPTS" != "" ]; then
	# Test BED files exist
	for T in $TRANSCRIPTS; do
		if [ ! -e $T ]; then
			echo "[ERROR] TRANSCRIPTS '$T' file DOES NOT exist!";
			exit 0;
		fi;
	done;
	TRANSCRIPTS_DEFAULT=$(echo $TRANSCRIPTS | cut -d" " -f1) #`date '+%Y%m%d-%H%M%S'`
fi;


# RUN
RUN_DEFAULT=`date '+%Y%m%d'`
if [ "$RUN" != "" ]; then
	RUN_DEFAULT=$(echo $RUN | cut -d" " -f1) #`date '+%Y%m%d-%H%M%S'`
fi;

# SAMPLE
# Auto Sample Name
#if [ "$SAMPLE" == "" ]; then

NB_SAMPLE=0
FASTQ_R2_ARRAY=($FASTQ_R2);
FASTQ_R2_CHECKED=""
SAMPLE_CHECKED=""
SAMPLE_ARRAY=($SAMPLE);
RUN_CHECKED=""
RUN_ARRAY=($RUN);
BED_CHECKED=""
BED_ARRAY=($BED);
BEDFILE_GENES_CHECKED=""
BEDFILE_GENES_ARRAY=($BEDFILE_GENES);
TRANSCRIPTS_CHECKED=""
TRANSCRIPTS_ARRAY=($TRANSCRIPTS);


for F in $FASTQ; do
	# Test FASTQ exists
	if [ ! -e $F ]; then
		echo "[ERROR] Input file '$F' does NOT exist! Please check input file --fastq (.fastq.gz|.fq.gz|.bam|.ubam|.cram|.ucram|.sam|.usam)";
		exit 0;
	fi;
	# Test FASTQ format
	if ! (($(echo "$F" | grep ".fastq.gz$\|.fq.gz$\|.bam$\|.ubam$\|.cram$\|.ucram$\|.sam$\|.usam$" -c))); then
		echo "[ERROR] Format of input file '$F' Unknown! Please check file format (.fastq.gz|.fq.gz|.bam|.ubam|.cram|.ucram|.sam|.usam)";
		exit 0;
	fi;
	# test FASTQ R2 exists
	if [ "${FASTQ_R2_ARRAY[$NB_SAMPLE]}" != "" ] && [ ! -e "${FASTQ_R2_ARRAY[$NB_SAMPLE]}" ]; then
		echo "[ERROR] Input file Read2 '${FASTQ_R2_ARRAY[$NB_SAMPLE]}' does NOT exist! Please check input file --fastq_R2 (.fastq.gz|.fq.gz)";
		exit 0;
	fi;
	# Test FASTQ R2 format if exists
	if [ "${FASTQ_R2_ARRAY[$NB_SAMPLE]}" != "" ]; then
		if ! (($(echo "${FASTQ_R2_ARRAY[$NB_SAMPLE]}" | grep ".fastq.gz$\|.fq.gz$" -c))); then
			echo "[ERROR] Format of input file Read2 '${FASTQ_R2_ARRAY[$NB_SAMPLE]}' Unknown! Please check file format (.fastq.gz|.fq.gz)";
			exit 0;
		fi;
	fi;
	# SAMPLE
	if [ "${SAMPLE_ARRAY[$NB_SAMPLE]}" == "" ]; then
		#F_PATTERN=`basename $F | tr "." "_" | sed 's/\.R[12]\.\|_R[12]_\|_[12]_\|\.[12]\./_/gi' `
		F_PATTERN=$(basename $F | sed "s/.fastq.gz$\|.fq.gz$\|.bam$\|.ubam$\|.cram$\|.ucram$\|.sam$\|.usam$//g" | tr "." "_" | sed "s/_S[0-9]*_R[12]_[0-9]*//gi")

		# IF SAMPLE name exist
		# try with the filename (removing FASTQ pattern filename) KEEP EXTENSION
		if (($(echo $SAMPLE_CHECKED | tr " " "\n" | grep $F_PATTERN -c))); then
			F_PATTERN=$(basename $F | tr "." "_" | sed "s/_S[0-9]*_R[12]_[0-9]*//gi")
		fi;
		# try with the filename (NOT removing FASTQ pattern filename) KEEP EXTENSION and fastq pattern
		if (($(echo $SAMPLE_CHECKED | tr " " "\n" | grep $F_PATTERN -c))); then
			F_PATTERN=$(basename $F | tr "." "_")
		fi;
		# try with the filename (NOT removing FASTQ pattern filename) KEEP EXTENSION and fastq pattern WITH number
		I=1
		while (($(echo $SAMPLE_CHECKED | tr " " "\n" | grep $F_PATTERN -c))); do
			F_PATTERN=$F_PATTERN"_"$I
			((I++))
		done
		SAMPLE_CHECKED="$SAMPLE_CHECKED$F_PATTERN "
		# | sed "s/.fastq.gz$\|.fq.gz$\|.bam$\|.ubam$\|.cram$\|.ucram$\|.sam$\|.usam$//g" | sed "s/_S[0-9]*_R[12]_[0-9]*//gi"
		# | tr "." "_" | sed "s/.fastq.gz$\|.fq.gz$\|.bam$\|.ubam$\|.cram$\|.ucram$\|.sam$\|.usam$//g" | sed "s/_S[0-9]*_R[12]_[0-9]*//gi"
		# echo "VBA-29_1_S1_R2_001.fastq.gz" | sed "s/.fastq.gz$\|.fq.gz$\|.bam$\|.ubam$\|.cram$\|.ucram$\|.sam$\|.usam$//g" | tr "." "_" | sed "s/_[0-9]\{3\}\|_R[12]\|_S[0-9]//gi"
		#SAMPLE_CHECKED="$SAMPLE_CHECKED$F_PATTERN "
	else
		SAMPLE_CHECKED="$SAMPLE_CHECKED${SAMPLE_ARRAY[$NB_SAMPLE]} "
	fi;

	# RUN
	if [ "${RUN_ARRAY[$NB_SAMPLE]}" == "" ]; then
		RUN_CHECKED="$RUN_CHECKED$RUN_DEFAULT "
	else
		RUN_CHECKED="$RUN_CHECKED${RUN_ARRAY[$NB_SAMPLE]} "
	fi;

	# BED
	if [ "${BED_ARRAY[$NB_SAMPLE]}" == "" ] || [ ! -e "${BED_ARRAY[$NB_SAMPLE]}" ]; then
		BED_CHECKED="$BED_CHECKED$BED_DEFAULT "
	else
		BED_CHECKED="$BED_CHECKED${BED_ARRAY[$NB_SAMPLE]} "
	fi;

	# BEDFILE_GENES
	if [ "${BEDFILE_GENES_ARRAY[$NB_SAMPLE]}" == "" ] || [ ! -e "${BEDFILE_GENES_ARRAY[$NB_SAMPLE]}" ]; then
		BEDFILE_GENES_CHECKED="$BEDFILE_GENES_CHECKED$BEDFILE_GENES_DEFAULT "
	else
		BEDFILE_GENES_CHECKED="$BEDFILE_GENES_CHECKED${BEDFILE_GENES_ARRAY[$NB_SAMPLE]} "
	fi;

	# TRANSCRIPTS
	if [ "${TRANSCRIPTS_ARRAY[$NB_SAMPLE]}" == "" ] || [ ! -e "${TRANSCRIPTS_ARRAY[$NB_SAMPLE]}" ]; then
		TRANSCRIPTS_CHECKED="$TRANSCRIPTS_CHECKED$TRANSCRIPTS_DEFAULT "
	else
		TRANSCRIPTS_CHECKED="$TRANSCRIPTS_CHECKED${TRANSCRIPTS_ARRAY[$NB_SAMPLE]} "
	fi;

	((NB_SAMPLE++))
done;


SAMPLE=$SAMPLE_CHECKED;
RUN=$RUN_CHECKED;
BED=$BED_CHECKED;
BEDFILE_GENES=$BEDFILE_GENES_CHECKED;
TRANSCRIPTS=$TRANSCRIPTS_CHECKED;
#echo -e  "RUN $RUN \n SAMPLE $SAMPLE \n BED $BED" | column -t; exit

#fi;
# Default Sample Name
#if [ "$SAMPLE" == "" ]; then
#	SAMPLE=`date '+%Y%m%d-%H%M%S'`
#fi;

# ARRAYS
FASTQ_ARRAY=($FASTQ);
FASTQ_R2_ARRAY=($FASTQ_R2);
SAMPLE_ARRAY=($SAMPLE);
RUN_ARRAY=($RUN);
BED_ARRAY=($BED);
BEDFILE_GENES_ARRAY=($BEDFILE_GENES);
TRANSCRIPTS_ARRAY=($TRANSCRIPTS);




# MULTIPLE BAM/CRAM/SAM/FASTQ
if (($BY_SAMPLE)); then
	if [ $(echo $FASTQ | wc -w ) -gt 1 ]; then
		echo "# Multiple SAMPLES: $SAMPLE";
		#FASTQ_ARRAY=($FASTQ);
		#FASTQ_R2_ARRAY=($FASTQ_R2);
		#SAMPLE_ARRAY=($SAMPLE);
		PARAM_MULTI=$(echo $PARAM | sed "s/--by_sample[^ |$]*//gi"  | sed "s/-g[^ |$]*//gi" | sed s/--$//gi)
		i=0
		for f in ${FASTQ_ARRAY[@]}; do
			#for i in $RUNS_ARRAY; do
			q_var=""
			if [ ! -z ${FASTQ_R2_ARRAY[$i]} ]; then q_var=" -q ${FASTQ_R2_ARRAY[$i]}"; fi;
			#echo "$SCRIPT_DIR/$(basename $0) $PARAM_MULTI --fastq=$f --fastq_R2="${FASTQ_R2_ARRAY[$i]}" --sample="${SAMPLE_ARRAY[$i]}" --run="${RUN_ARRAY[$i]}" --bed="${BED_ARRAY[$i]}" ;"
			#$SCRIPT_DIR/$(basename $0) $PARAM_MULTI --fastq=$f --fastq_R2="${FASTQ_R2_ARRAY[$i]}" --sample="${SAMPLE_ARRAY[$i]}" --run="${RUN_ARRAY[$i]}" --bed="${BED_ARRAY[$i]}" --bedfile_genes="${BEDFILE_GENES_ARRAY[$i]}"  --transcripts="${TRANSCRIPTS_ARRAY[$i]}" ;
			$SCRIPT_DIR/$(basename $0) $PARAM_MULTI --fastq=$f --fastq_R2="${FASTQ_R2_ARRAY[$i]}" --sample="${SAMPLE_ARRAY[$i]}" --run="${RUN_ARRAY[$i]}" --bed="${BED_ARRAY[$i]}" --genes="${BEDFILE_GENES_ARRAY[$i]}"  --transcripts="${TRANSCRIPTS_ARRAY[$i]}" ;
			((i++))
		done;
		exit 0;
	fi;
fi;
#echo "# Analysis of '$FASTQ'/'$FASTQ_R2'"






# PIPELIENS PARAMETERS
#ALIGNERS="bwamem"
#CALLERS="gatkHC"
#ANNOTATORS="vap"
if [ "$PIPELINES_INPUT" != "" ]; then
	PIPELINES=$PIPELINES_INPUT
fi;
if [ "$PIPELINES" == "" ]; then
	PIPELINES="bwamem.gatkHC.howard"
fi;

if ((0)); then
if [ "$PIPELINES" == "" ]; then
	if [ "$PIPELINES_INPUT" == "" ]; then
		#PIPELINES="bwamem.gatkHC.vap bwamem.gatkUG.vap"
		#PIPELINES="bwamem.gatkUG.howard"
		PIPELINES="bwamem.gatkHC.howard"
		echo "#[WARNING] NO PIPELINES defined. Default PIPELINES '$PIPELINES' used."
	else
		PIPELINES=$PIPELINES_INPUT
	fi;
fi;
fi;
#DATE
DATE_DAY=`date '+%Y%m%d'`
DATE_MIN=`date '+%Y%m%d-%H%M%S'`
ANALYSIS_REF=$DATE_MIN

#CORES=$(echo `ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w`" - 1 " | bc)
#THREADS=$CORES
#THREADS=1
re='^[0-9]+$'
CORES=$(ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w)
if ! [[ $THREADS =~ $re ]] || [ -z "$THREADS" ] || [ "$THREADS" == "" ] || [ $THREADS -gt $CORES ] ; then
	CORES_FREE=1
	THREADS=$(($CORES-$CORES_FREE))
fi;

# THREADS
if [[ $THREADS_INPUT =~ $re ]] && [ "$THREADS_INPUT" != "" ]; then
	THREADS=$THREADS_INPUT;
fi;

# TMP_SYS_FOLDER
if [ ! -d $TMP_SYS_FOLDER ]; then
	TMP_SYS_FOLDER=/tmp
fi;


#echo "THREADS_INPUT=$THREADS_INPUT"
#echo "THREADS=$THREADS"
#exit 0;

# INPUT
# Input is the folder of the first FASTQ by default
for F in $(echo $FASTQ | cut -d " " -f1); do
	INPUT=`dirname "$F"`
done


# DEFAULT OUTPUT
# OUTPUT is 1/ the OUTPUT folder in option, OR 2/ FOLDER_RESULT in ENV, OR 3/ folder of the first FASTQ by default
if [ -z "$OUTPUT" ]; then
	#OUTPUT=$FOLDER_RESULTS
	OUTPUT=$RESULTS_FOLDER
fi;
if [ -z "$OUTPUT" ]; then
	OUTPUT=$INPUT
fi;

#echo "$FASTQ | $INPUT | $OUTPUT"; exit 0;


# CREATE OUTPUT if necessary
if [ ! -d "$OUTPUT" ]; then
	mkdir -p $OUTPUT;
fi

# JAVA OPTIONS
if [ "$JAVA_FLAGS" == "" ]; then
	JAVA_FLAGS=" -Xmx"$JAVA_MEMORY"g  -Dsamjdk.try_use_intel_deflater=false -Dsnappy.disable=true -Dorg.xerial.snappy.tempdir="$TMP_FOLDER_TMP" -Djava.io.tmpdir="$TMP_FOLDER_TMP;
fi;




#echo "ENV=$ENV";
#APP_NAME_DEF=$($SCRIPT_DIR/extract_variable_from_env.sh $ENV "APP_NAME")
APP_NAME_DEF=$(source $ENV; echo $APP_NAME)
if [ -z "$APP_NAME" ]; then
	APP_NAME=$(echo $(basename $ENV) | sed "s/^env.//gi" | sed "s/.sh$//gi" | sed "s/sh$//gi")
fi;
if [ ! -z "$GROUP" ]; then
	SAMPLE_GROUP=$GROUP
else
	SAMPLE_GROUP=$(echo $APP_NAME | awk -F- '{print $1}');
fi;
if [ "$SAMPLE_GROUP" == "" ]; then SAMPLE_GROUP="UNKNOWN"; fi;
if [ ! -z "$PROJECT" ]; then
	SAMPLE_PROJECT=$PROJECT
else
	SAMPLE_PROJECT=$(echo $APP_NAME | awk -F- '{print $2}');
fi;
if [ "$SAMPLE_PROJECT" == "" ]; then SAMPLE_PROJECT="UNKNOWN"; fi;

#SAMPLE_PROJECT=$(echo $APP_NAME | awk -F- '{print $2}'); if [ "$SAMPLE_PROJECT" == "" ]; then SAMPLE_PROJECT="UNKNOWN"; fi;
SAMPLE_USER=$(echo $APP_NAME | awk -F- '{print $3}'); if [ "$SAMPLE_USER" == "" ]; then SAMPLE_USER="UNKNOWN"; fi;
if [ -z "$APP_NAME" ]; then APP_NAME="UNKNOWN"; fi;

#echo "ENV=$ENV";
#echo "APP=$APP_NAME";
#echo "APP_DEF=$APP_NAME_DEF";
#echo "GROUP=$GROUP";
#echo "PROJECT=$PROJECT";
#echo -e "APP $APP_NAME \n SAMPLE_GROUP $SAMPLE_GROUP \n SAMPLE_PROJECT $SAMPLE_PROJECT \n SAMPLE_USER $SAMPLE_USER " | column -t;
#exit 0;
#mkdir -p $RUN_DIR


#echo -e "FASTQ=$FASTQ\nFASTQ_R2=$FASTQ_R2\nSAMPLE=$SAMPLE\nRUN=$RUN\nINPUT$INPUT\nOUTPUT=$OUTPUT\n\n$PIPELINES\n$BED";
if (($VERBOSE)); then
	echo -e "FASTQ $FASTQ\nFASTQ_R2 $FASTQ_R2\nSAMPLE $SAMPLE\nRUN $RUN\nBED $BED\nGENES $BEDFILE_GENES\nTRANSCRIPTS $TRANSCRIPTS\nINPUT $INPUT\nOUTPUT $OUTPUT\n\nPIPELINES $PIPELINES" | column -t;
fi;

SAMPLE_ARRAY=($SAMPLE);
FASTQ_ARRAY=($FASTQ);
FASTQ_R2_ARRAY=($FASTQ_R2);
RUN_ARRAY=($RUN);
BED_ARRAY=($BED);
BEDFILE_GENES_ARRAY=($BEDFILE_GENES);
TRANSCRIPTS_ARRAY=($TRANSCRIPTS);

# FOREACH RUN
RUN_UNIQ=$(echo $RUN | tr " " "\n" | sort | uniq | tr "\n" " ");
#echo $RUN_UNIQ
#exit 0;


for RUU in $RUN_UNIQ; do


	#RUN_DIR=$OUTPUT/$RUN #$OUTPUT_RES_DIR/$RUN
	MAKEFILE_ANALYSIS_RUN=$OUTPUT/$RUU/analysis.V$ANALYSIS_REF.param.mk
	SHELL_ANALYSIS_RUN=$OUTPUT/$RUU/analysis.V$ANALYSIS_REF.param.sh
	LOGFILE_RES_RUN=$OUTPUT/$RUU/analysis.V$ANALYSIS_REF.log
	LOGFILE_RES_RUN_REPORT=$OUTPUT/$RUU/analysis.V$ANALYSIS_REF.report.log
	#LOGFILE_RES_RUN=$OUTPUT/$RUU/run_analysis.log
	#LOGFILE_RES_RUN_REPORT=$OUTPUT/$RUU/run_analysis.report.log
	RELEASE_RUN=$OUTPUT/$RUU/analysis.V$ANALYSIS_REF.release
	FINAL_REPORT_RUN=$OUTPUT/$RUU/analysis.V$ANALYSIS_REF.report
	FINAL_REPORT_FULL_RUN=$OUTPUT/$RUU/analysis.V$ANALYSIS_REF.full.report
	FINAL_REPORT_FULL_VCF=$OUTPUT/$RUU/$SAMPLE.V$ANALYSIS_REF.full.vcf

	# MKDIR & TOUCH
	mkdir -p $OUTPUT/$RUU
	touch $MAKEFILE_ANALYSIS_RUN
	touch $SHELL_ANALYSIS_RUN
	touch $LOGFILE_RES_RUN
	touch $LOGFILE_RES_RUN_REPORT


	echo "#[INFO] RUN '$RUU'"
	echo "## ANALYSIS" > $MAKEFILE_ANALYSIS_RUN
	echo "## ANALYSIS" > $SHELL_ANALYSIS_RUN
	RUNS_SAMPLES="";

	F_LIST=""
	Q_LIST=""
	S_LIST=""
	B_LIST=""
	G_LIST=""
	T_LIST=""


	I=0;
	NB_SAMPLE=0;
	for F in $FASTQ; do


		S=${SAMPLE_ARRAY[$I]};
		RU=${RUN_ARRAY[$I]};
		F_R2=${FASTQ_R2_ARRAY[$I]};
		B=${BED_ARRAY[$I]};
		G=${BEDFILE_GENES_ARRAY[$I]};
		T=${TRANSCRIPTS_ARRAY[$I]};

		#echo "$S $F $F_R2"; exit 0;

		# BY RUN
		if [ "$RU" != "$RUU" ]; then
			((I++))
			continue;
		fi;

		F_LIST="$F_LIST$F "
		Q_LIST="$Q_LIST$F_R2 "
		S_LIST="$S_LIST$S "
		B_LIST="$B_LIST$B "
		G_LIST="$G_LIST$G "
		T_LIST="$T_LIST$T "

		RUN_SAMPLE_DIR=$OUTPUT/$RU/$S

		if (($VERBOSE)); then
			echo "$F | $F_R2 | $S | $RU | $B | $PIPELINES | $INPUT | $OUTPUT | $RUN_SAMPLE_DIR ";
			#((I++)); continue;
		fi;

		mkdir -p $RUN_SAMPLE_DIR

		# INFOS
		echo "#[INFO] SAMPLE '$RU/$S' from file(s):"
		echo "#[INFO] $F $F_R2 $B $G"

		# Copy FASTQ
		PICARD_FLAGS="COMPRESSION_LEVEL=1 MAX_RECORDS_IN_RAM=500000"
		PICARD_UNALIGNED_FLAGS="SORT_ORDER=coordinate "
		PICARD_UNALIGNED_NAME_FLAGS="LIBRARY_NAME=001 PLATFORM=ILLUMINA PLATFORM_UNIT=PU READ_GROUP_NAME=A SAMPLE_NAME=$S "
		PICARD_VALIDATE_FLAGS=" MODE=SUMMARY "
		PICARD_VALIDATE_NAME_FLAGS="  "
		PICARD_ADDORREPLACEREADGROUPS_FLAGS="  "
		PICARD_ADDORREPLACEREADGROUPS_NAME_FLAGS=" RGLB=001 RGPL=ILLUMINA RGPU=PU RGSM=$S VALIDATION_STRINGENCY=SILENT "


		# DATA entry format is FASTQ/BAM/CRAM/SAM
		#if ((1)); then
		if [ ! -s $RUN_SAMPLE_DIR/$S.R1.fastq.gz ]; then
			echo "#[INFO] Input file '$RUN_SAMPLE_DIR/$S.R1.fastq.gz' does NOT exist"
			# FASTQ
			if (($(echo "$F" | grep ".fastq.gz$\|.fq.gz$" -c))); then
				echo "#[INFO] Create Input data from FASTQ file(s)"
				if [ "$ADAPTERS" != "" ]; then
					echo "#[INFO] Trim adapters."
					if [ -s $F ]; then
						if [ -s $F_R2 ]; then
							$JAVA -jar $TRIMMOMATIC PE -phred33 $F $F_R2 $F"_paired.fq.gz" $F"_unpaired.fq.gz" $F_R2"_paired.fq.gz" $F_R2"_unpaired.fq.gz" ILLUMINACLIP:$ADAPTERS:2:30:10 1>/dev/null 2>/dev/null #LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
							cat $F"_paired.fq.gz" $F"_unpaired.fq.gz" $F_R2"_paired.fq.gz" $F_R2"_unpaired.fq.gz" > $RUN_SAMPLE_DIR/$S.fastq.gz
							cat $F"_paired.fq.gz" $F"_unpaired.fq.gz" > $RUN_SAMPLE_DIR/$S.R1.fastq.gz
							cat $F_R2"_paired.fq.gz" $F_R2"_unpaired.fq.gz" > $RUN_SAMPLE_DIR/$S.R2.fastq.gz
						else
							$JAVA -jar $TRIMMOMATIC SE -phred33 $F $F"_paired.fq.gz" ILLUMINACLIP:$ADAPTERS:2:30:10 1>/dev/null 2>/dev/null #LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
							cat $F"_paired.fq.gz" > $RUN_SAMPLE_DIR/$S.R1.fastq.gz
						fi;

					fi;
				else
					# Create FASTQ
					#echo "#[INFO] FASTQ R1"
					#cat $F > $RUN_SAMPLE_DIR/$S.R1.fastq.gz;
					cp $F $RUN_SAMPLE_DIR/$S.R1.fastq.gz;

					#if [ -z $F_R2 ] && [ -e $F_R2 ]; then
					if [ -s "$F_R2" ]; then # && [ -e $F_R2 ]; then
						#echo "#[INFO] FASTQ R2 '$F_R2'"
						#exit 0;
						#cat $F_R2 >> $RUN_SAMPLE_DIR/$S.R2.fastq.gz;
						cp $F_R2 $RUN_SAMPLE_DIR/$S.R2.fastq.gz
					else
						touch $RUN_SAMPLE_DIR/$S.R2.fastq.gz
					fi;
				fi;

			# BAM/CRAM/SAM
			elif (($(echo $F | grep ".bam$\|.ubam$\|.cram$\|.ucram\|.sam$\|.usam$" -c))); then
				echo "#[INFO] Create Input data from BAM/CRAM/SAM file"
				# Other command with BEDTOOLS: bamToFastq [OPTIONS] -i <BAM> -fq <FASTQ> -fq2 <FASTQR2> -tags
				$SAMTOOLS bam2fq $F -1 $RUN_SAMPLE_DIR/$S.R1.fastq -2 $RUN_SAMPLE_DIR/$S.R2.fastq 1>$RUN_SAMPLE_DIR/$S.R0.fastq 2>/dev/null;
				if [ -s $RUN_SAMPLE_DIR/$S.R0.fastq ]; then
					cat $RUN_SAMPLE_DIR/$S.R2.fastq $RUN_SAMPLE_DIR/$S.R0.fastq >> $RUN_SAMPLE_DIR/$S.R1.fastq;
					> $RUN_SAMPLE_DIR/$S.R2.fastq;

				fi;
				$COMPRESS $RUN_SAMPLE_DIR/$S.R1.fastq --fast -f -q 1>/dev/null 2>/dev/null;
				$COMPRESS $RUN_SAMPLE_DIR/$S.R2.fastq --fast -f -q 1>/dev/null 2>/dev/null;
				rm -f $RUN_SAMPLE_DIR/$S.R0.fastq;
				if ((0)); then
				if (($($SAMTOOLS view -c -f 1 $F -@ $THREADS))); then # Paired-End
					#$SAMTOOLS bam2fq $F -1 $RUN_SAMPLE_DIR/$S.R1.fastq -2 $RUN_SAMPLE_DIR/$S.R2.fastq 1>$RUN_SAMPLE_DIR/$S.OTHERS.fastq 2>/dev/null;
					$SAMTOOLS bam2fq $F -2 $RUN_SAMPLE_DIR/$S.R2.fastq 1>$RUN_SAMPLE_DIR/$S.R1.fastq 2>/dev/null;
					#cat $RUN_SAMPLE_DIR/$S.OTHERS.fastq >> $RUN_SAMPLE_DIR/$S.R1.fastq
					$COMPRESS $RUN_SAMPLE_DIR/$S.R1.fastq --fast -f -q 1>/dev/null 2>/dev/null;
					$COMPRESS $RUN_SAMPLE_DIR/$S.R2.fastq --fast -f -q 1>/dev/null 2>/dev/null;
					rm -f $RUN_SAMPLE_DIR/$S.OTHERS.fastq
				else # Single-End
					$SAMTOOLS bam2fq $F > $RUN_SAMPLE_DIR/$S.R1.fastq 2>/dev/null;
					$COMPRESS $RUN_SAMPLE_DIR/$S.R1.fastq --fast -f -q 1>/dev/null 2>/dev/null;
				fi;
				fi;
			fi;
		else
			echo "#[INFO] Input file '$RUN_SAMPLE_DIR/$S.R1.fastq.gz' DOES exist"
		fi;



		# Copy BED
		if [[ $B =~ .bed$ ]]; then
			if [ -e $B ] && [ "$B" != "" ] && [ ! -e $RUN_SAMPLE_DIR/$S.bed ]; then
				echo "#[INFO] Copy BED file."
				cp -p $B $RUN_SAMPLE_DIR/$S.bed;
				cp -p $B $RUN_SAMPLE_DIR/$S.metrics.bed;
				touch $RUN_SAMPLE_DIR/$S.manifest;
				touch $RUN_SAMPLE_DIR/$S.manifest -r $B;
				touch $RUN_SAMPLE_DIR/$S.metrics.bed -r $B;
				echo -e $(basename $B)"\tfrom option - bed" > $RUN_SAMPLE_DIR/$S.bed_name
				echo -e $(basename $B)"\tfrom option - bed" > $RUN_SAMPLE_DIR/$S.manifest_name
				#export BED=$RUN_SAMPLE_DIR/$S.bed
			fi;
		else
			if [ -e $B ] && [ "$B" != "" ] && [ ! -e $RUN_SAMPLE_DIR/$S.manifest ]; then
				echo "#[INFO] Copy Manifest file."
				cp -p $B $RUN_SAMPLE_DIR/$S.manifest;
				touch $RUN_SAMPLE_DIR/$S.manifest -r $B;
				#echo -e "\tfrom manifest file" > $RUN_SAMPLE_DIR/$S.bed_name
				echo -e $(basename $B)"\tfrom option - manifest" > $RUN_SAMPLE_DIR/$S.bed_name
				echo -e $(basename $B)"\tfrom option - manifest" > $RUN_SAMPLE_DIR/$S.manifest_name
				#export MANIFEST=$RUN_SAMPLE_DIR/$S.manifest;
			fi;
		fi;

		# Copy BEDFILE_GENES
		if [ -e $G ] && [ "$G" != "" ] && [ ! -e $RUN_SAMPLE_DIR/$S.genes ]; then
			echo "#[INFO] Copy GENES file."
			cp -p $G $RUN_SAMPLE_DIR/$S.genes;
		fi;



		# Copy TRANSCRIPTS

		if [ -e $T ] && [ "$T" != "" ] && [ ! -e $RUN_SAMPLE_DIR/$S.transcripts ]; then
			echo "#[INFO] Copy TRANSCRIPTS file."
			cp -p $T $RUN_SAMPLE_DIR/$S.transcripts;
		fi;
		#echo "$T"; exit 0;

		# SampleSheet
		if ((1)); then
		if [ ! -e $RUN_SAMPLE_DIR/$S.SampleSheet.csv ]; then
			echo "#[INFO] Copy SampleSheet."
			touch $RUN_SAMPLE_DIR/$S.SampleSheet.csv;
			if [ -e $RUN_SAMPLE_DIR/$S.manifest ]; then
				touch -f $RUN_SAMPLE_DIR/$S.SampleSheet.csv -r $RUN_SAMPLE_DIR/$S.manifest;
			fi;
		fi;
		fi;


		# MAKEFILE_ANALYSIS_RUN
		RUNS_SAMPLES="$RUNS_SAMPLES$RU:$S "


		((I++))
		((NB_SAMPLE++))

	done;


	# THREADS
	############

	#NB_SAMPLE=1;
	#PIPELINES=$(echo $PIPELINES | tr " " "\n" | sort | uniq | tr "\n" " ")
	PIPELINES=$(echo $PIPELINES | tr "," "\n" | tr " " "\n" | awk '!x[$0]++' | tr "\n" " ")
	NB_PIPELINES=$(echo $PIPELINES | wc -w)
	NB_ALIGNERS=$(echo -e $(for P in $PIPELINES; do echo $P | cut -d. -f1; done) | tr " " "\n" | awk '!x[$0]++' | tr "\n" " " | wc -w)
	NB_CALLERS=$(echo -e $(for P in $PIPELINES; do echo $P | cut -d. -f2; done) | tr " " "\n" | awk '!x[$0]++' | tr "\n" " " | wc -w)

	# ADJUSTING threads parameters
	#ROUND="FLOOR" # CEIL or FLOOR;
	ROUND="CEIL" # CEIL or FLOOR;
	if [ "$ROUND" == "CEIL" ]; then
		THREADS_BY_SAMPLE=$( printf "%.0f" $(echo "scale=2;$THREADS/$NB_SAMPLE" | bc))			# Number of threads by sample
		THREADS_BY_PIPELINE=$( printf "%.0f" $(echo "scale=2;$THREADS_BY_SAMPLE/$NB_PIPELINES" | bc))	# Number of threads for a pipeline's command
		THREADS_BY_ALIGNER=$( printf "%.0f" $(echo "scale=2;$THREADS_BY_SAMPLE/$NB_ALIGNERS" | bc))	# Number of threads for a aligner's command
		THREADS_BY_CALLER=$( printf "%.0f" $(echo "scale=2;$THREADS_BY_SAMPLE/$NB_CALLERS" | bc))	# Number of threads for a caller's command
		THREADS_BWA=$THREADS_BY_ALIGNER									# NB of threads for BWA command
		THREADS_RTC=$THREADS_BY_ALIGNER									# NB of threads for GATK RealignerTargetCreator
	elif  [ "$ROUND" == "FLOOR" ]; then
		THREADS_BY_SAMPLE=$(($THREADS/$NB_SAMPLE))			# Number of threads by sample
		THREADS_BY_PIPELINE=$(($THREADS_BY_SAMPLE/$NB_PIPELINES))	# Number of threads for a pipeline's command
		THREADS_BY_ALIGNER=$(($THREADS_BY_SAMPLE/$NB_ALIGNERS))		# Number of threads for a aligner's command
		THREADS_BY_CALLER=$(($THREADS_BY_SAMPLE/$NB_CALLERS))		# Number of threads for a caller's command
		THREADS_BWA=$THREADS_BY_ALIGNER					# NB of threads for BWA command
		THREADS_RTC=$THREADS_BY_ALIGNER					# NB of threads for GATK RealignerTargetCreator
	else
		THREADS_BY_SAMPLE=1			# Number of threads by sample
		THREADS_BY_PIPELINE=1			# Number of threads for a pipeline's command
		THREADS_BY_ALIGNER=1			# Number of threads for a aligner's command
		THREADS_BY_CALLER=1			# Number of threads for a caller's command
		THREADS_BWA=$THREADS_BY_ALIGNER		# NB of threads for BWA command
		THREADS_RTC=$THREADS_BY_ALIGNER		# NB of threads for GATK RealignerTargetCreator
	fi;
	if [ "$THREADS_BY_SAMPLE" == "0" ]; then THREADS_BY_SAMPLE=1; fi;
	if [ "$THREADS_BY_ALIGNER" == "0" ]; then THREADS_BY_ALIGNER=1; fi;
	if [ "$THREADS_BY_CALLER" == "0" ]; then THREADS_BY_CALLER=1; fi;
	if [ "$THREADS_BY_PIPELINE" == "0" ]; then THREADS_BY_PIPELINE=1; fi;
	if [ "$THREADS_BWA" == "0" ]; then THREADS_BWA=1; fi;
	if [ "$THREADS_RTC" == "0" ]; then THREADS_RTC=1; fi;

	# THREAD_PARAMETERS
	THREAD_PARAMETERS=$THREAD_PARAMETERS" THREADS_BY_SAMPLE=$THREADS_BY_SAMPLE "
	THREAD_PARAMETERS=$THREAD_PARAMETERS" THREADS_BY_PIPELINE=$THREADS_BY_PIPELINE "
	THREAD_PARAMETERS=$THREAD_PARAMETERS" THREADS_BY_ALIGNER=$THREADS_BY_ALIGNER "
	THREAD_PARAMETERS=$THREAD_PARAMETERS" THREADS_BY_CALLER=$THREADS_BY_CALLER "
	THREAD_PARAMETERS=$THREAD_PARAMETERS" THREADS_BWA=$THREADS_BWA "
	THREAD_PARAMETERS=$THREAD_PARAMETERS" THREADS_RTC=$THREADS_RTC "


	# JAVA_MEMORY by sample and by caller
	###############

	# Default
	JAVA_MEMORY_BY_SAMPLE=$JAVA_MEMORY
	JAVA_MEMORY_BY_CALLER=$JAVA_MEMORY

	# Calcul
	JAVA_MEMORY_BY_SAMPLE=$(($MEMTOTAL/$NB_SAMPLE/1024/1024))				# Number of threads by sample
	JAVA_MEMORY_BY_PIPELINE=$(($JAVA_MEMORY_BY_SAMPLE/$NB_PIPELINES/1024/1024))	# Number of threads for a pipeline's command
	JAVA_MEMORY_BY_ALIGNER=$(($JAVA_MEMORY_BY_SAMPLE/$NB_ALIGNERS/1024/1024))	# Number of threads for a aligner's command
	JAVA_MEMORY_BY_CALLER=$(($JAVA_MEMORY_BY_SAMPLE/$NB_CALLERS/1024/1024))		# Number of threads for a caller's command

	# Test
	if [ "$JAVA_MEMORY_BY_SAMPLE" == "" ] || [ $JAVA_MEMORY_BY_SAMPLE -lt 1 ]; then JAVA_MEMORY_BY_SAMPLE=1; fi;
	if [ "$JAVA_MEMORY_BY_PIPELINE" == "" ] || [ $JAVA_MEMORY_BY_PIPELINE -lt 1 ]; then JAVA_MEMORY_BY_PIPELINE=1; fi;
	if [ "$JAVA_MEMORY_BY_ALIGNER" == "" ] || [ $JAVA_MEMORY_BY_ALIGNER -lt 1 ]; then JAVA_MEMORY_BY_ALIGNER=1; fi;
	if [ "$JAVA_MEMORY_BY_CALLER" == "" ] || [ $JAVA_MEMORY_BY_CALLER -lt 1 ]; then JAVA_MEMORY_BY_CALLER=1; fi;


	#JAVA_FLAGS_BY_SAMPLE=" -Xmx"$JAVA_MEMORY_BY_SAMPLE"g $JAVA_FLAGS_OTHER_PARAM $JAVA_FLAGS_TMP_FOLDER";


	echo "RUNS_SAMPLES=$RUNS_SAMPLES" >> $MAKEFILE_ANALYSIS_RUN
	echo "PIPELINES=$PIPELINES" >> $MAKEFILE_ANALYSIS_RUN
	echo "INTERSEC=2" >> $MAKEFILE_ANALYSIS_RUN

	echo "RUNS_SAMPLES=\"$RUNS_SAMPLES\"" >> $SHELL_ANALYSIS_RUN
	echo "PIPELINES=\"$PIPELINES\"" >> $SHELL_ANALYSIS_RUN
	echo "INTERSEC=2" >> $SHELL_ANALYSIS_RUN


	#echo "BEDFILE_GENES $BEDFILE_GENES"
	#echo "TRANSCRIPTS $TRANSCRIPTS"
	#echo "PARAMETERS $PARAMETERS"

	#PARAMETERS=" BEDFILE_GENES $BEDFILE_GENES"

	#exit 0;

	# OUTPUT
	echo "# "
	echo "# CONFIGURATION ($APP_NAME)"
	echo "################"
	echo "# NGS TOOLS Folder:      "$NGS_TOOLS
	#echo "# MISEQ Folder:          "$MISEQ_FOLDER
	#echo "# DEMULTIPLEXING Folder: "$DEMULTIPLEXING_FOLDER
	echo "# RESULTS Folder:        "$OUTPUT
	echo "# REPOSITORY Folder:     "$REPOSITORY
	#echo "# ANALYSIS Folder:       "$ANALYSIS_FOLDER
	echo "# ENV:                   "$ENV
	echo "# APP:                   "$APP_NAME
	echo "# APP RELEASE:           "$APP_RELEASE
	echo "# GROUP:                 "$GROUP
	echo "# PROJECT:               "$PROJECT
	echo "# ALIGNERS:              "$ALIGNERS
	echo "# CALLERS:               "$CALLERS
	echo "# ANNOTATORS:            "$ANNOTATORS
	echo "# PIPELINES:             "$PIPELINES
	echo "# "

	echo "[`date`] *** Start Analysis 	'$RELEASE_RUN'"
	echo "[`date`] MAKEFILE_ANALYSIS_RUN 	'$MAKEFILE_ANALYSIS_RUN'"
	echo "[`date`] SHELL_ANALYSIS_RUN 	'$SHELL_ANALYSIS_RUN'"
	echo "[`date`] LOGFILE_RES_RUN		'$LOGFILE_RES_RUN'"
	echo "[`date`] FASTQ/BAM/CRAM 		'$F_LIST'"
	echo "[`date`] FASTQ R2 		'$Q_LIST'"
	echo "[`date`] BED	 		'$B_LIST'"
	echo "[`date`] GENES	 		'$G_LIST'"
	echo "[`date`] TRANSCRIPTS	 	'$T_LIST'"
	echo "[`date`] SAMPLE	 		'$S_LIST'"
	echo "[`date`] RUN	 		'$RUU'"
	echo "[`date`] APP	 	        '$APP_NAME'"
	echo "[`date`] APP RELEASE	        '$APP_RELEASE'"
	echo "[`date`] GROUP	 	        '$SAMPLE_GROUP'"
	echo "[`date`] PROJECT	 	        '$SAMPLE_PROJECT'"
	echo "[`date`] PIPELINES	 	'$PIPELINES'"
	echo "[`date`] POST_ALIGNMENT	 	'"$(echo $POST_ALIGNMENT | tr "." "\n" | tac | tr "\n" " " )"'"


	#STARK_QUEUED=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "STARK_QUEUED")
	#STARK_RUNNING=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "STARK_RUNNING")
	#STARK_COMPLETE=$($SCRIPT_DIR/extract_variable_from_env.sh $RUN_ENV "STARK_COMPLETE")
	STARK_QUEUED=$(source $ENV; echo $STARK_QUEUED)
	STARK_RUNNING=$(source $ENV; echo $STARK_RUNNING)
	STARK_COMPLETE=$(source $ENV; echo $STARK_COMPLETE)
	if [ "$STARK_QUEUED" == "" ]; then STARK_QUEUED=STARKQueued.txt; fi;
	if [ "$STARK_RUNNING" == "" ]; then  STARK_RUNNING=STARKRunning.txt; fi;
	if [ "$STARK_COMPLETE" == "" ]; then  STARK_COMPLETE=STARKComplete.txt; fi;
	STARK_QUEUED_FILE=$OUTPUT/$RUU/$STARK_QUEUED
	STARK_RUNNING_FILE=$OUTPUT/$RUU/$STARK_RUNNING
	STARK_COMPLETE_FILE=$OUTPUT/$RUU/$STARK_COMPLETE

	# RUNNING
	echo "#["`date '+%Y%m%d-%H%M%S'`"] RUN $RUN running with STARK ($STARK_VERSION)" > $STARK_RUNNING_FILE

	# TREADS
	#JAVA_MEMORY=4
	THREADS_BY_SAMPLE=$THREADS; # Allocate all thread to the sample because no other sample analysed in parallele



	echo "["`date '+%Y%m%d-%H%M%S'`"] Main Analysis Process for Analysis '$RELEASE_RUN' START" >>$LOGFILE_RES_RUN
	make -k -j $THREADS -e ENV=$ENV PARAM=$MAKEFILE_ANALYSIS_RUN $PARAMETERS $THREAD_PARAMETERS JAVA_MEMORY=$JAVA_MEMORY SNAPSHOT=0 VALIDATION=1 INPUT=$INPUT OUTDIR=$OUTPUT RELEASE=$RELEASE_RUN FINAL_REPORT=$FINAL_REPORT_RUN ANALYSIS_REF=$ANALYSIS_REF -f $NGS_SCRIPTS/NGSWorkflow.mk 1>>$LOGFILE_RES_RUN 2>>$LOGFILE_RES_RUN
	echo "["`date '+%Y%m%d-%H%M%S'`"] Main Analysis Process for Analysis '$RELEASE_RUN' END" >>$LOGFILE_RES_RUN

	if (($(grep "\*\*\*" $LOGFILE_RES_RUN -c))); then
		echo "["`date '+%Y%m%d-%H%M%S'`"] Main Analysis Process for Analysis '$RELEASE_RUN' ERROR" >>$LOGFILE_RES_RUN
		echo "["`date`"] Main Analysis Process for Analysis '$RELEASE_RUN' finished with ERRORS"
		if (($VERBOSE)); then grep "\*\*\*" $LOGFILE_RES_RUN; fi;
		# RUNNING stop
		rm -f $STARK_RUNNING_FILE
		continue;
	fi;


	if ((0)); then
	# ANALYSIS

	echo "["`date '+%Y%m%d-%H%M%S'`"] Main Analysis Process for Analysis '$RELEASE_RUN' START" >>$LOGFILE_RES_RUN
	make -k -j $THREADS -e ENV=$ENV PARAM=$MAKEFILE_ANALYSIS_RUN $PARAMETERS $THREAD_PARAMETERS JAVA_MEMORY=$JAVA_MEMORY SNAPSHOT=0 VALIDATION=1 INPUT=$INPUT OUTDIR=$OUTPUT RELEASE=$RELEASE_RUN -f $NGS_SCRIPTS/NGSWorkflow.mk 1>>$LOGFILE_RES_RUN 2>>$LOGFILE_RES_RUN
	echo "["`date '+%Y%m%d-%H%M%S'`"] Main Analysis Process for Analysis '$RELEASE_RUN' END" >>$LOGFILE_RES_RUN

	if (($(grep "\*\*\*" $LOGFILE_RES_RUN -c))); then
		echo "["`date '+%Y%m%d-%H%M%S'`"] Main Analysis Process for Analysis '$RELEASE_RUN' ERROR" >>$LOGFILE_RES_RUN
		echo "["`date`"] Main Analysis Process for Analysis '$RELEASE_RUN' finished with ERRORS"
		if (($VERBOSE)); then grep "\*\*\*" $LOGFILE_RES_RUN; fi;
		# RUNNING stop
		rm -f $STARK_RUNNING_FILE
		continue;
	fi;


	# REPORT
	#echo "ANALYSIS_REF=$ANALYSIS_REF"
	echo "["`date '+%Y%m%d-%H%M%S'`"] Main Report Process for Analysis '$RELEASE_RUN' START" >>$LOGFILE_RES_RUN_REPORT
	#make -k -j $THREADS -e ENV=$ENV PARAM=$MAKEFILE_ANALYSIS_RUN $PARAMETERS $THREAD_PARAMETERS JAVA_MEMORY=$JAVA_MEMORY SNAPSHOT=0 VALIDATION=1 INPUT=$INPUT OUTDIR=$OUTPUT RELEASE=$RELEASE_RUN FINAL_REPORT=$FINAL_REPORT_RUN ANALYSIS_REF=$ANALYSIS_REF FINAL_REPORT_FULL=$FINAL_REPORT_FULL_RUN -f $NGS_SCRIPTS/NGSWorkflow.mk  1>>$LOGFILE_RES_RUN_REPORT 2>>$LOGFILE_RES_RUN_REPORT
	make -k -j $THREADS -e ENV=$ENV PARAM=$MAKEFILE_ANALYSIS_RUN $PARAMETERS $THREAD_PARAMETERS JAVA_MEMORY=$JAVA_MEMORY SNAPSHOT=0 VALIDATION=1 INPUT=$INPUT OUTDIR=$OUTPUT RELEASE=$RELEASE_RUN FINAL_REPORT=$FINAL_REPORT_RUN ANALYSIS_REF=$ANALYSIS_REF -f $NGS_SCRIPTS/NGSWorkflow.mk  1>>$LOGFILE_RES_RUN_REPORT 2>>$LOGFILE_RES_RUN_REPORT


	if (($(grep "\*\*\*" $LOGFILE_RES_RUN_REPORT -c))); then
		echo "["`date '+%Y%m%d-%H%M%S'`"] Main Report Process for Analysis '$RELEASE_RUN' ERROR" >>$LOGFILE_RES_RUN_REPORT
		echo "["`date`"] Main Report Process for Analysis '$RELEASE_RUN' finished with ERRORS"
		if (($VERBOSE)); then grep "\*\*\*" $LOGFILE_RES_RUN_REPORT; fi;
		# RUNNING stop
		rm -f $STARK_RUNNING_FILE
		continue;
	fi;
	fi;

	# STARK REport
	I=0
	RESULTS_FOLDER_COPY_FOLDER_RUU_STARKCopyComplete_list=""
	for S in $SAMPLE; do
		RU=${RUN_ARRAY[$I]};

		if [ "$RU" != "$RUU" ]; then
			((I++))
			continue;
		fi;

		# Reporting
		#echo "$NGS_SCRIPTS/stark_report.sh -r $OUTPUT -f $RUU -p $SAMPLE_PROJECT -g $SAMPLE_GROUP -u $SAMPLE_USER -s $S -e $ENV -d $ANALYSIS_REF"
		#echo "["`date`"] Reporting RUN/SAMPLE '$RUU/$S'..."
		#$NGS_SCRIPTS/stark_report.sh -r $OUTPUT -f $RUU -p $SAMPLE_PROJECT -g $SAMPLE_GROUP -u $SAMPLE_USER -s $S -e $ENV -i $(echo $PIPELINES | tr " " ",") -d $ANALYSIS_REF 1>>$LOGFILE_RES_RUN_REPORT 2>>$LOGFILE_RES_RUN_REPORT

		# REPOSITORY
		if ((1)); then
			# COPY of run/sample
			# List of folders
			if [ "$REPOSITORY" != "" ] ; then
				for RESULTS_FOLDER_COPY_FOLDER in $REPOSITORY; do
					#mkdir -p $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT;
					if [ ! -d $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT ]; then
						mkdir -p $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT;
					fi;
					if [ -d $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT ]; then
						RESULTS_FOLDER_COPY_ALL=$RESULTS_FOLDER_COPY_ALL" $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT";
					fi;

				done;
			fi;
			#RESULTS_SUBFOLDER_DATA="DATA"
			# Copy
			if [ "$RESULTS_FOLDER_COPY_ALL" != "$OUTPUT" ] && [ "$RESULTS_FOLDER_COPY_ALL" != "" ] ; then
				for RESULTS_FOLDER_COPY_FOLDER in $RESULTS_FOLDER_COPY_ALL;
				do
					echo "["`date`"] Copying RUN/SAMPLE '$RUU/$S' files from main result folder '$OUTPUT/$RUU/$S' to Copy folder '$RESULTS_FOLDER_COPY_FOLDER/$RUU/$S'..."

					# Copy SAMPLE files
					chmod $PERMS -R $OUTPUT/$RUU/$S 1>/dev/null 2>/dev/null
					mkdir -p $RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/$RESULTS_SUBFOLDER_DATA
					chmod $PERMS -R $RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/$RESULTS_SUBFOLDER_DATA 1>/dev/null 2>/dev/null
					#while ! $COMMAND_COPY $OUTPUT/$RUU/$S/* $RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/$RESULTS_SUBFOLDER_DATA 1>>$LOGFILE_RES_RUN_REPORT 2>>$LOGFILE_RES_RUN_REPORT ; do : ; done;
					$COMMAND_COPY $OUTPUT/$RUU/$S/* $RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/$RESULTS_SUBFOLDER_DATA 1>>$LOGFILE_RES_RUN_REPORT 2>>$LOGFILE_RES_RUN_REPORT
					# Copy RUN files
					#while ! nohup $COMMAND_COPY $(find $OUTPUT/$RUU -mindepth 1 -maxdepth 1 -type f) $RESULTS_FOLDER_COPY_FOLDER/$RUU 1>>$LOGFILE_RES_RUN_REPORT 2>>$LOGFILE_RES_RUN_REPORT ; do : ; done;
					$COMMAND_COPY $(find $OUTPUT/$RUU -mindepth 1 -maxdepth 1 -type f) $RESULTS_FOLDER_COPY_FOLDER/$RUU 1>>$LOGFILE_RES_RUN_REPORT 2>>$LOGFILE_RES_RUN_REPORT
					chmod $PERMS $RESULTS_FOLDER_COPY_FOLDER/$RUU/* 1>/dev/null 2>/dev/null
					# Copy ROOT FILE PATTERNS
					if [ $RESULTS_SUBFOLDER_DATA != "" ]; then
						for ROOT_FILE_PATTERN in $RESULTS_SUBFOLDER_ROOT_FILE_PATTERNS; do

							if [[ $ROOT_FILE_PATTERN =~ '$SAMPLE' ]]; then
								eval ROOT_FILE_PATTERN_VAR=$(echo $ROOT_FILE_PATTERN | sed s/\$SAMPLE/\$S/gi)
							else
								ROOT_FILE_PATTERN_VAR=$ROOT_FILE_PATTERN
							fi;

							for F_SOURCE in $(ls $RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/$RESULTS_SUBFOLDER_DATA/$ROOT_FILE_PATTERN_VAR); do
								L_SOURCE=$(echo $F_SOURCE | sed s#$RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/##g)
								F_TARGET="$RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/"$(basename $F_SOURCE)
								#D_TARGET="$RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/"
								#echo "F_SOURCE: $F_SOURCE";
								#echo "L_SOURCE: $L_SOURCE";
								#echo "F_TARGET: $F_TARGET";
								#echo "D_TARGET: $D_TARGET";
								if [ ! -e $F_TARGET ]; then
									$COMMAND_LINK $L_SOURCE $F_TARGET 1>>$LOGFILE_RES_RUN_REPORT 2>>$LOGFILE_RES_RUN_REPORT
								fi;
								if [ ! -e $F_TARGET ]; then
									rm -f $F_TARGET
									$COMMAND_COPY $F_SOURCE $F_TARGET 1>>$LOGFILE_RES_RUN_REPORT 2>>$LOGFILE_RES_RUN_REPORT
								fi;
							done;

						done;
					fi;

					# CREATE CopyComplete file
					echo "["`date '+%Y%m%d-%H%M%S'`"] Copy complete" >> $RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/STARKCopyComplete.txt
					chmod $PERMS $RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/STARKCopyComplete.txt 1>/dev/null 2>/dev/null
					RESULTS_FOLDER_COPY_FOLDER_RUU_STARKCopyComplete_list="$RESULTS_FOLDER_COPY_FOLDER_RUU_STARKCopyComplete_list $RESULTS_FOLDER_COPY_FOLDER/$RUU"

				done;
			fi;
		fi;

		((I++))
	done
	#echo "["`date '+%Y%m%d-%H%M%S'`"] Main Report Process for Analysis '$RELEASE_RUN' END" >>$LOGFILE_RES_RUN_REPORT


	echo "[`date`] *** Stop Analysis '$RELEASE_RUN'"

	# STARK COMPLETE
	echo "#["`date '+%Y%m%d-%H%M%S'`"] RUN $RUN analyzed by STARK ($STARK_VERSION)" >> $STARK_COMPLETE_FILE

	# CREATE CopyComplete file
	for RESULTS_FOLDER_COPY_FOLDER_RUU_STARKCopyComplete in $RESULTS_FOLDER_COPY_FOLDER_RUU_STARKCopyComplete_list; do
		echo "["`date '+%Y%m%d-%H%M%S'`"] Copy complete" >> $RESULTS_FOLDER_COPY_FOLDER_RUU_STARKCopyComplete/STARKCopyComplete.txt
		chmod $PERMS $RESULTS_FOLDER_COPY_FOLDER_RUU_STARKCopyComplete/STARKCopyComplete.txt 1>/dev/null 2>/dev/null
	done;

	# RUNNING stop
	rm -f $STARK_RUNNING_FILE

	#echo "# DONE @"`date '+%Y%m%d-%H%M%S'`

done;
