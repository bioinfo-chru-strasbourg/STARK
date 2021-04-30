#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="LaunchSample"
SCRIPT_DESCRIPTION="Launch a Sample Analysis"
SCRIPT_RELEASE="0.9.6.1b"
SCRIPT_DATE="19/03/2019"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-AGPL"

# Release note
#RELEASE_NOTES="#\n"
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.3b-09/10/2015: Add RUN configuration depending on Group and Project defined in RUN SampleSheet\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.4b-13/10/2015: Add Complete and Running flag files\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.5.1b-11/04/2016: Optimize Threads and modifiy options order\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.5.2b-13/04/2016: Adding reports process\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.6b-27/10/2016: Multi input FASTQ/BAM/CRAM/SAM and associated parameters\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.6.1b-19/03/2019: Minor change on FASTQ/BAM/CRAM/SAM input\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Configuration
ENV_CONFIG=$(find -L $SCRIPT_DIR/.. -name config.app)
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
	echo "#                                             Default: default.app if not defined";
	echo "# -f|--reads|--fastq_R1=<FILE1,FILE2,...>     List of FASTQ|BAM|SAM|CRAM (mandatory)";
	echo "#                                             Formats: *fastq.gz|*fq.gz|*bam|*ubam|*cram|*ucram|*sam|*usam";
	echo "# -q|--reads2|--fastq_R2=<FILE1,FILE2,...>    List of corresponding FASTQ Read2 file (beware of correspondance).";
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
	echo "# -o|--results=<FOLDER>                       RESULTS directory to generate RESULTS|RUN|SAMPLE|* files";
	echo "#                                             Default: Defined in APP, or first --reads file folder";
	echo "# -u|--repository=<FOLDER>                    Repository directory to generate GROUP|PROJECT|RUN|SAMPLE|* specific files";
	echo "#                                             Default: no copy in a repository";
	echo "# --archives=<FOLDER>                         Archives directory to generate GROUP|PROJECT|RUN|SAMPLE|* specific files";
	echo "#                                             Default: no copy in a archives";
	echo "# --databases=<FOLDER>                        Databases folder (requires STARK databases folder structure)";
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
# /home1/IRC/TOOLS/tools/stark/dev/STARK -f "$DATA/SAMPLE/HORIZON/HORIZON.R1.fastq.gz,$DATA/SAMPLE/HORIZON/HORIZON.bam,$DATA/SAMPLE/TEST/TEST.cram" -q "$DATA/SAMPLE/HORIZON/HORIZON.R2.fastq.gz" -s "H1,H2" -r "testH,testH,testT" -b "/home1/TOOLS/data/SAMPLE/HORIZON/HORIZON.manifest,/home1/TOOLS/data/SAMPLE/HORIZON/HORIZON.bed,/home1/TOOLS/data/SAMPLE/TEST/TEST.bed" -e env.HUSTUMSOL.sh --output="$DATA/RESULTS" --repository="$DATA/REPOSITORY"
# /home1/IRC/TOOLS/tools/stark/dev/STARK -f $DATA/SAMPLE/TEST/TEST.cram -b /home1/TOOLS/data/SAMPLE/TEST/TEST.bed --output="$DATA/RESULTS" -r testT

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "e:f:q:b:j:m:s:r:o:u:p:t:ga:vdnh" --long "env:,app:,application:,reads:,fastq:,reads1:,reads:,fastq_R1:,fastq_R2:,reads2:,index1:,index2:,analysis_tag:,sample_tag:,other_files:,design:,bed:,manifest:,genes:,transcripts:,sample:,sample_list:,sample_filter:,runs:,demultiplexing:,demultiplexing_only,samplesheet:,analysis:,analysis_name:,output:,analysis_dir:,results:,repository:,archives:,databases:,input:,output:,tmp:,log:,pipelines:,threads:,no_header,by_sample,adapters:,verbose,debug,release,help" -- "$@" 2> /dev/null)
[ $? -ne 0 ] && \
	echo "#[ERROR] Error in the argument list." "Use -h or --help to display the help." >&2 && usage && \
	exit 1
eval set -- "$ARGS"

#echo $PARAM
PARAM=$@

while true
do
	#echo "$1=$2"
	#echo "Eval opts";
	case "$1" in
		-e|--env|--app|--application)
			APP="$2"
			shift 2
			;;
		-f|--reads|--reads1|--fastq|--fastq_R1)
			FASTQ="$2"
			# transform RUNS list
			FASTQ=$(echo $FASTQ | tr "," " " | xargs)
			shift 2
			;;
		-q|--reads2|--fastq_R2)
			FASTQ_R2="$2"
			# transform RUNS list
			FASTQ_R2=$(echo $FASTQ_R2 | tr "," " " | xargs)
			shift 2
			;;
		--index1)
			INDEX1=$2 
			INDEX1=$(echo $INDEX1 | tr "," " " | xargs)
			shift 2
			;;
		--index2)
			INDEX2=$2 
			INDEX2=$(echo $INDEX2 | tr "," " " | xargs)
			shift 2
			;;
		--other_files)
			OTHER_FILES=$2 
			OTHER_FILES=$(echo $OTHER_FILES | tr "," " " | xargs)
			shift 2
			;;
		-b|--design|--bed|--manifest)
			BED_INPUT="$2"
			BED_INPUT=$(echo $BED_INPUT | tr "," " " | xargs)
			shift 2
			;;
		-j|--genes)
			BEDFILE_GENES_INPUT="$2"
			BEDFILE_GENES_INPUT=$(echo $BEDFILE_GENES_INPUT | tr "," " " | xargs)
			shift 2
			;;
		-m|--transcripts)
			TRANSCRIPTS_INPUT="$2"
			TRANSCRIPTS_INPUT=$(echo $TRANSCRIPTS_INPUT | tr "," " " | xargs)
			shift 2
			;;
		-s|--sample)
			SAMPLE="$2"
			SAMPLE=$(echo $SAMPLE | tr "," " " | xargs)
			shift 2
			;;
		--sample_tag)
			SAMPLE_TAG="$2"
			SAMPLE_TAG=$(echo $SAMPLE_TAG | tr "," " " | xargs)
			shift 2
			;;
		--sample_list|--sample_filter)
			SAMPLE_LIST="$2"
			SAMPLE_LIST=$(echo $SAMPLE_LIST | tr "," " " | xargs)
			shift 2
			;;
		--analysis_name)
			RUN="$2"
			RUN=$(echo $RUN | tr "," " " | xargs)
			shift 2
			;;
		--analysis_tag)
			ANALYSIS_TAG="$2"
			shift 2
			;;
		-l|--samplesheet)
			SAMPLESHEET_INPUT="$2" 
			shift 2
			;;
		--analysis_dir)
			ANALYSIS_DIR="$2"
			shift 2
			;;
		-u|--repository)
			REPOSITORY="$2"
			shift 2
			;;
		--archives)
			ARCHIVES="$2"
			shift 2
			;;
		--databases)
			DATABASES="$2"
			shift 2
			;;
		--input)
			INPUT="$2"
			shift 2
			;;
		--output)
			OUTPUT="$2"
			shift 2
			;;
		-o|--results)
			RESULTS="$2"
			shift 2
			;;
		--demultiplexing)
			DEMULTIPLEXING="$2" 
			shift 2
			;;
		--tmp)
			TMP="$2"
			shift 2
			;;
		--log)
			LOG="$2"
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
		--runs)
			#RUNS="$2"
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
		--no_header)
			NO_HEADER=1
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
		# UNUSED
		--demultiplexing)
			DEMULTIPLEXING="$2" 
			shift 2
			;;
		--demultiplexing_only)
			DEMULTIPLEXING_ONLY=1
			shift 1
			;;
		--) shift
			break
			;;
		*) 	echo "#[WARNING] Option $1 is not recognized. Use -h or --help to display the help." && usage && exit 1
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


#COMMAND_COPY="rsync -aucqpAXoghi --no-links --no-perms --no-owner --no-group" # "cp -auv" or "rsync -auv" # auvpAXog
COMMAND_COPY="rsync -auczqAXhi --no-links --no-perms --no-owner --no-group "
COMMAND_COPY_NO_COMPRESS="rsync -aucqAXhi --no-links --no-perms --no-owner --no-group "
COMMAND_LINK="ln " # "cp -auv" or "rsync -auv" # auvpAXog
PERMS="a+rwx"


# FUNCTIONS
#############

# function in_array
# input: $element $array
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

VERBOSE=$(echo $VERBOSE | awk '{print $0+0}')
DEBUG=$(echo $DEBUG | awk '{print $0+0}')



### APPLICATION

ENV=$(find_app "$APP" "$STARK_FOLDER_APPS")
source_app "$APP" "$STARK_FOLDER_APPS" 1

if (($DEBUG)); then
	echo "#[INFO] APP=$APP"
	echo "#[INFO] ASSEMBLY=$ASSEMBLY"
	echo "#[INFO] PIPELINES=$PIPELINES"
	echo "#[INFO] ALIGNERS=$ALIGNERS"
	echo "#[INFO] CALLERS=$CALLERS"
	echo "#[INFO] ANNOTATORS=$ANNOTATORS"
	echo "#[INFO] ASSEMBLY=$ASSEMBLY"
	#exit 0
fi;


(($DEBUG)) && echo "#[INFO] SAMPLE=$SAMPLE"



# FASTQ_R1
FASTQ_R1_RELOCATED=""
if [ "$FASTQ" != "" ]; then
	for F_INPUT in $FASTQ; do
		# Find file
		F=$(echo $(ls $F_INPUT 2>/dev/null | xargs realpath 2>/dev/null) $(cd $ANALYSIS_DIR && ls $F_INPUT 2>/dev/null | xargs realpath 2>/dev/null) | tr " " "\n" | sort -u)
		if [ "$F" == "" ] || [ $(echo $F | wc -w) -gt 1 ] || [ ! -f "$F" ]; then
			echo "[ERROR] No input FASTQ/BAM/CRAM/SAM '$F' file!";
			exit 0;
		fi
		(($DEBUG)) && echo "F=$F";
		# Relocate
		if ! (($(echo "$F" | grep ".fastq.gz$\|.fq.gz$\|.bam$\|.ubam$\|.cram$\|.ucram$" -c))); then
			echo "[ERROR]! Format of input file '$F' Unknown! Please check file format (.fastq.gz|.fq.gz|.bam|.ubam|.cram|.ucram)";
			exit 0;
		fi
		# Relocation
		FASTQ_R1_RELOCATED="$FASTQ_R1_RELOCATED $F"
	done;
	FASTQ=$FASTQ_R1_RELOCATED
fi;


# FASTQ_R2
FASTQ_R2_RELOCATED=""
if [ "$FASTQ_R2" != "" ]; then
	for F_INPUT in $FASTQ_R2; do
		# Find file
		F=$(echo $(ls $F_INPUT 2>/dev/null | xargs realpath 2>/dev/null) $(cd $ANALYSIS_DIR && ls $F_INPUT 2>/dev/null | xargs realpath 2>/dev/null) | tr " " "\n" | sort -u)
		if [ "$F" == "" ] || [ $(echo $F | wc -w) -gt 1 ] || [ ! -f "$F" ]; then
			echo "[ERROR] No input FASTQ '$F' file!";
			exit 0;
		fi
		(($DEBUG)) && echo "Q=$F";
		# Relocated
		if ! (($(echo "$F" | grep ".fastq.gz$\|.fq.gz$" -c))); then
			echo "[ERROR]! Format of input file '$F' Unknown! Please check file format (.fastq.gz|.fq.gz)";
			exit 0;
		fi
		# Relocation
		FASTQ_R2_RELOCATED="$FASTQ_R2_RELOCATED $F"
	done;
	FASTQ_R2=$FASTQ_R2_RELOCATED
fi;


# INDEX1
INDEX1_RELOCATED=""
if [ "$INDEX1" != "" ]; then
	for F_INPUT in $INDEX1; do
		# Find file
		F=$(echo $(ls $F_INPUT 2>/dev/null | xargs realpath 2>/dev/null) $(cd $ANALYSIS_DIR && ls $F_INPUT 2>/dev/null | xargs realpath 2>/dev/null) | tr " " "\n" | sort -u)
		if [ "$F" == "" ] || [ $(echo $F | wc -w) -gt 1 ] || [ ! -f "$F" ]; then
			echo "[ERROR] No input FASTQ '$F' file!";
			exit 0;
		fi
		(($DEBUG)) && echo "I1=$F";
		if ! (($(echo "$F" | grep ".fastq.gz$\|.fq.gz$" -c))); then
			echo "[ERROR]! Format of input file '$F' Unknown! Please check file format (.fastq.gz|.fq.gz)";
			exit 0;
		fi
		# Relocation
		INDEX1_RELOCATED="$INDEX1_RELOCATED $F"
	done;
	INDEX1=$INDEX1_RELOCATED
fi;


# INDEX2
INDEX2_RELOCATED=""
if [ "$INDEX2" != "" ]; then
	for F_INPUT in $INDEX2; do
		# Find file
		F=$(echo $(ls $F_INPUT 2>/dev/null | xargs realpath 2>/dev/null) $(cd $ANALYSIS_DIR && ls $F_INPUT 2>/dev/null | xargs realpath 2>/dev/null) | tr " " "\n" | sort -u)
		if [ "$F" == "" ] || [ $(echo $F | wc -w) -gt 1 ] || [ ! -f "$F" ]; then
			echo "[ERROR] No input FASTQ '$F' file!";
			exit 0;
		fi
		(($DEBUG)) && echo "I2=$F";
		if ! (($(echo "$F" | grep ".fastq.gz$\|.fq.gz$" -c))); then
			echo "[ERROR]! Format of input file '$F' Unknown! Please check file format (.fastq.gz|.fq.gz)";
			exit 0;
		fi
		# Relocation
		INDEX2_RELOCATED="$INDEX2_RELOCATED $F"
	done;
	INDEX2=$INDEX2_RELOCATED
fi;


# OTHER_FILES
OTHER_FILES_RELOCATED=""
if [ "$OTHER_FILES" != "" ]; then
	#for OFCF in $(echo $(ls $OTHER_FILES 2>/dev/null | xargs realpath 2>/dev/null) $(cd $ANALYSIS_DIR && ls $OTHER_FILES 2>/dev/null | xargs realpath 2>/dev/null) | tr " " "\n" | sort -u); do
	for OFCF in $OTHER_FILES; do
		(($DEBUG)) && echo "OFCF=$OFCF"
		OTHER_FILES_RELOCATED_ONE=""
		for OFCF_ONE in $(echo $OFCF | tr "+" " "); do
			# Original other file
			OF_ONE=$(echo $OFCF_ONE | awk -F: '{print $1}')
			# Copy other file
			CF_ONE=$(echo $OFCF_ONE | awk -F: '{print $2}')
			(($DEBUG)) && echo "OFCF_ONE=$OFCF_ONE"
			if [ "$OF_ONE" == "" ] || [ ! -e "$OF_ONE" ]; then
				echo "[ERROR] No input '$OF_ONE' file/folder!";
				exit 0;
			fi
			# Relocation
			OTHER_FILES_RELOCATED_ONE="$OTHER_FILES_RELOCATED_ONE+$OF_ONE"
			[ "$CF_ONE" != "" ] && OTHER_FILES_RELOCATED_ONE="$OTHER_FILES_RELOCATED_ONE:$CF_ONE"
		done;
		OTHER_FILES_RELOCATED_ONE=$(echo $OTHER_FILES_RELOCATED_ONE | sed "s/^+//" )
		OTHER_FILES_RELOCATED="$OTHER_FILES_RELOCATED $OTHER_FILES_RELOCATED_ONE"
	done;
	OTHER_FILES_RELOCATED=$(echo $OTHER_FILES_RELOCATED | sed "s/^,//" )
	OTHER_FILES=$OTHER_FILES_RELOCATED
fi;


# DETECT_ADAPTER_FOR_PE

# BED
BED_DEFAULT="";
BED_RELOCATED=""
if [ "$BED_INPUT" != "" ]; then
	# Test BED files exist
	for B in $BED_INPUT; do
		if [ -f "$B" ] || [ -f "$ANALYSIS_DIR/$B" ]; then
			(($DEBUG)) && echo "B=$B";
			# Relocation
			if [ -f "$B" ]; then
				BED_RELOCATED="$BED_RELOCATED $B"
		 	elif [ -f "$ANALYSIS_DIR/$B" ]; then
				BED_RELOCATED="$BED_RELOCATED $ANALYSIS_DIR/$B"
			fi;
		else
			echo "[ERROR] BED '$B' file DOES NOT exist!";
			exit 0;
		fi;
	done;


fi;
BED_INPUT=$BED_RELOCATED
BED_DEFAULT=$(echo $BED_INPUT | cut -d" " -f1) #`date '+%Y%m%d-%H%M%S'`


# BEDFILE_GENES
if [ "$BEDFILE_GENES_INPUT" != "" ]; then 
	BEDFILE_GENES=$BEDFILE_GENES_INPUT;
fi;

BEDFILE_GENES_DEFAULT="";
BEDFILE_GENES_RELOCATED=""
if [ "$BEDFILE_GENES" != "" ]; then
	# Test BED files exist
	for G_LIST_ONE in $(echo $BEDFILE_GENES); do
		BEDFILE_GENES_LIST=""
		for G in $(echo $G_LIST_ONE | tr "+" " "); do
			if [ -f "$G" ] || [ -f "$ANALYSIS_DIR/$G" ]; then
				(($DEBUG)) && echo "G=$G";
				# Relocation
				if [ -f "$G" ]; then
					BEDFILE_GENES_LIST="$BEDFILE_GENES_LIST+$G"
			 	elif [ -f "$ANALYSIS_DIR/$G" ]; then
					BEDFILE_GENES_LIST="$BEDFILE_GENES_LIST+$ANALYSIS_DIR/$G"
				fi;
			else
				echo "[ERROR] BEDFILE_GENES '$G' file DOES NOT exist!";
				exit 0;
			fi;
		done;
		BEDFILE_GENES_RELOCATED="$BEDFILE_GENES_RELOCATED "$(echo $BEDFILE_GENES_LIST | sed "s/+$//" | sed "s/^+//")
	done;
fi;
BEDFILE_GENES=$BEDFILE_GENES_RELOCATED
BEDFILE_GENES_DEFAULT=$(echo $BEDFILE_GENES | cut -d" " -f1) #`date '+%Y%m%d-%H%M%S'`


# TRANSCRIPTS
if [ "$TRANSCRIPTS_INPUT" != "" ]; then
	TRANSCRIPTS=$TRANSCRIPTS_INPUT;
fi;

TRANSCRIPTS_DEFAULT="";
TRANSCRIPTS_RELOCATED=""
if [ "$TRANSCRIPTS" != "" ]; then
	# Test BED files exist
	for T_LIST_ONE in $TRANSCRIPTS; do
		TRANSCRIPTS_LIST=""
		for T in $(echo $T_LIST_ONE | tr "+" " "); do
			if [ -f "$T" ] || [ -f "$ANALYSIS_DIR/$T" ]; then
				(($DEBUG)) && echo "T=$T";
				# Relocation
				if [ -f "$T" ]; then
					TRANSCRIPTS_LIST="$TRANSCRIPTS_LIST+$T"
				elif [ -f "$ANALYSIS_DIR/$G" ]; then
					TRANSCRIPTS_LIST="$TRANSCRIPTS_LIST+$ANALYSIS_DIR/$T"
				fi;
			else
				echo "[ERROR] TRANSCRIPTS '$T' file DOES NOT exist!";
				exit 0;
			fi;
		done;
		TRANSCRIPTS_RELOCATED="$TRANSCRIPTS_RELOCATED "$(echo $TRANSCRIPTS_LIST | sed "s/+$//" | sed "s/^+//")
	done;

fi;
TRANSCRIPTS=$TRANSCRIPTS_RELOCATED
TRANSCRIPTS_DEFAULT=$(echo $TRANSCRIPTS | cut -d" " -f1) #`date '+%Y%m%d-%H%M%S'`


# RUN
RUN_DEFAULT=`date '+%Y%m%d'`
if [ "$RUN" != "" ]; then
	RUN_DEFAULT=$(echo $RUN | cut -d" " -f1) #`date '+%Y%m%d-%H%M%S'`
fi;


NB_SAMPLE=0
FASTQ_R2_ARRAY=($FASTQ_R2);
FASTQ_R2_CHECKED=""
INDEX1_ARRAY=($INDEX1);
INDEX1_CHECKED=""
INDEX2_ARRAY=($INDEX2);
INDEX2_CHECKED=""
OTHER_FILES_ARRAY=($OTHER_FILES);
OTHER_FILES_CHECKED=""
SAMPLE_CHECKED=""
SAMPLE_ARRAY=($SAMPLE);
RUN_CHECKED=""
RUN_ARRAY=($RUN);
BED_CHECKED=""
BED_ARRAY=($BED_INPUT);
BEDFILE_GENES_CHECKED=""
BEDFILE_GENES_ARRAY=($BEDFILE_GENES);
TRANSCRIPTS_CHECKED=""
TRANSCRIPTS_ARRAY=($TRANSCRIPTS);
SAMPLE_TAG_CHECKED=""
SAMPLE_TAG_ARRAY=($SAMPLE_TAG);


(($DEBUG)) && echo "#[INFO] FASTQ=$FASTQ";

for F in $FASTQ; do

	# Test FASTQ format
	if ! (($(echo "$F" | grep ".fastq.gz$\|.fq.gz$\|.bam$\|.ubam$\|.cram$\|.ucram$\|.sam$\|.usam$" -c))); then
		echo "[ERROR] Format of input file '$F' Unknown! Please check file format (.fastq.gz|.fq.gz|.bam|.ubam|.cram|.ucram|.sam|.usam)";
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
	if [ "${BED_ARRAY[$NB_SAMPLE]}" == "" ] || [ ! -f "${BED_ARRAY[$NB_SAMPLE]}" ]; then
		BED_CHECKED="$BED_CHECKED$BED_DEFAULT "
	else
		BED_CHECKED="$BED_CHECKED${BED_ARRAY[$NB_SAMPLE]} "
	fi;

	# BEDFILE_GENES
	if [ "${BEDFILE_GENES_ARRAY[$NB_SAMPLE]}" == "" ]; then
		BEDFILE_GENES_CHECKED="$BEDFILE_GENES_CHECKED$BEDFILE_GENES_DEFAULT "
	else
		BEDFILE_GENES_CHECKED="$BEDFILE_GENES_CHECKED${BEDFILE_GENES_ARRAY[$NB_SAMPLE]} "
	fi;

	# TRANSCRIPTS
	if [ "${TRANSCRIPTS_ARRAY[$NB_SAMPLE]}" == "" ] || [ ! -f "${TRANSCRIPTS_ARRAY[$NB_SAMPLE]}" ]; then
		TRANSCRIPTS_CHECKED="$TRANSCRIPTS_CHECKED$TRANSCRIPTS_DEFAULT "
	else
		TRANSCRIPTS_CHECKED="$TRANSCRIPTS_CHECKED${TRANSCRIPTS_ARRAY[$NB_SAMPLE]} "
	fi;

	# SAMPLE TAG
	if [ "${SAMPLE_TAG_ARRAY[$NB_SAMPLE]}" == "" ] || [ ! -f "${SAMPLE_TAG_ARRAY[$NB_SAMPLE]}" ]; then
		SAMPLE_TAG_CHECKED="$SAMPLE_TAG_CHECKED"
	else
		SAMPLE_TAG_CHECKED="$SAMPLE_TAG_CHECKED${SAMPLE_TAG_ARRAY[$NB_SAMPLE]} "
	fi;


	((NB_SAMPLE++))
done;

SAMPLE=$SAMPLE_CHECKED;
RUN=$RUN_CHECKED;
BED=$BED_CHECKED;
BEDFILE_GENES=$BEDFILE_GENES_CHECKED;
TRANSCRIPTS=$TRANSCRIPTS_CHECKED;


# ARRAYS
FASTQ_ARRAY=($FASTQ);
FASTQ_R2_ARRAY=($FASTQ_R2);
INDEX1_ARRAY=($INDEX1);
INDEX2_ARRAY=($INDEX2);
OTHER_FILES_ARRAY=($OTHER_FILES);
SAMPLE_ARRAY=($SAMPLE);
RUN_ARRAY=($RUN);
BED_ARRAY=($BED);
BEDFILE_GENES_ARRAY=($BEDFILE_GENES);
TRANSCRIPTS_ARRAY=($TRANSCRIPTS);



# MULTIPLE BAM/CRAM/SAM/FASTQ
if (($BY_SAMPLE)); then
	if [ $(echo $FASTQ | wc -w ) -gt 1 ]; then
		echo "# Multiple SAMPLES: $SAMPLE";
		PARAM_MULTI=$(echo $PARAM | sed "s/--by_sample[^ |$]*//gi"  | sed "s/-g[^ |$]*//gi" | sed s/--$//gi)
		i=0
		for f in ${FASTQ_ARRAY[@]}; do
			q_var=""
			if [ ! -z ${FASTQ_R2_ARRAY[$i]} ]; then q_var=" -q ${FASTQ_R2_ARRAY[$i]}"; fi;
			$SCRIPT_DIR/$(basename $0) $PARAM_MULTI --reads=$f --reads2="${FASTQ_R2_ARRAY[$i]}" --index1="${INDEX1_ARRAY[$i]}" --index2="${INDEX2_ARRAY[$i]}" --other_files="${OTHER_FILES_ARRAY[$i]}" --sample="${SAMPLE_ARRAY[$i]}" --run="${RUN_ARRAY[$i]}" --bed="${BED_ARRAY[$i]}" --genes="${BEDFILE_GENES_ARRAY[$i]}"  --transcripts="${TRANSCRIPTS_ARRAY[$i]}" ;
			((i++))
		done;
		exit 0;
	fi;
fi;


# PIPELIENS PARAMETERS
if [ "$PIPELINES_INPUT" != "" ]; then
	PIPELINES=$PIPELINES_INPUT
fi;
if [ "$PIPELINES" == "" ]; then
	PIPELINES="bwamem.gatkHC.howard"
fi;


#DATE
DATE_DAY=`date '+%Y%m%d'`
DATE_MIN=`date '+%Y%m%d-%H%M%S'`
ANALYSIS_REF=$DATE_MIN


re='^[0-9]+$'
CORES=$(nproc)
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


# INPUT
# Input is the folder of the first FASTQ by default
for F in $(echo $FASTQ | cut -d " " -f1); do
	INPUT=`dirname "$F"`
done


# DEFAULT RESULTS
# RESULTS is 1/ the RESULTS folder in option, OR 2/ FOLDER_RESULT in ENV, OR 3/ folder of the first FASTQ by default
if [ -z "$RESULTS" ]; then
	#RESULTS=$FOLDER_RESULTS
	RESULTS=$RESULTS_FOLDER
fi;
if [ -z "$RESULTS" ]; then
	RESULTS=$INPUT
fi;

# DEFAULT RESULTS
# RESULTS is 1/ the RESULTS folder in option, OR 2/ FOLDER_RESULT in ENV, OR 3/ folder of the first FASTQ by default
if [ -z "$DEMULTIPLEXING" ]; then
	#RESULTS=$FOLDER_RESULTS
	DEMULTIPLEXING=$FOLDER_DEMULTIPLEXING
fi;
# if [ -z "$RESULTS" ]; then
# 	RESULTS=$INPUT
# fi;



# CREATE RESULTS if necessary
if [ ! -d "$RESULTS" ]; then
	mkdir -p $RESULTS;
fi

# JAVA OPTIONS
if [ "$JAVA_FLAGS" == "" ]; then
	JAVA_FLAGS=" -Xmx"$JAVA_MEMORY"g  -Dsamjdk.try_use_intel_deflater=false -Dsnappy.disable=true -Dorg.xerial.snappy.tempdir="$TMP_FOLDER_TMP" -Djava.io.tmpdir="$TMP_FOLDER_TMP;
fi;


APP_NAME=$(name_app "$APP" "$STARK_FOLDER_APPS");
SAMPLE_GROUP=$APP_GROUP
SAMPLE_PROJECT=$APP_PROJECT
if [ -z "$APP_NAME" ]; then APP_NAME="UNKNOWN"; fi;
if [ -z "$SAMPLE_GROUP" ]; then SAMPLE_GROUP="UNKNOWN"; fi;
if [ -z "$SAMPLE_PROJECT" ]; then SAMPLE_PROJECT="UNKNOWN"; fi;

if [ -z "$APP_NAME" ]; then APP_NAME="UNKNOWN"; fi;


if (($DEBUG)); then
	echo -e "FASTQ $FASTQ\nFASTQ_R2 $FASTQ_R2\nSAMPLE $SAMPLE\nRUN $RUN\nBED $BED\nGENES $BEDFILE_GENES\nTRANSCRIPTS $TRANSCRIPTS\nINPUT $INPUT\nOUTPUT $OUTPUT\nRESULTS $RESULTS\n\nPIPELINES $PIPELINES" | column -t;
fi;

SAMPLE_ARRAY=($SAMPLE);
FASTQ_ARRAY=($FASTQ);
FASTQ_R2_ARRAY=($FASTQ_R2);
INDEX1_ARRAY=($INDEX1);
INDEX2_ARRAY=($INDEX2);
OTHER_FILES_ARRAY=($OTHER_FILES);
RUN_ARRAY=($RUN);
BED_ARRAY=($BED);
BEDFILE_GENES_ARRAY=($BEDFILE_GENES);
TRANSCRIPTS_ARRAY=($TRANSCRIPTS);

# FOREACH RUN
RUN_UNIQ=$(echo $RUN | tr " " "\n" | sort | uniq | tr "\n" " ");



(($VERBOSE)) && echo "#[INFO] *** Start Analysis        ["$(date)"]"
START_ANALYSIS=$(date +%s)

echo "#[INFO] *** Input"

echo "#[INFO] Check Input Files..."


RUN_ANALYZED=0

for RUU in $RUN_UNIQ; do

	((RUN_ANALYZED++))

	ANALYSIS_PREFIX="STARK.$ANALYSIS_REF"

	# Analysis log files
	ANALYSIS_PREFIX_PATH=$RESULTS/$RUU/$ANALYSIS_PREFIX
	LOGFILE_DEM_RUN=$ANALYSIS_PREFIX_PATH.demultiplexing.log
	MAKEFILE_ANALYSIS_RUN=$ANALYSIS_PREFIX_PATH.analysis.param.mk
	LOGFILE_RES_RUN=$ANALYSIS_PREFIX_PATH.analysis.log
	CONFIG_RES_RUN=$ANALYSIS_PREFIX_PATH.config
	RELEASE_RUN=$ANALYSIS_PREFIX_PATH.analysis.release
	FASTQ_MK=$ANALYSIS_PREFIX_PATH.analysis.fastq.mk
	FASTP_MK=$ANALYSIS_PREFIX_PATH.analysis.fastp.mk
	RES_MK=$ANALYSIS_PREFIX_PATH.analysis.copy.mk

	# Analysis report file
	FINAL_REPORT_RUN=$RESULTS/$RUU/$ANALYSIS_PREFIX

	# MKDIR & TOUCH
	mkdir -p $RESULTS/$RUU
	touch $MAKEFILE_ANALYSIS_RUN
	#touch $SHELL_ANALYSIS_RUN

	### Results log
	touch $LOGFILE_RES_RUN

	### Demultiplexing log
	touch $LOGFILE_DEM_RUN

	### Find demultiplexing log
	if [ -e $DEMULTIPLEXING/$RUU/demultiplexing.configuration.log ]; then
		cp -pf $DEMULTIPLEXING/$RUU/demultiplexing.configuration.log $LOGFILE_DEM_RUN
	fi;
	

	# FASTQ
	> $FASTQ_MK
	FASTQ_MK_ALL=""

	# FASTP
	> $FASTP_MK
	FASTP_MK_ALL=""



	(($VERBOSE)) && echo "#[INFO] RUN '$RUU'"
	echo "## ANALYSIS" > $MAKEFILE_ANALYSIS_RUN
	#echo "## ANALYSIS" > $SHELL_ANALYSIS_RUN
	RUNS_SAMPLES="";

	F_LIST=""
	Q_LIST=""
	I1_LIST=""
	I2_LIST=""
	OF_LIST=""
	S_LIST=""
	B_LIST=""
	G_LIST=""
	T_LIST=""
	TAG_LIST=""
	#ATAG_LIST=""


	I=0;
	NB_SAMPLE=0;


	for F in $FASTQ; do


		S=${SAMPLE_ARRAY[$I]};
		RU=${RUN_ARRAY[$I]};
		F_R2=${FASTQ_R2_ARRAY[$I]};
		I1=${INDEX1_ARRAY[$I]};
		I2=${INDEX2_ARRAY[$I]};
		OFCF=${OTHER_FILES_ARRAY[$I]};
		B=${BED_ARRAY[$I]};
		G=${BEDFILE_GENES_ARRAY[$I]};
		T=${TRANSCRIPTS_ARRAY[$I]};
		TAG=${SAMPLE_TAG_ARRAY[$I]};

		# BY RUN
		if [ "$RU" != "$RUU" ]; then
			((I++))
			continue;
		fi;

		RUN_SAMPLE_DIR=$RESULTS/$RU/$S

		if (($DEBUG)); then
			echo "$F | $F_R2 | $I1 | $I2 | $OFCF | $S | $RU | $B | $G | $T | $PIPELINES | $INPUT | $RESULTS | $RUN_SAMPLE_DIR ";
		fi;

		mkdir -p $RUN_SAMPLE_DIR

		# INFOS
		(($VERBOSE)) && echo "#[INFO] RUN '$RUU' - SAMPLE '$S'"
		(($VERBOSE)) && echo "#[INFO] SAMPLE '$RU/$S' from file(s):"
		(($VERBOSE)) && echo "#[INFO] $F $F_R2 $I1 $I2 $OFCF $B $G $T "

		# Copy FASTQ
		PICARD_FLAGS="COMPRESSION_LEVEL=1 MAX_RECORDS_IN_RAM=500000"
		PICARD_UNALIGNED_FLAGS="SORT_ORDER=coordinate "
		PICARD_UNALIGNED_NAME_FLAGS="LIBRARY_NAME=001 PLATFORM=ILLUMINA PLATFORM_UNIT=PU READ_GROUP_NAME=A SAMPLE_NAME=$S "
		PICARD_VALIDATE_FLAGS=" MODE=SUMMARY "
		PICARD_VALIDATE_NAME_FLAGS="  "
		PICARD_ADDORREPLACEREADGROUPS_FLAGS="  "
		PICARD_ADDORREPLACEREADGROUPS_NAME_FLAGS=" RGLB=001 RGPL=ILLUMINA RGPU=PU RGSM=$S VALIDATION_STRINGENCY=SILENT "

		# DATA entry format is FASTQ/BAM/CRAM/SAM
		if [ ! -s $RUN_SAMPLE_DIR/$S.R1.fastq.gz ]; then
			# FASTQ
			if (($(echo "$F" | grep ".fastq.gz$\|.fq.gz$" -c))); then
				(($VERBOSE)) && echo "#[INFO] Create Input data from FASTQ file(s)"

				# Create FASTQ
				#$COMMAND_COPY $F $RUN_SAMPLE_DIR/$S.R1.fastq.gz;
				#$COMMAND_COPY_NO_COMPRESS $F $RUN_SAMPLE_DIR/$S.R1.fastq.gz;
				# FASTQ MK
				echo "$RUN_SAMPLE_DIR/$S.R1.fastq.gz: $F
					$COMMAND_COPY_NO_COMPRESS $F $RUN_SAMPLE_DIR/$S.R1.fastq.gz;
				" >> $FASTQ_MK
				FASTQ_MK_ALL=$FASTQ_MK_ALL" $RUN_SAMPLE_DIR/$S.R1.fastq.gz"

				if [ -s "$F_R2" ]; then
					#$COMMAND_COPY $F_R2 $RUN_SAMPLE_DIR/$S.R2.fastq.gz
					#$COMMAND_COPY_NO_COMPRESS $F_R2 $RUN_SAMPLE_DIR/$S.R2.fastq.gz
					# FASTQ MK
					echo "$RUN_SAMPLE_DIR/$S.R2.fastq.gz: $F_R2
						$COMMAND_COPY_NO_COMPRESS $F_R2 $RUN_SAMPLE_DIR/$S.R2.fastq.gz
					" >> $FASTQ_MK
					FASTQ_MK_ALL=$FASTQ_MK_ALL" $RUN_SAMPLE_DIR/$S.R2.fastq.gz"
				else
					#touch $RUN_SAMPLE_DIR/$S.R2.fastq
					#$GZ $RUN_SAMPLE_DIR/$S.R2.fastq
					echo "$RUN_SAMPLE_DIR/$S.R2.fastq.gz: $F
						touch $RUN_SAMPLE_DIR/$S.R2.fastq
						$GZ $RUN_SAMPLE_DIR/$S.R2.fastq
					" >> $FASTQ_MK
					FASTQ_MK_ALL=$FASTQ_MK_ALL" $RUN_SAMPLE_DIR/$S.R2.fastq.gz"
				fi;


				# INDEX1
				if [ -s "$I1" ]; then
					#$COMMAND_COPY_NO_COMPRESS $I1 $RUN_SAMPLE_DIR/$S.I1.fastq.gz
					# FASTQ MK
					echo "$RUN_SAMPLE_DIR/$S.I1.fastq.gz: $I1
						$COMMAND_COPY_NO_COMPRESS $I1 $RUN_SAMPLE_DIR/$S.I1.fastq.gz
					" >> $FASTQ_MK
					FASTQ_MK_ALL=$FASTQ_MK_ALL" $RUN_SAMPLE_DIR/$S.I1.fastq.gz"
				else
					#touch $RUN_SAMPLE_DIR/$S.I1.fastq
					#$GZ $RUN_SAMPLE_DIR/$S.I1.fastq
					echo "$RUN_SAMPLE_DIR/$S.I1.fastq.gz: $I1
						touch $RUN_SAMPLE_DIR/$S.I1.fastq
						$GZ $RUN_SAMPLE_DIR/$S.I1.fastq
					" >> $FASTQ_MK
					FASTQ_MK_ALL=$FASTQ_MK_ALL" $RUN_SAMPLE_DIR/$S.I1.fastq.gz"
				fi;

				# INDEX2
				if [ -s "$I2" ]; then
					#$COMMAND_COPY_NO_COMPRESS $I2 $RUN_SAMPLE_DIR/$S.I2.fastq.gz
					# FASTQ MK
					echo "$RUN_SAMPLE_DIR/$S.I2.fastq.gz: $I2
						$COMMAND_COPY_NO_COMPRESS $I2 $RUN_SAMPLE_DIR/$S.I2.fastq.gz
					" >> $FASTQ_MK
					FASTQ_MK_ALL=$FASTQ_MK_ALL" $RUN_SAMPLE_DIR/$S.I2.fastq.gz"
				else
					#touch $RUN_SAMPLE_DIR/$S.I2.fastq
					#$GZ $RUN_SAMPLE_DIR/$S.I2.fastq
					echo "$RUN_SAMPLE_DIR/$S.I2.fastq.gz: $I2
						touch $RUN_SAMPLE_DIR/$S.I2.fastq
						$GZ $RUN_SAMPLE_DIR/$S.I2.fastq
					" >> $FASTQ_MK
					FASTQ_MK_ALL=$FASTQ_MK_ALL" $RUN_SAMPLE_DIR/$S.I2.fastq.gz"
				fi;


			# BAM/CRAM/SAM
			elif (($(echo $F | grep ".bam$\|.ubam$\|.cram$\|.ucram\|.sam$\|.usam$" -c))); then

				# Generate FASTQ from BAM/SAM/CRAM
				if ((0)); then
					(($VERBOSE)) && echo "#[INFO] Create Input data from BAM/CRAM/SAM file"
					# Sort and FASTQ generation
					TMP_INPUT_BAM=$TMP_FOLDER_TMP/INPUT_BAM_$RANDOM
					$SAMTOOLS sort -n --reference $REF -@ $THREADS $F -T $TMP_INPUT_BAM -O SAM 2>/dev/null | $SAMTOOLS bam2fq - -1 $RUN_SAMPLE_DIR/$S.R1.fastq.gz -2 $RUN_SAMPLE_DIR/$S.R2.fastq.gz -O -@ $THREADS -c 1 2>/dev/null | $GZ -c 1> $RUN_SAMPLE_DIR/$S.R0.fastq.gz 2>/dev/null;

					# Ff R0 not empty or number of reads differs between R1 and R2 > Single END
					if (( $($UNGZ -c $RUN_SAMPLE_DIR/$S.R0.fastq.gz | head -n 1 | wc -l) )) \
						|| (( $(diff <($UNGZ -c $RUN_SAMPLE_DIR/$S.R1.fastq.gz | paste - - - - | cut -f1 -d$'\t') <($UNGZ -c $RUN_SAMPLE_DIR/$S.R2.fastq.gz | paste - - - - | cut -f1 -d$'\t') | head -n 1 | wc -l) )); then
						echo "#[WARNING] Read names or order differ between R1 and R2 fastq files, or R0 file is not empty. Reads are supposed not to be paired (Single-End). All reads in R1.";
						#$UNGZ -c $RUN_SAMPLE_DIR/$S.R2.fastq.gz $RUN_SAMPLE_DIR/$S.R0.fastq.gz >> $RUN_SAMPLE_DIR/$S.R1.fastq.gz;
						cat $RUN_SAMPLE_DIR/$S.R2.fastq.gz $RUN_SAMPLE_DIR/$S.R0.fastq.gz >> $RUN_SAMPLE_DIR/$S.R1.fastq.gz;
						$GZ < /dev/null > $RUN_SAMPLE_DIR/$S.R2.fastq.gz;
					fi;
					rm -f $RUN_SAMPLE_DIR/$S.R0.fastq.gz;
				fi;

				if ((1)); then
					(($VERBOSE)) && echo "#[INFO] Create Input data from BAM/CRAM/SAM file"

					# Sort and FASTQ generation
					TMP_INPUT_BAM=$TMP_FOLDER_TMP/INPUT_BAM_$RANDOM

					echo "$RUN_SAMPLE_DIR/$S.R1.fastq.gz: $F
						$SAMTOOLS sort -n --reference $REF -@ $THREADS $F -T $TMP_INPUT_BAM -O SAM 2>/dev/null | $SAMTOOLS bam2fq - -1 $RUN_SAMPLE_DIR/$S.R1.fastq.gz -2 $RUN_SAMPLE_DIR/$S.R2.fastq.gz -O -@ $THREADS -c 1 2>/dev/null | $GZ -c 1> $RUN_SAMPLE_DIR/$S.R0.fastq.gz 2>/dev/null;
						if (( \$\$($UNGZ -c $RUN_SAMPLE_DIR/$S.R0.fastq.gz | head -n 1 | wc -l) )) \
							|| [ \"\$\$($UNGZ -c $RUN_SAMPLE_DIR/$S.R1.fastq.gz | paste - - - - | cut -f1 | sha1sum | cut -d' ' -f1)\" != \"\$\$($UNGZ -c $RUN_SAMPLE_DIR/$S.R2.fastq.gz | paste - - - - | cut -f1 | sha1sum | cut -d' ' -f1)\" ]; then \
							echo '#[WARNING] Read names or order differ between R1 and R2 fastq files, or R0 file is not empty. Reads are supposed not to be paired (Single-End). All reads in R1.'; \
							cat $RUN_SAMPLE_DIR/$S.R2.fastq.gz $RUN_SAMPLE_DIR/$S.R0.fastq.gz >> $RUN_SAMPLE_DIR/$S.R1.fastq.gz; \
							rm -f $RUN_SAMPLE_DIR/$S.R2.fastq.gz; \
							touch $RUN_SAMPLE_DIR/$S.R2.fastq; \
							$GZ $RUN_SAMPLE_DIR/$S.R2.fastq; \
						fi;
						rm -f $RUN_SAMPLE_DIR/$S.R0.fastq.gz;
					" >> $FASTQ_MK
					FASTQ_MK_ALL=$FASTQ_MK_ALL" $RUN_SAMPLE_DIR/$S.R1.fastq.gz"

					echo "$RUN_SAMPLE_DIR/$S.R2.fastq.gz: $F
						#touch $RUN_SAMPLE_DIR/$S.I2.fastq
						#$GZ $RUN_SAMPLE_DIR/$S.I2.fastq
					" >> $FASTQ_MK
					FASTQ_MK_ALL=$FASTQ_MK_ALL" $RUN_SAMPLE_DIR/$S.R2.fastq.gz"

					echo "$RUN_SAMPLE_DIR/$S.I1.fastq.gz: $F
						touch $RUN_SAMPLE_DIR/$S.I1.fastq
						$GZ $RUN_SAMPLE_DIR/$S.I1.fastq
					" >> $FASTQ_MK
					FASTQ_MK_ALL=$FASTQ_MK_ALL" $RUN_SAMPLE_DIR/$S.I1.fastq.gz"

					echo "$RUN_SAMPLE_DIR/$S.I2.fastq.gz: $F
						touch $RUN_SAMPLE_DIR/$S.I2.fastq
						$GZ $RUN_SAMPLE_DIR/$S.I2.fastq
					" >> $FASTQ_MK
					FASTQ_MK_ALL=$FASTQ_MK_ALL" $RUN_SAMPLE_DIR/$S.I2.fastq.gz"

				fi;
				# || (( \$\$(diff <($UNGZ -c $RUN_SAMPLE_DIR/$S.R1.fastq.gz | paste - - - - | cut -f1 -d\$\$'\t') <($UNGZ -c $RUN_SAMPLE_DIR/$S.R2.fastq.gz | paste - - - - | cut -f1 -d\$\$'\t') | head -n 1 | wc -l) )); then
			fi;

		else
			(($VERBOSE)) && echo "#[INFO] Input file '$RUN_SAMPLE_DIR/$S.R1.fastq.gz' DOES exist"
		fi;

		# FASTQ PROCESSING
		if ((1)); then

			RUN_SAMPLE_DIR_SEQUENCING=$RUN_SAMPLE_DIR/$S.sequencing
			FASTP_HTML=$RUN_SAMPLE_DIR_SEQUENCING/$S.fastp.html
			FASTP_JSON=$RUN_SAMPLE_DIR_SEQUENCING/$S.fastp.json
			FASTP_LOG=$RUN_SAMPLE_DIR_SEQUENCING/$S.fastp.log
			FASTP_ERR=$RUN_SAMPLE_DIR_SEQUENCING/$S.fastp.err

			if [ ! -e $FASTP_LOG ]; then

				(($VERBOSE)) && echo "#[INFO] FASTQ processing (Adaptors, UMIs, quality)"


				# FASTP MK
				echo "$RUN_SAMPLE_DIR_SEQUENCING/$S: $RUN_SAMPLE_DIR/$S.R1.fastq.gz $RUN_SAMPLE_DIR/$S.R2.fastq.gz $RUN_SAMPLE_DIR/$S.I1.fastq.gz $RUN_SAMPLE_DIR/$S.I2.fastq.gz
					# Create sample sequencing folder
					mkdir -p $RUN_SAMPLE_DIR/$S.sequencing;
					# FASTP parameters
					echo ' --disable_trim_poly_g --disable_length_filtering ' > \$@.fastp.param;
					echo ' --thread \$(FASTP_THREADS_BY_SAMPLE) ' >> \$@.fastp.param;
					echo ' --compression \$(FASTQ_COMPRESSION_LEVEL) ' >> \$@.fastp.param;
					# OUTPUT
					echo '-h $FASTP_HTML -j $FASTP_JSON' >> \$@.fastp.param;
					# Paired-End or Single-End
					if (( \$\$($UNGZ -c $RUN_SAMPLE_DIR/$S.R2.fastq.gz | head -n 1 | wc -l) )); then \
						(($VERBOSE)) && echo '#[INFO] Paired-End Processing' ; \
						echo '-i $RUN_SAMPLE_DIR/$S.R1.fastq.gz -I $RUN_SAMPLE_DIR/$S.R2.fastq.gz' >> \$@.fastp.param ; \
						echo '-o $RUN_SAMPLE_DIR/$S.sequencing/$S.R1.processed.fastq.gz -O $RUN_SAMPLE_DIR/$S.sequencing/$S.R2.processed.fastq.gz' >> \$@.fastp.param ; \
						if (($DETECT_ADAPTER_FOR_PE)); then echo $FASTP_PARAM' --detect_adapter_for_pe' >> \$@.fastp.param ; fi ; \
					else \
						(($VERBOSE)) && echo '#[INFO] Single-End Processing' ; \
						echo '-i $RUN_SAMPLE_DIR/$S.R1.fastq.gz' >> \$@.fastp.param ; \
						echo '-o $RUN_SAMPLE_DIR/$S.sequencing/$S.R1.processed.fastq.gz' >> \$@.fastp.param ; \
					fi;
					# UMI Extraction
					if [ '$UMI_BARCODE_PATTERN' != '' ] && ! (( \$\$($UNGZ -c $RUN_SAMPLE_DIR/$S.R1.fastq.gz | head -n1 | awk -F' ' '{n_read_name=split(\$\$1,read_name,\":\"); if (read_name[n_read_name] ~ /[A-Z]/) { print read_name[n_read_name]} }' | wc -l) )); then \
						echo '--umi --umi_loc per_read --umi_len '${#UMI_BARCODE_PATTERN} >> \$@.fastp.param; \
					fi;
					# Read quality filtering
					if [ '$FASTQ_QUALITY_FILTERING' != '' ]; then \
						echo ' --cut_mean_quality=$FASTQ_QUALITY_FILTERING ' >> \$@.fastp.param; \
					else \
						echo ' --disable_quality_filtering ' >> \$@.fastp.param; \
					fi;
					# Report title
					echo ' --report_title=$RUU/$S ' >> \$@.fastp.param;
					# FASTP Process
					$FASTP \$\$(cat \$@.fastp.param) 1>$FASTP_LOG 2>$FASTP_ERR
					if (($FASTQ_DEMULTIPLEXING_KEEP)); then \
						if [ -e $RUN_SAMPLE_DIR/$S.R1.fastq.gz ] ; then $COMMAND_COPY_NO_COMPRESS $RUN_SAMPLE_DIR/$S.R1.fastq.gz $RUN_SAMPLE_DIR/$S.sequencing/$S.R1.demultiplexing.fastq.gz; fi; \
						if [ -e $RUN_SAMPLE_DIR/$S.R2.fastq.gz ] ; then $COMMAND_COPY_NO_COMPRESS $RUN_SAMPLE_DIR/$S.R2.fastq.gz $RUN_SAMPLE_DIR/$S.sequencing/$S.R2.demultiplexing.fastq.gz; fi; \
						if [ -e $RUN_SAMPLE_DIR/$S.I1.fastq.gz ] ; then $COMMAND_COPY_NO_COMPRESS $RUN_SAMPLE_DIR/$S.I1.fastq.gz $RUN_SAMPLE_DIR/$S.sequencing/$S.I1.demultiplexing.fastq.gz; fi; \
						if [ -e $RUN_SAMPLE_DIR/$S.I2.fastq.gz ] ; then $COMMAND_COPY_NO_COMPRESS $RUN_SAMPLE_DIR/$S.I2.fastq.gz $RUN_SAMPLE_DIR/$S.sequencing/$S.I2.demultiplexing.fastq.gz; fi; \
					fi;
					if [ -e $RUN_SAMPLE_DIR/$S.sequencing/$S.R1.processed.fastq.gz ] ; then mv $RUN_SAMPLE_DIR/$S.sequencing/$S.R1.processed.fastq.gz $RUN_SAMPLE_DIR/$S.R1.fastq.gz; fi;
					if [ -e $RUN_SAMPLE_DIR/$S.sequencing/$S.R2.processed.fastq.gz ] ; then mv $RUN_SAMPLE_DIR/$S.sequencing/$S.R2.processed.fastq.gz $RUN_SAMPLE_DIR/$S.R2.fastq.gz; fi;

				" >> $FASTP_MK
				FASTP_MK_ALL=$FASTP_MK_ALL" $RUN_SAMPLE_DIR_SEQUENCING/$S"

				#cat $FASTP_MK

			fi;

		fi;


		# OTHER_FILES
		if [ "$OFCF" != "" ]; then

			(($VERBOSE)) && echo "#[INFO] Copy OTHER FILES."
			for OFCF_ONE in $(echo $OFCF | tr "+" " "); do
				# Original other file
				OF=$(echo $OFCF_ONE | awk -F: '{print $1}')
				# Copy other file
				CF=$(echo $OFCF_ONE | awk -F: '{print $2}')
				[ "$CF" == "" ] && CF=$(basename $OF)
				(($DEBUG)) && echo "#[INFO] Copy OTHER FILES '$OF' to '$CF'."
				#cp $OF_ONE $RUN_SAMPLE_DIR/$(basename $OF_ONE)
				if [ -d $OF ]; then
					mkdir -p $RUN_SAMPLE_DIR/$CF
					#cp -Rf $OF/* $RUN_SAMPLE_DIR/$CF/
					$COMMAND_COPY $OF/* $RUN_SAMPLE_DIR/$CF/
				else
					mkdir -p $RUN_SAMPLE_DIR/$(dirname $CF)
					#cp $OF $RUN_SAMPLE_DIR/$CF
					$COMMAND_COPY $OF/* $RUN_SAMPLE_DIR/$CF/
				fi;
			done;

		fi;

		# Copy BED

		if true; then
			if [ ! -e $RUN_SAMPLE_DIR/$S.bed ]; then
				if [ "$B" != "" ]; then
					if [[ $B =~ .bed$ ]]; then
						if [ -e $B ] && [ "$B" != "" ] && [ ! -e $RUN_SAMPLE_DIR/$S.bed ]; then

							(($VERBOSE)) && echo "#[INFO] Copy original BED file."
							$COMMAND_COPY -p $B $RUN_SAMPLE_DIR/$S.original.bed;

							(($VERBOSE)) && echo "#[INFO] Sort/Merge/Normalize BED file."
							$BEDTOOLS sort -i $B  | $STARK_BED_NORMALIZATION | $BEDTOOLS merge -i - -c 4 -o distinct > $RUN_SAMPLE_DIR/$S.bed;
							
							touch $RUN_SAMPLE_DIR/$S.bed -r $B;
							$COMMAND_COPY $B $RUN_SAMPLE_DIR/$S.metrics.bed;
							touch $RUN_SAMPLE_DIR/$S.manifest;
							touch $RUN_SAMPLE_DIR/$S.manifest -r $B;
							touch $RUN_SAMPLE_DIR/$S.metrics.bed -r $B;
							echo -e $(basename $B)"\tfrom option - bed" > $RUN_SAMPLE_DIR/$S.bed_name
							echo -e $(basename $B)"\tfrom option - bed" > $RUN_SAMPLE_DIR/$S.manifest_name
							#export BED=$RUN_SAMPLE_DIR/$S.bed
						fi;
					else
						if [ -e $B ] && [ "$B" != "" ] && [ ! -e $RUN_SAMPLE_DIR/$S.manifest ]; then
							(($VERBOSE)) && echo "#[INFO] Copy Manifest file."
							#cp -p $B $RUN_SAMPLE_DIR/$S.manifest;
							$COMMAND_COPY -p $B $RUN_SAMPLE_DIR/$S.manifest;
							touch $RUN_SAMPLE_DIR/$S.manifest -r $B;
							#echo -e "\tfrom manifest file" > $RUN_SAMPLE_DIR/$S.bed_name
							echo -e $(basename $B)"\tfrom option - manifest" > $RUN_SAMPLE_DIR/$S.bed_name
							echo -e $(basename $B)"\tfrom option - manifest" > $RUN_SAMPLE_DIR/$S.manifest_name
							#export MANIFEST=$RUN_SAMPLE_DIR/$S.manifest;
						fi;
					fi;
				fi;
			fi;
		fi;

		# Genes list

		if true; then
			# Name of LIST.GENES
			BEDFILE_GENES_LIST_ONE=$RUN_SAMPLE_DIR/$S.list.genes
			if [ ! -e $BEDFILE_GENES_LIST_ONE ]; then
				# Add Genes
				if [ "$G" != "" ]; then
					(($VERBOSE)) && echo "#[INFO] Create LIST.GENES file '$RUN_SAMPLE_DIR/$S.list.genes' and copy .genes files."
					#> $RUN_SAMPLE_DIR/$S.list.genes;
					for G_ONE in $(echo $G | tr "+" " "); do
						# add '.genes' extension if not exists
						[ "${G_ONE##*.}" != "genes" ] && G_ONE_TARGET=$(basename $G_ONE)".genes" || G_ONE_TARGET=$(basename $G_ONE);
						# remove Sample name (mainly in case of relauch)
						G_ONE_TARGET=$(echo $G_ONE_TARGET | sed "s/$S\.//gi")
						#$COMMAND_COPY $G_ONE $RUN_SAMPLE_DIR/$S.$G_ONE_TARGET;
						$BEDTOOLS sort -i $G_ONE | $STARK_BED_NORMALIZATION > $RUN_SAMPLE_DIR/$S.$G_ONE_TARGET;
						echo $S.$G_ONE_TARGET >> $BEDFILE_GENES_LIST_ONE
					done;
				elif [ -s $B ] && [ "$B" != "" ]; then
					# GEnerate BEDFILE_GENES from design
					#if [ -s $B ]; then
					(($VERBOSE)) && echo "#[INFO] Create LIST.GENES file '$RUN_SAMPLE_DIR/$S.list.genes' from Design file '$B'."
					# Generate intermediate bed from Design (BED or manifest)
					BEDFILE_GENES_INTERMEDIATE=$RUN_SAMPLE_DIR/$S.genes.tmp
					BEDFILE_GENES_ONE=$RUN_SAMPLE_DIR/$S.from_design.genes

					if [[ $B =~ .bed$ ]]; then
						#cut -f1,2,3 $B > $BEDFILE_GENES_INTERMEDIATE;
						$BEDTOOLS sort -i $B | $BEDTOOLS merge -i - | cut -f1,2,3 > $BEDFILE_GENES_INTERMEDIATE;
					else
						$CAP_ManifestToBED --input $B --output $BEDFILE_GENES_INTERMEDIATE.tmp --output_type "region" --type=PCR 1>/dev/null;
						#cut -f1,2,3 $BEDFILE_GENES_INTERMEDIATE.tmp > $BEDFILE_GENES_INTERMEDIATE;
						$BEDTOOLS sort -i $BEDFILE_GENES_INTERMEDIATE.tmp | $BEDTOOLS merge -i - | cut -f1,2,3 > $BEDFILE_GENES_INTERMEDIATE;
						rm -f $BEDFILE_GENES_INTERMEDIATE.tmp;
					fi;

					if [ -s $BEDFILE_GENES_INTERMEDIATE ] ; then
						(($VERBOSE)) && echo "#[INFO] Create LIST.GENES file '$RUN_SAMPLE_DIR/$S.list.genes' with intersection between Design file '$B' and RefSeq '$REFSEQ_GENES'."
						#echo "MANIFEST.bed to bedfile_genes_list : $bedfile_genes_list";
						$BEDTOOLS intersect -wb -a $BEDFILE_GENES_INTERMEDIATE -b $REFSEQ_GENES | cut -f7 | sort -u > $BEDFILE_GENES_INTERMEDIATE.intersect;
						sort -k4 $REFSEQ_GENES > $BEDFILE_GENES_INTERMEDIATE.refseq;
						# Join
						join -1 1 -2 4 $BEDFILE_GENES_INTERMEDIATE.intersect $BEDFILE_GENES_INTERMEDIATE.refseq -o 2.1,2.2,2.3,2.4,2.5,2.6 | sort -u -k1,2 | tr " " "\t" | $BEDTOOLS sort | $BEDTOOLS merge -c 4,5,6 -o distinct,collapse,first | awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t0\t"$6}' | $STARK_BED_NORMALIZATION > $BEDFILE_GENES_ONE;
						echo $(basename $BEDFILE_GENES_ONE) > $BEDFILE_GENES_LIST_ONE;
					else
						(($VERBOSE)) && echo "#[WARNING] Generating GENES failed";
					fi;

					rm -f $BEDFILE_GENES_INTERMEDIATE*
					# ADD GENES in input files
					G=$BEDFILE_GENES_ONE

				fi;

			else
				(($VERBOSE)) && echo "#[INFO] LIST.GENES already exists."
			fi;
		fi;


		# Create genes.bed from genes
		(($VERBOSE)) && echo "#[INFO] Generate genes.bed from .genes files within LIST.GENES file '$RUN_SAMPLE_DIR/$S.list.genes'."
		if [ -e $BEDFILE_GENES_LIST_ONE ]; then
			for BEDFILE_GENES_LIST_ONE_GENES_FILE in $(cat $BEDFILE_GENES_LIST_ONE); do
				cp -p $RUN_SAMPLE_DIR/$BEDFILE_GENES_LIST_ONE_GENES_FILE $RUN_SAMPLE_DIR/$BEDFILE_GENES_LIST_ONE_GENES_FILE.bed
			done;
		fi;


		# Copy TRANSCRIPTS
		if [ "$T" != "" ] && [ ! -e $RUN_SAMPLE_DIR/$T.list.transcripts ]; then
			(($VERBOSE)) && echo "#[INFO] Create LIST.TRANSCRIPTS file '$RUN_SAMPLE_DIR/$S.list.transcripts' and concatenated file .transcripts"
			#echo $T | tr "+" "\n" > $RUN_SAMPLE_DIR/$S.list.transcripts;
			echo $T | tr "+" "\n" | xargs -l basename | sed "s/^/$S./gi" > $RUN_SAMPLE_DIR/$S.list.transcripts;
			#cat $(echo $T | tr "+" " ") > $RUN_SAMPLE_DIR/$S.transcripts;
			> $RUN_SAMPLE_DIR/$S.transcripts.tmp;
			for T_ONE in $(echo $T | tr "+" "\n"); do
				cat $T_ONE >> $RUN_SAMPLE_DIR/$S.transcripts.tmp;
				echo "" >> $RUN_SAMPLE_DIR/$S.transcripts.tmp;
			done
			cat $RUN_SAMPLE_DIR/$S.transcripts.tmp | grep -v "^$" > $RUN_SAMPLE_DIR/$S.transcripts
			rm -f $RUN_SAMPLE_DIR/$S.transcripts.tmp
		fi;


		# SAMPLE TAG
		if [ ! -e $RUN_SAMPLE_DIR/$S.tag ]; then
			(($VERBOSE)) && echo "#[INFO] Create TAG file."
			echo $TAG > $RUN_SAMPLE_DIR/$S.tag;
		fi;


		# ANALYSIS TAG
		if [ ! -e $RUN_SAMPLE_DIR/$S.analysis.tag ]; then
			(($VERBOSE)) && echo "#[INFO] Create Analysis TAG file."
			echo $ANALYSIS_TAG > $RUN_SAMPLE_DIR/$S.analysis.tag;
		fi;


		# SampleSheet
		if ((1)); then
		if [ ! -e $RUN_SAMPLE_DIR/$S.SampleSheet.csv ]; then
			(($VERBOSE)) && echo "#[INFO] Copy SampleSheet."
			#[ -f "$SAMPLESHEET_INPUT" ] && cp $SAMPLESHEET_INPUT $RUN_SAMPLE_DIR/$S.SampleSheet.csv
			[ -f "$SAMPLESHEET_INPUT" ] && $COMMAND_COPY $SAMPLESHEET_INPUT $RUN_SAMPLE_DIR/$S.SampleSheet.csv
			touch $RUN_SAMPLE_DIR/$S.SampleSheet.csv;
			if [ -e $RUN_SAMPLE_DIR/$S.manifest ]; then
				touch -f $RUN_SAMPLE_DIR/$S.SampleSheet.csv -r $RUN_SAMPLE_DIR/$S.manifest;
			fi;
		fi;
		fi;


		F_LIST="$F_LIST$F "
		Q_LIST="$Q_LIST$F_R2 "
		I1_LIST="$I1_LIST$I1 "
		I2_LIST="$I2_LIST$I2 "
		OFCF_LIST="$OFCF_LIST$OFCF "
		S_LIST="$S_LIST$S "
		B_LIST="$B_LIST$B "
		G_LIST="$G_LIST$G "
		T_LIST="$T_LIST$T "
		TAG_LIST="$TAG_LIST$TAG "


		# MAKEFILE_ANALYSIS_RUN
		RUNS_SAMPLES="$RUNS_SAMPLES$RU:$S "


		((I++))
		((NB_SAMPLE++))

	done;


	# THREADS
	############

	#NB_SAMPLE=1;

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


	echo "RUNS_SAMPLES=$RUNS_SAMPLES" >> $MAKEFILE_ANALYSIS_RUN
	echo "PIPELINES=$PIPELINES" >> $MAKEFILE_ANALYSIS_RUN
	#echo "INTERSEC=2" >> $MAKEFILE_ANALYSIS_RUN

	#echo "RUNS_SAMPLES=\"$RUNS_SAMPLES\"" >> $SHELL_ANALYSIS_RUN
	#echo "PIPELINES=\"$PIPELINES\"" >> $SHELL_ANALYSIS_RUN
	#echo "INTERSEC=2" >> $SHELL_ANALYSIS_RUN


	# RESULTS

	echo "#[INFO] *** Configuration"

	(($VERBOSE)) && echo "#[INFO] * ANALYSIS                "
	echo "#[INFO] ANALYSIS NAME             $RUU"
	echo "#[INFO] ANALYSIS TAG              $ANALYSIS_TAG"
	(($VERBOSE)) && echo "#[INFO] * SAMPLES                 "
	echo "#[INFO] SAMPLE NAMES              $S_LIST"
	echo "#[INFO] SAMPLE TAG                $TAG_LIST"
	(($VERBOSE)) && echo "#[INFO] FASTQ/BAM/CRAM            $F_LIST"
	(($VERBOSE)) && echo "#[INFO] FASTQ R2                  $Q_LIST"
	(($VERBOSE)) && echo "#[INFO] INDEX1                    $I1_LIST"
	(($VERBOSE)) && echo "#[INFO] INDEX2                    $I2_LIST"
	(($VERBOSE)) && echo "#[INFO] OTHER_FILES               $OFCF_LIST"
	(($VERBOSE)) && echo "#[INFO] DESIGN                    $B_LIST"
	(($VERBOSE)) && echo "#[INFO] GENES                     $G_LIST"
	(($VERBOSE)) && echo "#[INFO] TRANSCRIPTS               $T_LIST"
	(($VERBOSE)) && echo "#[INFO] * APPLICATION               "
	echo "#[INFO] APPLICATION NAME          $APP_NAME"
	(($VERBOSE)) && echo "#[INFO] APPLICATION FILE          "$(echo $ENV | sed s#$STARK_FOLDER_APPS/##gi)
	echo "#[INFO] GROUP                     $SAMPLE_GROUP"
	echo "#[INFO] PROJECT                   $SAMPLE_PROJECT"
	echo "#[INFO] PIPELINES                 $PIPELINES"
	echo "#[INFO] POST SEQUENCING           "$(echo $POST_SEQUENCING | tr "." "\n" | tac | tr "\n" " " )""
	echo "#[INFO] POST ALIGNMENT            "$(echo $POST_ALIGNMENT | tr "." "\n" | tac | tr "\n" " " )""
	echo "#[INFO] POST CALLING              "$(echo $POST_CALLING | tr "." "\n" | tac | tr "\n" " " )""
	echo "#[INFO] POST ANNOTATION           "$(echo $POST_ANNOTATION | tr "." "\n" | tac | tr "\n" " " )""
	echo "#[INFO] RESULTS                   $RESULTS"
	echo "#[INFO] REPOSITORY                $REPOSITORY"
	echo "#[INFO] ARCHIVES                  $ARCHIVES"
	(($VERBOSE)) && echo "#[INFO] RELEASE INFOS             $RELEASE_RUN"
	(($VERBOSE)) && echo "#[INFO] MAKEFILE CONFIGURATION    $MAKEFILE_ANALYSIS_RUN"
	#(($VERBOSE)) && echo "#[INFO] SHELL CONFIGURATION       $SHELL_ANALYSIS_RUN"
	(($VERBOSE)) && echo "#[INFO] LOGFILE                   $LOGFILE_RES_RUN"
	(($VERBOSE)) && echo "#[INFO] THREADS                   $THREADS"
	(($VERBOSE)) && echo "#[INFO] THREADS_BY_SAMPLE         $THREADS_BY_SAMPLE"
	

	STARK_QUEUED=$(source_app "$APP" "$STARK_FOLDER_APPS"; echo $STARK_QUEUED)
	STARK_RUNNING=$(source_app "$APP" "$STARK_FOLDER_APPS"; echo $STARK_RUNNING)
	STARK_COMPLETE=$(source_app "$APP" "$STARK_FOLDER_APPS"; echo $STARK_COMPLETE)
	if [ "$STARK_QUEUED" == "" ]; then STARK_QUEUED=STARKQueued.txt; fi;
	if [ "$STARK_RUNNING" == "" ]; then  STARK_RUNNING=STARKRunning.txt; fi;
	if [ "$STARK_COMPLETE" == "" ]; then  STARK_COMPLETE=STARKComplete.txt; fi;
	STARK_QUEUED_FILE=$RESULTS/$RUU/$STARK_QUEUED
	STARK_RUNNING_FILE=$RESULTS/$RUU/$STARK_RUNNING
	STARK_COMPLETE_FILE=$RESULTS/$RUU/$STARK_COMPLETE

	# RUNNING
	echo "#["$(date '+%Y%m%d-%H%M%S')"] RUN $RUN running with STARK ($STARK_VERSION)" > $STARK_RUNNING_FILE

	# TREADS
	#THREADS_BY_SAMPLE=$THREADS; # Allocate all thread to the sample because no other sample analysed in parallele

	echo "#[INFO] *** Process"

	# FASTQ && FASTP process

	echo "#[INFO] STARK Input Processing..."
	
	if [ -s $FASTQ_MK ] && [ -s $FASTP_MK ]; then
		(($VERBOSE)) && echo "#[INFO] Process Input data from FASTQ file(s) - multithreading mode [$THREADS]"
		if make -j $THREADS -e FASTP_THREADS_BY_SAMPLE=$THREADS_BY_SAMPLE -f $FASTP_MK -f $FASTQ_MK $FASTQ_MK_ALL $FASTP_MK_ALL 1>$FASTP_MK.log 2>$FASTP_MK.err; then
			(($VERBOSE)) && echo "#[INFO] STARK Input Processing done."
		else
			echo "#[ERROR] STARK Input Processing failed"
			cat $FASTP_MK.err
			exit 1
		fi;
	fi;
	! (($DEBUG)) && rm -f $FASTQ_MK $FASTQ_MK.log $FASTQ_MK.err $FASTP_MK $FASTP_MK.log $FASTP_MK.err


	echo "#[INFO] STARK Analysis Processing..."

	echo "["$(date '+%Y%m%d-%H%M%S')"] Main Analysis Process for Analysis '$RELEASE_RUN' START" >>$LOGFILE_RES_RUN


	MAKE_ERROR=0

	if ! make -k -j $THREADS -e ENV="$ENV" PARAM=$MAKEFILE_ANALYSIS_RUN $PARAMETERS $THREAD_PARAMETERS JAVA_MEMORY=$JAVA_MEMORY SNAPSHOT=0 VALIDATION=1 INPUT=$INPUT OUTDIR=$RESULTS RELEASE=$RELEASE_RUN FINAL_REPORT=$FINAL_REPORT_RUN ANALYSIS_REF=$ANALYSIS_REF -f $NGS_SCRIPTS/NGSWorkflow.mk 1>>$LOGFILE_RES_RUN 2>>$LOGFILE_RES_RUN; then
		MAKE_ERROR=1
	fi;
	rm -f $MAKEFILE_ANALYSIS_RUN $FINAL_REPORT_RUN $RELEASE_RUN
	echo "["$(date '+%Y%m%d-%H%M%S')"] Main Analysis Process for Analysis '$RELEASE_RUN' END" >>$LOGFILE_RES_RUN

	# LOG and END files

	rm -f $STARK_RUNNING_FILE

	I=0
	for S in $SAMPLE; do
		RU=${RUN_ARRAY[$I]};

		if [ "$RU" != "$RUU" ]; then
			((I++))
			continue;
		fi;
		
		# REPOSITORY
		# COPY of run/sample
		# List of folders
		# Repositories
		if [ "$REPOSITORY" != "" ] ; then
			for RESULTS_FOLDER_COPY_FOLDER in $(echo $REPOSITORY | tr "," " " | tr " " "\n" | sort -u ); do #
				if [ -d $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT/$RU ] || \
					mkdir -p $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT/$RU; then
					#echo "Copy $LOGFILE_RES_RUN $LOGFILE_DEM_RUN in $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT/$RU/ "
					cp -pf $LOGFILE_RES_RUN $LOGFILE_DEM_RUN $CONFIG_RES_RUN $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT/$RU/;
				fi;
			done;
		fi;
		# Archives
		if [ "$ARCHIVES" != "" ] ; then
			for RESULTS_FOLDER_COPY_FOLDER in $(echo $ARCHIVES | tr "," " " | tr " " "\n" | sort -u ); do #
				if [ -d $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT/$RU ] || \
					mkdir -p $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT/$RU; then
					#echo "Copy $LOGFILE_RES_RUN $LOGFILE_DEM_RUN in $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT/$RU/ "
					cp -pf $LOGFILE_RES_RUN $LOGFILE_DEM_RUN $CONFIG_RES_RUN $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT/$RU/;
				fi;
			done;
		fi;

	done



	### TEST ERROR

	if (( $(grep "\*\*\*" $LOGFILE_RES_RUN -c) || $MAKE_ERROR )); then

		echo "["`date '+%Y%m%d-%H%M%S'`"] Main Analysis Process for Analysis '$RELEASE_RUN' ERROR" >>$LOGFILE_RES_RUN
		STOP_ANALYSIS=$(date +%s)
		ANALYSIS_TIME=$(convertsecs $((STOP_ANALYSIS - START_ANALYSIS)))
		echo "#[INFO] Main Analysis Process for Analysis '$RELEASE_RUN' finished with ERRORS [`date`] - $ANALYSIS_TIME"
		#if (($VERBOSE)); then grep "\*\*\*" $LOGFILE_RES_RUN; fi;
		(($VERBOSE)) && grep "\*\*\*" $LOGFILE_RES_RUN -B10 | awk '{print "#[ERROR] LOG: "$0}';
		# RUNNING stop
		
		
		
		# EXIT
		exit 1;
		

		
		#continue;
	#fi;

	else

		# STARK COMPLETE
		echo "#["`date '+%Y%m%d-%H%M%S'`"] RUN $RUN analyzed by STARK ($STARK_VERSION)" >> $STARK_COMPLETE_FILE

		# STARK RUNNING stop
		#rm -f $STARK_RUNNING_FILE


		# STARK Copy Reports Files
		#if ((1)); then

		echo "#[INFO] STARK Copy Processing...";

		#THREADS_COPY=$THREADS
		THREADS_COPY=1

		I=0
		RESULTS_FOLDER_COPY_FOLDER_RUU_STARKCopyComplete_list=""

		> $RES_MK
		RES_MK_ALL=""

		#RES_MK_RSYNC_PARAM=" --temp-dir=$TMP_FOLDER_TMP/rsync_$RANDOM "

		for S in $SAMPLE; do
			RU=${RUN_ARRAY[$I]};

			if [ "$RU" != "$RUU" ]; then
				((I++))
				continue;
			fi;
			
			# REPOSITORY
			if ((1)); then
				# COPY of run/sample
				# List of folders
				RESULTS_FOLDER_COPY_ALL=""
				# Repositories
				if [ "$REPOSITORY" != "" ] ; then
					for RESULTS_FOLDER_COPY_FOLDER in $(echo $REPOSITORY | tr "," " " | tr " " "\n" | sort -u ); do #
						if [ -d $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT ] || \
							mkdir -p $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT; then
							RESULTS_FOLDER_COPY_ALL=$RESULTS_FOLDER_COPY_ALL" $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT:repository";
						fi;
					done;
				fi;
				# Archives
				if ((1)); then
					if [ "$ARCHIVES" != "" ] ; then
						for RESULTS_FOLDER_COPY_FOLDER in $(echo $ARCHIVES | tr "," " " | tr " " "\n" | sort -u ); do #
							if [ -d $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT ] || \
								mkdir -p $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT; then
								RESULTS_FOLDER_COPY_ALL=$RESULTS_FOLDER_COPY_ALL" $RESULTS_FOLDER_COPY_FOLDER/$SAMPLE_GROUP/$SAMPLE_PROJECT:archives";
							fi;
						done;
					fi;
				fi;


				RES_MK_ALL_RUU_S=""
				#RES_MK_REP_COMPLETE=""

				# Copy
				if [ "$RESULTS_FOLDER_COPY_ALL" != "$RESULTS" ] && [ "$RESULTS_FOLDER_COPY_ALL" != "" ] ; then
					for RESULTS_FOLDER_COPY_FOLDER_INFO in $RESULTS_FOLDER_COPY_ALL;
					do

						RESULTS_FOLDER_COPY_FOLDER=$(echo $RESULTS_FOLDER_COPY_FOLDER_INFO | awk -F: '{print $1}')
						RESULTS_FOLDER_COPY_FOLDER_INFO_TYPE=$(echo $RESULTS_FOLDER_COPY_FOLDER_INFO | awk -F: '{print $2}')

						### Repository type
						REPOSITORY_TYPE_FILE_PATTERNS=""
						if [ "$RESULTS_FOLDER_COPY_FOLDER_INFO_TYPE" == "repository" ]; then
							REPOSITORY_TYPE_FILE_PATTERNS=$REPOSITORY_FILE_PATTERNS
							ROOT_FILE_SOURCE=$RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/$RESULTS_SUBFOLDER_DATA
						elif [ "$RESULTS_FOLDER_COPY_FOLDER_INFO_TYPE" == "archives" ]; then
							REPOSITORY_TYPE_FILE_PATTERNS=$ARCHIVES_FILE_PATTERNS
							ROOT_FILE_SOURCE=$RESULTS/$RUU/$S
						else
							echo "#[ERROR] Repository type '$RESULTS_FOLDER_COPY_FOLDER_INFO_TYPE' failed"
						fi;
 

						#echo "#[INFO] Copying '$RUU/$S' files from '$RESULTS' to '$RESULTS_FOLDER_COPY_FOLDER' [$RESULTS_FOLDER_COPY_FOLDER_INFO_TYPE]..."
						(($VERBOSE)) && echo "#[INFO] STARK Copy '$RUU/$S' Process... from '$RESULTS' to '$RESULTS_FOLDER_COPY_FOLDER' [$RESULTS_FOLDER_COPY_FOLDER_INFO_TYPE]"
						 

						# Copy RUN files
						#RES_MK_COPY_ROOT_TMP=$TMP_FOLDER_TMP/rsync_$RANDOM$RANDOM
						if ! in_array $RESULTS_FOLDER_COPY_FOLDER/$RUU $RES_MK_ALL; then
							echo "$RESULTS_FOLDER_COPY_FOLDER/$RUU:
								#[INFO] Copying '$RUU'
								mkdir -p $RESULTS_FOLDER_COPY_FOLDER/$RUU
								$COMMAND_COPY_NO_COMPRESS \$\$(find -L $RESULTS/$RUU -mindepth 1 -maxdepth 1 -type f) $RESULTS_FOLDER_COPY_FOLDER/$RUU
								rm $RESULTS_FOLDER_COPY_FOLDER/$RUU/*copy.mk* $RESULTS_FOLDER_COPY_FOLDER/$RUU/$(basename $FINAL_REPORT_RUN).metrics $RESULTS_FOLDER_COPY_FOLDER/$RUU/$(basename $FINAL_REPORT_RUN).variants
								chmod $PERMS -R $RESULTS_FOLDER_COPY_FOLDER/$RUU $RESULTS_FOLDER_COPY_FOLDER/$RUU/* 
							" >> $RES_MK
							RES_MK_ALL=$RES_MK_ALL" $RESULTS_FOLDER_COPY_FOLDER/$RUU"
							RES_MK_ALL_RUU=$RES_MK_ALL_RUU" $RESULTS_FOLDER_COPY_FOLDER/$RUU"
						fi;
						

						# Copy STARK RESULTS files
						#RES_MK_COPY_ROOT_TMP=$TMP_FOLDER_TMP/rsync_$RANDOM$RANDOM
						if [ "$RESULTS_FOLDER_COPY_FOLDER_INFO_TYPE" == "repository" ]; then
							echo "$ROOT_FILE_SOURCE:
								#[INFO] Copying '$RUU/$S' - STARK Results
								mkdir -p $ROOT_FILE_SOURCE
								$COMMAND_COPY_NO_COMPRESS $RESULTS/$RUU/$S/* $ROOT_FILE_SOURCE
								chmod $PERMS -R $RESULTS_FOLDER_COPY_FOLDER/$RUU/$S $RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/* 
							" >> $RES_MK
						elif [ "$RESULTS_FOLDER_COPY_FOLDER_INFO_TYPE" == "archives" ]; then
							echo "$ROOT_FILE_SOURCE:
								#[INFO] Copying '$RUU/$S' - STARK Results
								#[INFO] Copying '$RUU/$S' - NO Copy
							" >> $RES_MK
						else
							echo "#[ERROR] Repository type failed"
						fi;
						RES_MK_ALL=$RES_MK_ALL" $ROOT_FILE_SOURCE"
						RES_MK_ALL_RUU_S=$RES_MK_ALL_RUU_S" $ROOT_FILE_SOURCE"

						
						# Copy ROOT FILE PATTERNS
						#echo "#[INFO] Copying REPOSITORY files patterns to $RESULTS_FOLDER_COPY_FOLDER/$RUU/$S..."
						if [ $RESULTS_SUBFOLDER_DATA != "" ]; then

							
							for ROOT_FILE_PATTERN_INFO in $REPOSITORY_TYPE_FILE_PATTERNS; do

								#RES_MK_COPY_ROOT_TMP=$TMP_FOLDER_TMP/rsync_$RANDOM$RANDOM

								ROOT_FILE_PATTERN=$(echo $ROOT_FILE_PATTERN_INFO | awk -F: '{print $1}')
								ROOT_FILE_PATTERN_TYPE=$(echo $ROOT_FILE_PATTERN_INFO | awk -F: '{print $2}')

								FILE_PATTERN_TYPE_LS=" ";
								[ "$ROOT_FILE_PATTERN_TYPE" == "FOLDER" ] && FILE_PATTERN_TYPE_LS=" -d ";

								if [[ $ROOT_FILE_PATTERN =~ '$SAMPLE' ]]; then
									eval ROOT_FILE_PATTERN_VAR="/"$(echo $ROOT_FILE_PATTERN | sed s/\$SAMPLE/\$S/gi)
								else
									ROOT_FILE_PATTERN_VAR="/$ROOT_FILE_PATTERN"
								fi;
								ROOT_FILE_PATTERN_VAR_TARGET=$(echo $ROOT_FILE_PATTERN_VAR | sed -e 's/[][\\^*+.$-]/____/g')

								echo "$RESULTS_FOLDER_COPY_FOLDER/$RUU/$S$ROOT_FILE_PATTERN_VAR: $ROOT_FILE_SOURCE
									#[INFO] Making  '$RUU/$S' - STARK Main Results
									mkdir -p $RESULTS_FOLDER_COPY_FOLDER/$RUU/$S
									#[INFO] Copying '$RUU/$S' - STARK Main Results
									for F_SOURCE in \$\$(ls $FILE_PATTERN_TYPE_LS $ROOT_FILE_SOURCE$ROOT_FILE_PATTERN_VAR 2>/dev/null); do \
										F_TARGET=$RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/\$\$(basename \$\$F_SOURCE); \
										F_SOURCE_BASE=\$\$(echo \$\$F_SOURCE | sed \"s#$ROOT_FILE_SOURCE/##\"); \
										F_TARGET_BASE=\$\$(echo \$\$F_TARGET | sed \"s#$RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/##\"); \
										[ '$ROOT_FILE_PATTERN_TYPE' == "FOLDER" ] && F_SOURCE=\$\$F_SOURCE/; \
										if [ ! -e \$\$F_TARGET ] && [ '$ROOT_FILE_PATTERN_TYPE' != "FOLDER" ] ; then \
											$COMMAND_LINK \$\$F_SOURCE \$\$F_TARGET 1>/dev/null 2>/dev/null; \
										fi; \
										if [ ! -e \$\$F_TARGET ]; then \
											rm -f \$\$F_TARGET; \
											$COMMAND_COPY_NO_COMPRESS \$\$F_SOURCE \$\$F_TARGET 1>/dev/null 2>/dev/null; \
										fi; \
										if [ ! -e \$\$F_TARGET ]; then \
											(($DEBUG)) && echo '#[ERROR] Copy file $F_SOURCE_BASE to $F_TARGET_BASE FAILED'; \
										fi; \
									done;
								" >> $RES_MK

								#RES_MK_ALL=$RES_MK_ALL" $ROOT_FILE_SOURCE$ROOT_FILE_PATTERN_VAR_TARGET"
								RES_MK_ALL_RUU_S=$RES_MK_ALL_RUU_S" $RESULTS_FOLDER_COPY_FOLDER/$RUU/$S$ROOT_FILE_PATTERN_VAR"
								
							done;
						fi;

						# CREATE CopyComplete file
						echo "$RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/STARKCopyComplete.txt: $RES_MK_ALL_RUU_S
							#(($VERBOSE)) && echo '#[INFO] STARK Copy $RUU/$S Complete'
							echo "["\$\$(date '+%Y%m%d-%H%M%S')"] Copy complete" >> $RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/STARKCopyComplete.txt
							chmod $PERMS $RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/STARKCopyComplete.txt
						" >> $RES_MK
						#RES_MK_REP_COMPLETE=$RES_MK_REP_COMPLETE" $RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/STARKCopyComplete.txt"
						RES_MK_ALL=$RES_MK_ALL" $RESULTS_FOLDER_COPY_FOLDER/$RUU/$S/STARKCopyComplete.txt"
						# CopyComplete file for RUN
						RESULTS_FOLDER_COPY_FOLDER_RUU_STARKCopyComplete_list="$RESULTS_FOLDER_COPY_FOLDER_RUU_STARKCopyComplete_list $RESULTS_FOLDER_COPY_FOLDER/$RUU"

					done;

				fi;
			fi;

			((I++))
		done

		echo "all: $RES_MK_ALL
		" >> $RES_MK

		#if make -j $THREADS --always-make -f $RES_MK all 1>$RES_MK.log 2>$RES_MK.err; then
		if [ -s $RES_MK ]; then
			(($VERBOSE)) && echo "#[INFO] Process Copy - multithreading mode [$THREADS_COPY]"
			if make -j $THREADS_COPY --always-make -f $RES_MK all 1>$RES_MK.log 2>$RES_MK.err; then
				(($VERBOSE)) && echo "#[INFO] STARK Copy Complete"
			else
				echo "#[ERROR] STARK Copy NOT Complete"
				cat $RES_MK.err;
				exit 1
			fi
		fi;
		! (($DEBUG)) && rm -f $RES_MK $RES_MK.log $RES_MK.err

	
		# CREATE CopyComplete file
		for RESULTS_FOLDER_COPY_FOLDER_RUU_STARKCopyComplete in $RESULTS_FOLDER_COPY_FOLDER_RUU_STARKCopyComplete_list; do
			echo "["$(date '+%Y%m%d-%H%M%S')"] Copy complete" >> $RESULTS_FOLDER_COPY_FOLDER_RUU_STARKCopyComplete/STARKCopyComplete.txt
			chmod $PERMS $RESULTS_FOLDER_COPY_FOLDER_RUU_STARKCopyComplete/STARKCopyComplete.txt 1>/dev/null 2>/dev/null
		done;


	fi;

	# Stop Analysis
	STOP_ANALYSIS=$(date +%s)
	ANALYSIS_TIME=$(convertsecs $((STOP_ANALYSIS - START_ANALYSIS)))

	echo "#[INFO] STARK Execution Time: $ANALYSIS_TIME"

	(($VERBOSE)) && echo "#[INFO] *** Stop Analysis         ["$(date)"] - $ANALYSIS_TIME"

	

done;


# Return
if (( $RUN_ANALYZED )); then
	exit 0;
else
	exit 1;
fi;
