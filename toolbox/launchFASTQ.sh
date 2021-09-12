#!/bin/bash
#################################
##
## NGS environment - Launch FASTQ
##
#################################

SCRIPT_NAME="LaunchFASTQ"
SCRIPT_DESCRIPTION="Launch FASTQ Analysis"
SCRIPT_RELEASE="0.9b"
SCRIPT_DATE="14/06/2015"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"

echo "#######################################";
echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
echo "# $SCRIPT_DESCRIPTION ";
echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © ";
echo "#######################################";

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
		    sh launchFASTQ.sh -f fastq -a aligners -c callers -n annotators -e environment [-h]

		Options:
		 	-f, --fastq
		 	This option is required. FASTQ, R1 +R2 optionaly, or BAM. NEEDED
		 	-a --aligners
		 	List of ALIGNERS to use.
		 	-c --callers
		 	List of CALLERS to use.
		 	-n, --annotators
		 	List of ANNOTATORS to use.
		 	-e --env
		 	Environnement file (define folders, tools...)
		  	-h, --help
		 	Print this message and exit the program.
		__EOF__
}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "f:a:c:n:e:h" --long "fastq:,aligners:,callers:,annotators:,env:,help" -- "$@" 2> /dev/null)
[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && \
	exit 1
eval set -- "$ARGS"
while true
do
	case "$1" in
		-f|--fastq)
			FASTQ="$2"
			shift 2 
			;;
		-a|--aligners)
			ALIGNERS="$2"
			shift 2 
			;;
		-c|--callers)
			CALLERS="$2"
			shift 2 
			;;
		-n|--annotators)
			ANNOTATORS="$2"
			shift 2 
			;;
		-e|--env)
			ENV="$2"
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

####################################################################################################################################
# Checking the input parameter
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if [ "$FASTQ" == "" ]
then
	echo "Option --fastq is required. " "Use -h or --help to display the help." && exit 1;
fi
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# OUTPUT
echo "# "
echo "# CONFIGURATION "
echo "################"
echo "# NGS BIN Folder:        "$NGS_BIN
echo "# MISEQ Folder:          "$MISEQ_FOLDER
echo "# DEMULTIPLEXING Folder: "$DEMULTIPLEXING_FOLDER
echo "# RESULTS Folder:        "$RESULTS_FOLDER
echo "# ANALYSIS Folder:       "$ANALYSIS_FOLDER
echo "# "

# MULTI_RUN_ANALYSIS_PROCESS
if [ "$MULTI_RUN_ANALYSIS_PROCESS" == "" ]; then
	MULTI_RUN_ANALYSIS_PROCESS=0 # 0 OR 1 OR 2 for only multireporting
fi;

# ALIGNERS
if [ "$ALIGNERS" == "" ]; then
	#ALIGNERS="bwamem bwasw" # bwasw bwamem
	#ALIGNERS="bwamem bwasw bwamemUnclipped" # bwasw
	ALIGNERS="bwamem bwamemUnclipped" # bwasw bwamemUnclipped
	#ALIGNERS="bwamemUnclipped" # bwasw bwamemUnclipped
	#ALIGNERS="bwamem" # bwasw
	echo "#[WARNING] NO ALIGNERS defined. Default ALIGNERS '$ALIGNERS' will be used."
fi;
# CALLERS
if [ "$CALLERS" == "" ]; then
	CALLERS="gatkHC gatkUG"
	#CALLERS="gatkHC"
	#CALLERS="gatkHC gatkUG MutaCaller"
	echo "#[WARNING] NO CALLERS defined. Default CALLERS '$CALLERS' will be used."
fi;
# ANNOTATORS
if [ "$ANNOTATORS" == "" ]; then
	#ANNOTATORS="trakxs"
	ANNOTATORS="vap"
fi;
# INTERSEC
if [ "$INTERSEC" == "" ]; then
	#INTERSEC="3"
	INTERSEC=$(echo "scale=0;  ( ( "$(echo $ALIGNERS | wc -w )" * "$(echo $CALLERS | wc -w )" * "$(echo $ANNOTATORS | wc -w )" ) + 1 ) / 2" | bc )
	if [ $INTERSEC -lt 2 ]; then INTERSEC=2; fi
	echo "#[WARNING] NO INTERSEC defined. Default INTERSEC '$INTERSEC' will be used."
fi;
# ENV
if [ "$ENV" == "" ] || [ ! -s $ENV ]; then
	ENV="";
	echo "#[WARNING] NO ENV defined. Default ENV used."
fi;
if [ -s $ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$ENV;
elif [ -s $SCRIPT_DIR/$ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$SCRIPT_DIR/$ENV;
elif [ "$ENV" == "" ] || [ ! -s $ENV ]; then
	if [ -s $SCRIPT_DIR/"env.sh" ]; then
		ENV=$SCRIPT_DIR/"env.sh";
	else
		ENV="";
	fi;
fi;


# ENV
source $SCRIPT_DIR/env.sh
#DEMULTIPLEXING=/media/IRC/RES2/demultiplexing
#DEMULTIPLEXING=$DEMULTIPLEXING_FOLDER
#$RESULTS_FOLDER=/media/miseq/RES
mkdir -p $DEMULTIPLEXING_FOLDER
mkdir -p $RESULTS_FOLDER
mkdir -p $ANALYSIS_FOLDER

# DATABASE_INTEGRATION_SCRIPT
#DATABASE_INTEGRATION_SCRIPT=$HOWARD/DBintegration.sh


# RUNS
if [ "$FASTQ" == "" ]; then
	echo "#[ERROR] NO FASTQ defined"
	exit 0;
fi;

# TEST RUNS validity
for F in $FASTQ;
do
	if [ ! -s "$F" ]; then
		echo "#[WARNING] FASTQ '$F' doesn't exist"
		#exit 0;
	fi;
done;

#RELEASE="V"`date '+%Y%m%d-%H%M%S'`
#echo "# RUNs '$RUNS' - Release '$RELEASE'"
#LOG=$ANALYSIS_FOLDER/run_analysis.$RELEASE.launch.log
#$SCRIPT_DIR/runs_analysis.sh $MULTI_RUN_ANALYSIS_PROCESS "$RUNS" "$ALIGNERS" "$CALLERS" "$ANNOTATORS" "$FILTER_SAMPLE" "$INTERSEC" $ENV "1" 1>>$LOG 2>>$LOG

FASTQ_R1=$(echo $FASTQ | awk '{print $1}')
FASTQ_R2=$(echo $FASTQ | awk '{print $2}')

OUTDIR=$(dirname $FASTQ_R1)

echo "FASTQ_R1=$FASTQ_R1 FASTQ_R2=$FASTQ_R2"

#touch /tmp/
#/media/IRCV2/V2/DEM/Somatic3Test
MAKEFILE_ANALYSIS_RUN="/tmp/$(basename $0).makafile.param.$$.tmp"
echo "##################################
# PARAMETERS   
##################################
#
#RUNS=
RUNS_SAMPLES=

#########################
# Additional parameters #
#########################

FASTQ_R1=$FASTQ_R1
FASTQ_R2=$FASTQ_R2
ALIGNERS=$ALIGNERS
CALLERS=$CALLERS
ANNOTATORS=$ANNOTATORS
INTERSEC=$INTERSEC
" > $MAKEFILE_ANALYSIS_RUN

RELEASE_RUN="/tmp/$(basename $0).release.$$.tmp"


make -k -j $THREADS -e FASTQ_R1=$FASTQ_R1 FASTQ_R2=$FASTQ_R2 NGSEnv=$NGS_FOLDER ENV=$ENV PARAM=$MAKEFILE_ANALYSIS_RUN OUTDIR=$OUTDIR RELEASE=$RELEASE_RUN -f $NGS_SCRIPTS/STARK.launch.analysis.mk


exit 0;


exit 1;



