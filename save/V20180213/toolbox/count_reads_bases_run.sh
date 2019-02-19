#!/bin/bash
#################################
##
## NGS environment
## Analysis of multiple RUNS
##
## author: Antony Le Bechec
##
#################################

SCRIPT_NAME="CountReadsBasesRun"
SCRIPT_DESCRIPTION="Count number of Reads and Bases of a RUN"
SCRIPT_RELEASE="0.9b"
SCRIPT_DATE="09/12/2015"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"

echo "#######################################";
echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
echo "# $SCRIPT_DESCRIPTION ";
echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © ";
echo "#######################################";

# Realse note
#RELEASE_NOTES="#\n"
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-09/12/2015: Script creation\n";

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
	echo "# USAGE:  $0 <ENV> <RUNS> ";
	echo "# ENV              Environnement file (define folders, tools...)";
	echo "# RUNS             List of RUN to launch. NEEDED";
	exit 0;
fi;


# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# INPUT
ENV=$1
RUNS=$2
ALIGNERS_INPUT=$3
CALLERS_INPUT=$4
ANNOTATORS_INPUT=$5
FILTER_SAMPLE=$6
SAMPLESHEET=$7
PARALLELIZATION=$8
REMOVE=$9

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


# RUNS
if [ "$RUNS" == "" ]; then
	echo "#[ERROR] NO RUN defined"
	exit 0;
fi;

for RUN in $RUNS; do
	if [ -d $RESULTS_FOLDER/$RUN ]; then
		echo "# RUN '$RUN'"
		nb=0;
		nbb=0;
		for f in $(ls $RESULTS_FOLDER/$RUN/*/*.fastqc/*fastqc/*_data.txt); do
			tot=$(grep "Total Sequence" $f | cut -f2);
			len=$(grep "Sequence length" $f | cut -f2);
			nb=$(($nb+$tot));
			nbb=$(($nbb+($tot*$len)));
		done;
		echo "# Reads   $nb";
		echo "# Bases   $nbb";
	fi;
done;


