#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARKDatabases"
SCRIPT_DESCRIPTION="STARK change contig ID to chr*"
SCRIPT_RELEASE="0.9b"
SCRIPT_DATE="24/12/2018"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-24/12/2018: Script creation\n";

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
	echo -e $RELEASE_NOTES
}

# Usage
function usage {
	echo "# USAGE: cat <VCF> | $(basename $0) [options...]";
	echo "# -c/--chr=<string>              Chromosome prefix (default 'chr').";

	echo "# -v/--verbose                   VERBOSE option";
	echo "# -d/--debug                     DEBUG option";
	echo "# -n/--release                   RELEASE option";
	echo "# -h/--help                      HELP option";
	echo "#";

}

# header
#header;

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "c:vdnh" --long "chr:,verbose,debug,release,help" -- "$@" 2> /dev/null)
# || [ -z $@ ]

PARAM=$@
#PARAM=$(echo $@ | tr "\n" " ")
#echo $PARAM;
#PARAM=$(echo $ARGS | sed s/--//gi);
#exit 0;

eval set -- "$ARGS"
while true
do
	#echo "$1=$2"
	#echo "Eval opts";
	case "$1" in
		-c|--chr)
			CHR="$2"
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
		*) 	echo "# Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done

# CHR prefix
if [ -z ${CHR+x} ]; then
	CHR="chr"
fi


awk '{
# table
chr[1]="1"
chr[2]="2"
chr[3]="3"
chr[4]="4"
chr[5]="5"
chr[6]="6"
chr[7]="7"
chr[8]="8"
chr[9]="9"
chr[10]="10"
chr[11]="11"
chr[12]="12"
chr[13]="13"
chr[14]="14"
chr[15]="15"
chr[16]="16"
chr[17]="17"
chr[18]="18"
chr[19]="19"
chr[20]="20"
chr[21]="21"
chr[22]="22"
chr[23]="X"
chr[24]="Y"
chr[12920]="M"
if($0 !~ /^#/) {
    if(match($0,/(NC_)(0*)([0-9]*)\.([0-9]*)(\t)(.*)/,m)) {
      print "'$CHR'"chr[m[3]]m[5]m[6]; 
    }
    else
    if(match($0,/([0-9XYM]*)(\t)(.*)/,m)) {
      print "'$CHR'"m[1]m[2]m[3]; 
    }
    else
      print "'$CHR'"$0
    }
else
  if(match($0,/(##contig=<ID=)(NC_)(0*)([0-9]*)\.([0-9]*)([,>])(.*)/,m)) {
    print m[1]"'$CHR'"chr[m[4]]m[6]m[7]m[8]; 
  }
else
  if(match($0,/(##contig=<ID=)(.*)/,m)) print m[1]"'$CHR'"m[2];
else
  print $0
}' 


