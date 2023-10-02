#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARKCheckVCF"
SCRIPT_DESCRIPTION="STARK Check 2 VCF"
SCRIPT_RELEASE="0.9b"
SCRIPT_DATE="17/10/2016"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-17/10/2016: Script creation\n";

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
	echo "# RELEASE NOTES:";
	echo -e $RELEASE_NOTES
}

# Usage
function usage {
	echo "# USAGE: $(basename $0) --vcf=<VCF>|--ref=<VCF_REF> [options...]";
	echo "# -f/--vcf            VCF to compare with the referecen VCF";
	echo "# -r/--ref            Refrence VCF";
	echo "# -v/--verbose          VERBOSE option";
	echo "# -d/--debug            DEBUG option";
	echo "# -n/--release          RELEASE option";
	echo "# -h/--help             HELP option";
	echo "#";
}

# header
header;

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "f:r:vdnh" --long "vcf:,ref:,verbose,debug,release,help" -- "$@" 2> /dev/null)
if [ $? -ne 0 ]; then
	:
fi;
PARAM=$@

 
eval set -- "$ARGS"
while true
do
	#echo "$1=$2"
	#echo "Eval opts";
	case "$1" in
		-f|--vcf)
			VCF="$2"
			shift 2 
			;;
		-r|--ref)
			REF=$2
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

if [ ! -e $VCF ]; then
	echo "#[ERROR] NO VCF";
	usage;
	exit;
fi
if [ ! -e $GENOME ]; then
	echo "#[ERROR] NO REF";
	usage;
	exit;
fi

TMP=/tmp


echo "# VCF: $VCF"
echo "# GENOME: $GENOME"

NBVARIANT_VCF=$(grep -cv ^# $VCF) 
NBVARIANT_REF=$(grep -cv ^# $GENOME) 

if [ "$NBVARIANT_VCF" -eq "$NBVARIANT_REF" ]; then
	echo "#[OK] VCF and REF have same number of variants: $NBVARIANT_VCF";
	
	# Variants VCF
	VARIANTS_VCF=$TMP/$RANDOM$RANDOM
	grep -v ^# $VCF | cut -f1-5 > $VARIANTS_VCF

	# Variants REF
	VARIANTS_REF=$TMP/$RANDOM$RANDOM
	grep -v ^# $GENOME | cut -f1-5 > $VARIANTS_REF
	
	DIFF=$(diff $VARIANTS_VCF $VARIANTS_REF);
	
	if [ "$DIFF" == "" ]; then
		echo "#[OK] VCF and REF have same $NBVARIANT_VCF variants";	
	else
		echo "#[ko] VCF and REF DO NOT have same variants";
		echo -e "$DIFF";
	fi;
	
else
	echo "#[ko] VCF and REF DO NOT have same number of variants: $NBVARIANT_VCF/$NBVARIANT_REF";
fi;



