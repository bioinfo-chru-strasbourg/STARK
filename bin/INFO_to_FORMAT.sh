#!/bin/bash
#################################
##
## STARK INFO to FORMAT transfert
##
#################################

SCRIPT_NAME="STARKINFOtoFORMAT"
SCRIPT_DESCRIPTION="STARK INFO to FORMAT transfert"
SCRIPT_RELEASE="1.0.0"
SCRIPT_DATE="17/01/2023"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="HUS/CPS"
SCRIPT_LICENCE="GNU GPLA V3"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 1.0.0-17/01/2023: Script creation\n";

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
	echo "# USAGE: $(basename $0) [--help] --vcf=<FILE> [options...]";
	echo "#";
	echo "### This script transfers INFO annotation to FORMAT annotation.";
	echo "#";
  	echo "# --vcf|input=<FILE>           Input VCF file ";
	echo "#                              Format: VCF with a uniq sample'";
  	echo "# --annotations=<STRING>       List of annotation in INFO";
	echo "#                              Format: 'annotation1,annotation2,...'";
	echo "#                              Default: ''";
  	echo "# --output=<FILE>              Output VCF file";
	echo "#                              Default: <VCF>.output.vcf.gz";
  	echo "# --threads=<INTEGER>          Number of threads for bcftools";
	echo "#                              Default: 1";
  	echo "# --bcftools=<BIN>             BCFTOOLS binary file";
  	echo "# --tabix=<BIN>                TABIX binary file";
	echo "# -v|--verbose                 Verbose mode";
	echo "# -d|--debug                   Debug mode";
	echo "# -n|--release                 Script Release";
	echo "# -h|--help                    Help message";
	echo "#";

}




####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "vdnh" --long "vcf:,input:,annotations:,output:,threads:,bcftools:,tabix:,verbose,debug,release,help" -- "$@" 2> /dev/null)
if [ $? -ne 0 ]; then
	:
	echo "#[ERROR] Error in the argument list:";
	echo "#[ERROR] $@"
	echo ""
	usage;
	exit;
fi;

PARAM=$@
DEBUG=0
VERBOSE=0

eval set -- "$ARGS"
while true
do
	case "$1" in
		--vcf|--input)
			VCF=$2
			shift 2
			;;
    	--annotations)
			ANNOT_LIST=$(echo "$2" | tr "," " ")
			shift 2
			;;
    	--output)
			VCF_OUTPUT=$2
			shift 2
			;;
  		--threads)
			THREADS=$2
			shift 2
			;;
  		--bcftools)
			BCFTOOLS=$2
			shift 2
			;;
  		--tabix)
			TABIX=$2
			shift 2
			;;
  		-h|--help)
			usage
			exit 0
			;;
		-v|--verbose)
			VERBOSE=1
			shift 1
			;;
		-d|--debug)
			DEBUG=1
			shift 1
			;;
		--) shift
			break
			;;
		*) 	echo "Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done

# header
(($NO_HEADER)) || header;



####################################################################################################################################
# Checking the input parameter
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if [ ! -e $VCF ] || [ "$VCF" == "" ]; then
	echo "#[ERROR] Required parameter: --vcf. Use --help to display the help." && echo "" && usage && exit 1;
fi
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


(($DEBUG)) && VERBOSE=1

# Input VCF
if [ ! -e $VCF ] || [ "$VCF" == "" ]; then
    echo "#[ERROR] NO input VCF '$VCF'"
    exit 1
fi;

# Output VCF
if [ "$VCF_OUTPUT" == "" ]; then
    VCF_OUTPUT=$VCF.output.vcf.gz
fi;

# threads
if [ "$THREADS" == "" ]; then
    THREADS=1
fi;

# bcftools
if [ "$BCFTOOLS" == "" ]; then
    BCFTOOLS="bcftools"
fi;

# tabix
if [ "$TABIX" == "" ]; then
    TABIX="tabix"
fi;


# Make output folder
mkdir -p $(dirname $VCF_OUTPUT)

# Verbose header
if (($VERBOSE)); then
    echo "#[INFO] Input VCF: '$VCF'"
    echo "#[INFO] Output VCF: '$VCF_OUTPUT'"
    echo "#[INFO] Annotations: '$ANNOT_LIST'"
    echo "#[INFO] Threads: $THREADS"
    echo "#[INFO] bcftools: $BCFTOOLS"
    echo "#[INFO] tabix: $TABIX"
fi

# Create Variables

# sample name
SAMPLE=$($BCFTOOLS query -l $VCF 2>/dev/null)
NB_SAMPLE=$(echo $SAMPLE | wc -w)
if [ $NB_SAMPLE -gt 1 ]; then
    echo "#[ERROR] Input VCF with more than 1 sample ($NB_SAMPLE)"
    exit 0
fi;
if (($DEBUG)); then
    echo "#[INFO] SAMPLE=$SAMPLE"
fi;

ANNOT_query=""
ANNOT_annotate=""
for ANNOT in $(echo $ANNOT_LIST | tr "," " "); do
    if (($($BCFTOOLS view --threads $THREADS -h $VCF 2>/dev/null | grep "##INFO=<ID=$ANNOT," -c))); then
        (($DEBUG)) && echo "#[INFO] Annotation '$ANNOT' in '$VCF'"
        ANNOT_query=$ANNOT_query"\t%$ANNOT"
        ANNOT_annotate=$ANNOT_annotate",FORMAT/$ANNOT"
    else
        (($VERBOSE)) && echo "#[WARNING] Annotation '$ANNOT' NOT in '$VCF'"
    fi;
done;


if (($DEBUG)); then
    echo "#[INFO] ANNOT_query=$ANNOT_query"
    echo "#[INFO] ANNOT_annotate=$ANNOT_annotate"
fi;

# Extract INFO/DP into a tab-delimited annotation file
$BCFTOOLS query -f "%CHROM\t%POS\t%REF\t%ALT$ANNOT_query\n" $VCF 2>/dev/null | bgzip -c > $VCF.tmp.annot.txt.gz
if (($DEBUG)); then
    echo "#[INFO] Annotation file:"
    bgzip -dc $VCF.tmp.annot.txt.gz | head
fi;

# Index the file with $TABIX
$TABIX -s1 -b2 -e2 $VCF.tmp.annot.txt.gz

# Create a header line for the new annotation
> $VCF.tmp.hdr.txt
for ANNOT in $(echo $ANNOT_LIST | tr "," " "); do
    $BCFTOOLS view --threads $THREADS -h $VCF 2>/dev/null | grep "##INFO=<ID=$ANNOT," | sed "s/^##INFO/##FORMAT/gi" >> $VCF.tmp.hdr.txt
done;

if (($DEBUG)); then
    echo "#[INFO] Annotation header file:"
    cat $VCF.tmp.hdr.txt | head
fi;

# Transfer the annotation to sample 'SAMPLE'
$BCFTOOLS annotate --threads $THREADS -s "$SAMPLE" -a $VCF.tmp.annot.txt.gz -h $VCF.tmp.hdr.txt -c CHROM,POS,REF,ALT$ANNOT_annotate $VCF -o $VCF_OUTPUT 2>/dev/null

$TABIX $VCF_OUTPUT 2>/dev/null

if (($DEBUG)); then
    echo "#[INFO] Output VCF file '$VCF_OUTPUT':"
    $BCFTOOLS view --threads $THREADS $VCF_OUTPUT 2>/dev/null | grep "^##" -v | head
fi;

# Check

if (($VERBOSE)); then
    echo "#[INFO] nb of variants in '$VCF': "$($BCFTOOLS view --threads $THREADS $VCF 2>/dev/null | grep "^#" -v | wc -l)
    echo "#[INFO] nb of variants in '$VCF_OUTPUT': "$($BCFTOOLS view --threads $THREADS $VCF_OUTPUT 2>/dev/null | grep "^#" -v | wc -l)
    
    if (($(diff <($BCFTOOLS view --threads $THREADS $VCF 2>/dev/null | grep "^#" -v | cut -f1,5) <($BCFTOOLS view --threads $THREADS $VCF_OUTPUT 2>/dev/null | grep "^#" -v | cut -f1,5) | wc -l))); then
        echo "#[INFO] diff:"
        diff <($BCFTOOLS view --threads $THREADS $VCF 2>/dev/null | grep "^#" -v | cut -f1,5) <($BCFTOOLS view --threads $THREADS $VCF_OUTPUT 2>/dev/null | grep "^#" -v | cut -f1,5)
        echo "#[ERROR] DIFF!!!"
        exit 1
    else
        echo "#[INFO] NO DIFF"
    fi;
fi;

# Clean
rm -f $VCF.tmp*

exit 0