#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="Launch.Validation"
SCRIPT_DESCRIPTION="Launch Method Validation"
SCRIPT_RELEASE="0.9b"
SCRIPT_DATE="03/06/2015"
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
	echo "# USAGE: $0 <VALIDATION> <GOLDSTANDARD> <DATA> <ENV>";
	echo "# VALIDATION           Validation release";
	echo "# GOLDSTANDARD         VCF with knwon variants. NEEDED";
	echo "# DATA                 VCF with found variants, NEEDED";
	exit 0;
fi;

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# INPUT
VALIDATION=$1
GOLDSTANDARD=$2
DATA=$3




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
if [ "$ENV" != "" ]; then
	source $ENV;
fi;

# VALIDATION
if [ "$VALIDATION" == "" ]; then
	VALIDATION="V"`date '+%Y%m%d-%H%M%S'`
fi;

# GOLDSTANDARD
if [ "$GOLDSTANDARD" == "" ] || [ ! -s $GOLDSTANDARD ]; then
	echo "#[ERROR] NO GOLDSTANDARD '$GOLDSTANDARD' valid."
	exit 1;
fi;

# DATA
if [ "$DATA" == "" ] || [ ! -s $DATA ]; then
	echo "#[ERROR] NO DATA '$DATA' valid."
	exit 1;
fi;



GENOME=$GENOMES/$ASSEMBLY/$ASSEMBLY.fa

#DEMULTIPLEXING=/media/IRC/RES2/demultiplexing
#DEMULTIPLEXING=$DEMULTIPLEXING_FOLDER
#$RESULTS_FOLDER=/media/miseq/RES
mkdir -p $DEMULTIPLEXING_FOLDER
mkdir -p $RESULTS_FOLDER
mkdir -p $VALIDATION_FOLDER

# OUTPUT
echo "# "
echo "# CONFIGURATION "
echo "################"
echo "# NGS BIN Folder:        "$NGS_BIN
echo "# MISEQ Folder:          "$MISEQ_FOLDER
echo "# DEMULTIPLEXING Folder: "$DEMULTIPLEXING_FOLDER
echo "# RESULTS Folder:        "$RESULTS_FOLDER
echo "# ANALYSIS Folder:       "$ANALYSIS_FOLDER
echo "# VALIDATION Folder:     "$VALIDATION_FOLDER
echo "# "


# TEST RUNS validity
for RUN in $RUNS;
do
	if [ ! -d "$MISEQ_FOLDER/$RUN" ]; then
		echo "#[WARNING] RUN '$RUN' doesn't exist in the MiSeq Folder '$MISEQ_FOLDER'"
		#exit 0;
	fi;
done;



echo "# VALIDATION      "$VALIDATION
echo "# "


# Preparing VCF files

# sorting

$JAVA -Xmx1g -jar $IGVTOOLS sort $DATA $DATA.sorted 1>/dev/null 2>/dev/null 
$JAVA -Xmx1g -jar $IGVTOOLS sort $GOLDSTANDARD $GOLDSTANDARD.sorted 1>/dev/null 2>/dev/null 


#export VCFDBSNP=$BDFOLDER/snp138.vcf.gz
#export VCFDBSNP137VCF=$BDFOLDER/dbsnp_137.hg19.vcf
#export VCF1000G=$BDFOLDER/1000G_phase1.indels.hg19.vcf
#export VCFMILLS1000G=$BDFOLDER/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
#export COSMIC=$DBFOLDER/cosmic.v54.hg19.vcf

if [ 0 == 1 ]; then

	$JAVA -jar $GATK \
		-T VariantEval \
		-R $GENOME \
		--eval:data $DATA.sorted \
		--eval:goldstandard $GOLDSTANDARD.sorted \
		--goldStandard $GOLDSTANDARD.sorted  \
		--comp $GOLDSTANDARD.sorted \
		--known_names goldstandard \
		-o $VALIDATION.VariantEval.txt \
		1>/dev/null 2>/dev/null 

	#--requireStrictAlleleMatch \
	# GOOD
	#	--eval:VCFDBSNP $VCFDBSNP \
	#	--eval:VCFMILLS1000G $VCFMILLS1000G \

	# BAD
	#	--eval:COSMIC $COSMIC \
	#	--eval:VCF1000G $VCF1000G \

	grep 'ValidationReport  CompRod  EvalRod       JexlExpression  Novelty\|ValidationReport  comp     data          none            all' $VALIDATION.VariantEval.txt | awk '{print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"}' | column -t
	
fi;


if [ 1 == 1 ]; then

	$JAVA -jar $GATK \
		-T GenotypeConcordance \
		-R $GENOME \
		-eval $DATA.sorted \
		-comp $GOLDSTANDARD.sorted \
		-o $VALIDATION.GenotypeConcordance.txt \
		--ignoreFilters \
		--printInterestingSites $VALIDATION.InterestingSites.txt \
		1>/dev/null 2>/dev/null 

	SENSITIVITY=$(grep "^ALL " $VALIDATION.GenotypeConcordance.txt  | perl -e 'print reverse <>' | head -n 1 | awk '{print $2}')
	echo "# SENSITIVITY="$SENSITIVITY

	#echo "Sample             Non-Reference Sensitivity  Non-Reference Discrepancy  Overall_Genotype_Concordance"
	#grep "ALL      " $VALIDATION.GenotypeConcordance.txt  | perl -e 'print reverse <>' | head -n 1

	#echo "InterestingSites: $VALIDATION.InterestingSites.txt"

fi;


rm -f $DATA.sorted $GOLDSTANDARD.sorted 

exit 0;



