#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARK_DEJAVU"
SCRIPT_DESCRIPTION="STARK DEJAVU ANNOVAR databases generation"
SCRIPT_RELEASE="0.9.2b"
SCRIPT_DATE="02/11/2017"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-05/09/2017: Script creation\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.1b-07/09/2017: Add generation of ANNOVAR generic database.\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.2b-02/11/2018: Use BCFTOOLS instead of VCFTOOLS.\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Configuration
ENV_CONFIG=$(find -L $SCRIPT_DIR/.. -name config.app)

source $ENV_CONFIG 1>/dev/null 2>/dev/null

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
	echo "# USAGE: $(basename $0) [options...]";
	#echo "# -u/--update                 Update annovar option";
	echo "# -h/--help                    HELP option";
	echo "#";
}

# header
header;


####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "e:r:uh" --long "update,app_folder:,application_folder:,repo_folder:,repo_folder:,tmp:,bcftools:,tabix:,bgzip:,annovar:,help" -- "$@" 2> /dev/null)

eval set -- "$ARGS"
while true
do
	case "$1" in
		-h|--help)
			usage
			exit 0
			;;
		-u|--update)
			UPDATE=1;
			shift 1
			;;
		-e|--app_folder|--application_folder)
			APP_FOLDER="$2";
			shift 2
			;;
		-r|--repo_folder|--repo_folder)
			REPO_FOLDER="$2";
			shift 2
			;;
		--tmp)
			TMP_FOLDER_TMP="$2";
			shift 2
			;;
		--bcftools)
			BCFTOOLS="$2";
			shift 2
			;;
		--tabix)
			TABIX="$2";
			shift 2
			;;
		--bgzip)
			BGZIP="$2";
			shift 2
			;;
		--annovar)
			ANNOVAR="$2";
			shift 2
			;;
		--) shift
			break
			;;
		*) 	echo "# Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done

# env
#######
#source $SCRIPT_DIR/env.sh


# FOLDERS
###########

#REP=/media/sf_NGS/DOCKER_RESULTS/REPOSITORY
#REP=/home1/L/Archives

# TEMOINS
###########

# TEMOIN="T_HORZON"
# #TEMOIN="P1335,P1439"
# if [ "$TEMOIN" != "" ]; then
# 	TEMOIN_FILTER="/"$(echo "$TEMOIN" | sed "s/,/.\*\\\|\//g")".*"
# else
# 	TEMOIN_FILTER="*"
# fi;
# echo "# TEMOINS: $TEMOIN"

# List GROUP/PROJECT folders

# VERBOSE
VERBOSE=0

if [ ! -z "$APP_FOLDER" ] && [ -d "$APP_FOLDER" ]; then
    STARK_FOLDER_APPS=$APP_FOLDER
fi;
if [ ! -z "$REPO_FOLDER" ] && [ -d "$REPO_FOLDER" ]; then
    FOLDER_REPOSITORY=$REPO_FOLDER
fi;

if [ -z "$FOLDER_ENV" ]; then
    FOLDER_ENV="."
fi;


echo "STARK_FOLDER_APPS=$STARK_FOLDER_APPS"
echo "APP_FOLDER=$APP_FOLDER"
echo "FOLDER_ENV=$FOLDER_ENV"
echo "FOLDER_REPOSITORY=$FOLDER_REPOSITORY"


#GP_FOLDER_LIST="$REPO_FOLDER/HUSDIAGGEN/TUBE"
GP_FOLDER_LIST=""
if [ -z $GP_FOLDER_LIST ]; then
    for ENV_DEF in $(find -L $STARK_FOLDER_APPS -name '*.app' -type f | sed s#$STARK_FOLDER_APPS/## | sort -f -t'/' -k2.3 -k2.2 -k2.1) $(find -L $STARK_FOLDER_APPS -name '*.plugapp' -type f | sed s#$STARK_FOLDER_APPS/## | sort -f -t'/' -k2.3 -k2.2 -k2.1); do
    #for ENV_DEF in "/home1/data/STARK/config/myapps/ONCO/HUSDIAGGEN.TUBE.app"; do
        #echo $ENV_DEF
        #exit 0
        #APP=$(source_app "$STARK_FOLDER_APPS/$ENV_DEF"  2>/dev/null; echo $APP_NAME)
        APP_FOLDER_REPOSITORY=$(source_app "$ENV_DEF"  2>/dev/null; echo $FOLDER_REPOSITORY)
        APP_GROUP=$(source_app "$ENV_DEF"  2>/dev/null; echo $APP_GROUP)
        APP_PROJECT=$(source_app "$ENV_DEF"  2>/dev/null; echo $APP_PROJECT)
        if [ ! -z "$REPO_FOLDER" ] && [ -d "$REPO_FOLDER" ] && [ -d "$REPO_FOLDER/$APP_GROUP/$APP_PROJECT" ]; then
            APP_FOLDER_REPOSITORY=$REPO_FOLDER
        fi;
        if [ -z "$APP_GROUP" ]; then
            APP_GROUP="UNKNOWN"
        fi;
        if [ -z "$APP_PROJECT" ]; then
            APP_PROJECT="UNKNOWN"
        fi;
        echo "$APP_FOLDER_REPOSITORY $APP_GROUP $APP_PROJECT"
        GP_FOLDER_LIST="$GP_FOLDER_LIST\n$APP_FOLDER_REPOSITORY/$APP_GROUP/$APP_PROJECT"
    done;
fi;
GP_FOLDER_LIST_UNIQ=$(echo -e $GP_FOLDER_LIST | sort -u)

echo $GP_FOLDER_LIST_UNIQ



if ((0)); then
    GP_FOLDER_LIST=""
    for ENV in $STARK/env*.sh; do
	    #echo $ENV
	    source $ENV
	    GP_FOLDER_LIST="$GP_FOLDER_LIST\n$FOLDER_REPOSITORY/$GROUP/$PROJECT"
	    #echo $FOLDER_REPOSITORY/$GROUP/$PROJECT
    done
    GP_FOLDER_LIST_UNIQ=$(echo -e $GP_FOLDER_LIST | sort -u)
fi;
#GP_FOLDER_LIST_UNIQ="/home1/L/Archives/HUSHEMATO/TSOMYELOID"
#GP_FOLDER_LIST_UNIQ="/home1/L/Archives/CPSGEN/HCSOP"


# DEJAVU
if [ "$FOLDER_ENV" == "" ]; then FOLDER_ENV=.; fi;
DEJAVU=$FOLDER_ENV/dejavu
RELEASE=$(date +%Y%m%d)
mkdir -p $DEJAVU/$RELEASE
echo "# RELEASE: $RELEASE"

# ASSEMBLY PREFIX
ASSEMBLY_PREFIX_DEFAULT="hg19"

# INFO SUFFIX (e.g. release): $PREFIX_dejavu.$GROUP.$PROJECT$SUFFIX.txt
#SUFFIX=".$RELEASE"
SUFFIX=""

# TMP
if [ "$TMP_FOLDER_TMP" == "" ]; then TMP_FOLDER_TMP=/tmp; fi;
TMP=$TMP_FOLDER_TMP/dejavu_$RANDOM$RANDOM$RANDOM$RANDOM
mkdir -p $TMP
echo "# TMP: $TMP"

# MK
MK=$DEJAVU/$RELEASE/$RELEASE.mk
> $MK

# LOG
LOG=$DEJAVU/$RELEASE/$RELEASE.log
> $LOG


for GP_FOLDER in $GP_FOLDER_LIST_UNIQ; do

	GROUP=$(basename $(dirname $GP_FOLDER))
	PROJECT=$(basename $GP_FOLDER)
	
    echo "# $GP_FOLDER / $GROUP / $PROJECT "

	#ls -l $DEJAVU/$RELEASE/dejavu.$GROUP.$PROJECT.txt

	if [ ! -s $DEJAVU/$RELEASE/dejavu.$GROUP.$PROJECT.done ]; then

		#echo $GP_FOLDER
		#NB_VCF=$(ls $GP_FOLDER/*/*/*final.vcf 2>/dev/null | grep -v "$TEMOIN_FILTER" | wc -w)
		NB_VCF=$(find $GP_FOLDER/*/*/ -maxdepth 1 -name '*.final.vcf' -a ! -name '*.*-*.final.vcf' 2>/dev/null | wc -l)
		#VCFS=$(find $GP_FOLDER/*/*/ -maxdepth 1 -name '*.final.vcf' -a ! -name '*.*-*.final.vcf')
		#echo "#vcf: "$NB_VCF
		if [ $NB_VCF -gt 0 ]; then
			echo ""
			echo "# NGS REPOSITORY FOLDER '$GROUP/$PROJECT' ($NB_VCF VCF)"
			> $MK.$GROUP.$PROJECT.log
			> $MK.$GROUP.$PROJECT.err

			# TMP folder creation
			mkdir -p $TMP/$GROUP/$PROJECT
			cp -f $(find $GP_FOLDER/*/*/ -maxdepth 1 -name '*.final.vcf' -a ! -name '*.*-*.final.vcf') $TMP/$GROUP/$PROJECT/ 2>/dev/null
			NB_VCF_FOUND=$(ls $TMP/$GROUP/$PROJECT/ | wc -w)
			echo "# $NB_VCF_FOUND copied files (some files may be found multiple times)"
			

			# MK files
			> $MK

			VCFGZ_LIST=""
			#for VCF in $(ls $TMP/$GROUP/$PROJECT/*); do
			for VCF in $TMP/$GROUP/$PROJECT/*; do
				#echo "VCF: "$VCF
				if [ -s $VCF ] && (($(grep ^# -cv $VCF))); then

					echo "$VCF.gz: $VCF
					mkdir $<.sort.
					$BCFTOOLS sort $< -T $<.sort. 2>/dev/null | $BCFTOOLS annotate -x FILTER,QUAL,^INFO/AN,^FORMAT/GT -o \$@ -O z 2>/dev/null;
					$TABIX \$@;
					rm -f $<.sort.*
					

					" >> $MK
					VCFGZ_LIST="$VCFGZ_LIST $VCF.gz"
					#| sed s/ID=PL,Number=G/ID=PL,Number=./gi
				fi;
			done

			#echo "$TMP/$GROUP/$PROJECT/dejavu.vcf: $TMP/$GROUP/$PROJECT/VCF_LIST
			#$BCFTOOLS merge $^ | $BCFTOOLS norm -m -any > \$@
			#$BCFTOOLS merge $(cat $TMP/$GROUP/$PROJECT/VCF_LIST) | $BCFTOOLS norm -m -any > \$@
			#$BCFTOOLS merge \$\$(find . -empty -ls $^) | $BCFTOOLS norm -m -any > \$@
			
			echo "$TMP/$GROUP/$PROJECT/dejavu.vcf: $VCFGZ_LIST
			$BCFTOOLS merge $TMP/$GROUP/$PROJECT/*.vcf.gz | $BCFTOOLS norm -m -any > \$@
			
			" >> $MK

			echo "$TMP/$GROUP/$PROJECT/VCF_LIST: $VCFGZ_LIST
			ls $^ > $TMP/$GROUP/$PROJECT/VCF_LIST

			" >> $MK

			echo "%.vcf.gz: %.vcf
				$BGZIP -c $< > \$@

			" >> $MK

			echo "%.gz.tbi: %.gz
				$TABIX $<

			" >> $MK

			#echo "awk -F'\t' '{AF=\$\$6/($NB_VCF_FOUND*2)} {print \$\$1\"\\t\"\$\$2\"\\t\"\$\$3\"\\t\"\$\$4\"\\t\"AF}'"
			#exit 0

			#echo "$DEJAVU/$RELEASE/dejavu.$GROUP.$PROJECT.txt: $TMP/$GROUP/$PROJECT/dejavu.vcf
			echo "$TMP/$GROUP/$PROJECT/dejavu.percent: $TMP/$GROUP/$PROJECT/dejavu.vcf
				$BCFTOOLS query -f'%CHROM\t%POS\t%REF\t%ALT\t%AN\n' $< | awk -F'\t' '{AF=\$\$5/($NB_VCF_FOUND*2)} {print \$\$1\"\\t\"\$\$2\"\\t\"\$\$3\"\\t\"\$\$4\"\\t\"AF}' > \$@

			" >> $MK

			echo "$TMP/$GROUP/$PROJECT/dejavu.annovar: $TMP/$GROUP/$PROJECT/dejavu.vcf
				perl $ANNOVAR/convert2annovar.pl --format vcf4old  --allallele --outfile \$@.tmp $<
				cat \$@.tmp | cut -f1-5 > \$@
				rm -f \$@.tmp

			" >> $MK

			echo "$TMP/$GROUP/$PROJECT/dejavu.annovar.percent: $TMP/$GROUP/$PROJECT/dejavu.annovar $TMP/$GROUP/$PROJECT/dejavu.percent
				paste $^ | cut -f1-5,10 > \$@

			" >> $MK
			
			echo "$TMP/$GROUP/$PROJECT/dejavu.txt: $TMP/$GROUP/$PROJECT/dejavu.annovar.percent
				perl $ANNOVAR/index_annovar.sh $< --outfile \$@

			" >> $MK

			#cat $MK
			echo "# DEJAVU..."
			make -j $THREADS -f $MK $TMP/$GROUP/$PROJECT/dejavu.txt $TMP/$GROUP/$PROJECT/dejavu.vcf $TMP/$GROUP/$PROJECT/dejavu.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.vcf.gz.tbi 1>>$MK.$GROUP.$PROJECT.log 2>>$MK.$GROUP.$PROJECT.err
			# grep "\*\*\*" $MK.err
			#grep "\*\*\*" $MK.log
			#grep "\*\*\*" $MK.err
			
			
			# echo
			if (($(cat $MK.$GROUP.$PROJECT.log $MK.$GROUP.$PROJECT.err | grep "\*\*\*" -c))); then
				echo "# ERROR File '$DEJAVU/$RELEASE/dejavu.$GROUP.$PROJECT.txt' generation..."
				(($VERBOSE)) && cat $MK.$GROUP.$PROJECT.log $MK.$GROUP.$PROJECT.err | grep "\*\*\*" -B 80
			else
			
				# find assembly prefix
				ASSEMBLY_PREFIX_VCF=$(grep "^##reference=" $TMP/$GROUP/$PROJECT/dejavu.vcf | cut -d= -f2 | xargs basename | sed 's/.fa$//' | sed 's/.fasta$//')
				if [ "$ASSEMBLY_PREFIX_VCF" == "" ]; then
					ASSEMBLY_PREFIX=$ASSEMBLY_PREFIX_DEFAULT
				else 
					ASSEMBLY_PREFIX=$ASSEMBLY_PREFIX_VCF
				fi;
				ASSEMBLY_PREFIX_ANNOVAR=$ASSEMBLY_PREFIX"_"
				
				# SUFFIX
				SUFFIX_ANNOVAR=$SUFFIX
				
				PATTERN=$ASSEMBLY_PREFIX_ANNOVAR"dejavu."$GROUP.$PROJECT$SUFFIX_ANNOVAR
				cp $TMP/$GROUP/$PROJECT/dejavu.txt $DEJAVU/$RELEASE/$PATTERN.txt
				cp $TMP/$GROUP/$PROJECT/dejavu.txt.idx $DEJAVU/$RELEASE/$PATTERN.txt.idx
				cp $TMP/$GROUP/$PROJECT/dejavu.vcf $DEJAVU/$RELEASE/$PATTERN.vcf
				cp $TMP/$GROUP/$PROJECT/dejavu.vcf.gz $DEJAVU/$RELEASE/$PATTERN.vcf.gz
				cp $TMP/$GROUP/$PROJECT/dejavu.vcf.gz.tbi $DEJAVU/$RELEASE/$PATTERN.vcf.gz.tbi
				NB_SAMPLES=$(echo $(grep "^#CHROM" $DEJAVU/$RELEASE/$PATTERN.vcf | wc -w)" - 9" | bc)
				NB_VARIANTS=$(grep -cv "^#" $DEJAVU/$RELEASE/$PATTERN.vcf)
				echo "# File '$DEJAVU/$RELEASE/$PATTERN.txt' generated with $NB_SAMPLES samples for $NB_VARIANTS variants"
				
				 for p in 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do awk -F"\t" -v p=$p '$6>p{SUM++} {TOT++} END {PERC=((SUM+0)/TOT)*100; print "# "SUM+0" variants out of "TOT" ("PERC"%) found more than "(p*100)"% in the set"}' $DEJAVU/$RELEASE/$PATTERN.txt; done;
				
				if ((0)); then
					awk -F"\t" '$6>0.1{SUM++} {TOT++} END {PERC=(SUM/TOT)*100; print "# "SUM" variants out of "TOT" ("PERC"%) found more than 10% "}' $DEJAVU/$RELEASE/$ASSEMBLY_PREFIX_ANNOVARdejavu.$GROUP.$PROJECT$SUFFIX_ANNOVAR.txt
					awk -F"\t" '$6>0.05{SUM++} {TOT++} END {PERC=(SUM/TOT)*100; print "# "SUM" variants out of "TOT" ("PERC"%) found more than 5% "}' $DEJAVU/$RELEASE/$ASSEMBLY_PREFIX_ANNOVARdejavu.$GROUP.$PROJECT$SUFFIX_ANNOVAR.txt
					awk -F"\t" '$6>0.01{SUM++} {TOT++} END {PERC=(SUM/TOT)*100; print "# "SUM" variants out of "TOT" ("PERC"%) found more than 1% "}' $DEJAVU/$RELEASE/$ASSEMBLY_PREFIX_ANNOVARdejavu.$GROUP.$PROJECT$SUFFIX_ANNOVAR.txt
				fi;
				
				echo "[INFO] DEJAVU ANNOVAR database generated for '$GROUP/$PROJECT' release $RELEASE" > $DEJAVU/$RELEASE/dejavu.$GROUP.$PROJECT.done
				
				# CLeaning
				rm -rf $TMP/$GROUP/$PROJECT

				# UPDATE ANNOVAR DATABASE ???
				#cp dejavu.CPSGEN.HCSOP.txt /home1/TOOLS/annovar_sources/hg19_dejavu.CPSGEN.HCSOP.txt

				# LOG
				#echo "# Database $GROUP/$PROJECT [release ''$RELEASE', file 'dejavu.$GROUP.$PROJECT.txt'] generated with $NB_SAMPLES samples for $NB_VARIANTS variants" >> $LOG

				#if (($UPDATE)); then
				#		echo "# Update annovar database"
				#fi;


			fi;
			#exit 0;
		fi;
	fi;
done

rm -rf $TMP

exit 0;
