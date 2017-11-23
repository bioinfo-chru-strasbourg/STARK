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
		    sh launch.validation.reproducibility.sh -v validation -r runs -l validation_folder -p replicat -e env [-h]

		Description:
		    This script defines all the parameters to use and launch runs_analysis.sh

		Options:
		 	-v, --validation
		 	This option is optionnal. Validation release. If not exists, creation of a new Validation release
		 	-r, --runs
		 	/!\ This option is required /!\  List of RUN to launch.
		 	-l, --validation_folder
		 	Validation Folder.
		 	-p, --replicat
		 	Number of analysis replication.
		 	-e, --env
		 	Environnement file (define folders, tools...).
		  	-h, --help
		 	Print this message and exit the program.
		__EOF__
}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "v:r:l:p:e:h" --long "validation:,runs:,validation_folder:,replicat:,env:,help" -- "$@" 2> /dev/null)
[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && \
	exit 1
eval set -- "$ARGS"
while true
do
	case "$1" in
		-v|--validation)
			VALIDATION="$2"
			shift 2 
			;;
		-r|--runs)
			RUNS="$2"
			shift 2 
			;;
		-l|--validation_folder)
			FOLDER="$2"
			shift 2 
			;;
		-p|--replicat)
			REPLICAT="$2"
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
if [ "$RUNS" == "" ]
then
	echo "Option --runs is required. " "Use -h or --help to display the help." && exit 1;
fi
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

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

# Validation FOLDER
if [ -d $FOLDER ] && [ "$FOLDER" != "" ]; then
	VALIDATION_FOLDER=$FOLDER
else
	if [ -d $VALIDATION_FOLDER ] && [ "$VALIDATION_FOLDER" != "" ]; then
		echo "#[WARNING] NO VALIDATION_FOLDER defined. Default VALIDATION_FOLDER '$VALIDATION_FOLDER' used."
	else
		echo "#[ERROR] NO VALIDATION_FOLDER defined. NO Default VALIDATION_FOLDER defined."
		exit 1;
	fi;
fi;

# VALIDATION RELEASE
if [ "$VALIDATION" == "" ]; then
	VALIDATION_RELEASE="V"`date '+%Y%m%d-%H%M%S'`
else
		VALIDATION_RELEASE=$VALIDATION
#elif [ -d $VALIDATION_FOLDER/$VALIDATION ]; then
#	VALIDATION_RELEASE=$VALIDATION
fi;



# Replicat
re='^[0-9]+$'
REPLICAT_DEFAULT=3;
if ! [[ $REPLICAT =~ $re ]] ; then
	echo "#[WARNING] REPLICAT '$REPLICAT' is NOT a number. Default REPLICAT '$REPLICAT_DEFAULT' will be used.";
	REPLICAT=$REPLICAT_DEFAULT
fi;

#if [ "$REPLICAT" =="" ] || [ $REPLICAT gt 0 ]; then
#	echo $REPLICAT;
#fi;


GENOME=$GENOMES/$ASSEMBLY/$ASSEMBLY.fa

DATABASE_INTEGRATION_SCRIPT=""

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
echo "# REPLICAT:              "$REPLICAT
echo "# "


# RUNS
if [ "$RUNS" == "" ]; then
	echo "#[ERROR] NO RUN defined"
	exit 0;
fi;

# TEST RUNS validity
for RUN in $RUNS;
do
	if [ ! -d "$MISEQ_FOLDER/$RUN" ]; then
		echo "#[WARNING] RUN '$RUN' doesn't exist in the MiSeq Folder '$MISEQ_FOLDER'"
		#exit 0;
	fi;
done;



echo "# VALIDATION Release     "$VALIDATION_RELEASE
echo "# "


#if [ ! -d $VALIDATION_FOLDER/$VALIDATION_RELEASE ] || [ ! -s $VALIDATION_FOLDER/$VALIDATION_RELEASE/validation.done ]; then
#	VALIDATION_RELEASE=$VALIDATION

	echo "#"
	echo "# GENERATION of VALIDATION RELEASE '$VALIDATION_RELEASE'"
	echo "#"


	mkdir -p $VALIDATION_FOLDER/$VALIDATION_RELEASE


	#echo "# RUNs '$RUNS' - Release '$RELEASE'"
	#for i in {1..$REPLICAT }
	for (( REP=1; REP<=$REPLICAT; REP++ ))
	do
		echo "# REPLICATION '$REP'"
		REP_FOLDER=$VALIDATION_FOLDER/$VALIDATION_RELEASE/$REP
		REP_FOLDER_DEM=$REP_FOLDER/$DEM
		REP_FOLDER_RES=$REP_FOLDER$RESULTS_SUBFOLDER
		REP_ENV=$REP_FOLDER/env.R$REP.sh

		mkdir -p $REP_FOLDER
		mkdir -p $REP_FOLDER_DEM
		mkdir -p $REP_FOLDER_RES
	
		# TEST
		#MAIN_FOLDER=$REP_FOLDER
		#if [ ! -z $MAIN_FOLDER ]; then
		#	echo "MAIN_FOLDER '$MAIN_FOLDER' exists";
		#else
		#	echo "MAIN_FOLDER DOESNOT '$MAIN_FOLDER' exists";
		#fi;
		#exit 1;

		echo "## STARK environment for VALIDATION Release '$VALIDATION_RELEASE' - REPLICAT '$REP'" >  $REP_ENV
		echo "source $ENV" >>  $REP_ENV


		echo "export DEMULTIPLEXING_FOLDER=$REP_FOLDER_DEM		# DEMULTIPLEXING folder" >>  $REP_ENV
		echo "export RESULTS_FOLDER=$REP_FOLDER_RES			# RESULTS folder ALL" >>  $REP_ENV
		echo "export ANALYSIS_FOLDER=$REP_FOLDER/$ANA			# ANALYSIS folder" >>  $REP_ENV

		touch $REP_ENV

		#continue;
		#exit 1;

		for RUN in $RUNS;
		do
			echo "#"
			#cat $VALIDATION_FOLDER/$VALIDATION_RELEASE/$REP/done.log; exit 0;
			if [ -d "$MISEQ_FOLDER/$RUN" ]; then
				#if [ "`cat $VALIDATION_FOLDER/$VALIDATION_RELEASE/$REP/done.log 2>/dev/null`" != "done" ]; then
				#echo "`cat $REP_FOLDER_RES/$RUN/done.log`";
				#continue
				if [ "`cat $REP_FOLDER_RES/$RUN/done.log 2>/dev/null`" != "done" ]; then
					RELEASE="V"`date '+%Y%m%d-%H%M%S'`
					echo "# RUN '$RUN' - Release '$RELEASE'"
					LOG=$VALIDATION_FOLDER/$VALIDATION_RELEASE/run_analysis.$RELEASE.launch.log
					#echo "#$SCRIPT_DIR/runs_analysis.sh "$MULTI_RUN_ANALYSIS_PROCESS" "$RUN" "$ALIGNERS" "$CALLERS" "$ANNOTATORS" "$FILTER_SAMPLE" "$INTERSEC" "$REP_ENV" "1" "$DATABASE_INTEGRATION_SCRIPT" 1>>$LOG 2>>$LOG"
					#echo "# Just for test, copy files..."
					#cp -fR $DEMULTIPLEXING_FOLDER/$RUN  $REP_FOLDER_DEM 
					#cp -fR $RESULTS_FOLDER/$RUN  $REP_FOLDER_RES
					$SCRIPT_DIR/runs_analysis.sh -m "$MULTI_RUN_ANALYSIS_PROCESS" -r "$RUN" -a "$ALIGNERS" -c "$CALLERS" -n "$ANNOTATORS" -f "$FILTER_SAMPLE" -i "$INTERSEC" -e "$REP_ENV" -u "1" -d "$DATABASE_INTEGRATION_SCRIPT" 1>>$LOG 2>>$LOG
					if [ "`grep '***' $LOG -c `" == "0" ]; then
					#	echo "NO error";
					echo "done" > $REP_FOLDER_RES/$RUN/done.log
					#else 
					#	echo "error";
					fi;


					echo "# DONE @"`date '+%Y%m%d-%H%M%S'`
				else
					echo "# RUN '$RUN' already analyzed"
				fi;
			else
				echo "# RUN '$RUN' doesn't exist in the MiSeq Folder '$MISEQ_FOLDER'"
			fi;
		done;
		echo "#"
	done
	echo "#"

	echo "# DATA Generation DONE @"`date '+%Y%m%d-%H%M%S'` > $VALIDATION_FOLDER/$VALIDATION_RELEASE/validation.done 

#else

#	echo "#"
#	echo "# VALIDATION RELEASE '$VALIDATION_RELEASE' Exists"
#	echo "#"

#fi;


echo "#"
echo "# VALIDATION check"
echo "#"

VALIDATION_FOLDER_LOG=$VALIDATION_FOLDER/$VALIDATION_RELEASE/LOG
mkdir -p $VALIDATION_FOLDER_LOG
VAL_LOG=$VALIDATION_FOLDER_LOG/$VALIDATION_RELEASE.log
VAL_LOG_CP=$VALIDATION_FOLDER/$VALIDATION_RELEASE/$(basename VAL_LOG)

for RUN in $RUNS;
do
	echo "#"
	echo "# RUN '$RUN'"
	echo "#"

	VAL_LOG_RUN=$VALIDATION_FOLDER_LOG/$RUN

	for (( REP=1; REP<=$REPLICAT; REP++ ))
	do

		#echo "# REPLICATION '$REP'"
		REP_FOLDER=$VALIDATION_FOLDER/$VALIDATION_RELEASE/$REP
		REP_FOLDER_DEM=$REP_FOLDER/$DEM
		REP_FOLDER_RES=$REP_FOLDER$RESULTS_SUBFOLDER
		#REP_FOLDER_VAL=$REP_FOLDER/$VAL
		REP_ENV=$REP_FOLDER/env.R$REP.sh

		

		# LOOP ON REP
		REP2=$(( $REP )) #$($REP+1)
		for (( REP2=$(( $REP )); REP2<=$REPLICAT; REP2++ ))
		do
			if [ ! $REP2 -eq $REP ]; then

				REP2_FOLDER=$VALIDATION_FOLDER/$VALIDATION_RELEASE/$REP2
				REP2_FOLDER_DEM=$REP2_FOLDER/$DEM
				REP2_FOLDER_RES=$REP2_FOLDER$RESULTS_SUBFOLDER
				REP2_ENV=$REP2_FOLDER/env.R$REP2.sh


				COMP="$REP-$REP2";
				echo "# COMPARISON $COMP"
				VAL_LOG_RUN_COMP=$VALIDATION_FOLDER_LOG/$RUN.$COMP
				rm $VAL_LOG_RUN_COMP # remove
				#touch $VAL_LOG_RUN_COMP

				if [ 1 == 0 ]; then

				BAMs=$(find $REP_FOLDER_RES/$RUN  -mindepth 2 -maxdepth 2  -type f -name *.bam)
				for BAM in $BAMs
				do
					T=f
					BAM2=$(echo $BAM | sed "s|$REP_FOLDER_RES/$RUN|$REP2_FOLDER_RES/$RUN|gi")
					VAL_LOG_RUN_COMP_BAM=$VAL_LOG_RUN_COMP.$(basename $BAM)

					# DIFF TEST

					diff $BAM $BAM2 > $VAL_LOG_RUN_COMP_BAM
					if [ "$(grep ^ -c $VAL_LOG_RUN_COMP_BAM)" != "0" ]; then
						echo "#[ERROR] File '$(basename $BAM)' error validation"
					fi;
					cat $VAL_LOG_RUN_COMP_BAM >> $VAL_LOG_RUN_COMP

					# OTHER TEST

					# TODO
				done;
				fi;



				VCFs=$(find $REP_FOLDER_RES/$RUN  -mindepth 2 -maxdepth 2  -type f -name *.vcf)
				for VCF in $VCFs
				do
					T=f
					VCF2=$(echo $VCF | sed "s|$REP_FOLDER_RES/$RUN|$REP2_FOLDER_RES/$RUN|gi")
					#VCF2=`echo $VCF | sed -i "s/^$REP_FOLDER_RES\/$RUN/$REP2_FOLDER_RES\/$RUN/gi"`
					#echo "# VCF file $(basename $VCF)"
					#echo "# VCF1: $VCF"
					#echo "# VCF2: $VCF2"
					VAL_LOG_RUN_COMP_VCF=$VAL_LOG_RUN_COMP.$(basename $VCF)
					rm -f $VAL_LOG_RUN_COMP_VCF # remove
					#touch $VAL_LOG_RUN_COMP_VCF

					# DIFF TEST
					if [ 0 == 1 ]; then
					diff $VCF $VCF2 > $VAL_LOG_RUN_COMP_VCF
					if [ "$(grep ^ -c $VAL_LOG_RUN_COMP_VCF)" != "0" ]; then
						echo "#[ERROR] File '$(basename $VCF)' error validation"
					fi;
					cat $VAL_LOG_RUN_COMP_VCF >> $VAL_LOG_RUN_COMP
					#ls $REP_FOLDER_RES/$RUN/15C0015/15C0015.bwamem.gatkHC.vap.vcf
					#ls $REP2_FOLDER_RES/$RUN/15C0015/15C0015.bwamem.gatkHC.vap.vcf
					#diff $REP_FOLDER_RES/$RUN/15C0015/15C0015.bwamem.gatkHC.vap.vcf $REP2_FOLDER_RES/$RUN/15C0015/15C0015.bwamem.gatkHC.vap.vcf
					fi;

					# VARIANT TEST

					# GATK
					if [ 0 == 1 ]; then
					java -jar $GATK \
					      -T CombineVariants \
					      -R $GENOME \
					      -V:VCF1 /media/IRCV2/V2/VAL/VALIDATION_TEST1/1/RES/ALL/150626_M01656_0049_000000000-D0JH7/15C0015/15C0015.bwamem.gatkHC.vap.vcf \
					      -V:VCF2 /media/IRCV2/V2/VAL/VALIDATION_TEST1/2/RES/ALL/150626_M01656_0049_000000000-D0JH7/15C0015/15C0015.bwamem.gatkHC.vap.vcf \
					      -V:VCF3 /media/IRCV2/V2/VAL/VALIDATION_TEST1/3/RES/ALL/150626_M01656_0049_000000000-D0JH7/15C0015/15C0015.bwamem.gatkUG.vap.vcf \
					      --genotypemergeoption UNIQUIFY \
					      -o  /media/IRCV2/V2/VAL/VALIDATION_TEST1/LOG/merge.vcf; \
					java -jar /media/IRCV2/NGSEnv/tools2/gatk/current/bin/GenomeAnalysisTK.jar \
					     -T VariantEval \
					     -R /media/IRCV2/NGSEnv/genomes/hg19/hg19.fa \
					     -select 'set=="Intersection"' -selectName Intersection \
					     -o /media/IRCV2/V2/VAL/VALIDATION_TEST1/LOG/merge.gatkreport \
					     -eval  /media/IRCV2/V2/VAL/VALIDATION_TEST1/LOG/merge.vcf  \
					     -l INFO;
					fi;

					#GREP AWK SED 
					if [ 1 == 1 ]; then
					grep -v ^# $VCF | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$7}' > $VCF.reduced
					grep -v ^# $VCF2 | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$7}' > $VCF2.reduced
					diff $VCF.reduced $VCF2.reduced > $VAL_LOG_RUN_COMP_VCF
					rm -f $VCF.reduced $VCF2.reduced
					if [ "$(grep ^ -c $VAL_LOG_RUN_COMP_VCF)" != "0" ]; then
						echo "#[ERROR] File '$(basename $VCF)' error validation "
					fi;
					cat $VAL_LOG_RUN_COMP_VCF >> $VAL_LOG_RUN_COMP

					fi;

				done;

				

			fi;
		done;
		
		cat $VAL_LOG_RUN_COMP >> $VAL_LOG_RUN
		
	done;
	echo "#"

	cat $VAL_LOG_RUN >> $VAL_LOG

done
echo "#"

cp $VAL_LOG $VAL_LOG_CP


echo "# VALIDATION LOG:"
if [ "$(grep ^ -c $VAL_LOG_RUN_COMP_VCF)" == "0" ]; then
	echo "# NO ERROR in validation"
else
	echo "# ERROR in validation"
	cat $VAL_LOG	
fi;
echo "#"

exit 1;



