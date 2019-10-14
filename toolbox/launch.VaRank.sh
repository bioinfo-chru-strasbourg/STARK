#!/bin/bash
#################################
##
## VARANK ANALYSIS
## need alamut-batch license
## author: Sinthuja Pachchek
##
#################################

####################################################################################################################################
# Define the function to print the usage of the script
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
			sh launch.VaRank.sh -a VARANK_FOLDER -r STARK RUN [-h]

		Options:
			-a, --varank_folder 	OUTPUT_VARANK_FOLDER path
			-r, --run             STARK RUN path
			-h, --help          Print this message and exit the program.

		__EOF__
}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "r:a:h" --long "run,varank_folder,help" -- "$@" 2> /dev/null)


[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && \
	exit 1
eval set -- "$ARGS"
while true  
do
	case "$1" in
		-a|--varank_folder)
			OUTPUT_VARANK_FOLDER="$2"
			shift 2 
			;;
		-r|--run)
			RUN="$2"
			shift 2 
			;;	
		-h|--help)
			usage
			exit 0
			;;
		--) shift
			break 
			;;
		*)  echo "Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### if it's a unitary Varank launch code for a specific varank folder
### mist to have .vcf file to analyses in the varank older
STARK_FOLDER="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $STARK_FOLDER/env.DIAG.sh

if [ ! -z $OUTPUT_VARANK_FOLDER ]; then
	[ -f $OUTPUT_VARANK_FOLDER/VaRank.log ] || touch $OUTPUT_VARANK_FOLDER/VaRank.log
	exec 3>&1
	{
		[ -f $OUTPUT_VARANK_FOLDER/configfile ] || rsync -r $VARANK/configfile $OUTPUT_VARANK_FOLDER/
		[ -f $OUTPUT_VARANK_FOLDER/*tsv ] && rm $OUTPUT_VARANK_FOLDER/*tsv
		unset GZIP; $VARANK/bin/VaRank -vcfdir "$OUTPUT_VARANK_FOLDER" -alamutHumanDB hg19 -SamVa "yes" -AlamutProcesses $THREADS
		sleep 1m
	} 1>$OUTPUT_VARANK_FOLDER/VaRank.log 2>&1 
fi;

#### if it's a RUN analysed by STARK we have samplesheet in the folder
## split samplesheet by sample information
if [ ! -z $RUN ]; then
	SAMPLESHEET=$(ls $RUN/*/*.SampleSheet.csv | head -1)
	echo SAMPLESHEET=$SAMPLESHEET
	i=0;
	while read ligne
	do 	
		### for each samples
		for line in $(sed -n -e '/\[Data\]/,$p' | tail -n +2)
		do 
			for line2 in $(echo $line | grep 'Sample_ID,')
			do
				SAMPLE_ID_num=$(echo $line | grep 'Sample_ID,' | tr ',' '\n' | grep -n Sample_ID | awk -F ':' '{print $1}')
				SAMPLE_PROJECT_num=$(echo $line | grep 'Sample_Project,' | tr ',' '\n' | grep -n Sample_Project | awk -F ':' '{print $1}')
				SAMPLE_DESCRIPTION_num=$(echo $line | grep 'Description,' | tr ',' '\n' | grep -n Description | awk -F ':' '{print $1}')
			done

			SAMPLE_NAME=$(echo $line | awk -F ',' '{ print $'$SAMPLE_ID_num'}')
			SAMPLE_PROJECT=$(echo $line | awk -F ',' '{ print $'$SAMPLE_PROJECT_num'}')
			SAMPLE_DESCRIPTION=$(echo $line | awk -F ',' '{ print $'$SAMPLE_DESCRIPTION_num'}')
			
			if [ $SAMPLE_NAME != "Sample_ID" ]; then
				ENV=$(echo $SAMPLE_PROJECT | awk -F '-' '{ print "env."$1"."$2".sh"}')
				source $STARK_FOLDER/$ENV


				### find POOL
				SAMPLE_DESCRIPTION_POOL=$(echo $SAMPLE_DESCRIPTION | tr '-' '\n' | grep "^POOL" | tr '\n' ' ')
				
				#### if VARANK_FOLDER and configfile doesn't exist
				if [ ! -f "${VARANK_FOLDER}/${GROUP}/${PROJECT}/configfile" ]
				then
					rsync -r $VARANK/configfile ${VARANK_FOLDER}/${GROUP}/${PROJECT}/
				fi;
				CONFIG_FILE=${VARANK_FOLDER}/${GROUP}/${PROJECT}/configfile

				## copy final.vcf.gz to $VARANK_FOLDER
				[ -f $RUN/$SAMPLE_NAME/*.reports/$SAMPLE_NAME.final.vcf.gz ] && rsync $RUN/$SAMPLE_NAME/*.reports/$SAMPLE_NAME.final.vcf.gz $VARANK_FOLDER/$GROUP/$PROJECT/ || rsync $RUN/$SAMPLE_NAME/*.reports/$SAMPLE_NAME.final.vcf $VARANK_FOLDER/$GROUP/$PROJECT/
			
				### add sample in configfil
				last_lign_config=$(tail -1 $CONFIG_FILE | awk '{print $1}')
				if [ $last_lign_config = "#" ]; then
					echo "fam1: $SAMPLE_NAME $SAMPLE_DESCRIPTION_POOL" >> $CONFIG_FILE
				else 
					tail -n+333 $CONFIG_FILE | awk -F " " '{print $0}' > ${VARANK_FOLDER}/${GROUP}/${PROJECT}/config_sample_list.txt
					last_number_config=$(echo $last_lign_config | sed "s/fam//gi" | sed "s/://gi");
					find=0;

					if grep -q $SAMPLE_NAME ${VARANK_FOLDER}/${GROUP}/${PROJECT}/config_sample_list.txt; then
    					find=1;
    					echo "[error 1 ] : sample $SAMPLE_NAME is already present in VaRank configfile"
    					# break;
					fi
				
					rm ${VARANK_FOLDER}/${GROUP}/${PROJECT}/config_sample_list.txt
					if [ $find == 0 ];
						then {
						let "last_number_config++"
						echo "fam$last_number_config: ${SAMPLE_NAME} $SAMPLE_DESCRIPTION_POOL" >> $CONFIG_FILE;
						}
					fi;
				fi;

				#### detected varank folder to launch one time varank for all same project sample
				ARRAY[$i]="$STARK_FOLDER/$ENV"
				let i++	;

			fi;
		done
	done <$SAMPLESHEET

	### launch varank commande for each project
	uniq=($(printf "%s\n" "${ARRAY[@]}" | sort -u)); 
	for i in ${uniq[@]}; do
		source ${i}
		[ -f ${VARANK_FOLDER}/${GROUP}/${PROJECT}/VaRank.log ] || touch ${VARANK_FOLDER}/${GROUP}/${PROJECT}/VaRank.log
		exec 3>&1
		{
			VARANK_FOLDER_GROUP_PROJECT=${VARANK_FOLDER}/${GROUP}/${PROJECT}
			FOLDER_REPOSITORY_VARANK_FOLDER_GROUP_PROJECT=${FOLDER_REPOSITORY}/${GROUP}/${PROJECT}/VARANK
			[ -f $VARANK_FOLDER_GROUP_PROJECT/*tsv ] && rm $VARANK_FOLDER_GROUP_PROJECT/*tsv
			unset GZIP; $VARANK/bin/VaRank -vcfdir "$VARANK_FOLDER_GROUP_PROJECT" -alamutHumanDB hg19 -SamVa "yes" -AlamutProcesses $THREADS
			sleep 1m


			### Copy VaRank results in FOLDER_REPOSITORY
			if [ -n "${FOLDER_REPOSITORY}" ]
			then
				[ -d $FOLDER_REPOSITORY_VARANK_FOLDER_GROUP_PROJECT ] && rm $FOLDER_REPOSITORY_VARANK_FOLDER_GROUP_PROJECT/*tsv
				[ -d $FOLDER_REPOSITORY_VARANK_FOLDER_GROUP_PROJECT ] || mkdir -p $FOLDER_REPOSITORY_VARANK_FOLDER_GROUP_PROJECT
				rsync $VARANK_FOLDER_GROUP_PROJECT/*.tsv $FOLDER_REPOSITORY_VARANK_FOLDER_GROUP_PROJECT/
			else
				echo "FOLDER_REPOSITORY is not specify "
			fi;
		} 1>${VARANK_FOLDER}/${GROUP}/${PROJECT}/VaRank.log 2>&1 
	done

fi; ### end of varank run stark analysis
