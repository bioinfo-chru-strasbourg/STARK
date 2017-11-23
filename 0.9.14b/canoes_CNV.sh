#!/bin/bash
#################################
##
## CANOES
##
#################################

####################################################################################################################################
# Define the function to print the usage of the script
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
			sh canoes_CNV.sh -i sample_name -o output_path -s SampleSheet -b bed -r R -t BEDTOOLS -u GATK -g hg19 -c bams_list -n blank -a AnnotSV[-h]

		Description:
			CANOES is an algorithm for the detection of rare copy number variants from exome sequencing data.

		Options:
			-i, --sample_name	SAMPLE Name
			-o, --output_path	The output path
			-s, --SampleSheet	SampleSheet path
			-b, --bed 			bed use to analyse the sample
			-r, --r 			Path to R
			-t, --BEDTOOLS 		Path to bedtools
			-u, --GATK			Path to GATK : GenomeAnlysisTK.jar
			-g, --hg19			Reference genome
			-c, --bams_list		file with all bam path
			-n, --blank 		list of blank sample
			-a, --AnnotCNV       Path to ANNOTCNV tools
		 	-h, --help			Print this message and exit the program.

		__EOF__
}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
# ARGS=$(getopt -o "i:o:s::b:r:t:u:g:c:n:a:h" --long "sample_name:,output_path:,SampleSheet:,bed:,r:,bedtools:,gatk:,hg19:,bams_list:,blank:,AnnotSV:,help" -- "$@" 2> /dev/null)
ARGS=$(getopt -o "i:o:s:b:r:t:u:g:c:a:h" --long "sample_name:,output_path:,SampleSheet:,bed:,r:,bedtools:,gatk:,hg19:,bams_list:,AnnotCNV:,help" -- "$@" 2> /dev/null)

[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && \
	exit 1
eval set -- "$ARGS"
while true	
do
	case "$1" in
		-i|--sample_name)
			SAMPLE_NAME="$2"
			shift 2 
			;;
		-o|--output_path)
			OUTPUT="$2"
			shift 2 
			;;
		-s|--SampleSheet)
			SAMPLESHEET="$2"
			shift 2 
			;;
		-b|--bed)
			BED="$2"
			shift 2 
			;;
		-r|--r)
			R="$2"
			shift 2 
			;;
		-t|--bedtools)
			BEDTOOLS="$2"
			shift 2 
			;;
		-u|--gatk)
			GATK="$2"
			shift 2 
			;;
		-g|--hg19)
			HG19_GENOME="$2"
			shift 2 
			;;
		-c|--bams_list)
			BAM_LIST="$2"
			shift 2 
			;;
		# -n|--blank)
		# 	BLANK="$2"
		# 	shift 2 
		# 	;;
		-a|--AnnotCNV)
			ANNOTCNV="$2"
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

canoe_log=$OUTPUT/$SAMPLE_NAME.canoes.log
touch $canoe_log

exec 3>&1
{

[ -f $OUTPUT/Final_CANOES_Anlaysis_$SAMPLE_NAME.tsv ] || touch $OUTPUT/Final_CANOES_Anlaysis_$SAMPLE_NAME.tsv
Final_output=$OUTPUT/Final_CANOES_Anlaysis_$SAMPLE_NAME.tsv

echo -e "### CNV CANOE ANALYSIS RESULTS\n" | tee -a $Final_output

echo "SAMPLE_NAME : $SAMPLE_NAME" | tee -a $Final_output

### detection of samplesheet
### if we have no samplesheet CANOE analysis all patients witout sex chromosome
if [ ! -s "$SAMPLESHEET" ]; then
   	echo -e "\n##################  WARNING  ###########################" | tee -a $Final_output
   	echo -e "############### No SampleSheet #########################" | tee -a $Final_output
   	echo -e "##### CANOES analysis all samples without Sexuel Chromosome ########" | tee -a $Final_output
   	echo -e "########################################################" | tee -a $Final_output
   	echo -e "\n" | tee -a $Final_output
    sample_sheet_info=0
else {
 	sample_sheet_info=1
 	echo -e "\n##################  With SampleSheet  ###########################" | tee -a $Final_output
}
fi;

### check if we have blank information and created a list for BLANK_space
### default blank  
BLANK="BlcADN,blanc,BlcPCR,blcPCR,T_NTC"
BLANK_space=$(echo $BLANK | tr "," " ")



 ### TEST
 echo -e "SAMPLE_NAME=$SAMPLE_NAME\n"
 echo -e "OUTPUT=$OUTPUT\n"
 echo -e "BED=$BED\n"
 echo -e "SAMPLESHEET=$SAMPLESHEET\n"
 echo -e "R=$R\n"
 echo -e "BEDTOOLS=$BEDTOOLS\n"
 echo -e "GATK=$GATK\n"
 echo -e "HG19_GENOME=$HG19_GENOME\n"
 echo -e "BAM_LIST=$BAM_LIST\n"
 echo -e "BLANK=$BLANK\n"
 echo -e "ANNOTCNV=$ANNOTCNV\n"

 ### 

#### created a correct bed for canoes with 3 colums : chr start end 
[ -f $OUTPUT/$SAMPLE_NAME.canoes.bed ] || awk '{print $1"\t"$2"\t"$3}' $BED >> $OUTPUT/$SAMPLE_NAME.canoes.bed
cp $BAM_LIST $OUTPUT/bamlist.txt

SEX=0
if [ $sample_sheet_info == 1 ];
then {
 		########## created a bamlist_sample.txt corresponding for each specifique manifest and sex(bed)
 		SAMPLE_SECTION_FIRST_LINE=$(grep "^Sample_ID," $SAMPLESHEET -n | awk -F: '{print $1}')
 		Correct_line_sample=$(echo $(($SAMPLE_SECTION_FIRST_LINE+1)))
 		tail -n+$Correct_line_sample $SAMPLESHEET >> $OUTPUT/tmp1_$SAMPLE_NAME.txt
 		SAMPLE_SAMPLE_SHEET_MISEQ=$(grep -i ^Sample_ID $SAMPLESHEET | tr -d '\r\n' | sed -e 's/,/\n/g' | grep -i -n Sample_ID | cut -d \: -f 1)
 		MANIFEST_SAMPLE_SHEET_MISEQ=$(grep -i ^Sample_ID $SAMPLESHEET | tr -d '\r\n' | sed -e 's/,/\n/g' | grep -i -n Manifest | cut -d \: -f 1)
 		SEX_SAMPLE_SHEET_MISEQ=$(grep -i ^Sample_ID $SAMPLESHEET | tr -d '\r\n' | sed -e 's/,/\n/g' | grep -i -n Description | cut -d \: -f 1)
	
 		## to maintain the correct order of manifest and sex
 		cat $OUTPUT/tmp1_$SAMPLE_NAME.txt| cut -d, -f$SAMPLE_SAMPLE_SHEET_MISEQ >> $OUTPUT/SAMPLE_SAMPLE_SHEET_MISEQ.txt
 		cat $OUTPUT/tmp1_$SAMPLE_NAME.txt| cut -d, -f$MANIFEST_SAMPLE_SHEET_MISEQ >> $OUTPUT/MANIFEST_SAMPLE_SHEET_MISEQ.txt
 		cat $OUTPUT/tmp1_$SAMPLE_NAME.txt| cut -d, -f$SEX_SAMPLE_SHEET_MISEQ | sed -e 's/-/\t/g' >> $OUTPUT/SEX_SAMPLE_SHEET_MISEQ.txt
 		paste $OUTPUT/SAMPLE_SAMPLE_SHEET_MISEQ.txt $OUTPUT/MANIFEST_SAMPLE_SHEET_MISEQ.txt $OUTPUT/SEX_SAMPLE_SHEET_MISEQ.txt >> $OUTPUT/tmp2_$SAMPLE_NAME.txt

 	while read sample manifest sex
	do
		if [ $sample == $SAMPLE_NAME ];
		then {
			sample_manifest_correct=$manifest
			sample_correct_sex=$sex
			echo "sample_manifest_correct=$sample_manifest_correct"
			echo "sample_correct_sex=$sample_correct_sex"
		}
		fi;
	done < $OUTPUT/tmp2_$SAMPLE_NAME.txt

	### if we have no sex information
	SEX=1
	if [ -z "$sample_correct_sex" ]; then
    	echo -e "##################  WARNING  ###########################" | tee -a $Final_output
    	echo -e "##### No sex information for this patient ##############" | tee -a $Final_output
    	echo -e "##### CANOES analysis without Sexuel Chromosome ########" | tee -a $Final_output
    	echo -e "########################################################" | tee -a $Final_output
    	echo -e "\nManifest for sample $SAMPLE_NAME is : $sample_manifest_correct" | tee -a $Final_output
    	echo -e "\n" | tee -a $Final_output
    	SEX=0
    else {
 	echo -e "\nManifest for sample $SAMPLE_NAME is : $sample_manifest_correct" | tee -a $Final_output
	echo -e "Sex for sample $SAMPLE_NAME is : $sample_correct_sex\n" | tee -a $Final_output
    }
	fi


 		if [ $SEX == 1 ];
 			#### if we have sex information
 		then {
 			while read sample manifest sex
 			do
 				if [ "$manifest" == "$sample_manifest_correct" ] && [ "$sex" == "$sample_correct_sex" ];
 					then {
						sample_same_manifest_sex=$sample
 						grep -w $sample_same_manifest_sex $OUTPUT/bamlist.txt >> $OUTPUT/manifest_specific_bamlist.txt
 					}
 					else {
 						sample_not_manifest_sex=$sample
 						echo "sample remove from CANOE analysis because sex are different : $sample_not_manifest_sex" | tee -a $Final_output
 						# grep -w $sample_not_manifest_sex $OUTPUT/bamlist.txt >> $OUTPUT/no_use_bamlist.txt
 					}
 				fi;
 			done < $OUTPUT/tmp2_$SAMPLE_NAME.txt
 		}
 		else {
 			### if we have no sex information
 			while read sample manifest
 			do
 				if [ "$manifest" == "$sample_manifest_correct" ];
 					then {
 						sample_same_manifest=$sample
 						cat $OUTPUT/bamlist.txt >> $OUTPUT/manifest_specific_bamlist.txt
 					}
 				fi;
 			done < $OUTPUT/tmp2_$SAMPLE_NAME.txt
 		}
 		fi;
 	echo -e "\n\n"
 	CANOE_BAM_LIST=$OUTPUT/manifest_specific_bamlist.txt
}
else {
	CANOE_BAM_LIST=$OUTPUT/bamlist.txt
}
fi;


#### remove POOL samples and Blank samples from CANOE_BAM_LIST 
while read line_canoe_bam_list
do 
	cpt=0
	if [ $(grep POOL <<<$line_canoe_bam_list) ]; then { patient_POOL=1; }; else { patient_POOL=0; }; fi;
 	if [ "$patient_POOL" == 1 ];
 		then {
 				echo "POOL samples are remove from CANOE analysis" | tee -a $Final_output
 				cpt=1;
 		}
 	fi;

	for blank in $BLANK_space
	do
		if [ $(grep $blank <<<$line_canoe_bam_list) ]; then { sample_blank=1; }; else { sample_blank=0; }; fi;
		if [ "$sample_blank" == 1 ];
		then {
			echo "Blank samples are remove from CANOE analysis : $blank" | tee -a $Final_output
			cpt=1;
		}
		fi;
	done;

	if [ "$cpt" == 0 ];
	then {
	grep -w $line_canoe_bam_list $CANOE_BAM_LIST >> $OUTPUT/CANOE_BAM_LIST.txt
	}
	fi; 	
done < $CANOE_BAM_LIST



### read counts with bedtools
$BEDTOOLS/bedtools multicov -bams $(cat $OUTPUT/CANOE_BAM_LIST.txt) -bed $OUTPUT/$SAMPLE_NAME.canoes.bed -q 20 > $OUTPUT/canoes.reads_$SAMPLE_NAME.txt
 	
#### calculate the GC content for each exome capture region
java -Xmx2000m -jar $GATK -T GCContentByInterval -L $OUTPUT/$SAMPLE_NAME.canoes.bed -R $HG19_GENOME -o $OUTPUT/gc.$SAMPLE_NAME.txt
	

 #### with sexual chromosome for canoe analysis
 if [ $SEX == 1 ];
 	then {
 		sed -e "s/^chrX/chr23/" $OUTPUT/canoes.reads_$SAMPLE_NAME.txt | sed -e "s/^chrY/chr24/" > $OUTPUT/newCount.reads_$SAMPLE_NAME.txt
 		sed -e "s/^chrX/chr23/" $OUTPUT/gc.$SAMPLE_NAME.txt | sed -e "s/^chrY/chr24/" > $OUTPUT/Newgc.$SAMPLE_NAME.txt
 	}
 fi;

samplenames_R=$(cat $OUTPUT/CANOE_BAM_LIST.txt | rev | cut -d'/' -f 1 | rev | cut -d"." -f1 | sed 's/^/\"/g' | sed 's/$/\"/g' | sed ':a;N;$!ba;s/\n/,/g')


#### created the R script for CANOES analysis
echo -e "\nSamples use for CANOES Analysis (same Manifest and Sex) : $samplenames_R" | tee -a $Final_output
echo -e "\nCANOES recommended at least 30 samples" | tee -a $Final_output

### create the script.R
echo "setwd(\"$OUTPUT\")" >> $OUTPUT/script_create_$SAMPLE_NAME.R
if [ $SEX == 1 ];
then {
	echo "gc <- read.table(\"Newgc.$SAMPLE_NAME.txt\")\$V2" >> $OUTPUT/script_create_$SAMPLE_NAME.R
	echo "canoes.reads <- read.table(\"newCount.reads_$SAMPLE_NAME.txt\")" >> $OUTPUT/script_create_$SAMPLE_NAME.R
}
else {
	echo "gc <- read.table(\"gc.$SAMPLE_NAME.txt\")\$V2" >> $OUTPUT/script_create_$SAMPLE_NAME.R
	echo "canoes.reads <- read.table(\"canoes.reads_$SAMPLE_NAME.txt\")" >> $OUTPUT/script_create_$SAMPLE_NAME.R
}
fi;
echo "sample.names <- c($samplenames_R)" >> $OUTPUT/script_create_$SAMPLE_NAME.R
echo "names(canoes.reads) <- c(\"chromosome\", \"start\", \"end\", sample.names)" >> $OUTPUT/script_create_$SAMPLE_NAME.R
echo "target <- seq(1, nrow(canoes.reads))" >> $OUTPUT/script_create_$SAMPLE_NAME.R
echo "canoes.reads <- cbind(target, gc, canoes.reads)" >> $OUTPUT/script_create_$SAMPLE_NAME.R
if [ $SEX == 1 ];
	then {
 		echo "source(\"$CANOE_DIR/CANOES_with_SEX.R\")" >> $OUTPUT/script_create_$SAMPLE_NAME.R
	}
 else {
 	echo "source(\"$CANOE_DIR/CANOES_without_SEX.R\")" >> $OUTPUT/script_create_$SAMPLE_NAME.R
 }
 fi;
 echo "xcnv.list <- vector('list', length(sample.names))" >> $OUTPUT/script_create_$SAMPLE_NAME.R
 echo "for (i in 1:length(sample.names)){" >> $OUTPUT/script_create_$SAMPLE_NAME.R
 echo "xcnv.list[[i]] <- CallCNVs(sample.names[i], canoes.reads)" >> $OUTPUT/script_create_$SAMPLE_NAME.R
 echo "}" >> $OUTPUT/script_create_$SAMPLE_NAME.R
 echo "xcnvs <- do.call('rbind', xcnv.list)" >> $OUTPUT/script_create_$SAMPLE_NAME.R
 echo "write.csv(xcnvs, \"results_$SAMPLE_NAME.csv\", row.names=FALSE)" >> $OUTPUT/script_create_$SAMPLE_NAME.R
 
 cd $OUTPUT/
 $R CMD BATCH $OUTPUT/script_create_$SAMPLE_NAME.R # || exit 1 >&1 ; 

entete=$(echo "$samplenames_R" | sed 's/,/\t/g')
sed -i "1iSamples\t\t$entete" $OUTPUT/canoes.reads_$SAMPLE_NAME.txt


if [ -f  $OUTPUT/results_$SAMPLE_NAME.csv ];
	then {
		cat $OUTPUT/results_$SAMPLE_NAME.csv | sed 's/"//g' | sed 's/,/\t/g' | sed '1d' >> $OUTPUT/results_CNV_all_Sample.txt
    	
    	### AnnotCNV annotation for canoe results ###
    	## transforme to AnnotCNV format 
		while IFS=, read ligne
		do
			set $(echo $ligne)
			interval=$3
			chrom=$(echo $interval | cut -d: -f1)
			if [ $chrom == "23" ] ;
				then {
					chrom="X"
				}
			elif [ $chrom == 24 ]
				then {
					chrom="Y"
				}
			fi ;
			start=$(echo $interval | cut -d: -f2 | cut -d- -f1)
			end=$(echo $interval | cut -d: -f2 | cut -d- -f2)
			CNV=$2
			patient=$1
			echo -e "$chrom\t$start\t$end\t$CNV\t$patient" >> $OUTPUT/CNV_Results_all_Sample.bed				
		done < $OUTPUT/results_CNV_all_Sample.txt

		if [ -f  $OUTPUT/CNV_Results_all_Sample.bed ];
			then {
    			################################
    			### AnnotSV annotation for caoes results ###
				### cd $OUTPUT
				source /home1/TOOLS/tools/stark/dev/env.sh	
				$ANNOTCNV/Sources/AnnotCNV-main.tcl -bedFile $OUTPUT/CNV_Results_all_Sample.bed
			}
		fi;


		### seperated results by sample in Final file
		if [ -f  $OUTPUT/CNV_Results_all_Sample.annotated.tsv ];
			then {
				while IFS=, read ligne
				do
					set $(echo $ligne)
					sample_CNV=$(eval echo $5)

					if [ $sample_CNV == "$SAMPLE_NAME" ];
						then {
							echo -e "$ligne" >> $OUTPUT/CNV_Results_$SAMPLE_NAME.annotated.tsv
						}
					fi;
				done < $OUTPUT/CNV_Results_all_Sample.annotated.tsv
			}
		fi;

		if [ -f  $OUTPUT/CNV_Results_$SAMPLE_NAME.annotated.tsv ];
			then {
				### copy AnnotSV result on Final file
				echo -e "\n\n" >> $Final_output
				echo -e "##########################" >> $Final_output
				echo -e "##########################" >> $Final_output
				echo -e "##### CNV Annotation #####" >> $Final_output
				echo -e "##########################" >> $Final_output
				echo -e "##########################" >> $Final_output
				echo -e "\n" >> $Final_output
				head -1 CNV_Results_all_Sample.annotated.tsv >> $Final_output
				#cat $Final_output $OUTPUT/CNV_Results_$SAMPLE_NAME.annotated.tsv >> $Final_output
				cat $OUTPUT/CNV_Results_$SAMPLE_NAME.annotated.tsv >> $Final_output
				clean=YES
			}
			else {
				echo -e "\n\n" >> $Final_output
				echo -e "#########################################" >> $Final_output
				echo -e "#########################################" >> $Final_output
				echo -e "##### No CNV found for this SAMPLE ######" >> $Final_output
				echo -e "#########################################" >> $Final_output
				echo -e "#########################################" >> $Final_output
				echo -e "\n" >> $Final_output	
				clean=YES	
			}
		fi;
	}
	else {
		echo -e "\n\n" >> $Final_output
		echo -e "####################################" >> $Final_output
		echo -e "##### Error on CANOES Analysis #####" >> $Final_output
		echo -e "####################################" >> $Final_output
		grep "Error" $OUTPUT/script_create_$SAMPLE_NAME.Rout -A 1 >> $Final_output
		clean=NO
	}
fi;

if [ "$clean" == "YES" ] ;
	then {
		#### clean	
		rm -rf $OUTPUT/bamlist.txt
		rm -rf $OUTPUT/manifest_specific_bamlist.txt
		rm -rf $OUTPUT/SAMPLE_SAMPLE_SHEET_MISEQ.txt
		rm -rf $OUTPUT/$SAMPLE_NAME.canoes.bed
		rm -rf $OUTPUT/Newgc.$SAMPLE_NAME.txt
		rm -rf $OUTPUT/script_create_$SAMPLE_NAME.R
		rm -rf $OUTPUT/entete.txt
		rm -rf $OUTPUT/MANIFEST_SAMPLE_SHEET_MISEQ.txt
		#rm -rf $OUTPUT/script_create_$SAMPLE_NAME.Rout
		rm -rf $OUTPUT/tmp1_$SAMPLE_NAME.txt
		rm -rf $OUTPUT/results_CNV_all_Sample.txt
		rm -rf $OUTPUT/SEX_SAMPLE_SHEET_MISEQ.txt
		rm -rf $OUTPUT/tmp2_$SAMPLE_NAME.txt
		rm -rf $OUTPUT/newCount.reads_$SAMPLE_NAME.txt
		rm -rf $OUTPUT/results_$SAMPLE_NAME.csv
		rm -rf $OUTPUT/CNV_Results_all_Sample.bed
		rm -rf $OUTPUT/CNV_Results_$SAMPLE_NAME.annotated.tsv
	}
fi ;

echo -e "\n#### END of CNV CANOE Analysis" | tee -a $Final_output
echo -e "\n\n#####################################################\n\n"
	
} 1>$canoe_log 2>&1 

if [ ! -f $OUTPUT/script_create_$SAMPLE_NAME.Rout ] ;
then {
	exit 1;
}
else {
	### to have the canoe results on log file
	cat $OUTPUT/script_create_$SAMPLE_NAME.Rout >> $canoe_log
	rm -rf $OUTPUT/script_create_$SAMPLE_NAME.Rout
}
fi;
