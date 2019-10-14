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
			sh canoes_CNV.sh -o output_path -b bed -r R -t BEDTOOLS -u GATK -g hg19 -c bams_list -a AnnotSV -p sample_name [-h]

		Description:
			CANOES is an algorithm for the detection of rare copy number variants from exome sequencing data.

		Options:
			-o, --output_path   The output path
			-p, --sample        The sample
			-s, --SampleSheet   SAMPLESHEET_LIST path
			-l, --ana			part to analyse ( "ALL,SEX_F,SEX_M" )
			-b, --bed           bed use to analyse the sample
			-r, --r             Path to R
			-t, --BEDTOOLS      Path to bedtools
			-u, --GATK          Path to GATK : GenomeAnlysisTK.jar
			-g, --hg19          Reference genome
			-c, --bams_list     file with all bam path (or multicov_list)
			-a, --AnnotCNV      Path to ANNOTCNV tools
			-h, --help          Print this message and exit the program.

		__EOF__
}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
# ARGS=$(getopt -o "i:o:s::b:r:t:u:g:c:n:a:h" --long "output_path:,SampleSheet:,bed:,r:,bedtools:,gatk:,hg19:,bams_list:,blank:,AnnotSV:,help" -- "$@" 2> /dev/null)
ARGS=$(getopt -o "o:p:s:l:b:r:t:u:g:c:a:h" --long "output_path:,sample:,SampleSheet:,ana:,bed:,r:,bedtools:,gatk:,hg19:,bams_list:,AnnotCNV:,help" -- "$@" 2> /dev/null)

[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && \
	exit 1
eval set -- "$ARGS"
while true  
do
	case "$1" in
		-o|--output_path)
			OUTPUT="$2"
			shift 2 
			;;
		-p|--sample)
			SAMPLENAME="$2"
			shift 2 
			;;			
		-s|--SampleSheet)
			SAMPLESHEET_LIST="$2"
			shift 2 
			;;
		-l|--ana)
			list="$2"
			shift 2 
			;;
		-b|--beds_list)
			BED_LIST="$2"
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
		-c|--multicov_list)
			MULTICOV_LIST="$2"
			shift 2 
			;;
		# -n|--blank)
		#   BLANK="$2"
		#   shift 2 
		#   ;;
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
		*)  echo "Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


ANA=$(date +"%Y%m%d")
clean="YES"

exec &> $OUTPUT/CANOES.log

echo -e "### CNV CANOES ANALYSIS RESULTS\n"

### BED
BED_CONTROL=$(cat $BED_LIST | awk '{print NF}' | sort -nu | tail -n 1);
if [ $BED_CONTROL == 1 ]
then
	BED=$(head -1 $BED_LIST);
else
	BED=$BED_LIST;
fi;

### Sample
[ ! -z $SAMPLENAME ] || SAMPLENAME=$(basename $MULTICOV_LIST | cut -d"." -f1)

### informations about paths for Bioinformaticians
echo -e "OUTPUT=$OUTPUT\n"
echo -e "BED=$BED\n"
echo -e "SAMPLESHEET=$SAMPLESHEET_LIST\n"
echo -e "R=$R\n"
echo -e "BEDTOOLS=$BEDTOOLS\n"
echo -e "GATK=$GATK\n"
echo -e "HG19_GENOME=$HG19_GENOME\n"
echo -e "MULTICOV_LIST=$MULTICOV_LIST\n"
echo -e "CANOES=$CANOE\n"
echo -e "ANNOTCNV=$ANNOTCNV\n"
echo -e "SAMPLE=$SAMPLENAME\n"


### created a correct bed for canoes
### format:
### chr start end
[ -f $OUTPUT/$ANA.canoes.bed ] || awk '{print $1"\t"$2"\t"$3}' $BED > $OUTPUT/$ANA.canoes.bed


### copy bamlist in canoes folder
### format:
### FOLDER/SAMPLE.caller.bam
cp $MULTICOV_LIST $OUTPUT/multicovlist.txt
if [[ -s $SAMPLESHEET_LIST ]]; then
	sed -i 's/\r//g' $SAMPLESHEET_LIST
fi;


### Select Manifest associated SampleSheet
### if no SampleSheet, Generate some
if [ -s $SAMPLESHEET_LIST ]
then
	SAMPLESHEET_CONTROL=$(grep "\[Header\]" $SAMPLESHEET_LIST);
	if [ "$SAMPLESHEET_CONTROL" == "[Header]" ]
	then
		SAMPLESHEET=$SAMPLESHEET_LIST;
	else
		SAMPLESHEET=$(head -1 $SAMPLESHEET_LIST);
	fi;
	cat $SAMPLESHEET > $OUTPUT/SampleSheet.csv
	[ -z "$SAMPLENAME" ] || SAMLEMANIFEST=$(cat $SAMPLESHEET | sed -n -e '/\[Data\]/,$p' | tail -n +2 | awk -F ',' 'NR==1 { for (i=1; i<=NF; i++) { f[$i] = i } } { print $(f["Sample_ID"]), $(f["Manifest"]) }' | grep $SAMPLENAME | awk '{print $2}');
	if [ ! -z "$SAMLEMANIFEST" ]
	then
		cat $SAMPLESHEET | grep -A 1 -B 1000 '\[Data\]' > $OUTPUT/SampleSheet.csv;
		cat $SAMPLESHEET | sed -n -e '/\[Data\]/,$p' | tail -n +3 | grep ",$SAMLEMANIFEST," >> $OUTPUT/SampleSheet.csv;
		#grep ",$SAMLEMANIFEST," $SAMPLESHEET >> $OUTPUT/SampleSheet.csv;
	fi;
	if [ ! -s $OUTPUT/SampleSheet.csv ]
	then
		echo "Auto generated SampleSheet" > $OUTPUT/SampleSheet.csv
		echo "[Data]" >> $OUTPUT/SampleSheet.csv
		echo "Sample_ID,Description" >> $OUTPUT/SampleSheet.csv
		while read MULTISAMPLE
		do 
			filename=$(basename $MULTISAMPLE)
			echo "${filename%.*.*}," >> $OUTPUT/SampleSheet.csv;
		done < $OUTPUT/multicovlist.txt
	fi;
else
	echo "Generic SampleSheet" > $OUTPUT/SampleSheet.csv
	echo "[Data]" >> $OUTPUT/SampleSheet.csv
	echo "Sample_ID,Description" >> $OUTPUT/SampleSheet.csv
	while read MULTISAMPLE
	do 
		filename=$(basename $MULTISAMPLE)
		echo "${filename%.*.*}," >> $OUTPUT/SampleSheet.csv;
	done < $OUTPUT/multicovlist.txt
fi;


sed -i '/^$/d' $OUTPUT/SampleSheet.csv;
### remove POOL and BLANCK samples / define a specific order of patients in multicov & canoes reads
### format:
### Sample Sex
#[ -s $OUTPUT/SampleSheet.SampleSex.txt ] || cat $OUTPUT/SampleSheet.csv | grep -vE 'POOL|BlcADN|blanc|BlcPCR|blcPCR|_NTC' | sed -n -e '/\[Data\]/,$p' | tail -n +2 | awk -F ',' 'NR==1 { for (i=1; i<=NF; i++) { f[$i] = i } } { print $(f["Sample_ID"]), $(f["Description"]) }' | sed '1d' > $OUTPUT/SampleSheet.SampleSex.txt
#[ -s $OUTPUT/SampleSheet.SampleSex.Complet.txt ] || cat $OUTPUT/SampleSheet.csv | sed -n -e '/\[Data\]/,$p' | tail -n +2 | awk -F ',' 'NR==1 { for (i=1; i<=NF; i++) { f[$i] = i } } { print $(f["Sample_ID"]), $(f["Description"]) }' | sed '1d' > $OUTPUT/SampleSheet.SampleSex.Complet.txt
SAMP_ID=$(cat $OUTPUT/SampleSheet.csv | grep 'Sample_ID,' | tr ',' '\n' | grep -n Sample_ID | awk -F ':' '{print $1}')
DESC=$(cat $OUTPUT/SampleSheet.csv | grep 'Sample_ID,' | tr ',' '\n' | grep -n Description | awk -F ':' '{print $1}')
[ -s $OUTPUT/SampleSheet.SampleSex.txt ] || cat $OUTPUT/SampleSheet.csv | grep -vE 'POOL|BlcADN|blanc|BlcPCR|blcPCR|_NTC' | sed -n -e '/\[Data\]/,$p' | tail -n +2 | awk -F ',' '{ print $'$SAMP_ID', $'$DESC'}' | sed '1d' > $OUTPUT/SampleSheet.SampleSex.txt
[ -s $OUTPUT/SampleSheet.SampleSex.Complet.txt ] || cat $OUTPUT/SampleSheet.csv | sed -n -e '/\[Data\]/,$p' | tail -n +2 | awk -F ',' '{ print $'$SAMP_ID', $'$DESC'}' | sed '1d' > $OUTPUT/SampleSheet.SampleSex.Complet.txt


### bamlist for CANOES analysis
### format:
### FOLDER/SAMPLE.caller.bam
[ -s $OUTPUT/SampleSheet.ALL.txt ] || cat $OUTPUT/SampleSheet.SampleSex.txt | awk '{ print $1}' > $OUTPUT/SampleSheet.ALL.txt
[ -s $OUTPUT/SampleSheet.SEX_M.txt ] || cat $OUTPUT/SampleSheet.SampleSex.txt | grep SEX_M | awk '{ print $1}' > $OUTPUT/SampleSheet.SEX_M.txt
[ -s $OUTPUT/SampleSheet.SEX_F.txt ] || cat $OUTPUT/SampleSheet.SampleSex.txt | grep SEX_F | awk '{ print $1}' > $OUTPUT/SampleSheet.SEX_F.txt

### %GC
### format:
### chrom:start-end %GC
[ -s $OUTPUT/$ANA.CANOES_GC.tsv ] || java -Xmx2000m -jar $GATK -T GCContentByInterval -L $OUTPUT/$ANA.canoes.bed -R $HG19_GENOME -o $OUTPUT/$ANA.CANOES_GC.tsv;
sed -e "s/^chrX/chr23/" $OUTPUT/$ANA.CANOES_GC.tsv | sed -e "s/^chrY/chr24/" > $OUTPUT/new_$ANA.CANOES_GC.txt;
cat $OUTPUT/new_$ANA.CANOES_GC.txt | grep -vE 'chr23|chr24' > $OUTPUT/$ANA.CANOES_GC.ALL.txt;
[ -s $OUTPUT/SampleSheet.SEX_M.txt ] && cp $OUTPUT/new_$ANA.CANOES_GC.txt $OUTPUT/$ANA.CANOES_GC.SEX_M.txt;
[ -s $OUTPUT/SampleSheet.SEX_F.txt ] && cp $OUTPUT/new_$ANA.CANOES_GC.txt $OUTPUT/$ANA.CANOES_GC.SEX_F.txt;


### merge multicov files from multicov.list
### format:
### chrom start end coverage
while read SAMPLE SEX
do 
	multicovfile=$(grep -F "$SAMPLE." $OUTPUT/multicovlist.txt);
	[ "$multicovfile" == "" ] && echo "Error: No multicov file for sample: $SAMPLE";
	[ "$multicovfile" == "" ] && exit
	EXTENSION=$(head -1 $OUTPUT/multicovlist.txt | awk -F'.' '{print $NF}')
	[ "$EXTENSION" == "bam" ] && $BEDTOOLS/bedtools multicov -bams $multicovfile -bed $BED -q 20 > $OUTPUT/tmp.multicov;
	[ "$EXTENSION" == "bam" ] && multicovfile=$OUTPUT/tmp.multicov;
	[ -f $OUTPUT/$ANA.CANOES_Reads.tsv ] || cat $multicovfile | awk '{print $1"\t"$2"\t"$3}' > $OUTPUT/$ANA.CANOES_Reads.tsv
	cat $multicovfile | awk '{print $NF}' > $OUTPUT/$SAMPLE.multicov.txt;
	paste -d'\t' $OUTPUT/$ANA.CANOES_Reads.tsv $OUTPUT/$SAMPLE.multicov.txt > $OUTPUT/tmp.txt;
	cat $OUTPUT/tmp.txt > $OUTPUT/$ANA.CANOES_Reads.tsv;
	rm -rf $OUTPUT/tmp.multicov;
	rm -rf $OUTPUT/$SAMPLE.multicov.txt;
done < $OUTPUT/SampleSheet.SampleSex.txt
cp $OUTPUT/$ANA.CANOES_Reads.tsv $OUTPUT/all.CANOES.coverage.tsv
entete=$(cat $OUTPUT/SampleSheet.SampleSex.txt | awk '{print $1}' | paste -sd "\t" | sed 's/\r//g');
sed -i 1i"Chr\tStart\tEnd\t$entete" $OUTPUT/all.CANOES.coverage.tsv;
sed -e "s/^chrX/chr23/" $OUTPUT/$ANA.CANOES_Reads.tsv | sed -e "s/^chrY/chr24/" > $OUTPUT/newCount.reads_$ANA.txt;


### AUTOSOMAL CHROM for ALL
### format:
### chrom start end coverage
[ -s $OUTPUT/$ANA.CANOES_Reads.ALL.tsv ] || cat $OUTPUT/newCount.reads_$ANA.txt | grep -vE 'chr23|chr24' > $OUTPUT/$ANA.CANOES_Reads.ALL.tsv


### SEXUAL CHROM
### catch chrom X & Y for the sex specific analysis
### format:
### chrom start end coverage
cp $OUTPUT/newCount.reads_$ANA.txt $OUTPUT/$ANA.CANOES_Reads.withSEX.tsv;
sed -i 1i"Chr\tStart\tEnd\t$entete" $OUTPUT/$ANA.CANOES_Reads.withSEX.tsv;
while read SAMPLE SEX
do
	if [[ $SEX == "SEX_M" || $SEX == "SEX_F" ]]
	then
		[ -s $OUTPUT/$ANA.CANOES_Reads.$SEX.tsv ] || cat $OUTPUT/$ANA.CANOES_Reads.withSEX.tsv | awk '{print $1"\t"$2"\t"$3}' > $OUTPUT/$ANA.CANOES_Reads.$SEX.tsv
		[ ! -s $OUTPUT/$ANA.CANOES_Reads.withSEX.tsv ] || cat $OUTPUT/$ANA.CANOES_Reads.withSEX.tsv | awk -v header=$SAMPLE ' BEGIN { FS="\t"; c=0 } NR == 1 { for (i=1;i<=NF;i++) { if ($i==header) { c=i }} } NR > 1 && c>0 { print $c } ' > $OUTPUT/$SAMPLE.reads.txt;
		[ ! -s $OUTPUT/$ANA.CANOES_Reads.withSEX.tsv ] || sed -i 1i" $SAMPLE" $OUTPUT/$SAMPLE.reads.txt;
		[ ! -s $OUTPUT/$ANA.CANOES_Reads.withSEX.tsv ] || paste -d'\t' $OUTPUT/$ANA.CANOES_Reads.$SEX.tsv $OUTPUT/$SAMPLE.reads.txt > $OUTPUT/tmp.txt;
		[ ! -s $OUTPUT/$ANA.CANOES_Reads.withSEX.tsv ] || cat $OUTPUT/tmp.txt > $OUTPUT/$ANA.CANOES_Reads.$SEX.tsv
		rm -rf $OUTPUT/$SAMPLE.reads.txt;
	fi;
done < $OUTPUT/SampleSheet.SampleSex.txt
# delete header
[ ! -s $OUTPUT/$ANA.CANOES_Reads.SEX_M.tsv ] || sed -i '1d' $OUTPUT/$ANA.CANOES_Reads.SEX_M.tsv;
[ ! -s $OUTPUT/$ANA.CANOES_Reads.SEX_F.tsv ] || sed -i '1d' $OUTPUT/$ANA.CANOES_Reads.SEX_F.tsv;


###coverage boxplot
grep -vE "chrX|chrY" $OUTPUT/all.CANOES.coverage.tsv > $OUTPUT/all.CANOES.stats.tsv.tmp;
head -1 $OUTPUT/all.CANOES.stats.tsv.tmp | awk '{printf $1":"$2"-"$3"\t"} {for(i=4;i<=NF;i++){printf "%s\t", $i}; printf "\n"}' > $OUTPUT/all.CANOES.stats.tsv;
#sed 's/    0/      1/g'
cat $OUTPUT/all.CANOES.stats.tsv.tmp | sed -e "1d"| awk '{printf $1":"$2"-"$3"\t"} {for(i=4;i<=NF;i++){printf "%s\t", $i/($3-$2)}; printf "\n"}' | sed 's/\,/./g' >> $OUTPUT/all.CANOES.stats.tsv;
sed -i 's/.$//' $OUTPUT/all.CANOES.stats.tsv;

###coverage barplot
head -1 $OUTPUT/all.CANOES.coverage.tsv | awk '{printf $1":"$2"-"$3"\t"} {for(i=4;i<=NF;i++){printf "%s\t", $i}; printf "\n"}' > $OUTPUT/all.CANOES.stats2.tsv;
cat $OUTPUT/all.CANOES.coverage.tsv | sed -e "1d"| awk '{printf $1":"$2"-"$3"\t"} {for(i=4;i<=NF;i++){printf "%s\t", $i/($3-$2)}; printf "\n"}' | sed 's/\,/./g' >> $OUTPUT/all.CANOES.stats2.tsv;
sed -i 's/.$//' $OUTPUT/all.CANOES.stats2.tsv;
#cat $OUTPUT/all.CANOES.stats.tsv.tmp1 | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$16}' > $OUTPUT/all.CANOES.stats.tsv
###attention dernière colonne vide!!


###write script
echo "setwd(\"$OUTPUT\")" > $OUTPUT/script_create_coverage.$ANA.R;
echo "source(\"$CANOE/plotCoverage.R\")" >> $OUTPUT/script_create_coverage.$ANA.R;
echo "PlotCoverage(\"$SAMPLENAME\",\"$OUTPUT/all.CANOES.stats.tsv\",\"$OUTPUT/$SAMPLENAME.CANOES.boxplot.coverage.png\")" >> $OUTPUT/script_create_coverage.$ANA.R;
echo "PlotCoverage2(\"$SAMPLENAME\",\"$OUTPUT/all.CANOES.stats2.tsv\",\"$OUTPUT/$SAMPLENAME.CANOES.barplot.coverage.png\")" >> $OUTPUT/script_create_coverage.$ANA.R;
cd $OUTPUT;
$R CMD BATCH $OUTPUT/script_create_coverage.$ANA.R;


### R analysis
### CANOES script
#[[ -z $list ]] && list=( "ALL" "SEX_M" "SEX_F" )
[[ -z $list ]] && list="ALL,SEX_M,SEX_F"
LISTS=$(echo $list | tr , "\n")

for SEX in $LISTS    
do
	if [ -s $OUTPUT/$ANA.CANOES_Reads.$SEX.tsv ]
	then
		SEX_R=$(cat $OUTPUT/SampleSheet.$SEX.txt | sed 's/^/\"/g' | sed 's/$/\"/g' | sed ':a;N;$!ba;s/\n/,/g');
		echo "setwd(\"$OUTPUT\")" > $OUTPUT/script_create.$SEX.$ANA.R;
		echo "gc <- read.table(\"$OUTPUT/$ANA.CANOES_GC.$SEX.txt\")\$V2" >> $OUTPUT/script_create.$SEX.$ANA.R;
		echo "canoes.reads <- read.table(\"$OUTPUT/$ANA.CANOES_Reads.$SEX.tsv\")" >> $OUTPUT/script_create.$SEX.$ANA.R;
		echo "sample.names <- c($SEX_R)" >> $OUTPUT/script_create.$SEX.$ANA.R;
		echo "names(canoes.reads) <- c(\"chromosome\", \"start\", \"end\", sample.names)" >> $OUTPUT/script_create.$SEX.$ANA.R;
		echo "target <- seq(1, nrow(canoes.reads))" >> $OUTPUT/script_create.$SEX.$ANA.R;
		echo "canoes.reads <- cbind(target, gc, canoes.reads)" >> $OUTPUT/script_create.$SEX.$ANA.R;
		echo "source(\"$CANOE/CANOES_with_SEX.R\")" >> $OUTPUT/script_create.$SEX.$ANA.R;
		echo "xcnv.list <- vector('list', length(sample.names))" >> $OUTPUT/script_create.$SEX.$ANA.R;
		echo "for (i in 1:length(sample.names)){" >> $OUTPUT/script_create.$SEX.$ANA.R;
		echo "xcnv.list[[i]] <- CallCNVs(sample.names[i], canoes.reads)" >> $OUTPUT/script_create.$SEX.$ANA.R;
		echo "}" >> $OUTPUT/script_create.$SEX.$ANA.R;
		echo "xcnvs <- do.call('rbind', xcnv.list)" >> $OUTPUT/script_create.$SEX.$ANA.R;
		touch $OUTPUT/results_$ANA.$SEX.csv
		echo "write.csv(xcnvs, \"$OUTPUT/results_$ANA.$SEX.csv\", row.names=FALSE)" >> $OUTPUT/script_create.$SEX.$ANA.R;
		cd $OUTPUT;
		$R CMD BATCH $OUTPUT/script_create.$SEX.$ANA.R;
		entete=$(echo "$SEX_R" | sed 's/,/\t/g');
		sed -i "1iChr\tStart\tEnd\t$entete" $OUTPUT/$ANA.CANOES_Reads.$SEX.tsv;
		if [[ $SEX == "ALL" ]];
		then
			cat $OUTPUT/results_$ANA.$SEX.csv | sed 's/"//g' | sed 's/,/\t/g' | sed '1d' >> $OUTPUT/final_results_$ANA.unsorted.csv;
		else
			cat $OUTPUT/results_$ANA.$SEX.csv | grep -E '"23:|"24:' | sed 's/"//g' | sed 's/,/\t/g' >> $OUTPUT/final_results_$ANA.unsorted.csv;
		fi;
	fi;
done

[[ -f $OUTPUT/script_create.ALL.$ANA.Rout ]] && RSCRIPT_ERROR_ALL=$(grep -E "Error|Warning" $OUTPUT/script_create.ALL.Rout -A 1);
[[ -f $OUTPUT/script_create.SEX_M.$ANA.Rout ]] && RSCRIPT_ERROR_SEX_M=$(grep -E "Error|Warning" $OUTPUT/script_create.SEX_M.Rout -A 1);
[[ -f $OUTPUT/script_create.SEX_F.$ANA.Rout ]] && RSCRIPT_ERROR_SEX_F=$(grep -E "Error|Warning" $OUTPUT/script_create.SEX_F.Rout -A 1);
RSCRIPT_ERROR=$(grep "Error" $OUTPUT/script_create.*.Rout -A 1);
[[ -z $RSCRIPT_ERROR ]] || { echo $RSCRIPT_ERROR; exit 0; }


sort -u $OUTPUT/final_results_$ANA.unsorted.csv > $OUTPUT/final_results_$ANA.csv;
### AnnotCNV script
### Annot canoes results
if [ -f $OUTPUT/final_results_$ANA.csv ];
then {
	### AnnotCNV annotation for canoe results ###
	## transform to AnnotCNV format 
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
		echo -e "$chrom\t$start\t$end\t$CNV\t$patient" >> $OUTPUT/$ANA.CANOES.bed;
	done < $OUTPUT/final_results_$ANA.csv;

	if [ -f  $OUTPUT/$ANA.CANOES.bed ];
	then {
		################################
		### AnnotCNV annotation for canoes results
		#source /home1/TOOLS/tools/stark/dev/env.sh 
		$ANNOTCNV/Sources/AnnotCNV-main.tcl -bedFile $OUTPUT/$ANA.CANOES.bed;
		### output fiel : $OUTPUT/all.CANOES.annotated.tsv
	}
	fi ;
	cat $OUTPUT/$ANA.CANOES.annotated.tsv | head -1 > $OUTPUT/all.CANOES.annotated.tsv
	sed -i '1d' $OUTPUT/$ANA.CANOES.annotated.tsv
	sort -g $OUTPUT/$ANA.CANOES.annotated.tsv >> $OUTPUT/all.CANOES.annotated.tsv
}
else {
	echo -e "#########  WARNING  ################" >> $OUTPUT/all.CANOES.annotated.tsv;
	echo -e "##### Error on CANOES Analysis #####" >> $OUTPUT/all.CANOES.annotated.tsv;
	echo -e "####################################" >> $OUTPUT/all.CANOES.annotated.tsv;
	grep "Error" $OUTPUT/script_create.$SEX.$ANA.Rout -A 1 >> $OUTPUT/all.CANOES.annotated.tsv;
	clean="NO";
}
fi ;


### Create Sample specific canoes_results-AnnotCNV files
### format:
### SEX, SAMPLES USED, CNV results by samples
if [ ! -z "$SAMPLENAME" ]
then
	SAMPLEINFO=$(cat $OUTPUT/SampleSheet.SampleSex.Complet.txt | grep "$SAMPLENAME\ ");
	echo $SAMPLEINFO > $OUTPUT/SampleSheet.SampleSex.sample.txt;
else
	cat $OUTPUT/SampleSheet.SampleSex.Complet.txt > $OUTPUT/SampleSheet.SampleSex.sample.txt;
fi;

echo -e "##\n#$(head -1 $OUTPUT/all.CANOES.annotated.tsv)" > $OUTPUT/header.txt;
while read SAMPLE SEX
do 
	CATCH_ERROR=$(grep "Error" $OUTPUT/CANOES.log -A 1);
	[ -z $CATCH_ERROR ] || echo "WARNING: $CATCH_ERROR" >> $OUTPUT/$SAMPLENAME.CANOES.Annotated.tsv;
	SAMPLETYPE=""
	SAMPLETYPE=$(echo "$SAMPLE" | grep -E 'POOL|BlcADN|blanc|BlcPCR|blcPCR|_NTC')
	if [[ -z $CATCH_ERROR ]]
	then
		if [ -n "$SAMPLETYPE" ]
		then
			echo "Sample $SAMPLE removed from CANOES analysis because POOL of sample or BLANK sample" > $OUTPUT/$SAMPLE.CANOES.Annotated.tsv;
		else
			if [[ $SEX == "SEX_M" || $SEX == "SEX_F" ]]
			then
				echo "##Sex for sample $SAMPLE is : $SEX" >> $OUTPUT/$SAMPLE.CANOES.Annotated.tsv;
				echo "##Autosomal chromosomes: Samples used for CANOES Analysis : $(cat $OUTPUT/SampleSheet.SampleSex.txt | wc -l)/$(cat $OUTPUT/SampleSheet.SampleSex.Complet.txt | wc -l) (recommended 30 samples minimum) : $(cat $OUTPUT/SampleSheet.SampleSex.txt | awk -F ' ' '{print $1}' | tr '\n' ' ')" >> $OUTPUT/$SAMPLE.CANOES.Annotated.tsv;
				echo "##Sexual chromosomes - Samples used for CANOES Analysis with $SEX : $(cat $OUTPUT/SampleSheet.SampleSex.txt | grep $SEX | wc -l) (recommended 30 samples minimum) : $(cat $OUTPUT/SampleSheet.SampleSex.txt | grep $SEX | awk -F ' ' '{print $1}' | tr '\n' ' ')"$'' >> $OUTPUT/$SAMPLE.CANOES.Annotated.tsv;
			else
				echo "##WARNING: /!\/!\/!\ Sex is UNKNOWN for Sample $SAMPLE" >> $OUTPUT/$SAMPLE.CANOES.Annotated.tsv;
				echo "##WARNING: /!\/!\/!\ CANOES launched without sexual chromosomes" >> $OUTPUT/$SAMPLE.CANOES.Annotated.tsv;
				echo "##Autosomal chromosomes - Samples used for CANOES Analysis : $(cat $OUTPUT/SampleSheet.SampleSex.txt | wc -l)/$(cat $OUTPUT/SampleSheet.SampleSex.Complet.txt | wc -l) (recommended 30 samples minimum) : $(cat $OUTPUT/SampleSheet.SampleSex.txt | awk -F ' ' '{print $1}' | tr '\n' ' ')"$'' >> $OUTPUT/$SAMPLE.CANOES.Annotated.tsv;
			fi;
			RSCRIPT_ERROR_ALL=$(grep -E "Error|Warning" $OUTPUT/script_create.ALL.$ANA.Rout -A 2);
			RSCRIPT_ERROR_SEX=$(grep -E "Error|Warning" $OUTPUT/script_create.SEX*.$ANA.Rout -A 2);
			if [[ ! -z $RSCRIPT_ERROR_ALL ]]
			then
				echo '##Autosomal chromosomes - /!\/!\/!\ A WARNING message appears, please check your results CAREFULLY: $RSCRIPT_ERROR_ALL' >> $OUTPUT/$SAMPLE.CANOES.Annotated.tsv;
				echo -e "##Autosomal chromosomes - /!\/!\/!\ A WARNING message appears, please check your results CAREFULLY: $RSCRIPT_ERROR_ALL\n$(cat $OUTPUT/all.CANOES.annotated.tsv)" > $OUTPUT/all.CANOES.annotated.tsv;
			else
				echo "##Autosomal chromosomes - The analysis ENDED WELL!" >> $OUTPUT/$SAMPLE.CANOES.Annotated.tsv;
				echo -e "##Autosomal chromosomes - The analysis ENDED WELL!\n$(cat $OUTPUT/all.CANOES.annotated.tsv)" > $OUTPUT/all.CANOES.annotated.tsv;
			fi;
			if [[ ! -z $RSCRIPT_ERROR_SEX ]]
			then
				echo "##Sexual chromosomes - /!\/!\/!\ A WARNING message appears, please check your results CAREFULLY: $RSCRIPT_ERROR_SEX" >> $OUTPUT/$SAMPLE.CANOES.Annotated.tsv;
				echo -e "##Sexual chromosomes - /!\/!\/!\ A WARNING message appears, please check your results CAREFULLY: $RSCRIPT_ERROR_SEX\n$(cat $OUTPUT/all.CANOES.annotated.tsv)" > $OUTPUT/all.CANOES.annotated.tsv;
			else
				echo "##Sexual chromosomes - The analysis ENDED WELL!" >> $OUTPUT/$SAMPLE.CANOES.Annotated.tsv;
				echo -e "##Sexual chromosomes - The analysis ENDED WELL!\n$(cat $OUTPUT/all.CANOES.annotated.tsv)" > $OUTPUT/all.CANOES.annotated.tsv;
			fi;
			head -2 $OUTPUT/header.txt >> $OUTPUT/$SAMPLE.CANOES.Annotated.tsv;
			cat $OUTPUT/all.CANOES.annotated.tsv | grep $SAMPLE >> $OUTPUT/$SAMPLE.CANOES.Annotated.tsv;
		fi;
	else
		echo "WARNING: $CATCH_ERROR" >> $OUTPUT/$SAMPLENAME.CANOES.Annotated.tsv;
		echo "WARNING: $CATCH_ERROR" >> $OUTPUT/all.CANOES.annotated.tsv;
	fi;
done < $OUTPUT/SampleSheet.SampleSex.sample.txt

rm -rf $OUTPUT/.RData
rm -rf $OUTPUT/header.txt


### Clean tempory files
if [ "$clean" == "YES" ] ;
	then 
	{
		#### clean
		rm -rf $OUTPUT/bamlist.txt
		rm -rf $OUTPUT/*.bed
		rm -rf $OUTPUT/*.CANOES_GC.txt
		rm -rf $OUTPUT/*.R
		rm -rf $OUTPUT/*.Rout
		rm -rf $OUTPUT/all.CANOES.stats*
		rm -rf $OUTPUT/SAMPLE.*
		rm -rf $OUTPUT/SampleSheet*
		rm -rf $OUTPUT/newCount.reads*
		rm -rf $OUTPUT/results_*
		rm -rf $OUTPUT/CNV_Results*
		rm -rf $OUTPUT/final_results*
		rm -rf $OUTPUT/*CANOES_Reads*
		rm -rf $OUTPUT/tmp.txt
		rm -rf $OUTPUT/tmp.ALL.txt
		rm -rf $OUTPUT/*.tmp
		rm -rf $OUTPUT/*.tmp1
		rm -rf $OUTPUT/*.CANOES_GC.*
		rm -rf $OUTPUT/*multicovlist.txt
		rm -rf $OUTPUT/$ANA.CANOES.annotated.tsv
	}
fi ;

exit 0
