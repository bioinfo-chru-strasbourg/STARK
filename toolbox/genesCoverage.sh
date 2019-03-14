#! /bin/sh
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Authors: Amandine VELT, Antony Le Béchec
# Date: 04/04/2017
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Define the function to print the usage of the script
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    sh genesCoverage.sh  -b bedfile-genes -c coverage-criteria -n nb-bases-arounds -t bedtools -u bedtools2 -s samtools -o output [-h]

		Description:
		    This script allows to calculate the coverage at $coverage_criteria X for each gene of the given bedfile. The output is a file with a line per gene with
		    the corresponding percent of bases cover by more than $coverage_criteria

		Options:
		 	-f, --bam-file the bam file of the sample
		 	This option is required.
		 	-b, --bedfile-genes the bed file containing the coordinates of all the exons (5'UTR/3'UTR and CDS)
		 	This option is required. 
		 	-c, --coverage-criteria the coverage criteria to calculte the percentage of bases cover more than this value (eg 100 for 100 X)
		 	This option is required.
		 	-n, --nb-bases-arounds The number of bases to look around the exons coordinates (eg 10)
		 	This option is required.
		 	-t, --bedtools Path to bedtools
		 	This option is required.
		 	-u, --bedtools2 Path to more recent bedtools
		 	This option is required.
		 	-s, --samtools Path to samtools
		 	This option is required.
		 	--threads number of threads to use for the coverage calculation
		 	This option is not required (default=1).
		 	-o, --output output file
		 	This option is required.
		 	-h, --help
		        	Print this message and exit the program.
		__EOF__
}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ARGS=$(getopt -o "f:b:c:n:t:u:s:o:h" --long "bam-file:,bedfile-genes:,coverage-criteria:,nb-bases-arounds:,bedtools:,bedtools2:,samtools:,output:,threads:,help" -- "$@" 2> /dev/null)
[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && \
	exit 1
eval set -- "$ARGS"
while true
do
	case "$1" in
		-f|--bam-file)
			BAM_FILE="$2"
			shift 2 
			;;
		-b|--bedfile-genes)
			BEDFILE_GENES="$2"
			shift 2 
			;;
		-c|--coverage-criteria)
			COVERAGE_CRITERIA="$2"
			shift 2 
			;;
		-n|--nb-bases-arounds)
			NB_BASES_AROUND="$2"
			shift 2 
			;;
		-t|--bedtools)
			BEDTOOLS="$2"
			shift 2 
			;;
		-u|--bedtools2)
			BEDTOOLS2="$2"
			shift 2 
			;;
		-s|--samtools)
			SAMTOOLS="$2"
			shift 2 
			;;
		-o|--output)
			OUTPUT="$2"
			shift 2 
			;;
		--threads)
			THREADS="$2"
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
[ "$BAM_FILE" == "" ] || [ "$BEDFILE_GENES" == "" ] || [ "$COVERAGE_CRITERIA" == "" ] || [ "$NB_BASES_AROUND" == "" ] || [ "$BEDTOOLS" == "" ] || [ "$BEDTOOLS2" == "" ] || [ "$SAMTOOLS" == "" ] || [ "$OUTPUT" == "" ] && \
	echo "Options --bam_file, --bedfile-genes, --coverage-criteria, --nb-bases-arounds, --output, --samtools --bedtools2 and --bedtools are required. " "Use -h or --help to display the help." && exit 1;
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## INPUT PARAMETERS
if [ -z "$THREADS" ]; then
	THREADS=1
fi;


BEDFILE_GENES_CUT=$BEDFILE_GENES.$RANDOM.cut

if (($NB_BASES_AROUND)); then
	# add of number of bases around the bed
	awk -F"\t" -v x=$NB_BASES_AROUND '{ print $1"\t"$2-x"\t"$3+x"\t"$4 }' $BEDFILE_GENES > $BEDFILE_GENES.bases_arounds.bed; mv $BEDFILE_GENES.bases_arounds.bed $BEDFILE_GENES_CUT
	# sort the bed file
	cat  $BEDFILE_GENES_CUT | sort -k1,1 -k2,2n > $BEDFILE_GENES.sort.bed; mv $BEDFILE_GENES.sort.bed $BEDFILE_GENES_CUT
	# merge the overlapping coordinates of the bed file
	cat $BEDFILE_GENES_CUT | $BEDTOOLS2/mergeBed -i - -c 4 -o distinct -delim "|" > $BEDFILE_GENES.merge.bed
	mv $BEDFILE_GENES.merge.bed $BEDFILE_GENES_CUT
	#rm $BEDFILE_GENES.merge.bed
else
	cp $BEDFILE_GENES $BEDFILE_GENES_CUT
fi;
# we recover the gene list from our bed file
#list_genes=$( cut -f4 $BEDFILE_GENES | tr "\n" "" | sort | uniq )
list_genes=$( cut -f4 $BEDFILE_GENES_CUT | sort | uniq  | tr "\n" " ")
#echo $list_genes; exit 0;
# we recover the samplename of the current sample
samplename=$( basename $BAM_FILE | cut -d. -f1 )
# we did the coverage along all the genes. "$SAMTOOLS view -uF 0x400" allows to ignore the duplicate reads
#echo "$SAMTOOLS view -uF 0x400 $BAM_FILE | $BEDTOOLS/coverageBed -abam - -b $BEDFILE_GENES -d > $BEDFILE_GENES.coverage_bases" 


## COVERAZGE BASES calculation
################################

#head $BEDFILE_GENES_CUT
#grep -P "^$chr\t" $BEDFILE_GENES_CUT | head
#exit 0;

MK=$BEDFILE_GENES_CUT.mk

echo "all: $BEDFILE_GENES_CUT.coverage_bases" > $MK

BEDFILE_GENES_CHR_COVERAGE_BASES_LIST=""

#if (($($SAMTOOLS idxstats $BAM_FILE | awk '{SUM+=$3+$4} END {print SUM}'))); then
#	for chr in $($SAMTOOLS idxstats $BAM_FILE | grep -v "\*" | awk '{ if ($3+$4>0) print $1 }'); do
if (($(cat $BEDFILE_GENES_CUT | cut -f1 | sort -u | wc -l))); then
	for chr in $(cat $BEDFILE_GENES_CUT | cut -f1 | sort -u); do
		#echo $chr;
		BEDFILE_GENES_CHR_COVERAGE_BASES_LIST=$BEDFILE_GENES_CHR_COVERAGE_BASES_LIST" $BEDFILE_GENES_CUT.$chr.coverage_bases"
		echo "

$BEDFILE_GENES_CUT.$chr.bed: $BEDFILE_GENES_CUT
	-grep -P \"^$chr\\t\" $BEDFILE_GENES_CUT > $BEDFILE_GENES_CUT.$chr.bed
	if [ ! -s $BEDFILE_GENES_CUT.$chr.bed ]; then touch $BEDFILE_GENES_CUT.$chr.bed; fi;

$BEDFILE_GENES_CUT.$chr.coverage_bases: $BAM_FILE $BEDFILE_GENES_CUT.$chr.bed
	-$SAMTOOLS view -uF 0x400 $BAM_FILE $chr | $BEDTOOLS/coverageBed -abam - -b $BEDFILE_GENES_CUT.$chr.bed -d > $BEDFILE_GENES_CUT.$chr.coverage_bases;
	if [ ! -s $BEDFILE_GENES_CUT.$chr.coverage_bases ]; then touch $BEDFILE_GENES_CUT.$chr.coverage_bases; fi;

" >> $MK
	#head $BEDFILE_GENES_CUT.$chr.coverage_bases
done; fi;

echo "$BEDFILE_GENES_CUT.coverage_bases: $BEDFILE_GENES_CHR_COVERAGE_BASES_LIST
	cat $BEDFILE_GENES_CHR_COVERAGE_BASES_LIST > $BEDFILE_GENES_CUT.coverage_bases
" >> $MK

#cat $MK


#THREADS=1
make -j $THREADS -f $MK $BEDFILE_GENES_CUT.coverage_bases #1>/dev/null 2>/dev/null


## OLD but OK
if ((0)); then
	echo "$SAMTOOLS view -uF 0x400 $BAM_FILE | $BEDTOOLS/coverageBed -abam - -b $BEDFILE_GENES.cut -d > $BEDFILE_GENES.coverage_bases"
	$SAMTOOLS view -uF 0x400 $BAM_FILE | $BEDTOOLS/coverageBed -abam - -b $BEDFILE_GENES.cut -d > $BEDFILE_GENES.coverage_bases
fi;



#cat $BEDFILE_GENES.coverage_bases
#echo "test"; exit 0;



#################################
# coverage table to open in excel
#################################



TMPDIR=${OUTPUT}.$RANDOM.split_genes_DP
mkdir -p $TMPDIR

DP_LIST_FILE=""
DP_LIST_FILE_LATEX=""
#DP_LIST_FILE_LATEX_WARNING=""
DP_HEADER=$TMPDIR/header.genes_stats;
#for DP in $COVERAGE_CRITERIA; do
for DP in $(echo $COVERAGE_CRITERIA | tr "," " "); do
	HEADER1="#Gene\tNbases";
	HEADER2="%bases >"$DP"X\tNbases <"$DP"X";
	echo -e $HEADER2 > $TMPDIR/$DP.genes_stats;
	#awk -v DIR=$TMPDIR -v DP=$DP -F"\t" '{g=$4; gc=g; gsub(/\|/,"_",gc)} {a[g]++} {b[g]=b[g]+0} $6>=DP{b[g]++} {c[g]=(b[g]/a[g])*100} {c2[g]=sprintf("%.2f", c[g])} {d[g]=a[g]-b[g]} {print g"\t"a[g]"\t"c2[g]"\t"d[g] > DIR"/file_"gc".csv"}' $BEDFILE_GENES_CUT.coverage_bases ;
	#echo "awk -v DIR=$TMPDIR -v DP=$DP -F"\t" '{g=$4; gc=g; gsub(/\|/,"_",gc)} {a[g]++} {b[g]=b[g]+0} $6>=DP{b[g]++} {c[g]=(b[g]/a[g])*100} {c2[g]=sprintf("%.2f", c[g])} {d[g]=a[g]-b[g]} END  {print g"\t"a[g]"\t"c2[g]"\t"d[g] > DIR"/file_"gc".csv"}' $BEDFILE_GENES_CUT.coverage_bases"
	#awk -v DIR=$TMPDIR -v DP=$DP -F"\t" '{g=$4; gc=g; gsub(/\|/,"_",gc)} {a[g]++} {b[g]=b[g]+0} $6>=DP{b[g]++} {c[g]=(b[g]/a[g])*100} {c2[g]=sprintf("%.2f", c[g])} {d[g]=a[g]-b[g]} END  {print g"\t"a[g]"\t"c2[g]"\t"d[g] > DIR"/file_"gc".csv"}' $BEDFILE_GENES_CUT.coverage_bases ;
	awk -v DIR=$TMPDIR -v DP=$DP -F"\t" '{g=$4; gc=g; gsub(/\|/,"_",gc)} {a[g]++} {b[g]=b[g]+0} $6>=DP{b[g]++} {c[g]=(b[g]/a[g])*100} {c2[g]=sprintf("%.2f", c[g])} {d[g]=a[g]-b[g]} END { for (gene in a)  {print gene"\t"a[gene]"\t"c2[gene]"\t"d[gene] } }' $BEDFILE_GENES_CUT.coverage_bases | sort > $TMPDIR/$DP.genes_stats.tmp # OK new
	cat $TMPDIR/$DP.genes_stats.tmp | cut -f3,4 >> $TMPDIR/$DP.genes_stats;
	#awk -v DIR=$TMPDIR -v DP=$DP -F"\t" '{g=$4; gc=g; gsub(/\|/,"_",gc)} {a[g]++} {b[g]=b[g]+0} $6>=DP{b[g]++} {c[g]=(b[g]/a[g])*100} {c2[g]=sprintf("%.2f", c[g])} {d[g]=a[g]-b[g]} END { for (gene in a)  {print gene"\t"a[gene]"\t"c2[gene]"\t"d[gene] } }' $BEDFILE_GENES_CUT.coverage_bases | sort | cut -f3,4 >> $TMPDIR/$DP.genes_stats
	#tail -n 1 -q $TMPDIR/file_*.csv | cut -f3,4 >> $TMPDIR/$DP.genes_stats;
	#find $TMPDIR -type f -name file_\*.csv | xargs tail -n 1 -q | cut -f3,4 >> $TMPDIR/$DP.genes_stats;
	echo -e $HEADER2 > $TMPDIR/$DP.genes_stats.latex;
	#echo -e $HEADER2 > $TMPDIR/$DP.genes_stats.warning.latex;
	#tail -n 1 -q $TMPDIR/file_*.csv | cut -f3,4 | awk -F"\t" '{C="gray"} $1<95{C="orange"} $1<90{C="red"} {print "\\\\color{"C"}{\\\\databar{"$1"}}\t\\\\color{"C"}{"$2"}"}' >> $TMPDIR/$DP.genes_stats.latex;
	cat $TMPDIR/$DP.genes_stats.tmp | cut -f3,4 | awk -F"\t" '{C="gray"} $1<100{C="yellow"} $1<95{C="orange"} $1<90{C="red"} {print "\\\\color{gray}{\\\\databar{"$1"}}\t\\\\color{"C"}{"$2"}"}' >> $TMPDIR/$DP.genes_stats.latex;
	#tail -n 1 -q $TMPDIR/file_*.csv | cut -f3,4 | awk -F"\t" '{C="gray"} $1<95{C="orange"} $1<90{C="red"} $1<100{print "\\\\color{"C"}{\\\\databar{"$1"}}\t\\\\color{gray}{"$2"}"}' >> $TMPDIR/$DP.genes_stats.warning.latex;
	#if [ ! -s $DP_HEADER ]; then echo -e $HEADER1 > $DP_HEADER; tail -n 1 -q $TMPDIR/*csv | cut -f1,2 >> $DP_HEADER; fi;
	if [ ! -s $DP_HEADER ]; then echo -e $HEADER1 > $DP_HEADER; cat $TMPDIR/$DP.genes_stats.tmp | cut -f1,2 >> $DP_HEADER; fi;
	DP_LIST_FILE="$DP_LIST_FILE $TMPDIR/$DP.genes_stats";
	DP_LIST_FILE_LATEX="$DP_LIST_FILE_LATEX $TMPDIR/$DP.genes_stats.latex";
	#DP_LIST_FILE_LATEX_WARNING="$DP_LIST_FILE_LATEX_WARNING $TMPDIR/$DP.genes_stats.warning.latex";
	#rm $TMPDIR/*csv
	#find $TMPDIR -type f -name \*.csv | xargs rm
	
done;
paste $DP_HEADER $DP_LIST_FILE > ${OUTPUT}-all.txt
paste $DP_HEADER $DP_LIST_FILE_LATEX  | sed 's/\t/ \& /gi' | sed 's/%/\\\\%/gi' | sed 's/>/\\\\textgreater/gi' | sed 's/</\\\\textless/gi' > ${OUTPUT}-all-latex.txt
#paste $DP_HEADER $DP_LIST_FILE_LATEX_WARNING  | sed 's/\t/ \& /gi' | sed 's/%/\\\\%/gi' | sed 's/>/\\\\textgreater/gi' | sed 's/</\\\\textless/gi' > ${OUTPUT}-warning-latex.txt
rm -rf $TMPDIR $BEDFILE_GENES.coverage_bases


rm $BEDFILE_GENES_CUT*


exit 0;








### OLD PART####
################



# GOOD except BEDFILE_GENES_CUT with "_CUT"
if ((0)); then

TMPDIR=${OUTPUT}.$RANDOM.split_genes_DP
mkdir -p $TMPDIR

DP_LIST_FILE=""
DP_LIST_FILE_LATEX=""
#DP_LIST_FILE_LATEX_WARNING=""
DP_HEADER=$TMPDIR/header.genes_stats;
#for DP in $COVERAGE_CRITERIA; do
for DP in $(echo $COVERAGE_CRITERIA | tr "," " "); do
	HEADER1="#Gene\tNbases";
	HEADER2="%bases >$DP\tNbases <$DP";
	#awk -v DIR=$TMPDIR -v DP=$DP -F"\t" '{g=$4; gc=g; gsub(/\|/,"_",gc)} {a[g]++} {b[g]=b[g]+0} $6>=DP{b[g]++} {c[g]=(b[g]/a[g])*100} {c2[g]=sprintf("%.2f", c[g])} {d[g]=a[g]-b[g]} {print g"\t"a[g]"\t"c2[g]"\t"d[g] > DIR"/file_"gc".csv"}' $BEDFILE_GENES.coverage_bases ;
	awk -v DIR=$TMPDIR -v DP=$DP -F"\t" '{g=$4; gc=g; gsub(/\|/,"_",gc)} {a[g]++} {b[g]=b[g]+0} $6>=DP{b[g]++} {c[g]=(b[g]/a[g])*100} {c2[g]=sprintf("%.2f", c[g])} {d[g]=a[g]-b[g]} {print g"\t"a[g]"\t"c2[g]"\t"d[g] > DIR"/file_"gc".csv"}' $BEDFILE_GENES.coverage_bases ;
	echo -e $HEADER2 > $TMPDIR/$DP.genes_stats;
	#tail -n 1 -q $TMPDIR/file_*.csv | cut -f3,4 >> $TMPDIR/$DP.genes_stats;
	find $TMPDIR -type f -name file_\*.csv | xargs tail -n 1 -q | cut -f3,4 >> $TMPDIR/$DP.genes_stats;
	echo -e $HEADER2 > $TMPDIR/$DP.genes_stats.latex;
	#echo -e $HEADER2 > $TMPDIR/$DP.genes_stats.warning.latex;
	#tail -n 1 -q $TMPDIR/file_*.csv | cut -f3,4 | awk -F"\t" '{C="gray"} $1<95{C="orange"} $1<90{C="red"} {print "\\\\color{"C"}{\\\\databar{"$1"}}\t\\\\color{"C"}{"$2"}"}' >> $TMPDIR/$DP.genes_stats.latex;
	find $TMPDIR -type f -name file_\*.csv | xargs tail -n 1 -q | cut -f3,4 | awk -F"\t" '{C="gray"} $1<100{C="yellow"} $1<95{C="orange"} $1<90{C="red"} {print "\\\\color{gray}{\\\\databar{"$1"}}\t\\\\color{"C"}{"$2"}"}' >> $TMPDIR/$DP.genes_stats.latex;
	#tail -n 1 -q $TMPDIR/file_*.csv | cut -f3,4 | awk -F"\t" '{C="gray"} $1<95{C="orange"} $1<90{C="red"} $1<100{print "\\\\color{"C"}{\\\\databar{"$1"}}\t\\\\color{gray}{"$2"}"}' >> $TMPDIR/$DP.genes_stats.warning.latex;
	#if [ ! -s $DP_HEADER ]; then echo -e $HEADER1 > $DP_HEADER; tail -n 1 -q $TMPDIR/*csv | cut -f1,2 >> $DP_HEADER; fi;
	if [ ! -s $DP_HEADER ]; then echo -e $HEADER1 > $DP_HEADER; find $TMPDIR -type f -name file_\*.csv | xargs tail -n 1 -q | cut -f1,2 >> $DP_HEADER; fi;
	DP_LIST_FILE="$DP_LIST_FILE $TMPDIR/$DP.genes_stats";
	DP_LIST_FILE_LATEX="$DP_LIST_FILE_LATEX $TMPDIR/$DP.genes_stats.latex";
	#DP_LIST_FILE_LATEX_WARNING="$DP_LIST_FILE_LATEX_WARNING $TMPDIR/$DP.genes_stats.warning.latex";
	#rm $TMPDIR/*csv
	find $TMPDIR -type f -name \*.csv | xargs rm
done;
paste $DP_HEADER $DP_LIST_FILE > ${OUTPUT}-all.txt
paste $DP_HEADER $DP_LIST_FILE_LATEX  | sed 's/\t/ \& /gi' | sed 's/%/\\\\%/gi' | sed 's/>/\\\\textgreater/gi' | sed 's/</\\\\textless/gi' > ${OUTPUT}-all-latex.txt
#paste $DP_HEADER $DP_LIST_FILE_LATEX_WARNING  | sed 's/\t/ \& /gi' | sed 's/%/\\\\%/gi' | sed 's/>/\\\\textgreater/gi' | sed 's/</\\\\textless/gi' > ${OUTPUT}-warning-latex.txt
rm -rf $TMPDIR $BEDFILE_GENES.coverage_bases

fi;

#################################
# coverage table to open in excel
#################################
if ((0)); then
	first=1
	# for each coverage criteria we did the calculations.
	# there is a single table for all the coverage criteria, so for the first entry in the for loop, we write the first header columns (gene name, gene length, ...) and then for the other entries, we put only the calculations for the current coverage criteria
	#for cov in $COVERAGE_CRITERIA
	for cov in $(echo $COVERAGE_CRITERIA | tr "," " ")
	do
		if [[ "$first" == "1" ]]
		then
			# first header of the table
			echo "#Gene name	Nb of studied genes bases	% bases > ${cov}	Nb bases < ${cov}" > ${OUTPUT}_${cov}
			# coverage calculation fo each gene, put in the table
			for gene in `echo $list_genes`
			do
				grep "\<$gene\>" $BEDFILE_GENES.coverage_bases > $BEDFILE_GENES.coverage_bases_tmp
				total_number_bases=$( wc -l  $BEDFILE_GENES.coverage_bases_tmp | cut -d" " -f1 )
				number_bases_over_criteria=$( awk -v x=$cov ' $6 >= x ' $BEDFILE_GENES.coverage_bases_tmp | wc -l | cut -d" " -f1 )
				percent_bases_over_criteria=$( echo "scale=2; $number_bases_over_criteria*100/$total_number_bases" | bc ); nb_bases_under_criteria=$( echo "scale=2; ${total_number_bases}-${number_bases_over_criteria}" | bc ); echo "${gene}	${total_number_bases}	${percent_bases_over_criteria}	${nb_bases_under_criteria}" >> ${OUTPUT}_${cov}
				rm $BEDFILE_GENES.coverage_bases_tmp
			done
		else
			# after the first entry, we write only the percent bases > $cov ... for each coverage criteria
			echo "% bases > ${cov}	Nb bases < ${cov}" > ${OUTPUT}_${cov}
			# coverage calculation fo each gene, put in the table
			for gene in `echo $list_genes`
			do
				grep "\<$gene\>" $BEDFILE_GENES.coverage_bases > $BEDFILE_GENES.coverage_bases_tmp
				total_number_bases=$( wc -l  $BEDFILE_GENES.coverage_bases_tmp | cut -d" " -f1 )
				number_bases_over_criteria=$( awk -v x=$cov ' $6 >= x ' $BEDFILE_GENES.coverage_bases_tmp | wc -l | cut -d" " -f1 )
				percent_bases_over_criteria=$( echo "scale=2; $number_bases_over_criteria*100/$total_number_bases" | bc ); nb_bases_under_criteria=$( echo "scale=2; ${total_number_bases}-${number_bases_over_criteria}" | bc ); echo "${percent_bases_over_criteria}	${nb_bases_under_criteria}" >> ${OUTPUT}_${cov}
				rm $BEDFILE_GENES.coverage_bases_tmp
			done
		fi
		let "first=first+1"
	done
	# we paste all the tables of all the coverage criteria to have a single table per bed at the end
	paste ${OUTPUT}_* > ${OUTPUT}-all.txt
fi;
#################################
# same table but write in latex language
#################################

#cat $BEDFILE_GENES.coverage_bases; exit 0;


if ((0)); then

	echo -n "#Gene name & Nb of studied bases" > ${OUTPUT}_latex;
	for cov in $(echo $COVERAGE_CRITERIA | tr "," " "); do
		echo -n " & \\\\% bases \\\\textgreater ${cov} & Nb bases \\\\textless ${cov}" >> ${OUTPUT}_latex;
	done;
	echo "" >> ${OUTPUT}_latex;
	for gene in `echo $list_genes`; do
		grep "\<$gene\>" $BEDFILE_GENES.coverage_bases > $BEDFILE_GENES.coverage_bases_tmp
		echo -n "${gene} & ${total_number_bases}" >> ${OUTPUT}_latex
		for cov in $(echo $COVERAGE_CRITERIA | tr "," " "); do
			#echo ${OUTPUT}_${cov}_latex
			total_number_bases=$( wc -l  $BEDFILE_GENES.coverage_bases_tmp | cut -d" " -f1 )
			number_bases_over_criteria=$( awk -v x=$cov ' $6 >= x ' $BEDFILE_GENES.coverage_bases_tmp | wc -l | cut -d" " -f1 )
			percent_bases_over_criteria=$( echo "scale=2; $number_bases_over_criteria*100/$total_number_bases" | bc ); 
			nb_bases_under_criteria=$( echo "scale=2; ${total_number_bases}-${number_bases_over_criteria}" | bc );
			echo -n " & \\\\databar{${percent_bases_over_criteria}} & ${nb_bases_under_criteria}" >> ${OUTPUT}_latex
			#rm $BEDFILE_GENES.coverage_bases_tmp
		done;
		echo " " >> ${OUTPUT}_latex
	done;
	#echo ${OUTPUT}_latex;
	#cat ${OUTPUT}_latex;
	
	cp ${OUTPUT}_latex ${OUTPUT}-all-latex.txt
	sed -i "s/\t//" ${OUTPUT}-all-latex.txt
	rm ${OUTPUT}_* $BEDFILE_GENES.coverage_bases $BEDFILE_GENES.coverage_bases_tmp

fi;


if ((0)); then

	first=1
	#for cov in $COVERAGE_CRITERIA
	for cov in $(echo $COVERAGE_CRITERIA | tr "," " ")
	do
		if [[ "$first" == "1" ]]
		then
			# first header of the table
			echo "#Gene name & Nb of studied bases & \\\\% bases \\\\textgreater ${cov} & Nb bases \\\\textless ${cov}" > ${OUTPUT}_${cov}_latex
			# coverage calculation fo each gene, put in the table
			for gene in `echo $list_genes`
			do
				grep "\<$gene\>" $BEDFILE_GENES.coverage_bases > $BEDFILE_GENES.coverage_bases_tmp
				total_number_bases=$( wc -l  $BEDFILE_GENES.coverage_bases_tmp | cut -d" " -f1 )
				number_bases_over_criteria=$( awk -v x=$cov ' $6 >= x ' $BEDFILE_GENES.coverage_bases_tmp | wc -l | cut -d" " -f1 )
				percent_bases_over_criteria=$( echo "scale=2; $number_bases_over_criteria*100/$total_number_bases" | bc ); nb_bases_under_criteria=$( echo "scale=2; ${total_number_bases}-${number_bases_over_criteria}" | bc ); echo "${gene} & ${total_number_bases} & \\\\databar{${percent_bases_over_criteria}} & ${nb_bases_under_criteria}" >> ${OUTPUT}_${cov}_latex
				rm $BEDFILE_GENES.coverage_bases_tmp
			done
		else
			# after the first entry, we write only the percent bases > $cov ... for each coverage criteria
			echo " & \\\\% bases \\\\textgreater ${cov} & Nb bases \\\\textless ${cov}" > ${OUTPUT}_${cov}_latex
			for gene in `echo $list_genes`
			do
				grep "\<$gene\>" $BEDFILE_GENES.coverage_bases > $BEDFILE_GENES.coverage_bases_tmp
				total_number_bases=$( wc -l  $BEDFILE_GENES.coverage_bases_tmp | cut -d" " -f1 )
				number_bases_over_criteria=$( awk -v x=$cov ' $6 >= x ' $BEDFILE_GENES.coverage_bases_tmp | wc -l | cut -d" " -f1 )
				percent_bases_over_criteria=$( echo "scale=2; $number_bases_over_criteria*100/$total_number_bases" | bc ); nb_bases_under_criteria=$( echo "scale=2; ${total_number_bases}-${number_bases_over_criteria}" | bc ); echo " & \\\\databar{${percent_bases_over_criteria}} & ${nb_bases_under_criteria}" >> ${OUTPUT}_${cov}_latex
				rm $BEDFILE_GENES.coverage_bases_tmp
			done
		fi
		let "first=first+1"
	done



	# as before, we paste the latex tables of all coverage criteria and we have juste to print this table in the report
	paste ${OUTPUT}_*_latex > ${OUTPUT}-all-latex.txt
	sed -i "s/\t//" ${OUTPUT}-all-latex.txt
	rm ${OUTPUT}_* $BEDFILE_GENES.coverage_bases

fi;



