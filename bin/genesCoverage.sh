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
		 	This option is optional (default="1,30,100").
			-n, --nb-bases-arounds The number of bases to look around the exons coordinates (eg 10)
		 	This option is optional (default=0).
			--dp_fail
		 	This option is optional (default=30).
			--dp_warn
		 	This option is optional (default=100).
			--dp_threshold
		 	This option is optional (default=1).
		 	-t, --bedtools Path to bedtools
		 	This option is required.
		 	-s, --samtools Path to samtools
		 	This option is required.
		 	--threads number of threads to use for the coverage calculation
		 	This option is optional (default=1).
		 	-o, --output output file
		 	This option is required.
		 	-h, --help
		    Print this message and exit the program.
			-v, --verbose
			Print verbose information.
		__EOF__
}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ARGS=$(getopt -o "f:b:c:n:t:u:s:o:vhd" --long "bam-file:,bedfile-genes:,coverage-criteria:,nb-bases-arounds:,dp_fail:,dp_warn:,dp_threshold:,bedtools:,samtools:,output:,threads:,verbose,debug,help" -- "$@" 2> /dev/null)
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
		--dp_fail)
			DP_FAIL="$2"
			shift 2
			;;
		--dp_warn)
			DP_WARN="$2"
			shift 2
			;;
		--dp_threshold)
			DP_THRESHOLD="$2"
			shift 2
			;;
		-t|--bedtools)
			BEDTOOLS="$2"
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
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Checking the input parameter
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
[ "$BAM_FILE" == "" ] || [ "$BEDFILE_GENES" == "" ] || [ "$BEDTOOLS" == "" ] || [ "$SAMTOOLS" == "" ] || [ "$OUTPUT" == "" ] && \
	echo "Options --bam_file, --bedfile-genes, --output, --samtools and --bedtools are required. " "Use -h or --help to display the help." && exit 1;
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## INPUT PARAMETERS
if [ -z "$THREADS" ]; then
	THREADS=1
fi;
if [ -z "$COVERAGE_CRITERIA" ]; then
	COVERAGE_CRITERIA="1,30,100"
fi;
if [ -z "$NB_BASES_AROUND" ]; then
	NB_BASES_AROUND=0
fi;
if [ -z "$DP_FAIL" ]; then
	DP_FAIL=30
fi;
if [ -z "$DP_WARN" ]; then
	DP_WARN=100
fi;
if [ -z "$DP_THRESHOLD" ]; then
	DP_THRESHOLD=1
fi;
if [ -z "$VERBOSE" ]; then
	VERBOSE=0
fi;

(($DEBUG)) && echo $BEDFILE_GENES && head $BEDFILE_GENES



BEDFILE_GENES_CHECKED=$BEDFILE_GENES.$RANDOM.checked

# Add strand in bed file if only 4 column
if [ "$(grep '^#' -v $BEDFILE_GENES | head -n 1  | awk -F'\t' '{print NF}')" == 4 ] || ! [[ "$(grep '^#' -v $BEDFILE_GENES | head -n 1 | awk -F'\t' '{print $4}')" =~ [+-] ]]; then
	awk -F"\t" '{ print $1"\t"$2"\t"$3"\t+\t"$4 }' $BEDFILE_GENES > $BEDFILE_GENES_CHECKED
else
	cp $BEDFILE_GENES $BEDFILE_GENES_CHECKED
fi

(($DEBUG)) && echo $BEDFILE_GENES_CHECKED && head $BEDFILE_GENES_CHECKED

BEDFILE_GENES_CUT=$BEDFILE_GENES_CHECKED.$RANDOM.cut

if (($NB_BASES_AROUND)); then
	# add of number of bases around the bed
	awk -F"\t" -v x=$NB_BASES_AROUND '{ print $1"\t"$2-x"\t"$3+x"\t"$4"\t"$5 }' $BEDFILE_GENES_CHECKED > $BEDFILE_GENES_CHECKED.bases_arounds.bed; mv $BEDFILE_GENES_CHECKED.bases_arounds.bed $BEDFILE_GENES_CUT
	# sort the bed file
	cat  $BEDFILE_GENES_CUT | sort -k1,1 -k2,2n > $BEDFILE_GENES_CHECKED.sort.bed;
	mv $BEDFILE_GENES_CHECKED.sort.bed $BEDFILE_GENES_CUT
	#echo ""; head $BEDFILE_GENES_CUT
	# merge the overlapping coordinates of the bed file
	#cat $BEDFILE_GENES_CUT | $BEDTOOLS2/mergeBed -i - -c 4 -o distinct -delim "|" > $BEDFILE_GENES.merge.bed
	cat $BEDFILE_GENES_CUT | $BEDTOOLS merge -i - -c 4,5 -o distinct -delim "|" > $BEDFILE_GENES_CHECKED.merge.bed
	mv $BEDFILE_GENES_CHECKED.merge.bed $BEDFILE_GENES_CUT
	#echo ""; head $BEDFILE_GENES_CUT
	#rm $BEDFILE_GENES.merge.bed
else
	cp $BEDFILE_GENES_CHECKED $BEDFILE_GENES_CUT
fi;

(($DEBUG)) && echo $BEDFILE_GENES_CUT && head $BEDFILE_GENES_CUT

#echo ""; head $BEDFILE_GENES_CUT

#exit 0
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
	#-$SAMTOOLS view -uF 0x400 $BAM_FILE $chr | $BEDTOOLS/coverageBed -abam - -b $BEDFILE_GENES_CUT.$chr.bed -d > $BEDFILE_GENES_CUT.$chr.coverage_bases;
	$SAMTOOLS depth $BAM_FILE -b $BEDFILE_GENES_CUT.$chr.bed -a -r $chr -d 0 | awk -F\"\t\" '{print \$\$1\"\t\"\$\$2\"\t\"\$\$2\"\t+\t\"\$\$3}' | $BEDTOOLS intersect -b stdin -a $BEDFILE_GENES_CUT.$chr.bed -wb | awk -F\"\t\" '{print \$\$1\"\t\"\$\$7\"\t\"\$\$7\"\t+\t\"\$\$10\"\t\"\$\$5}' > $BEDFILE_GENES_CUT.$chr.coverage_bases;
	if [ ! -s $BEDFILE_GENES_CUT.$chr.coverage_bases ]; then touch $BEDFILE_GENES_CUT.$chr.coverage_bases; fi;

" >> $MK
	#head $BEDFILE_GENES_CUT.$chr.coverage_bases
done; fi;

echo "$BEDFILE_GENES_CUT.coverage_bases: $BEDFILE_GENES_CHR_COVERAGE_BASES_LIST
	cat $BEDFILE_GENES_CHR_COVERAGE_BASES_LIST > $BEDFILE_GENES_CUT.coverage_bases
" >> $MK

#cat $MK


#THREADS=1
if [ ! -e ${OUTPUT}.coverage_bases ]; then
	echo "#[INFO] Coverage bases generation"
	make -j $THREADS -f $MK $BEDFILE_GENES_CUT.coverage_bases #1>/dev/null 2>/dev/null
	cp $BEDFILE_GENES_CUT.coverage_bases ${OUTPUT}.coverage_bases
else
	echo "#[INFO] Coverage stats already generated"
fi;
#echo ""; head $BEDFILE_GENES_CUT.coverage_bases
#exit 0


## OLD but OK
#if ((0)); then
#	echo "$SAMTOOLS view -uF 0x400 $BAM_FILE | $BEDTOOLS/coverageBed -abam - -b $BEDFILE_GENES.cut -d > $BEDFILE_GENES.coverage_bases"
#	$SAMTOOLS view -uF 0x400 $BAM_FILE | $BEDTOOLS/coverageBed -abam - -b $BEDFILE_GENES.cut -d > $BEDFILE_GENES.coverage_bases
#fi;



(($VERBOSE)) && echo "" && head -n 20 ${OUTPUT}.coverage_bases
#echo "test"; exit 0;


#################################
# coverage table to open in excel
#################################

echo "#[INFO] Genes Coverage stats calculation"

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

	awk -v DIR=$TMPDIR -v DP=$DP -F"\t" '
	{g=$6; gc=g; gsub(/\|/,"_",gc)}
	{a[g]++}
	{b[g]=b[g]+0}
	$5>=DP {b[g]++}
	{c[g]=(b[g]/a[g])*100}
	{c2[g]=sprintf("%.2f", c[g])}
	{d[g]=a[g]-b[g]}
	END {
		for (gene in a) {
			print gene"\t"a[gene]"\t"c2[gene]"\t"d[gene]
		}
	}' ${OUTPUT}.coverage_bases | sort > $TMPDIR/$DP.genes_stats.tmp # OK new
	cat $TMPDIR/$DP.genes_stats.tmp | cut -f3,4 >> $TMPDIR/$DP.genes_stats;

	echo -e $HEADER2 > $TMPDIR/$DP.genes_stats.latex;

	cat $TMPDIR/$DP.genes_stats.tmp | cut -f3,4 | awk -F"\t" '{C="gray"} $1<100{C="yellow"} $1<95{C="orange"} $1<90{C="red"} {print "\\\\color{gray}{\\\\databar{"$1"}}\t\\\\color{"C"}{"$2"}"}' >> $TMPDIR/$DP.genes_stats.latex;

	if [ ! -s $DP_HEADER ]; then echo -e $HEADER1 > $DP_HEADER; cat $TMPDIR/$DP.genes_stats.tmp | cut -f1,2 >> $DP_HEADER; fi;
	DP_LIST_FILE="$DP_LIST_FILE $TMPDIR/$DP.genes_stats";
	DP_LIST_FILE_LATEX="$DP_LIST_FILE_LATEX $TMPDIR/$DP.genes_stats.latex";

done;
paste $DP_HEADER $DP_LIST_FILE > ${OUTPUT}-all.txt
paste $DP_HEADER $DP_LIST_FILE_LATEX  | sed 's/\t/ \& /gi' | sed 's/%/\\\\%/gi' | sed 's/>/\\\\textgreater/gi' | sed 's/</\\\\textless/gi' > ${OUTPUT}-all-latex.txt
#paste $DP_HEADER $DP_LIST_FILE_LATEX_WARNING  | sed 's/\t/ \& /gi' | sed 's/%/\\\\%/gi' | sed 's/>/\\\\textgreater/gi' | sed 's/</\\\\textless/gi' > ${OUTPUT}-warning-latex.txt

(($VERBOSE)) && echo "" && head -n 20 ${OUTPUT}-all.txt

echo "#[INFO] Genes Coverage information messages"

awk -v DIR=$TMPDIR -v FAIL=$DP_FAIL -v WARN=$DP_WARN -v THRESHOLD=$DP_THRESHOLD -F"\t" '
length(THRESHOLD)==0 {THRESHOLD=1}
{g=$6; gc=g; gsub(/\|/,"_",gc)}
{a[g]++}
{bFAIL[g]=bFAIL[g]+0}
$5>=FAIL {bFAIL[g]++}
{cFAIL[g]=bFAIL[g]/a[g]}
{bWARN[g]=bWARN[g]+0}
$5>=WARN {bWARN[g]++}
{cWARN[g]=bWARN[g]/a[g]}
{M[g]="PASS"; D[g]="100% bases with DP >="WARN}
cWARN[g]<THRESHOLD {M[g]="WARN"; D[g]="only "sprintf("%.2f", cWARN[g]*100)"% bases with DP >="WARN}
cFAIL[g]<THRESHOLD {M[g]="FAIL"; D[g]="only "sprintf("%.2f", cFAIL[g]*100)"% bases with DP >="FAIL}
END {
	for (gene in a) {
		print gene"\t"M[gene]"\t"D[gene]
	}
}' ${OUTPUT}.coverage_bases | sort >> $TMPDIR/$DP.genes_message;

cp $TMPDIR/$DP.genes_message ${OUTPUT}.msg


#cat $TMPDIR/$DP.genes_message.tmp | cut -f3,4 >> $TMPDIR/$DP.genes_message;
#" (accepted threshold "sprintf("%.2f", THRESHOLD)"%)"

(($VERBOSE)) && echo "" && head -n 20 ${OUTPUT}.msg;


(($VERBOSE)) && echo "Genes with coverage passing: "$(awk -F"\t" '$2=="PASS" {print $1}' ${OUTPUT}.msg)
(($VERBOSE)) && echo "Genes with coverage warning: "$(awk -F"\t" '$2=="WARN" {print $1" "}' ${OUTPUT}.msg)
(($VERBOSE)) && echo "Genes with coverage failing: "$(awk -F"\t" '$2=="FAIL" {print $1}' ${OUTPUT}.msg)

rm -rf $TMPDIR $BEDFILE_GENES.coverage_bases #${OUTPUT}.coverage_bases

rm  $BEDFILE_GENES_CHECKED* $BEDFILE_GENES_CUT*


exit 0;
