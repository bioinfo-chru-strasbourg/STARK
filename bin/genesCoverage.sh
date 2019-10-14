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
			--coverage-bases the coverage bases files in format "chr pos depth"
		 	This option is optional (default generated with samtools depth from bam file).
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
ARGS=$(getopt -o "f:b:c:n:t:u:s:o:vhd" --long "bam-file:,bedfile-genes:,coverage-criteria:,coverage-bases:,nb-bases-arounds:,dp_fail:,dp_warn:,dp_threshold:,bedtools:,samtools:,output:,threads:,verbose,debug,help" -- "$@" 2> /dev/null)
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
		--coverage-bases)
			COVERAGE_BASES="$2"
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

if [ -d $TMP_FOLDER_TMP ] && [ "$TMP_FOLDER_TMP" != "" ]; then
	TMP_GENESCOVERAGE=$TMP_FOLDER_TMP/genesCoverage.$RANDOM
else
	TMP_GENESCOVERAGE=/tmp/genesCoverage.$RANDOM
fi;
mkdir -p $TMP_GENESCOVERAGE
(($DEBUG)) && echo $TMP_GENESCOVERAGE

PRECISION=4

(($DEBUG)) && echo $BEDFILE_GENES && head $BEDFILE_GENES

SORT_ORDER="-k1,1 -k2,2n"

#BEDFILE_GENES_CHECKED=$BEDFILE_GENES.$RANDOM.checked
BEDFILE_GENES_CHECKED=$TMP_GENESCOVERAGE/BEDFILE_GENES.$RANDOM.checked

# Normalize bed
awk -F"\t" '
{chr=$1}
{start=$2}
{stop=$3}
{strand=$4}
{gene=$5}
strand !~ /[+-]/ {strand="+"}
gene == "" { if ($4 !~ /[+-]/ && $4 != "") {gene=$4} else {gene=chr"_"start"_"stop} }
{print chr"\t"start"\t"stop"\t"strand"\t"gene}
' $BEDFILE_GENES | sort $SORT_ORDER > $BEDFILE_GENES_CHECKED


(($DEBUG)) && echo $BEDFILE_GENES_CHECKED && head $BEDFILE_GENES_CHECKED

#BEDFILE_GENES_CUT=$BEDFILE_GENES_CHECKED.$RANDOM.cut
BEDFILE_GENES_CUT=$TMP_GENESCOVERAGE/BEDFILE_GENES_CHECKED.$RANDOM.checked


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

if ((0)); then

	MK=$BEDFILE_GENES_CUT.mk

	echo "all: $BEDFILE_GENES_CUT.coverage_bases" > $MK

	BEDFILE_GENES_CHR_COVERAGE_BASES_LIST=""


	if (($(cat $BEDFILE_GENES_CUT | cut -f1 | sort -u | wc -l))); then
		for chr in $(cat $BEDFILE_GENES_CUT | cut -f1 | sort -u); do
			#echo $chr;
			BEDFILE_GENES_CHR_COVERAGE_BASES_LIST=$BEDFILE_GENES_CHR_COVERAGE_BASES_LIST" $BEDFILE_GENES_CUT.$chr.coverage_bases"
			echo "

	$BEDFILE_GENES_CUT.$chr.bed: $BEDFILE_GENES_CUT
		-grep -P \"^$chr\\t\" $BEDFILE_GENES_CUT > $BEDFILE_GENES_CUT.$chr.bed
		if [ ! -s $BEDFILE_GENES_CUT.$chr.bed ]; then touch $BEDFILE_GENES_CUT.$chr.bed; fi;

			" >> $MK


			if  ((0)) && [ "$COVERAGE_BASES" != "" ] && [ -s $COVERAGE_BASES ]; then
			# Coverage bases generation from Coverage Bases file

			if [ "${COVERAGE_BASES##*.}" = "gz" ]; then
				COVERAGE_BASES_VIEW=" $GZ -d -c "
			else
				COVERAGE_BASES_VIEW=" cat "
			fi

			echo "

	$BEDFILE_GENES_CUT.$chr.coverage_bases: $BAM_FILE $BEDFILE_GENES_CUT.$chr.bed
		$COVERAGE_BASES_VIEW $COVERAGE_BASES | grep -P \"^$chr\\t\"  | sort $SORT_ORDER | awk -F\"\t\" '{print \$\$1\"\t\"\$\$2\"\t\"\$\$2\"\t+\t\"\$\$3}' | $BEDTOOLS intersect -b stdin -a $BEDFILE_GENES_CUT.$chr.bed -wb -sorted | awk -F\"\t\" '{print \$\$10\"\t\"\$\$5}' > $BEDFILE_GENES_CUT.$chr.coverage_bases;
		if [ ! -s $BEDFILE_GENES_CUT.$chr.coverage_bases ]; then touch $BEDFILE_GENES_CUT.$chr.coverage_bases; fi;

			" >> $MK

			else
			# Coverage bases generation from samtools depth on BAM
			echo "

	$BEDFILE_GENES_CUT.$chr.coverage_bases: $BAM_FILE $BEDFILE_GENES_CUT.$chr.bed
		cat $BEDFILE_GENES_CUT.$chr.bed | awk -F"\t" '{print \$\$1\"\t\"\$\$2-1\"\t\"\$\$3\"\t4\t\"\$\$5}' > $BEDFILE_GENES_CUT.$chr.bed.0_based;
		$SAMTOOLS depth $BAM_FILE -b $BEDFILE_GENES_CUT.$chr.bed.0_based -a -r $chr -d 0 | sort $SORT_ORDER | awk -F\"\t\" '{print \$\$1\"\t\"\$\$2\"\t\"\$\$2\"\t+\t\"\$\$3}' | $BEDTOOLS intersect -b stdin -a $BEDFILE_GENES_CUT.$chr.bed -w -sorted | awk -F\"\t\" '{print \$\$10\"\t\"\$\$5}' > $BEDFILE_GENES_CUT.$chr.coverage_bases;
		rm $BEDFILE_GENES_CUT.$chr.bed.0_based;
		if [ ! -s $BEDFILE_GENES_CUT.$chr.coverage_bases ]; then touch $BEDFILE_GENES_CUT.$chr.coverage_bases; fi;

			" >> $MK
			fi;
		#head $BEDFILE_GENES_CUT.$chr.coverage_bases
	done; fi;

	echo "$BEDFILE_GENES_CUT.coverage_bases: $BEDFILE_GENES_CHR_COVERAGE_BASES_LIST
		cat $BEDFILE_GENES_CHR_COVERAGE_BASES_LIST > $BEDFILE_GENES_CUT.coverage_bases
	" >> $MK

	(($DEBUG)) && cat $MK



	#cat $$one_bed | awk -F"\t" '{print $$1"\t"$$2-1"\t"$$3"\t"$$4"\t"$$5}' > $$one_bed.0_based; \
	#$(SAMTOOLS) mpileup $(SAMTOOLS_METRICS_MPILEUP_PARAM) -aa -l $$one_bed.0_based $< | cut -f1,2,4

	#THREADS=1
	time if [ ! -e ${OUTPUT}.coverage_bases ]; then
		echo "#[INFO] Coverage bases generation"
		make -j $THREADS -f $MK $BEDFILE_GENES_CUT.coverage_bases 1>/dev/null 2>/dev/null
	else
		echo "#[INFO] Coverage stats already generated"
	fi;

fi;

COVERAGE_BASES_VIEW=" cat "
if  [ "$COVERAGE_BASES" != "" ] && [ -s $COVERAGE_BASES ]; then

	if [ "${COVERAGE_BASES##*.}" = "gz" ]; then
		COVERAGE_BASES_VIEW=" $GZ -d -c "
	fi

	echo "#[INFO] Coverage bases provided"
	#ln -s $COVERAGE_BASES $BEDFILE_GENES_CUT.coverage_bases
	cp $COVERAGE_BASES $BEDFILE_GENES_CUT.coverage_bases


else

	echo "#[INFO] Coverage bases generation"
	make -j $THREADS -f $MK $BEDFILE_GENES_CUT.coverage_bases 1>/dev/null 2>/dev/null

fi;

#echo $BEDFILE_GENES_CUT.coverage_bases
#ls -l $BEDFILE_GENES_CUT.coverage_bases
#echo $COVERAGE_BASES_VIEW $BEDFILE_GENES_CUT.coverage_bases
#$COVERAGE_BASES_VIEW $BEDFILE_GENES_CUT.coverage_bases | head
#zcat $BEDFILE_GENES_CUT.coverage_bases | head
#exit 0

#$COVERAGE_BASES_VIEW $BEDFILE_GENES_CUT.coverage_bases | head

(($VERBOSE)) && $COVERAGE_BASES_VIEW $BEDFILE_GENES_CUT.coverage_bases | head -n 20
#exit 0

#################################
# coverage table to open in excel
#################################



if ((1)); then

	echo "#[INFO] Genes Coverage stats calculation"
	#$COVERAGE_BASES_VIEW $BEDFILE_GENES_CUT.coverage_bases | awk -v FAIL=$DP_FAIL -v WARN=$DP_WARN -v THRESHOLD=$DP_THRESHOLD -v COVERAGE_CRITERIA=$COVERAGE_CRITERIA -v PRECISION=$PRECISION -F"\t" '
	$COVERAGE_BASES_VIEW $COVERAGE_BASES | awk -v FAIL=$DP_FAIL -v WARN=$DP_WARN -v THRESHOLD=$DP_THRESHOLD -v COVERAGE_CRITERIA=$COVERAGE_CRITERIA -v PRECISION=$PRECISION -F"\t" '
	BEGIN {
		{n=split(COVERAGE_CRITERIA,CC,",")}
	}
	{
		{split($2,genes,"|")}
		{for (j in genes) {nb_base[genes[j]]++}}
		{dp=$1}
		if (dp>0)
		{ for (i = 1; i <= n; i++) {
				COV=CC[i]
				if (dp>=COV) {
					for (j in genes) {bases[genes[j]][COV]++}
				} else {
					break
				}
			}
		}
	}
	END {
		#printf "%f","#Gene	Nbases"/NR
		#printf "%s","#Gene\tNbases"
		printf "#Gene\tNbases\tThreshold"
		for (i = 1; i <= n; i++) {
			COV=CC[i]
			printf "\t%Coverage "COV"X\t#Bases <"COV"X"
		}
		print ""
		#for (gene in nb_base) {
		for (gene in bases) {
			GENE_HEAD=gene"\t"nb_base[gene]
			GENE_COV=""
			GENE_MSG="PASS"
			for (i = 1; i <= n; i++) {
				COV=CC[i]
				{percent[gene][COV]=((bases[gene][COV]+0)/nb_base[gene])*100}
				{percent_precision[gene][COV]=sprintf("%."PRECISION"f", percent[gene][COV])}
				{nb_bases_failed[gene][COV]=nb_base[gene]-(bases[gene][COV]+0)}
				GENE_COV=GENE_COV"\t"percent_precision[gene][COV]"\t"nb_bases_failed[gene][COV]
				if ( COV==FAIL && percent[gene][COV]<(THRESHOLD*100) ) { GENE_MSG="FAIL" }
				if ( COV==WARN && GENE_MSG!="FAIL" && percent[gene][COV]<(THRESHOLD*100) ) { GENE_MSG="WARN" }
			}
			#print " "
			print GENE_HEAD "\t"GENE_MSG GENE_COV
		}
	}' | sort > ${OUTPUT}.genes.txt

	#(($VERBOSE)) && echo "genes.txt" && head -n 20 ${OUTPUT}.genes.txt
	(($VERBOSE)) && head -n 20 ${OUTPUT}.genes.txt


	# FAIL file
	grep "^#" ${OUTPUT}.genes.txt > ${OUTPUT}.genes.FAIL.txt
	cat ${OUTPUT}.genes.txt | awk -F'\t' '{ if ($3 == "FAIL") {print $0} }' >> ${OUTPUT}.genes.FAIL.txt

	# WARN file
	grep "^#" ${OUTPUT}.genes.txt > ${OUTPUT}.genes.WARN.txt
	cat ${OUTPUT}.genes.txt | awk -F'\t' '{ if ($3 == "WARN") {print $0} }' >> ${OUTPUT}.genes.WARN.txt

	# FAIL WARN file
	grep "^#" ${OUTPUT}.genes.txt > ${OUTPUT}.genes.FAIL.WARN.txt
	cat ${OUTPUT}.genes.txt | awk -F'\t' '{ if  ( ($3 == "FAIL") || ($3 == "WARN") ) {print $0} }' >> ${OUTPUT}.genes.FAIL.WARN.txt

	# WARN file
	grep "^#" ${OUTPUT}.genes.txt > ${OUTPUT}.genes.PASS.txt
	cat ${OUTPUT}.genes.txt | awk -F'\t' '{ if ($3 == "PASS") {print $0} }' >> ${OUTPUT}.genes.PASS.txt

	# MISS file
	grep -v "^#" ${OUTPUT}.genes.txt | cut -f1 | sort -u > ${OUTPUT}.genes.FOUND.txt
	grep -v "^#" $BEDFILE_GENES_CUT  | cut -f5 | sort -u > ${OUTPUT}.genes.CHECKED.txt
	echo "#Gene" > ${OUTPUT}.genes.MISS.txt
	diff ${OUTPUT}.genes.FOUND.txt ${OUTPUT}.genes.CHECKED.txt | grep "^>" | cut -d" " -f2 >> ${OUTPUT}.genes.MISS.txt


fi;


(($VERBOSE)) && echo "Genes with coverage passing (max 20): "$(awk -F"\t" '$3=="PASS" {print $1}' ${OUTPUT}.genes.txt | head -n 20)
(($VERBOSE)) && echo "Genes with coverage warning (max 20): "$(awk -F"\t" '$3=="WARN" {print $1" "}' ${OUTPUT}.genes.txt | head -n 20)
(($VERBOSE)) && echo "Genes with coverage failing (max 20): "$(awk -F"\t" '$3=="FAIL" {print $1}' ${OUTPUT}.genes.txt | head -n 20)

#rm -rf $TMPDIR $BEDFILE_GENES.coverage_bases #${OUTPUT}.coverage_bases $BEDFILE_GENES_CUT.coverage_bases

#rm -f $BEDFILE_GENES_CHECKED* $BEDFILE_GENES_CUT*

if ((1)); then
	if [ -d $TMP_GENESCOVERAGE ] && [ $TMP_GENESCOVERAGE != "" ]; then
		rm -rf $TMP_GENESCOVERAGE
	fi;
fi;


exit 0;
