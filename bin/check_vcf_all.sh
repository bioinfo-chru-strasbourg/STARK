#!/bin/bash
# Launch ALL analysis with data test

RUNS=$1
#SAMPLES=$2
#VCFS=$3

if [ "$RUNS" == "" ]; then RUNS="*"; fi;
#if [ "$SAMPLES" == "" ]; then SAMPLES="*"; fi;

# INPUT
echo "# RUNS     $RUNS"
#echo "# SAMPLES  $SAMPLES"
echo "# VCFS     "$(ls *.vcf)

for VCF in *.vcf; do
	echo -e "#\n###### VCF: $VCF ###"
	VCF_SAMPLE=$(echo $(basename "$VCF") | sed "s/.vcf$//gi")
	./check_vcf.sh "$RUNS" "*$VCF_SAMPLE*" "$VCF"
done;

exit 0;
for RUN in "$RUNS"; do
	for SAMPLE in "$SAMPLES"; do

		for VCF_i in $RUN/$SAMPLE/*vcf; do
			VCF_i_NBVARIANT=$(grep -cv ^# $VCF_i) #| sed "s/:/ /gi" | column -t 2>/dev/null
			#echo "$VCF_i nbv:"$(grep -cv ^# $VCF_i);
			for VCF in "$VCFS"; do
			 	# compare VCFs
				#echo "$VCF nbv:"$(grep -cv ^# $VCF);
				VCF_NBVARIANT=$(grep -cv ^# $VCF) #| sed "s/:/ /gi" | column -t 2>/dev/null
				#echo "$VCF_i_NBVARIANT" -eq "$VCFNB_VARIANT"
				if [ "$VCF_i_NBVARIANT" -eq "$VCF_NBVARIANT" ]; then
					echo "# OK";
				else
					echo "# ko";
				fi;

			done;
		done;
	done;
done;





