RUN1=130730_M01658_0013_000000000-A3TH8
#RUN1=131028_M01658_0017_000000000-A541T
RUN2=131028_M01658_0017_000000000-A541T
SAMPLES="CPS0117I04 CPS0433III24 CPS2195II04 CPS2721III04"
RESULTS_FOLDER=/mnt/IRC/RES/ALL
#RUN1_LA=Alignment13
RUN1_LA=Alignment16
RUN2_LA=Alignment3
ISEC=/NGS/bin/vcftools/vcf-isec

echo "########################################################"
echo "# RUN1=$RUN1 (Alignment $RUN1_LA)"
echo "# RUN2=$RUN2 (Alignment $RUN2_LA)"
echo "# SAMPLES=$SAMPLES"
echo "#"

for SAMPLE in $SAMPLES; do
	VCF1=$RESULTS_FOLDER/$RUN1/$SAMPLE/DATA/MiSeq-Aligner-$RUN1_LA/MiSeq-Caller/trakxs/$SAMPLE.annotated.vcf
	VCF2=$RESULTS_FOLDER/$RUN2/$SAMPLE/DATA/MiSeq-Aligner-$RUN2_LA/MiSeq-Caller/trakxs/$SAMPLE.annotated.vcf
	
	bgzip -f $VCF1
	tabix $VCF1.gz
	bgzip -f $VCF2
	tabix $VCF2.gz
	mkdir -p $SAMPLE
	perl $ISEC -p $SAMPLE/$SAMPLE"_" $VCF1.gz $VCF2.gz
	perl $ISEC --nfiles =2 $VCF1.gz $VCF2.gz > $SAMPLE/$SAMPLE"_isecannotation.vcf"
	
	
	NBVAR_VCF1=`gunzip $VCF1.gz -c | grep -v "^#" -c`
	NBVAR_VCF2=`gunzip $VCF2.gz -c | grep -v "^#" -c`
	NBVAR_INTERSECTION_VCF1_VCF2=`gunzip $SAMPLE/$SAMPLE"_0_1.vcf.gz" -c | grep -v "^#" -c`
	NBVAR_UNIQUE_VCF1=`gunzip $SAMPLE/$SAMPLE"_0.vcf.gz" -c | grep -v "^#" -c`
	NBVAR_UNIQUE_VCF2=`gunzip $SAMPLE/$SAMPLE"_1.vcf.gz" -c | grep -v "^#" -c`
	NBVAR_TOTAL=`echo "scale=2; $NBVAR_UNIQUE_VCF1+$NBVAR_UNIQUE_VCF2+$NBVAR_INTERSECTION_VCF1_VCF2" | bc`
	NBVAR_PERCENT=`echo "scale=2; ($NBVAR_INTERSECTION_VCF1_VCF2/$NBVAR_TOTAL)*100" | bc`
	
	NBSNP_VCF1=`gunzip $VCF1.gz -c | grep -v "^#" | grep "dbSNP=rs" -c`
	NBSNP_VCF2=`gunzip $VCF2.gz -c | grep -v "^#" | grep "dbSNP=rs" -c`
	NBSNP_INTERSECTION_VCF1_VCF2=`gunzip $SAMPLE/$SAMPLE"_0_1.vcf.gz" -c | grep -v "^#" | grep "dbSNP=rs" -c`
	NBSNP_UNIQUE_VCF1=`gunzip $SAMPLE/$SAMPLE"_0.vcf.gz" -c | grep -v "^#" | grep "dbSNP=rs" -c`
	NBSNP_UNIQUE_VCF2=`gunzip $SAMPLE/$SAMPLE"_1.vcf.gz" -c | grep -v "^#" | grep "dbSNP=rs" -c`
	#echo "echo scale=2; $NBSNP_UNIQUE_VCF1+$NBSNP_UNIQUE_VCF2+$NBSNP_INTERSECTION_VCF1_VCF2 | bc"
	NBSNP_TOTAL=`echo "scale=2; $NBSNP_UNIQUE_VCF1+$NBSNP_UNIQUE_VCF2+$NBSNP_INTERSECTION_VCF1_VCF2" | bc`
	NBSNP_PERCENT=`echo "scale=2; ($NBSNP_INTERSECTION_VCF1_VCF2/$NBSNP_TOTAL)*100" | bc`
	
	
	echo "# $SAMPLE"
	echo "#    Variant distribution:"
	echo "#       NBVAR_VCF1:                     $NBVAR_VCF1"
	echo "#       NBVAR_VCF2:                     $NBVAR_VCF2"
	echo "#       NBVAR_INTERSECTION_VCF1_VCF2:   $NBVAR_INTERSECTION_VCF1_VCF2 ($NBVAR_PERCENT%)"
	echo "#       NBVAR_UNIQUE_VCF1:              $NBVAR_UNIQUE_VCF1"
	echo "#       NBVAR_UNIQUE_VCF2:              $NBVAR_UNIQUE_VCF2"
	echo "#    SNP distribution:"
	echo "#       NBSNP_VCF1:                     $NBSNP_VCF1"
	echo "#       NBSNP_VCF2:                     $NBSNP_VCF2"
	echo "#       NBSNP_INTERSECTION_VCF1_VCF2:   $NBSNP_INTERSECTION_VCF1_VCF2 ($NBSNP_PERCENT%)"
	echo "#       NBSNP_UNIQUE_VCF1:              $NBSNP_UNIQUE_VCF1"
	echo "#       NBSNP_UNIQUE_VCF2:              $NBSNP_UNIQUE_VCF2"
	echo "#"
	
done;

