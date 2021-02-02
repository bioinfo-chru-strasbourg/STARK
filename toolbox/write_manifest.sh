#!/bin/bash

BED_FILE=$1
DIRECTORY=$(dirname $BED_FILE)
TMP_FILE="$DIRECTORY/tmp.bed_to_manifest.bed"
MANIFEST="${BED_FILE%.*}.manifest"
MANIFEST_GENES="$MANIFEST.genes"

### create manifest:
# remove header
sed '1,2d' $BED_FILE > $TMP_FILE
# select specific columns
awk '{print $4"_"$1"_"$2"_"$3"\t"$1"\t"$2"\t"$3}' $TMP_FILE > $TMP_FILE.sort
# create manifest header
echo "#Some of the genomic variants, genes, nucleic acid sequences, or genomic regions on this list, and their use in specific applications, may be protected by patents.  Customers are advised to determine whether they are required to obtain licenses from the party that owns or controls such patents in order to use the product in Customer's specific application.  Unless expressly stated otherwise in writing by Illumina, Customer, and not Illumina, is solely responsible for obtaining such licenses." > $MANIFEST
echo "[Header]" >> $MANIFEST
echo "ReferenceGenome	Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFASTA" >> $MANIFEST
echo "[Regions]" >> $MANIFEST
#echo "Name	Chromosome	Amplicon Start	Amplicon End" >> $MANIFEST
echo "Name	Chromosome	Start	End" >> $MANIFEST
# copy in manifest file
cat $TMP_FILE.sort >> $MANIFEST
# remove tempory files
rm -rf $TMP_FILE
rm -rf $TMP_FILE.sort
[ -s $MANIFEST ] || echo "error"


### create manfest.genes:
# cut using RefSeq
/home1/TOOLS/tools/bedtools/current/bin/intersectBed -a /home1/DIAG/TOOLS/annovar_sources/genetrack-byExon.withoutUTR.RefSeq.GRCh37.bed -b $BED_FILE -wa > $DIRECTORY/intersectBed_genes.txt
# select specific columns
awk '{print $1"\t"$2"\t"$3"\t"$5}'  $DIRECTORY/intersectBed_genes.txt  >  $DIRECTORY/intersectBed_genes_no_sort.txt
# sort file
cat $DIRECTORY/intersectBed_genes_no_sort.txt | uniq > $MANIFEST_GENES
#remove tempory files
rm -rf $DIRECTORY/intersectBed_genes_no_sort.txt
rm -rf $DIRECTORY/intersectBed_genes.txt
[ -s $MANIFEST_GENES ] || echo "error"

echo "done"






