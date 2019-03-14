



BAM=/home1/IRC/DATA/DEV/RES/170519_M01656_0172_000000000-B4W4H/Z_Horizon/Z_Horizon.bwamem.bam
cp /home1/IRC/DATA/RAW/Manifests/IFU447_TumorHotspotMASTRPlus_Manifest_v150804.AmpliconManifest.genes /home1/IRC/DATA/DEV/RES/170519_M01656_0172_000000000-B4W4H/IFU447_TumorHotspotMASTRPlus_Manifest_v150804.AmpliconManifest.genes
BEDFILE_GENES=/home1/IRC/DATA/DEV/RES/170519_M01656_0172_000000000-B4W4H/IFU447_TumorHotspotMASTRPlus_Manifest_v150804.AmpliconManifest.genes

MK=test.mk
echo "all: $BEDFILE_GENES.coverage_bases" > $MK

BEDFILE_GENES_CHR_COVERAGE_BASES_LIST=""

if (($($SAMTOOLS idxstats $BAM | awk '{SUM+=$3+$4} END {print SUM}'))); then
	for chr in $($SAMTOOLS idxstats $BAM | grep -v "\*" | awk '{ if ($3+$4>0) print $1 }'); do
		echo $chr;
		
		BEDFILE_GENES_CHR_COVERAGE_BASES_LIST=$BEDFILE_GENES_CHR_COVERAGE_BASES_LIST" $BEDFILE_GENES.$chr.coverage_bases"
		echo "

$BEDFILE_GENES.$chr.bed: $BEDFILE_GENES
	grep -P \"^$chr\\t\" $BEDFILE_GENES > $BEDFILE_GENES.$chr.bed

$BEDFILE_GENES.$chr.coverage_bases: $BAM $BEDFILE_GENES.$chr.bed
	$SAMTOOLS view -uF 0x400 $BAM $chr | $BEDTOOLS/coverageBed -abam - -b $BEDFILE_GENES.$chr.bed -d > $BEDFILE_GENES.$chr.coverage_bases

" >> $MK
	head $BEDFILE_GENES.$chr.coverage_bases
done; fi;

echo "$BEDFILE_GENES.coverage_bases: $BEDFILE_GENES_CHR_COVERAGE_BASES_LIST
	cat $BEDFILE_GENES_CHR_COVERAGE_BASES_LIST > $BEDFILE_GENES.coverage_bases
" >> $MK

cat $MK

make -j 1 -f $MK all





exit 0;



cp /home1/DIAG/DATA/NGS/Manifests/EXOME_DI_Agilent.manifest.genes /home1/DIAG/DATA/NGS/DEV/RES/170331_NB551027_0105_AHVGCYBGXY/
bedfile_genes=/home1/DIAG/DATA/NGS/DEV/RES/170331_NB551027_0105_AHVGCYBGXY/EXOME_DI_Agilent.manifest.genes
bedfile_name=TEST

/home1/TOOLS/tools/stark/0.9.14d/genesCoverage.sh -f /home1/DIAG/DATA/NGS/DEV/RES/170331_NB551027_0105_AHVGCYBGXY/SGT171378/SGT171378.bwamem.bam -b $bedfile_genes -c "1,30" -n 20 -t /home1/TOOLS/tools/bedtools/2.17.0/bin -u /home1/TOOLS/tools/bedtools/2.25.0/bin -s /home1/TOOLS/tools/samtools/1.3.1/bin/samtools -o /home1/DIAG/DATA/NGS/DEV/RES/170331_NB551027_0105_AHVGCYBGXY/SGT171378/SGT171378.bwamem.bam.metrics/$bedfile_name;




 grep -P '^chr1\t' $BEDFILE_GENES > $BEDFILE_GENES.$chr







