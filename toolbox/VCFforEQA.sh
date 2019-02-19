#!/bin/bash
#################################
##
## NGS environment
## VCF manipulation for EQA
##
## version: 0.9
## date: 18/12/2014
## author: Antony Le Bechec
##
#################################

BEDTOOLS=/media/IRCV2/NGSEnv/bin/bedtools
VCFTOOLS=/media/IRCV2/NGSEnv/bin/vcftools
SCRIPTS=/media/IRCV2/NGSEnv/scripts
SAMTOOLS=samtools
BED_GENELIST=bed_genelist.pl

VCF_INPUT=$1 	# EQA_NEQAS.bwamem.gatkHC.vap.vcf
BED_LAB=$2 	# Lab0036_Covered_bed.bed OR Lab0036_Covered_bed.bed.overlap.bed
VCF_BEDFILTERED=$VCF_INPUT.bedfiltered.vcf
VCF_EXONIC=$VCF_INPUT.bedfiltered.exonic.vcf
VCF_EXONIC_TXT=$VCF_INPUT.bedfiltered.exonic.vcf.txt
VCF_SIMPLIFIED=$VCF_INPUT.bedfiltered.exonic.simplified.vcf

# PARAMS
ANNOTATIONS="PZScore,PZFlag,PZComment,Symbol,hgvs,location,outcome,AlleleFrequency,AD,dbSNP,dbSNPNonFlagged,GGCPolymorphismsBRCA,DP,AF,VF,GQ,Ensembl,TI,FC,GWASCatalog,COSMIC,LocalDB,1000genomesALL,1000genomesEUR,6500NHLBIALL,6500NHLBIEUR,PolyPhen2HumanVarPred,PolyPhen2HumanDivPred,MutationTasterPred,MutationAssessorPred,LTRPred,IARCTP53,SIFT,phastCons,PhyloP,SiPhy,FATHMM,LRT,GERP,PolyPhen2HumanVar,PolyPhen2HumanDiv,MutationTaster,MutationAssessor,TFBS,FilterComment,ALL"
SORT_BY="CHROM,POS"
ORDER_BY="ASC,ASC"

# Hearder
grep ^# $VCF_INPUT > $VCF_INPUT.header

# Bedfilter
cat $VCF_INPUT.header > $VCF_BEDFILTERED
$BEDTOOLS/intersectBed -a $VCF_INPUT -b $BED_LAB >> $VCF_BEDFILTERED

# Exonic
cat $VCF_INPUT.header > $VCF_EXONIC
grep "location=exonic" $VCF_BEDFILTERED >> $VCF_EXONIC
grep "location=splic" $VCF_BEDFILTERED >> $VCF_EXONIC
$SCRIPTS/VCFtranslation.pl --input_file=$VCF_EXONIC --output_file=$VCF_EXONIC_TXT  --annotation="$ANNOTATIONS" --sort_by="$SORT_BY" --order_by="$ORDER_BY" 1>/dev/null 2>/dev/null

# Simplified and cleaned VCF
FieldsToRevome="Symbol,1000genomesALL,Ensembl,outcome,dbSNPNonFlagged,location,UCSC,knownGene,RefSeq,1000genomesEUR";
FieldsToRevome=$FieldsToRevome",LRT,SiPhy,MutationTaster,MutationAssessor,PolyPhen2HumanDivPred,MutationTasterPred,6500NHLBIALL,PhyloP";
FieldsToRevome=$FieldsToRevome",MutationAssessorPred,PolyPhen2HumanVar,6500NHLBIEUR,PolyPhen2HumanVarPred,LRTPred,FATHMM,PolyPhen2HumanDiv";
FieldsToRevome=$FieldsToRevome",SIFT,PolyPhen2,GERP,CLINVAR,COSMIC,GGCPolymorphismsBRCA,IARCTP53,LocalDB,phastCons,snpid,hgvs,dbSNP";
FieldsToRevome=$FieldsToRevome",AC,FS,QD,ReadPosRankSum,MLEAF,MLEAC,BaseQRankSum,DP,MQ0,AF,AN,MQ,MQRankSum,ClippingRankSum,AD,GQ,PL";
FieldsToRevome=$FieldsToRevome",DB,DS,HaplotypeScore,InbreedingCoeff,GATKCommandLine/HaplotypeCaller,FILTER/LowQual,contig/*";
#AC=1;FS=0.000;QD=14.85;ReadPosRankSum=1.244;MLEAF=0.500;MLEAC=1;BaseQRankSum=0.659;DP=28;MQ0=0;AF=0.500;AN=2;MQ=60.00;MQRankSum=-0.073;ClippingRankSum=-1.537
cat EQA_NEQAS.bwamem.gatkHC.vap.vcf.bedfiltered.exonic.vcf | /media/IRCV2/NGSEnv/bin/vcftools/vcf-annotate -r "$FieldsToRevome" | grep -v "^##source_" | sed 's/EQA_NEQAS/Lab_0036/g'  > $VCF_SIMPLIFIED


# OUTPUT
echo "# VCF_INPUT:          $VCF_INPUT";
echo "# BED_LAB:            $BED_LAB";
echo "# VCF_BEDFILTERED:    $VCF_BEDFILTERED";
echo "# VCF_EXONIC:         $VCF_EXONIC";
echo "# VCF_EXONIC_TXT:     $VCF_EXONIC_TXT";
echo "# VCF_SIMPLIFIED:     $VCF_SIMPLIFIED";

# Validation
echo "# VCF Validation:"
VCF_VALIDATION=`/media/IRCV2/NGSEnv/bin/vcftools/vcf-validator $VCF_SIMPLIFIED`




