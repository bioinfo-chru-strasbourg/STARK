
FILE=$1

#for sample in `$BCFTOOLS view -h $FILE | grep "^#CHROM" | cut -f10-`; do
#	$BCFTOOLS view -c1 -Oz -s $sample -o ${FILE/.vcf*/.$sample.vcf.gz} $FILE;
#done;
java -jar $GATK -R /home1/IRC/TOOLS/genomes/hg19/hg19.fa -T CombineVariants $(for sample in $($BCFTOOLS view -h $FILE | grep "^#CHROM" | cut -f10-); do $BCFTOOLS view -c1 -Oz -s $sample -o ${FILE/.vcf*/.$sample.vcf.gz} $FILE; $TABIX ${FILE/.vcf*/.$sample.vcf.gz}; echo " --variant:"$sample" "${FILE/.vcf*/.$sample.vcf.gz}; done) -o ${FILE/.vcf*/.merged.vcf} -genotypeMergeOptions PRIORITIZE --assumeIdenticalSamples -priority $(for sample in $($BCFTOOLS view -h $FILE | grep "^#CHROM" | cut -f10-); do echo $sample"," | tr -d '[[:space:]]'; done;)

# PRIORITIZE UNIQUIFY
