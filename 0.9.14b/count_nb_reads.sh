FILES="P974.unaligned.bam P974.bwamem.bam P974.bwasw.bam"

for F in $FILES
do
	echo "BAM File '$F'"
	echo "/NGS/bin/samtools view $F | cut -f1 | sort | uniq | wc -l"
	/NGS/bin/samtools view $F | cut -f1 | sort | uniq | wc -l
	echo "/NGS/bin/samtools view -c $F"
	/NGS/bin/samtools view -c $F
	echo "/NGS/bin/samtools flagstat $F"
	/NGS/bin/samtools flagstat $F
done;

