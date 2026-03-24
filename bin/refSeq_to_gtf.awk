# AWK script to convert RefSeq gene annotation from NCBI into GTF format, removing stop codons from the ends of coding sequences and appending a running number to duplicate uses of the same transcript ID. The output is sorted by chromosome, strand, start and end position, and feature type.
# Usage: cat refSeq_genes.txt | awk -F '\t' -v OFS='\t' -f refSeq_to_gtf.awk | sort -k1,1V -k4,4n -k5,5n -k3,3 -S4G > refSeq_genes.gtf
function min(x, y) { return (x>y) ? y : x }
function max(x, y) { return (x<y) ? y : x }
{
	split($10, start, ",")
	split($11, end, ",")
	split($16, frame, ",")
	# remove stop codon from left end for coding genes on the minus strand
	if ($4=="-" && $14=="cmpl" && (start[1]!=$7 || (min(end[1],$8)-start[1]+frame[1])%3==0)) {
		$7+=3
		for (i in end)
			if ($7>=end[i] && $7<=end[i]+2)
				$7+=start[i+1]-end[i]
	}
	# remove stop codon from right end for coding genes on the plus strand
	if ($4=="+" && $15=="cmpl" && (end[$9]!=$8 || (end[$9]-max(start[$9],$7)+frame[$9])%3==0)) {
		$8-=3
		for (i in start)
			if ($8<=start[i] && $8>=start[i]-2)
				$8-=start[i]-end[i-1]
	}
	# append running number to duplicate uses of the same transcript ID
	gene_id=$13
	if (transcripts[$2]++) {
		gene_id=$13"_"transcripts[$2]
		$2=$2"_"transcripts[$2]
	}
	# print one line for each exon
	for (i=1; i<=$9; i++) {
		exon=($4=="+") ? i : $9-i+1
		attributes="gene_id \""gene_id"\"; transcript_id \""$2"\"; exon_number \""exon"\"; exon_id \""$2"."exon"\"; gene_name \""$13"\";"
		print $3,"RefSeq","exon",start[i]+1,end[i],".",$4,".",attributes
		# print one line for each coding region
		if ($14~/cmpl/ && $7<=end[i] && $8>=start[i])
			print $3,"RefSeq","CDS",max($7,start[i])+1,min($8,end[i]),".",$4,frame[i],attributes
	}
}