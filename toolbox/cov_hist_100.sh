/home1/TOOLS/tools/samtools/1.3.1/bin/samtools depth HORIZON_FASTQ.bwamem.bam | awk ' {DPT++} { for (i=10; i<=100; i++10) { if(>=i) {DP[i]++} } } END { for (k in DP) { print k DP[k] DP[k]/DPT } }' | sort