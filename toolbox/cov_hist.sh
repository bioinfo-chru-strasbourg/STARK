/home1/TOOLS/tools/samtools/1.3.1/bin/samtools depth HORIZON_FASTQ.bwamem.bam | awk '{DP[]++} END { for (k in DP) { print k DP[k] } }'
