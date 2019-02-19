BAM=$1 #B00HLN0.bwamem.bam
MDP=100
CHR_F="" #"-r chr22:1-10000000"
BED=$2 #B00HLN0.bed
BED_F="-b $BED"
M=" -d 100 "
#time /home1/TOOLS/tools/samtools/1.3.1/bin/samtools depth $BAM $M $BED_F $CHR_F | head -n 50 #awk -v MDP=$MDP '{ if($3>=MDP) {DP[MDP]++} } { if($3<MDP) {DP[$3]++} } END { for (k in DP) { print k" "DP[k] } }' | sort -g -r | awk ' {SUM+=$2} {a[$1]=SUM} END { for (k in a) { print k" "a[k]/SUM} }' | sort -g
#time /home1/TOOLS/tools/samtools/1.3.1/bin/samtools depth $BAM $M $BED_F $CHR_F | awk -v MDP=$MDP '{ if($3>=MDP) {DP[MDP]++} } { if($3<MDP) {DP[$3]++} } END { for (k in DP) { print k" "DP[k] } }' | sort -g -r | awk ' {SUMM+=$2} {a[$1]=SUMM} END { for (k in a) { print k" "a[k]/SUMM" "a[k]"/"SUMM } }' | sort -g | head -n 100
#time /home1/TOOLS/tools/samtools/1.3.1/bin/samtools depth $BAM $BED_F $M | awk -v MDP=$MDP ' {DP[$3]++} END { for (k in DP) { print k" "DP[k] } }' | sort -g -r #| awk ' {SUM+=$2} {a[$1]=SUM} END { for (k in a) { print k" "a[k]/SUM}" "SUM }' | sort -g | head -n 101
#time $BEDTOOLS/bedtools genomecov -i $BAM   | more #awk -v MDP=$MDP '{ if($3>=MDP) {DP[MDP]++} } { if($3<MDP) {DP[$3]++} } END { for (k in DP) { print k" "DP[k] } }' | sort -g -r | awk ' {SUMM+=$2} {a[$1]=SUMM} END { for (k in a) { print k" "a[k]/SUMM" "a[k]"/"SUMM } }' | sort -g | head -n 100
#time /home1/TOOLS/tools/samtools/1.3.1/bin/samtools depth $BAM $M $BED_F $CHR_F | awk -v MDP=$MDP '{ if($3>=MDP) {DP[MDP]++} else {DP[$3]++} } END { for (i=0; i<=MDP; i++) { if (DP[i]!="") {previous=DP[i] }; print i" "previous } }' | sort -g -r | awk ' {SUMM+=$2} {a[$1]=SUMM} END { for (k in a) { print k" "a[k]/SUMM" "a[k]"/"SUMM } }' | sort -g #| head -n 100
#time $SAMTOOLS stats $BAM -t $BED -c 1,100,1 | grep ^COV | tac | awk '{SUM+=$4} {CUM[$3]=SUM} END { for (D in CUM) { print D" "CUM[D]" "CUM[D]/SUM} }' | sort -g 
time $SAMTOOLS depth $BAM -b $BED | awk -v MDP=$MDP '{SUM++} { if ($3>MDP) {DP[MDP]++} else {DP[$3]++} } END { for (i in DP) {print i" "DP[i]" SUM"SUM}}' | sort -g -r | awk '{SUM+=$2} {CUM[$1]=SUM}  END { for (j in CUM) {print j" "CUM[j]" "SUM" "(CUM[j]/SUM)} }' | sort -g | column -t
