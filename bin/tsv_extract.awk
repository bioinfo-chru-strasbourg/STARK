#!/usr/bin/awk
# Use tsv_extract.awk -v cols=cal1,col2,col3...
#cat FILE | awk -f tsv_extract.awk -F"\t" -v COLS=col1,col2,col3...



BEGIN {
	if (length(DEFAULT_NA) == 0) DEFAULT_NA=""
	if (length(DEFAULT_SEP) == 0) DEFAULT_SEP="\t"
	if (length(COLS) == 0) COLS="#CHROM,POS,ID,REF,ALT,QUAL,FILTER"
	split(COLS,out,",")
}

NR==1 {
    for (i=1; i<=NF; i++) {
        ix[$i] = i
	}
	SEP=""
	for (i=1; i<=length(out); i++) {
		printf "%s%s", SEP, out[i]
		SEP=DEFAULT_SEP
	}
	print ""
}
NR>1 {
	SEP=""
    for (i=1; i<=length(out); i++) {
		if (ix[out[i]]>0) {
			printf "%s%s", SEP, $ix[out[i]]
		} else {
			printf "%s%s", SEP, DEFAULT_NA
		}
		SEP=DEFAULT_SEP
	}
    print ""
}
