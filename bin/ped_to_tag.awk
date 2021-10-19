BEGIN {
    if (TAGS_sep=="") TAGS_sep="!"
    if (TAGS_header=="") TAGS_header="PEDIGREE"
    TAGS=""
}
{
    if (NR==1) {
        #print "HEADER " $0
        split($0, HEADER)
    } else {
        #print SAMPLE " - "  " | " $0
        split($0, IND)
        if (SAMPLE=="" || SAMPLE==IND[2]) {
            # Fixed columns
            if (IND[1]!="") TAGS = TAGS TAGS_sep TAGS_header "#family:" IND[1]
            if (IND[2]!="") TAGS = TAGS TAGS_sep TAGS_header "#individual:" IND[2]
            if (IND[3]!="") TAGS = TAGS TAGS_sep TAGS_header "#maternal:" IND[3]
            if (IND[4]!="") TAGS = TAGS TAGS_sep TAGS_header "#paternal:" IND[4]
            if (IND[5]!="") TAGS = TAGS TAGS_sep TAGS_header "#sex:" IND[5]
            if (IND[6]!="") TAGS = TAGS TAGS_sep TAGS_header "#phenotype:" IND[6]
            # Additional columns
            for (i in HEADER) { #print HEADER[i] "=" IND[i]
                if (i>6) {
                    if (IND[i]!="") {
                        # TAG column
                        if (tolower(HEADER[i])=="tag" || tolower(HEADER[i])=="#tag" ) {
                            gsub(" ", TAGS_sep, IND[i])
                            TAGS = TAGS TAGS_sep IND[i]
                        # Custom columns
                        } else {
                            TAGS = TAGS TAGS_sep TAGS_header "#" HEADER[i] ":" IND[i]
                        }
                    }
                }
            }
        }
    }
}
END {
    if (TAGS=="") {
        print TAGS_sep
    } else {
        print TAGS
    }
}