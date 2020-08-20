#!/bin/awk -f

# Normalize BED/GENES file
# chr	start	stop	gene/target score   strand	...

function join(array, sep, start, end, result, i)
{
    # sep
    if (sep == "")
       sep = "\t"
    else if (sep == SUBSEP) # magic value
       sep = ""
    
    # start
    if (length(start) == 0)
        start = 1

    # start
    if (length(end) == 0)
        end = 1

    

    # result init
    result = array[start]

    for (i = start + 1; i <= end; i++)
        #print i "=>" array[i]
        result = result sep array[i]
    return result
}

{
    # original number of columns
    nb_col=split($0,line,"\t");
    
    # fields
    chr=line[1]
    start=line[2]
    stop=line[3]
    name=line[4]
    score=line[5]
    strand=line[6]

    # Normalizatino
    if (name == "") {line[4]=chr"_"start"_"stop}
    #if (score == "") {line[5]="0"}
    #if (strand !~ /[+-]/) {line[6]="+"}

    # Adjust number of fields
    if (nb_col<4) {nb_col=4}

    # Prin new line
    print join(line,"\t",1,nb_col)

}
