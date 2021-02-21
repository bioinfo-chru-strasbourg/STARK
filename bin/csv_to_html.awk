#!/bin/awk -f

function str2map(str,fs1,fs2,map,   n,tmp) {
   n=split(str,map,fs1)
   for (;n>0;n--) {
     split(map[n],tmp,fs2);
     map[tmp[1]]=tmp[2]; map[tmp[1],"full"]=map[n]
     delete map[n]
   }
}

# Set field separator as comma for csv and print the HTML header line
BEGIN {
    FS="\t";
    nb_line=0;
    if (limit=="") limit=2000;
    if (limit_string=="") limit_string=50;
    #header_replace_by_space="[_]";
    str2map(special_value_color,"|","=",special_value_color_array)
}

# Function to print a row with one argument to handle either a 'th' tag or 'td' tag
function printRow(tag,tag_row_plus,tag_col_plus) {
    if (tag=="th") print "<thead>";
    cells=""
    for(i=1; i<=NF; i++) {
        cell_value=$i;
        cell_title=cell_value;
        cell_show=cell_value;
        gsub(/<[^>]+>/,"",cell_title);
        if (length(cell_show)>limit_string && tag=="td") {
            cell_show=substr(cell_show,1,limit_string)"...";
        };
        if (tag=="th" && header_replace_by_space!="") {
             gsub(header_replace_by_space, " ", cell_show)
        }
        cells=cells "<"tag" "tag_col_plus" title='"cell_title"' style='color:"special_value_color_array[cell_value]"'>"cell_show"</"tag"> ";
    };
    print "<tr "tag_row_plus">";
    print cells;
    print "</tr>"
    if (tag=="th") print "</thead>";
}
# If CSV file line number (NR variable) is 1, call printRow fucntion with 'th' as argument
NR==1 {
    printRow("th",tag_row_head_plus,tag_col_head_plus)
}
# If CSV file line number (NR variable) is greater than 1, call printRow fucntion with 'td' as argument
NR==2 {
    print "<tbody>"
}
# If CSV file line number (NR variable) is greater than 1, call printRow fucntion with 'td' as argument
NR>1 {
        nb_line++;
    if (nb_line<=limit) printRow("td",tag_row_body_plus,tag_col_body_plus);
}
# Print HTML footer
END {
        if (nb_line>1) {print "</tbody>"}
}