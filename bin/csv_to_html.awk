#!/bin/awk -f

# Set field separator as comma for csv and print the HTML header line
BEGIN {
    FS="\t";
	nb_line=0;
    if (limit=="") limit=2000;
    #print "<table>";
	#print "tag_row_head_plus="tag_row_head_plus;
	#print "tag_col_head_plus="tag_col_head_plus;
	#print "tag_row_body_plus="tag_row_body_plus;
	#print "tag_col_body_plus="tag_col_body_plus;
}
# Function to print a row with one argument to handle either a 'th' tag or 'td' tag
function printRow(tag,tag_row_plus,tag_col_plus) {
	if (tag=="th") print "<thead>";
    print "<tr "tag_row_plus">";
    for(i=1; i<=NF; i++) print "<"tag" "tag_col_plus">"$i"</"tag">";
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
