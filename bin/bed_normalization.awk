#!/bin/awk -f

# Normalize BED/GENES file
# chr	start	stop	strand	gene/target

{chr=$1}
{start=$2}
{stop=$3}
{strand=$4}
{gene=$5}
strand !~ /[+-]/ {strand="+"}
gene == "" { if ($4 !~ /[+-]/ && $4 != "") {gene=$4} else {gene=chr"_"start"_"stop} }
{print chr"\t"start"\t"stop"\t"strand"\t"gene}
