#!/bin/awk
#############################
# Normalize OutLyzer Output #
#############################
# Input: OuLyzer VCF output
# Output: VCF normalized
{
	if (match($1, /^##/)) {
		# Print header
		print $0
	} else if (match($1, /^#CHROM/)) {
		# Add header
		print "##FORMAT=<ID=GT,Number=.,Type=String,Description=\"Genotype\">"
		print "##FORMAT=<ID=DP,Number=.,Type=Integer,Description=\"Raw Depth\">"
		print "##FORMAT=<ID=FREQ,Number=.,Type=String,Description=\"Allele Frequency Percent\">"
		print "##FORMAT=<ID=VAF,Number=.,Type=Float,Description=\"Allele Frequency\">"
		# CHROM line
		print $0
	} else {
		# Split INFO field
		n_INFO=split($8,INFO,";");
		# Calculation GT
		GT="./.";
		FREQ=".";
		VAF=".";
		DP=".";
		for (i = 1; i <= n_INFO; ++i) {
			# Split Annotations
			n_ANN=split(INFO[i],ANN,"=");
			# Find AF annotation value
			if (ANN[1]=="AF") {
				FREQ=ANN[2];
				VAF=FREQ/100;
				# CalculatiON GT
				if (ANN[2]>80) {
					GT="1/1";
				} else {
					GT="0/1";
				}
			}
			# Find DP
			if (ANN[1]=="DP") {
				DP=ANN[2];
			};
		}
		# Print VCF first fields
		for (i = 1; i <= 8; ++i) {
			printf $i"\t"
		}
		# Print FORMAT and Sample format
		print "GT:DP:FREQ:VAF:"$9"\t"GT":"DP":"FREQ":"VAF":"$10
	}
}
