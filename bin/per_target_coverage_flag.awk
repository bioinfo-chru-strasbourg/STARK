#!/bin/awk
############################
# FLAG per target coverage #
############################
# Input: HsMetrics file from Picard
# EXPECTED_DEPTH
# MINIMUM_DEPTH
BEGIN {
	if (EXPECTED_DEPTH=="") {EXPECTED_DEPTH=100}
	if (MINIMUM_DEPTH=="") {MINIMUM_DEPTH=30}
}

{

	# INIT
	split($0,columns,"\t")
	before=""
	after=""
	flag=""
	average=columns[7]
	minimum=columns[11]
	average_flag="PASS"
	minimum_flag="PASS"

	# BEFORE
	for (i=1; i<6; i++) {
		before=before""columns[i]"\t"
	}
	for (i=6; i<=NF; i++) {
		after=after"\t"columns[i]
	}

	# FLAG
	if (NR==1) {
		flag="average_threshold\tminimum_threshold"
	} else {
		if (average<EXPECTED_DEPTH) {average_flag="WARN"}
		if (minimum<EXPECTED_DEPTH) {minimum_flag="WARN"}
		if (average<MINIMUM_DEPTH) {average_flag="FAIL"}
		if (minimum<MINIMUM_DEPTH) {minimum_flag="FAIL"}
		if (average==0) {average_flag="MISS"}
		if (minimum==0) {minimum_flag="MISS"}
		flag=average_flag"\t"minimum_flag
	}

	# PRINT
	print before""flag""after

}
