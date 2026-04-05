# AWK script to convert RefSeq gene annotation from NCBI into GTF format.
# For each gene, emits one gene line spanning all its transcripts (min txStart to
# max txEnd), followed by transcript, exon, and CDS lines for each transcript.
# Stop codons are removed from coding sequences.
# Duplicate transcript IDs are handled by appending a counter suffix.
#
# Usage: awk -F '\t' -v OFS='\t' -f refSeq_to_gtf.awk ncbiRefSeq.txt > refSeq_genes.gtf

function min(x, y) { return (x > y) ? y : x }
function max(x, y) { return (x < y) ? y : x }

# Return transcript_type based on transcript ID prefix or gene name
function transcript_type(tx_id, gene_name) {
	if (tx_id ~ /^NM_/)            return "protein_coding"
	if (gene_name ~ /^MIR[0-9]/)   return "miRNA"
	return "pseudogene"
}

{
	split($10, exon_starts, ",")
	split($11, exon_ends,   ",")
	split($16, exon_frames, ",")

	gene_name  = $13
	chrom      = $3
	strand     = $4
	txStart    = int($5)
	txEnd      = int($6)
	cdsS       = int($7)
	cdsE       = int($8)
	exon_cnt   = int($9)
	cds_s_stat = $14
	cds_e_stat = $15
	orig_tx    = $2

	# Remove stop codon from left end for coding genes on the minus strand
	if (strand == "-" && cds_s_stat == "cmpl" &&
	    (exon_starts[1] != cdsS ||
	     (min(exon_ends[1], cdsE) - exon_starts[1] + exon_frames[1]) % 3 == 0)) {
		cdsS += 3
		for (i in exon_ends)
			if (cdsS >= exon_ends[i] && cdsS <= exon_ends[i] + 2)
				cdsS += exon_starts[i+1] - exon_ends[i]
	}
	# Remove stop codon from right end for coding genes on the plus strand
	if (strand == "+" && cds_e_stat == "cmpl" &&
	    (exon_ends[exon_cnt] != cdsE ||
	     (exon_ends[exon_cnt] - max(exon_starts[exon_cnt], cdsS) + exon_frames[exon_cnt]) % 3 == 0)) {
		cdsE -= 3
		for (i in exon_starts)
			if (cdsE <= exon_starts[i] && cdsE >= exon_starts[i] - 2)
				cdsE -= exon_starts[i] - exon_ends[i-1]
	}

	# Track gene extent (min txStart, max txEnd) across all transcripts of the same gene
	if (!(gene_name in gene_min_start)) {
		gene_order[gene_count++] = gene_name
		gene_chrom[gene_name]     = chrom
		gene_strand[gene_name]    = strand
		gene_min_start[gene_name] = txStart
		gene_max_end[gene_name]   = txEnd
		gene_first_tx[gene_name]  = orig_tx
	} else {
		if (txStart < gene_min_start[gene_name]) gene_min_start[gene_name] = txStart
		if (txEnd   > gene_max_end[gene_name])   gene_max_end[gene_name]   = txEnd
	}

	# Handle duplicate transcript IDs by appending a counter suffix
	tx_seen[orig_tx]++
	tx_id   = orig_tx
	gid_eff = gene_name
	if (tx_seen[orig_tx] > 1) {
		suffix  = tx_seen[orig_tx] - 1
		tx_id   = orig_tx   "_" suffix
		gid_eff = gene_name "_" suffix
	}

	tx_type  = transcript_type(tx_id, gene_name)
	tx_extra = " transcript_type \"" tx_type "\"; transcript_status \"KNOWN\"; transcript_name \"" tx_id "\";"

	# Build transcript line
	out = chrom "\tRefSeq\ttranscript\t" txStart "\t" txEnd "\t.\t" strand "\t.\t" \
	      "gene_id \"" gid_eff "\"; transcript_id \"" tx_id "\"; gene_name \"" gene_name "\";" tx_extra "\n"

	# Build exon and (optionally) CDS lines
	for (i = 1; i <= exon_cnt; i++) {
		exon_num = (strand == "+") ? i : exon_cnt - i + 1
		attr = "gene_id \"" gid_eff "\"; transcript_id \"" tx_id \
		       "\"; exon_number \"" exon_num "\"; exon_id \"" tx_id "." exon_num \
		       "\"; gene_name \"" gene_name "\";" tx_extra
		out = out chrom "\tRefSeq\texon\t" (exon_starts[i] + 1) "\t" exon_ends[i] \
		          "\t.\t" strand "\t.\t" attr "\n"
		if (cds_s_stat ~ /cmpl/ && cdsS <= exon_ends[i] && cdsE >= exon_starts[i])
			out = out chrom "\tRefSeq\tCDS\t" (max(cdsS, exon_starts[i]) + 1) \
			          "\t" min(cdsE, exon_ends[i]) "\t.\t" strand "\t" exon_frames[i] "\t" attr "\n"
	}

	gene_lines[gene_name] = gene_lines[gene_name] out
}

END {
	# Print genes in the order they were first encountered
	for (g = 0; g < gene_count; g++) {
		gname = gene_order[g]
		print gene_chrom[gname], "RefSeq", "gene", \
		      gene_min_start[gname], gene_max_end[gname], \
		      ".", gene_strand[gname], ".", \
		      "gene_id \"" gname "\"; transcript_id \"" gene_first_tx[gname] "\"; gene_name \"" gname "\";"
		printf "%s", gene_lines[gname]
	}
}