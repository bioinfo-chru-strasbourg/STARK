#!/usr/bin/awk -f
###############################################################################
# outlyzer_norm.awk - Normalize a single-sample VCF produced by outLyzer
###############################################################################
#
# WHAT THE SCRIPT DOES
# --------------------
# outLyzer stores variant quality information in the INFO field (AF, DP, ...)
# while the sample FORMAT only contains BAL, GT, WTbal, MSS, and its GT is
# always 0/1. This script cleans and completes the VCF so that downstream
# tools (bcftools filters, annotation, STARK...) can use it:
#
#   1. GT is RECOMPUTED from INFO/AF (a percentage, e.g. 54.884):
#        AF >  HOM_AF  -> 1/1 (homozygous)
#        AF <= HOM_AF  -> 0/1 (heterozygous)
#      HOM_AF is an input variable (default 80, see USAGE). GT is replaced in
#      place. If AF is missing or multi-allelic, the original GT is kept.
#
#   2. Values equal to "None" are replaced by "." (VCF missing value), in the
#      INFO values and in the FORMAT sample values. This also applies to
#      elements of lists such as "None,1.2" or "None/None".
#
#   3. New FORMAT keys are appended at the END of FORMAT:
#        DP     <- INFO/DP   (raw depth)
#        FREQ   <- INFO/AF   (raw value, in %, e.g. 11.429)
#        VAF    <- INFO/AF / 100 (e.g. 0.11429)
#        QUAL   <- QUAL column
#        FILTER <- FILTER column
#      (the QUAL and FILTER columns themselves are not modified)
#
#   4. (deprecated) A new INFO field VARTYPE (one value per ALT allele) is appended at the
#      end of INFO, computed from REF and ALT:
#        SNV   single base substitution (A>G)
#        MNV   substitution of several bases, same length (AT>GC)
#        INS   pure insertion, REF is a prefix/suffix of ALT (A>ATG)
#        DEL   pure deletion, ALT is a prefix/suffix of REF (ATG>A)
#        OTHER anything else (complex delins, symbolic ALT <DEL>, "*", ".")
#
# The other FORMAT keys and sample values (BAL, WTbal, MSS...) are kept in
# their original order. If DP, FREQ, VAF, QUAL or FILTER already exist in the
# input FORMAT, the old key/value is removed and replaced by the new one
# (avoids duplicated keys, which make the VCF invalid). Likewise an existing
# INFO/VARTYPE is replaced (deprecated). Missing ##INFO / ##FORMAT header lines are added.
#
# Example (HOM_AF=80):
#   before: chr7 55249063 . G A 29.39 PASS AF=54.884;DP=1423;SDQ=None
#           GT:BAL:WTbal:MSS   0/1:391/390:320/322:None
#   after : chr7 55249063 . G A 29.39 PASS AF=54.884;DP=1423;SDQ=.
#           GT:BAL:WTbal:MSS:DP:FREQ:VAF:QUAL:FILTER
#           0/1:391/390:320/322:.:1423:54.884:0.54884:29.39:PASS
#
# NOTES
# -----
#  - Single-sample VCF only (10 columns): any other line aborts with an error.
#  - For multi-allelic sites, normalizing first (bcftools norm -m-) is
#    recommended: VAF is then computed per allele, but GT is not recomputed.
#
# USAGE
# -----
#   awk -f outlyzer_norm.awk in.vcf > out.vcf
#   awk -v HOM_AF=85 -f outlyzer_norm.awk in.vcf > out.vcf   # other threshold
#
# Variable: HOM_AF (default 80) = AF threshold in % above which GT = 1/1
###############################################################################

# Replace every "None" value by "." in a string. "None" is only replaced when it
# is a whole element, i.e. delimited by the string ends or by "," or "/".
function clean(v,    st) {
	while (match(v, /(^|[,\/])None($|[,\/])/)) {
		st = (substr(v, RSTART, 1) == "N") ? RSTART : RSTART + 1
		v = substr(v, 1, st - 1) "." substr(v, st + 4)
	}
	return v
}

# # Variant type of one ALT allele, given REF
# function vartype(ref, alt,    lr, la) {
# 	ref = toupper(ref); alt = toupper(alt)
# 	if (alt == "." || alt == "*" || alt ~ /[<>\[\]]/ || ref == alt) return "OTHER"
# 	lr = length(ref); la = length(alt)
# 	if (lr == 1 && la == 1) return "SNV"
# 	if (lr == la) return "MNV"
# 	if (lr < la && (substr(alt, 1, lr) == ref || substr(alt, la - lr + 1) == ref)) return "INS"
# 	if (lr > la && (substr(ref, 1, la) == alt || substr(ref, lr - la + 1) == alt)) return "DEL"
# 	return "OTHER"                          # complex event (e.g. delins)
# }

BEGIN {
	FS = OFS = "\t"                         # VCF fields are strictly tab-separated
	if (HOM_AF == "") HOM_AF = 80           # homozygous threshold (in %)

	# New FORMAT keys, in output order (appended at the end of FORMAT)
	n_new = split("DP FREQ VAF QUAL FILTER", NEWKEY, " ")
	for (k = 1; k <= n_new; k++) ISNEW[NEWKEY[k]] = 1

	# Header declarations, written only if not already present in the input
	FMT_HDR["GT"]     = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
	FMT_HDR["DP"]     = "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">"
	FMT_HDR["FREQ"]   = "##FORMAT=<ID=FREQ,Number=1,Type=String,Description=\"Allele Frequency Percent\">"
	FMT_HDR["VAF"]    = "##FORMAT=<ID=VAF,Number=1,Type=Float,Description=\"Allele Frequency\">"
	FMT_HDR["QUAL"]   = "##FORMAT=<ID=QUAL,Number=1,Type=Float,Description=\"Variant quality (copy of the QUAL column)\">"
	FMT_HDR["FILTER"] = "##FORMAT=<ID=FILTER,Number=1,Type=String,Description=\"Variant filter (copy of the FILTER column)\">"
	#INFO_HDR_VARTYPE  = "##INFO=<ID=VARTYPE,Number=A,Type=String,Description=\"Variant type for each ALT allele: SNV, INS, DEL, MNV or OTHER\">"
	n_fmt = split("GT DP FREQ VAF QUAL FILTER", FMT_ORDER, " ")
}

# --- 1. Meta-information lines: copied unchanged ------------------------------
/^##/ {
	# Remember which ##FORMAT / ##INFO IDs are already declared
	if (match($0, /^##FORMAT=<ID=[^,>]+/))
		FMT_DEFINED[substr($0, 14, RLENGTH - 13)] = 1
	else if (match($0, /^##INFO=<ID=[^,>]+/))
		INFO_DEFINED[substr($0, 12, RLENGTH - 11)] = 1
	print
	next
}

# --- 2. #CHROM line: first add the missing header declarations ----------------
/^#CHROM/ {
	#if (!("VARTYPE" in INFO_DEFINED)) print INFO_HDR_VARTYPE
	for (k = 1; k <= n_fmt; k++)
		if (!(FMT_ORDER[k] in FMT_DEFINED)) print FMT_HDR[FMT_ORDER[k]]
	print
	next
}

# --- 3. Variant lines ---------------------------------------------------------
{
	if (NF != 10) {
		printf "outlyzer_norm.awk: line %d: %d columns (10 expected, single-sample VCF)\n", NR, NF > "/dev/stderr"
		exit 1
	}

	# 3a. Parse INFO: replace "None" by ".", extract AF and DP,
	#     drop any pre-existing VARTYPE (recomputed below) - deprecated
	FREQ = "."; VAF = "."; DP = "."; GT_NEW = ""
	new_info = ""
	n_info = split($8, INFO, ";")
	for (i = 1; i <= n_info; i++) {
		if (INFO[i] == "" || INFO[i] == ".") continue
		eq = index(INFO[i], "=")
		if (eq == 0) {                              # flag without value: keep
			item = INFO[i]
		} else {
			key = substr(INFO[i], 1, eq - 1)
			val = clean(substr(INFO[i], eq + 1))
			# if (key == "VARTYPE") continue
			item = key "=" val

			if (key == "AF" && val != "." && val != "") {
				FREQ = val + "%"                    # raw value, in %
				n_af = split(val, AFS, ",")
				VAF = ""
				for (j = 1; j <= n_af; j++)         # one value per ALT allele
					VAF = VAF (j > 1 ? "," : "") sprintf("%.6g", AFS[j] / 100)
				# Genotype from AF (single-allelic only): > HOM_AF -> 1/1, else 0/1
				if (n_af == 1)
					GT_NEW = (AFS[1] + 0 > HOM_AF + 0) ? "1/1" : "0/1"
			} else if (key == "DP" && val != ".") {
				DP = val
			}
		}
		new_info = new_info (new_info == "" ? "" : ";") item
	}

	# # 3b. VARTYPE: one value per ALT allele, comma-separated (Number=A)
	# n_alt = split($5, ALTS, ",")
	# vt = ""
	# for (j = 1; j <= n_alt; j++)
	# 	vt = vt (j > 1 ? "," : "") vartype($4, ALTS[j])
	# $8 = new_info (new_info == "" ? "" : ";") "VARTYPE=" vt

	# 3c. Rebuild FORMAT and sample columns: original keys/values first
	#     (GT recomputed, "None" -> "."), then the new keys at the end
	n_key = split($9,  FK, ":")
	n_val = split($10, FV, ":")
	keys = ""; vals = ""
	for (j = 1; j <= n_key; j++) {
		if (FK[j] in ISNEW) continue                # avoid duplicated keys
		# trailing sample values may be omitted in a VCF: fill with "."
		v = (j <= n_val && FV[j] != "") ? clean(FV[j]) : "."
		if (FK[j] == "GT" && GT_NEW != "") v = GT_NEW
		keys = keys (keys == "" ? "" : ":") FK[j]
		vals = vals (vals == "" ? "" : ":") v
	}
	q = clean($6); if (q == "") q = "."
	f = clean($7); if (f == "") f = "."
	$9  = keys ":DP:FREQ:VAF:QUAL:FILTER"
	$10 = vals ":" DP ":" FREQ ":" VAF ":" q ":" f
	print                                           # $0 rebuilt with OFS = tab
}