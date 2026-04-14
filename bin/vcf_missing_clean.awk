#!/usr/bin/env awk
BEGIN { FS=OFS="\t" }

/^#/ { print; next }

{
    # Split FORMAT
    split($9, fmt, ":")
    nfmt = length(fmt)

    # Find the index of GT
    gt_idx = 0
    for (i = 1; i <= nfmt; i++)
        if (fmt[i] == "GT") { gt_idx = i; break }
    if (gt_idx == 0) { print; next }

    # For each sample
    for (s = 10; s <= NF; s++) {
        split($s, f, ":")
        gt = f[gt_idx]

        # Replace targeted genotypes with missing
        if (gt ~ /^0[\/|]0$/ || gt ~ /^0[\/|]\.$/ || gt ~ /^\.[\/|]0$/ || gt ~ /^\.[\/|]\.$/) {
            sep = (gt ~ /\|/) ? "|" : "/"
            gt = "." sep "."
            # keep only GT for this sample
            $s = gt
        }
        # Otherwise, the sample remains unchanged
    }

    # The FORMAT ($9) remains unchanged
    print
}
