/^#:GATKTable:RecalTable0:/ { table=1; next }
/^#:GATKTable:/ && table { exit }
table && /^ReadGroup[[:space:]]+EventType/ { next }
table && NF { found=1 }
END { exit !found }