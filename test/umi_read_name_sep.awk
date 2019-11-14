#!/bin/awk -f

# Modify read name header to change umi sep "_" to ":"
{
	if (NR%4==1) {
		n_head=split($0,head," ");
		n_read=split(head[1],read,"_");
		gsub("_",":", head[1])
		printf head[1];
		for (i = 2; i <= n_head; ++i) {
			printf " "head[i];
		}
		print "";
	}
	else
	{
		print $0
	}
}
