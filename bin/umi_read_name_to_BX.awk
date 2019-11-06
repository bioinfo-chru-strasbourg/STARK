#!/bin/awk -f

# Modify read name header to change umi sep "_" to ":"
{
	if (substr($0,1,1) != "@" ) {
		n_read=split($0,read,"\t");
		n_name=split(read[1],name,":");
		umi=name[n_name]
		print $0"\tBX:Z:"umi;
	}
	else
	{
		print $0
	}
}
