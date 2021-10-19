###
# Clean fastq header
# reformat UMI ("-" as normalized separator)
# reformat FASTQ comment TODO
# parameters:
#    NO_COMMENT: remove fastq comment. default 0
#    SAM_TAG: reformat FASTQ comment to SAM TAG. Will remove comment not in SAM TAG format. Move CASAVA TAG as BC TAG. default 0
#    UMI_REFORMAT: Reformat UMI on read nam column 8. default 0
#    UMI_TAG: Add/replace UMI SAM TAG "RX:Z:XXX". default 0
#    UMI_SEP: UMI separator to use. default "-"
#    UMI_LOC: defined UMI location, if is *index* (e.g. index1, index2, per_index), UMI will be removed from BC SAM TAG tag (CASAVA moved to BC if SAM_TAG=1). default ""
#    INPUT_FORMAT: either "FASTQ" or "BAM". default "FASTQ"
BEGIN {
	if (NO_COMMENT=="") {
		NO_COMMENT=0
	}
	if (SAM_TAG=="") {
		SAM_TAG=0
	}
	if (UMI_REFORMAT=="") {
		UMI_REFORMAT=0
	}
	if (UMI_TAG=="") {
		UMI_TAG=0
	}
	if (UMI_SEP=="") {
		UMI_SEP="-"
	}
	if (REMOVE_UMI_FROM_BC=="") {
		REMOVE_UMI_FROM_BC=0
	}
	#REMOVE_UMI_FROM_BC=0
	if (match(UMI_LOC"", /.*index.*/)) {
		REMOVE_UMI_FROM_BC=1
	}
	if (INPUT_FORMAT=="") {
		INPUT_FORMAT="FASTQ"
	}
	COMMENT_SEP=" "
	FIRST_SEP=" "
	if (INPUT_FORMAT=="BAM") {
		NO_COMMENT=0
		SAM_TAG=0
		COMMENT_SEP="\t"
		FIRST_SEP="\t"
	}
}
{
	#if ( ( INPUT_FORMAT=="FASTQ" && match($0, /^@/) ) \
	if ( ( INPUT_FORMAT=="FASTQ" && NR%4 == 1 ) \
		|| \
		( INPUT_FORMAT=="BAM" && ! match($0, /^@/) ) \
		) {
		nb_a=split($0,a,"[ \t]");
		nb_b=split(a[1],b,"[:]");
		UMI="";
		# UMI field in 7
		if (UMI_REFORMAT) {
			if (nb_b==7 && b[7] ~ /[a-zA-Z]/) {
				# check sep
				sep1_index=index(b[7],"_");
				sep2_index=index(b[7],"-");
				# if one sep is found
				if ((sep1_index+sep2_index)) {
					# find first sep
					if (sep1_index==0) {sep_first="-"}
					else if (sep2_index==0) {sep_first="_"}
					else if (sep1_index>sep2_index) {sep_first="-"} else {sep_first="_"}
					# split 
					nb_c=split(b[7],c,sep_first);
					# recreate column 7
					b[7]=c[1];
					# create column 8 UMI
					if (nb_c>1) {
						b[8] = sep = "";
						for (i=2; i in c; i++) {
							gsub("[^A-Za-z]", UMI_SEP, c[i]);
							b[8] = b[8] sep c[i];
							sep = UMI_SEP;
						};
					};
				}
			};
			# UMI field
			if (nb_b>7) {
				gsub("[^A-Za-z]", UMI_SEP, b[8]);
			};
		}
		if (nb_b>7) {
			UMI=b[8];
		};
		# Recreate read name
		bjoined = sep = "";
		for (i=1; i in b; i++) {
			bjoined = bjoined sep b[i];
			sep = ":";
		};
		# Recreate read header and Reformat FASTQ read comment
		a[1]=bjoined;
		ajoined = a[1];
		if (!NO_COMMENT) {
			sep = FIRST_SEP;
			UMI_already_exists_in_comment=0;
			BC_already_exists_in_comment=0;
			if (match($0, /BC:Z:[^:]*/)) {
				BC_already_exists_in_comment=1;
			}
			nb_sam_tag=0;
			if (nb_a>1) {
				for (i=2; i in a; i++) {
					# Check comment format to SAM specification
					if (match(a[i], /^[^:]*:N:0:[^:]*$/)) {
						# BC - reformat x:N:0:xxxx
						if (SAM_TAG) {
							nb_CASAVA=split(a[i],CASAVA,"[:]");
							BC=CASAVA[4]
							gsub("+","-",BC)
							# check if UMI
							if (UMI!="" && REMOVE_UMI_FROM_BC) {
								UMI_AS_BC=UMI
								#gsub("-","+",UMI_AS_BC)
								if (UMI_AS_BC == BC) { BC=""}
								gsub(UMI_AS_BC,"",BC)
								#gsub("^+","",BC)
								#gsub("+$","",BC)
								gsub("^-","",BC)
								gsub("-$","",BC)
							}
							if (BC!="" && ! BC_already_exists_in_comment) {
								a[i]="BC:Z:"BC
							} else {
								a[i]=""
							}
						}
					} else if (match(a[i], /^BC:Z:[^:]*$/)) {
						# BC: if UMI already exists and need to be removed (see UMI_LOC)
						nb_BC=split(a[i],BCZ,"[:]");
						BC=BCZ[3]
						gsub("+","-",BC)
						if (UMI!="" && REMOVE_UMI_FROM_BC) {
							UMI_AS_BC=UMI
							#gsub("-","+",UMI_AS_BC)
							if (UMI_AS_BC == BC) { BC=""}
							gsub(UMI_AS_BC,"",BC)
							#gsub("^+","",BC)
							#gsub("+$","",BC)
							gsub("^-","",BC)
							gsub("-$","",BC)
						}
						if (BC!="") {
							a[i]="BC:Z:"BC
						} else {
							a[i]=""
						}
					} else if (match(a[i], /^RX:Z:[^:]*$/)) {
						# RX: if UMI already exists and UMI_TAG needed, remove if another is previously detected, else keep it
						UMI_already_exists_in_comment=1;
						if (UMI!="" && UMI_TAG) {
							a[i]="";
						}
					} else if (match(a[i], /^[^:]*:[^:]*:[^:]*$/)) {
						# apparently in SAM format
					} else if (SAM_TAG) {
						a[i]=""
					} else  {
						# nothing to do
					};
					# join
					if (a[i]!="") {
						nb_sam_tag++;
						ajoined = ajoined sep a[i];
						if (SAM_TAG) {
							sep = "\t";
						} else {
							sep = COMMENT_SEP;
						}	
					}
				};
			};
			# Add UMI if exists
			if (UMI!="" && UMI_TAG) {
				a[nb_a]= sep "RX:Z:" UMI
				ajoined = ajoined sep "RX:Z:" UMI
			}
		}
		print ajoined
	} else {
		print $0
	}
}