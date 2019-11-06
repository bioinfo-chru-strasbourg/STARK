
#DIR=/STARK/data/SAMPLE/HORIZON
#FASTQ=$DIR/HORIZON.R1.fastq.gz

if ((0)); then
	DIR=/STARK/data/SAMPLE/TEST
	FASTQ=$DIR/TEST.R1b.fastq.gz
	FASTQ_R2=$DIR/TEST.R2.fastq.gz
	SAMPLE_NAME=TEST
	BC_PATTERN="NNNNNNNN"
fi;

if ((0)); then
	DIR=/STARK/data/SAMPLE/190830_M01656_0328_000000000-CLK2J
	FASTQ=$DIR/F10_S2_R1_001.fastq.gz
	FASTQ_R2=$DIR/F10_S2_R2_001.fastq.gz
	SAMPLE_NAME=F10
	BC_PATTERN="NNNNNNNN"
fi;

if ((0)); then
	DIR=/STARK/data/SAMPLE/190830_M01656_0328_000000000-CLK2J
	FASTQ=$DIR/F09_S1_R1_001.fastq.gz
	FASTQ_R2=$DIR/F09_S1_R2_001.fastq.gz
	SAMPLE_NAME=F09
	BC_PATTERN="NNNNNNNN"
fi;

if ((1)); then
	DIR=/STARK/data/SAMPLE/190830_M01656_0328_000000000-CLK2J
	FASTQ="$DIR/F09_S1_R1_001.fastq.gz,$DIR/F10_S2_R1_001.fastq.gz"
	FASTQ_R2="$DIR/F09_S1_R2_001.fastq.gz,$DIR/F10_S2_R2_001.fastq.gz"
	SAMPLE_NAME="F09,F10"
	BC_PATTERN="NNNNNNNN"
fi;




time if ((1)); then

	echo "# STARK analysis"
	#STARK --reads=$FASTQ --reads2=$FASTQ_R2 --analysis_name=ANALYSIS_UMI_2 --app=default+UMI+NO_ANNOTATION --verbose
	STARK --reads=$FASTQ --reads2=$FASTQ_R2 --analysis_name=ANALYSIS_UMI_9 --app=default+UMI+NO_ANNOTATION --verbose

	exit 0

fi;



DIR_OUTPUT=/STARK/data/umi_test.dir
FASTQ_UMI=$DIR_OUTPUT/umi.R1.fastq.gz
FASTQ_R2_UMI=$DIR_OUTPUT/umi.R2.fastq.gz

LOG_UMI=$DIR_OUTPUT/umi.log
SAM_OUT=$DIR_OUTPUT/out.sam
BAM_OUT=$DIR_OUTPUT/out.bam
BAM_GROUP_OUT=$DIR_OUTPUT/out.group.bam
BAM_GROUP_LOG=$DIR_OUTPUT/out.group.bam.log
BAM_GROUP_MX=$DIR_OUTPUT/out.group.bam.mx
BAM_DEDUP_OUT=$DIR_OUTPUT/out.dedup.bam
BAM_DEDUP_LOG=$DIR_OUTPUT/out.dedup.bam.log
BAM_DEDUP_MX=$DIR_OUTPUT/out.dedup.bam.mx
BAM_MARKDUP_OUT=$DIR_OUTPUT/out.markdup.bam
BAM_MARKDUP_LOG=$DIR_OUTPUT/out.markdup.bam.log
BAM_MARKDUP_MX=$DIR_OUTPUT/out.markdup.bam.mx

#FASTQ_OUT=$DIR_OUTPUT/out.fastq

PICARD=/STARK/tools/picard/current/bin/picard.jar
BWA=/STARK/tools/bwa/current/bin/bwa
SAMTOOLS=/STARK/tools/samtools/current/bin/samtools
GENOME=/STARK/databases/genomes/hg19/hg19.fa
GZ=gzip

TMP_DIR=$DIR_OUTPUT/tmp

mkdir -p $DIR_OUTPUT
mkdir -p $TMP_DIR



# FASTP
FASTP=/STARK/tools/fastp/current/bin/fastp
FASTQ_FASTP=$DIR_OUTPUT/fastp.R1.fastq.gz
FASTQ_R2_FASTP=$DIR_OUTPUT/fastp.R2.fastq.gz
FASTP_HTML=$DIR_OUTPUT/fastp.html
FASTP_JSON=$DIR_OUTPUT/fastp.json




if ((0)); then

	# Test if already UMI
	if ! (( $(zcat $FASTQ | head -n1 | awk -F" " '{n_read_name=split($1,read_name,":"); if (read_name[n_read_name] ~ /[A-Z]/) { print read_name[n_read_name]} }' | wc -l) )); then
		echo "NO UMI TRIM YET"
	else
		echo "UMI TRIM TODO"
	fi;

	if ! (( $(zcat $FASTQ_FASTP | head -n1 | awk -F" " '{n_read_name=split($1,read_name,":"); if (read_name[n_read_name] ~ /[A-Z]/) { print read_name[n_read_name]} }' | wc -l) )); then
		echo "NO UMI TRIM YET"
	else
		echo "UMI TRIM TODO"
	fi;
	#zcat $FASTQ_FASTP | head -n1 | awk -F" " '{n_read_name=split($1,read_name,":"); if (read_name[n_read_name] ~ /[A-Z]/) { print read_name[n_read_name]} }' | wc -l

	exit 0


fi;


if ((1)); then


	echo "## FASTP"

	time if ((1)); then

		# ./fastp -i /STARK/output/results/ANALYSIS_UMI_4/F09/F09.R1.fastq.gz -h /STARK/output/results/ANALYSIS_UMI_4/F09/F09.sequencing/fastp.single.html -j /STARK/output/results/ANALYSIS_UMI_4/F09/F09.sequencing/fastp.single.json --help
		FASTP_PARAM="--disable_trim_poly_g --disable_quality_filtering --disable_length_filtering --detect_adapter_for_pe "
		$FASTP -i $FASTQ -I $FASTQ_R2 \
			-o $FASTQ_FASTP -O $FASTQ_R2_FASTP \
			--umi --umi_loc per_read --umi_len ${#BC_PATTERN} \
			$FASTP_PARAM \
			--thread 3 \
			-h $FASTP_HTML -j $FASTP_JSON

		echo "### R1"
		zcat $FASTQ | wc -l
		zcat $FASTQ | head -n 8
		zcat $FASTQ_FASTP | wc -l
		zcat $FASTQ_FASTP | head -n 8
		echo "### R2"
		zcat $FASTQ_R2 | wc -l
		zcat $FASTQ_R2 | head -n 8
		zcat $FASTQ_R2_FASTP | wc -l
		zcat $FASTQ_R2_FASTP | head -n 8

		exit 0

	fi;



	# Alignment
	echo "## ALIGNMENT"

	time if ((0)); then

		$BWA mem -C -M -t 3 -R "@RG\tID:1\tPL:ILLUMINA\tPU:PU\tLB:001\tSM:$SAMPLE_NAME" $GENOME $FASTQ_FASTP $FASTQ_R2_FASTP | \
		sed 's/[1-2]\:N\:0\:\([A-Z]*\)/BC\:Z\:\1/' | \
		$SAMTOOLS view -b -1 -S -T $GENOME - -@ 3 | \
		$SAMTOOLS sort - -l 1 -O BAM -o $BAM_OUT -T $BAM_OUT.SAMTOOLS_PREFIX -@ 3
		$SAMTOOLS index $BAM_OUT

		echo $BAM_OUT
		$SAMTOOLS flagstat $BAM_OUT
		#$SAMTOOLS view -h $BAM_OUT | head -n 30

	fi;


	# UMI group
	echo "## UMI Group"

	time if ((0)); then

		#echo $BAM_OUT
		#$SAMTOOLS flagstat $BAM_OUT
		#$SAMTOOLS view -h $BAM_OUT | head -n 30

		umi_tools group -I $BAM_OUT --paired --unpaired-reads=use --unmapped-reads=use --chimeric-pairs=use \
		 --group-out=$GROUP_OUT \
		 --umi-separator=":" \
		 --output-bam -S $BAM_GROUP_OUT \
		 -L $BAM_GROUP_LOG \
		 --temp-dir=$TMP_DIR

		# --chimeric-pairs=discard --chimeric-pairs=use --ignore-umi

		echo $BAM_GROUP_OUT
		$SAMTOOLS flagstat $BAM_GROUP_OUT

	fi;

	time if ((0)); then

		#echo $BAM_OUT
		#$SAMTOOLS flagstat $BAM_OUT
		#$SAMTOOLS view -h $BAM_OUT | head -n 30

		#$SAMTOOLS view $BAM_OUT -c
		$SAMTOOLS view -h $BAM_OUT | /tool/bin/umi_read_name_to_BX.awk | $SAMTOOLS view -S -b > $BAM_GROUP_OUT

		# --chimeric-pairs=discard --chimeric-pairs=use --ignore-umi

		echo $BAM_GROUP_OUT
		$SAMTOOLS flagstat $BAM_GROUP_OUT

	fi;


	# MARKDuplicates
	echo "## MarkUMIDuplicates"

	time if ((0)); then
		java -jar $PICARD MarkDuplicates \
		      I=$BAM_GROUP_OUT \
		      O=$BAM_MARKDUP_OUT \
		      M=$BAM_MARKDUP_MX \
			  BARCODE_TAG=BX

		echo $BAM_MARKDUP_OUT
		$SAMTOOLS flagstat $BAM_MARKDUP_OUT

		$SAMTOOLS view $BAM_MARKDUP_OUT | head -n 10
		grep "^#" -v $BAM_MARKDUP_MX | column -t


	fi;




fi;


if ((1)); then

	echo "# MAPPED"
	$SAMTOOLS view $BAM_OUT | head -n 4
	echo "# GROUPED"
	$SAMTOOLS view $BAM_GROUP_OUT | grep "BX:Z" -c
	$SAMTOOLS view $BAM_GROUP_OUT | grep "^$" -v -c
	$SAMTOOLS view $BAM_GROUP_OUT | grep -o 'BX:Z:[^\t\n]*' | sed s/^BX:Z://gi | sort -u | wc -l
	$SAMTOOLS view $BAM_GROUP_OUT | head -n 4
	echo "# MARKDUP"
	$SAMTOOLS view $BAM_MARKDUP_OUT | head -n 4

	exit 0

fi;


exit 0

#$GZ -d -c $FASTQ | tr " " "|" | head -n 30
#$SAMTOOLS view -h $BAM_MARKDUP_OUT | head -n 30 | tr " " "|"
#exit 0

#java -jar $PICARD SamToFastq I=/STARK/data/SAMPLE/TEST/TEST.bam FASTQ=/STARK/data/SAMPLE/TEST/TEST.R1.fastq SECOND_END_FASTQ=/STARK/data/SAMPLE/TEST/TEST.R2.fastq

if ((0)); then

	# FASTQ
	echo $FASTQ
	echo "Number of reads R1:"
	$GZ -d -c $FASTQ | grep "^@" | wc -l
	$GZ -d -c $FASTQ | head -n 10

	# FASTQ R2
	echo $FASTQ_R2
	echo "Number of reads R2:"
	$GZ -d -c $FASTQ_R2 | grep "^@" | wc -l
	$GZ -d -c $FASTQ_R2 | head -n 10

fi;


if ((0)); then
	#samtools view -u -f 1 -F 12 $BAM_MARKDUP_OUT -c

	#bedtools bamtofastq -i $BAM_MARKDUP_OUT -fq $FASTQ.test.fastq -fq2 $FASTQ_R2.test.2.fastq
	#samtools bam2fq $BAM_MARKDUP_OUT > $FASTQ.test.fastq
	java -Xmx2g -jar $PICARD SamToFastq I=$BAM_MARKDUP_OUT F=$FASTQ.test.fastq F2=$FASTQ_R2.test.fastq UNPAIRED_FASTQ=$FASTQ_R2.U.test.fastq

	echo $FASTQ.test.fastq
	echo "Number of reads R1:"
	cat $FASTQ.test.fastq | grep "^@" | wc -l
	#cat $FASTQ.test.fastq | head -n 10

	echo $FASTQ_R2.test.fastq
	echo "Number of reads R2:"
	cat $FASTQ_R2.test.fastq | grep "^@" | wc -l
	#cat $FASTQ_R2.test.fastq | head -n 10

	echo $FASTQ_R2.U.test.fastq
	echo "Number of reads RU:"
	cat $FASTQ_R2.U.test.fastq | grep "^@" | wc -l
	#cat $FASTQ_R2.U.test.fastq | head -n 10


	exit 0

fi;

# Extract
echo "## EXTRACT UMI"

if ((0)); then # DO NOT USE (or for single read)

	# EXTRACT
	umi_tools extract --stdin=$FASTQ --bc-pattern=$BC_PATTERN --log=$LOG_UMI | awk '{
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
	}' | $GZ -c -1 > $FASTQ_UMI

	echo $FASTQ_UMI
	$GZ -d -c $FASTQ_UMI | head -n 10


	# EXTRACT
	umi_tools extract --stdin=$FASTQ_R2 --bc-pattern=$BC_PATTERN --log=$LOG_UMI | awk '{
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
	}' | $GZ -c -1 > $FASTQ_R2_UMI

	echo $FASTQ_R2_UMI
	$GZ -d -c $FASTQ_R2_UMI | head -n 10

fi;


time if ((1)); then # Paired reads


	if ((0)); then
		$GZ -d -c $FASTQ $FASTQ_R2 | $GZ -c > $FASTQ.R1R2.fastq.gz
		umi_tools extract --bc-pattern=$BC_PATTERN \
		-I $FASTQ.R1R2.fastq.gz  \
		--stdout=$FASTQ_UMI.tmp.gz \
		--reconcile-pairs --either-read-resolve=quality \
		--temp-dir=$TMP_DIR
	fi;

	if ((1)); then
		umi_tools extract --bc-pattern=$BC_PATTERN \
		-I $FASTQ --read2-in=$FASTQ_R2 \
		--stdout=$FASTQ_UMI.tmp.gz --read2-out=$FASTQ_R2_UMI.tmp.gz \
		--reconcile-pairs --either-read-resolve=quality -L $LOG_UMI \
		--temp-dir=$TMP_DIR
	fi;

	if ((0)); then # error
		umi_tools extract --bc-pattern=$BC_PATTERN \
		-I $FASTQ  \
		--stdout=$FASTQ_UMI.tmp.gz \
		--reconcile-pairs --either-read-resolve=quality \
		--temp-dir=$TMP_DIR
		umi_tools extract --bc-pattern=$BC_PATTERN \
		-I $FASTQ_R2 \
		--stdout=$FASTQ_R2_UMI.tmp.gz \
		--reconcile-pairs --either-read-resolve=quality \
		--temp-dir=$TMP_DIR
	fi;

	# --reconcile-pairs

	echo "Number of reads R1 & R2:"
	$GZ -d -c $FASTQ_UMI.tmp.gz | grep "^@" | wc -l
	$GZ -d -c $FASTQ_R2_UMI.tmp.gz | grep "^@" | wc -l


	zcat $FASTQ_UMI.tmp.gz | /tool/bin/umi_read_name_sep.awk | $GZ -c -1 > $FASTQ_UMI

	echo $FASTQ_UMI
	echo "Number of reads:"
	$GZ -d -c $FASTQ_UMI | grep "^@" | wc -l
	$GZ -d -c $FASTQ_UMI | head -n 10

	zcat $FASTQ_R2_UMI.tmp.gz | /tool/bin/umi_read_name_sep.awk | $GZ -c -1 > $FASTQ_R2_UMI

	echo $FASTQ_R2_UMI
	echo "Number of reads:"
	$GZ -d -c $FASTQ_R2_UMI | grep "^@" | wc -l
	$GZ -d -c $FASTQ_R2_UMI | head -n 10

fi;



# Alignment
echo "## ALIGNMENT"

time if ((1)); then

	$BWA mem -C -M -t 3 -R "@RG\tID:1\tPL:ILLUMINA\tPU:PU\tLB:001\tSM:$SAMPLE_NAME" $GENOME $FASTQ_UMI $FASTQ_R2_UMI | \
	sed 's/[1-2]\:N\:0\:\([A-Z]*\)/BC\:Z\:\1/' | \
	$SAMTOOLS view -b -1 -S -T $GENOME - -@ 3 | \
	$SAMTOOLS sort - -l 1 -O BAM -o $BAM_OUT -T $BAM_OUT.SAMTOOLS_PREFIX -@ 3
	$SAMTOOLS index $BAM_OUT

	echo $BAM_OUT
	$SAMTOOLS flagstat $BAM_OUT
	#$SAMTOOLS view -h $BAM_OUT | head -n 30

fi;


# UMI group
echo "## UMI Group"

time if ((1)); then

	#echo $BAM_OUT
	#$SAMTOOLS flagstat $BAM_OUT
	#$SAMTOOLS view -h $BAM_OUT | head -n 30

	umi_tools group -I $BAM_OUT --paired --unpaired-reads=use --unmapped-reads=use --chimeric-pairs=use \
	 --group-out=$GROUP_OUT \
	 --umi-separator=":" \
	 --output-bam -S $BAM_GROUP_OUT \
	 -L $BAM_GROUP_LOG \
	 --temp-dir=$TMP_DIR

	# --chimeric-pairs=discard --chimeric-pairs=use --ignore-umi

	echo $BAM_GROUP_OUT
	$SAMTOOLS flagstat $BAM_GROUP_OUT

fi;


# MARKDuplicates
echo "## MarkUMIDuplicates"

time if ((1)); then
	java -jar $PICARD MarkDuplicates \
	      I=$BAM_GROUP_OUT \
	      O=$BAM_MARKDUP_OUT \
	      M=$BAM_MARKDUP_MX \
		  BARCODE_TAG=BX

	echo $BAM_MARKDUP_OUT
	$SAMTOOLS flagstat $BAM_MARKDUP_OUT

	$SAMTOOLS view $BAM_MARKDUP_OUT | head -n 10
	grep "^#" -v $BAM_MARKDUP_MX | column -t


fi;



time if ((1)); then

	echo "## resume"

	$SAMTOOLS view $BAM_OUT -c
	$SAMTOOLS view -u -f 1 -F 12 $BAM_OUT -c
	$GZ -d -c $FASTQ | grep "^@" | wc -l
	$GZ -d -c $FASTQ_R2 | grep "^@" |wc -l
	$GZ -d -c $FASTQ_UMI | grep "^@" |wc -l
	$GZ -d -c $FASTQ_R2_UMI | grep "^@" |wc -l

	echo "# number of molecules"
	samtools view $BAM_MARKDUP_OUT | grep -o 'BX:Z:[^\t\n]*' | sed s/^BX:Z://gi | sort -u | wc -l

fi;






exit 0




#
