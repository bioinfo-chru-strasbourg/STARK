
#DIR=/STARK/data/SAMPLE/HORIZON
#FASTQ=$DIR/HORIZON.R1.fastq.gz

#DIR=/STARK/data/SAMPLE/TEST
#FASTQ=$DIR/TEST.R1b.fastq.gz
#FASTQ_R2=$DIR/TEST.R2.fastq.gz
#SAMPLE_NAME=TEST
#BC_PATTERN="NNNNNNNN"


DIR=/STARK/data/SAMPLE/190830_M01656_0328_000000000-CLK2J
FASTQ=$DIR/F10_S2_R1_001.fastq.gz
FASTQ_R2=$DIR/F10_S2_R2_001.fastq.gz
SAMPLE_NAME=F10
BC_PATTERN="NNNNNNNN"


DIR_OUTPUT=/STARK/data/umi_test.dir
FASTQ_UMI=$DIR_OUTPUT/umi.R1.fastq.gz
FASTQ_R2_UMI=$DIR_OUTPUT/umi.R2.fastq.gz

LOG_UMI=$DIR_OUTPUT/umi.log
SAM_OUT=$DIR_OUTPUT/out.sam
BAM_OUT=$DIR_OUTPUT/out.bam
BAM_GROUP_OUT=$DIR_OUTPUT/out.group.bam
BAM_DEDUP_OUT=$DIR_OUTPUT/out.dedup.bam
BAM_MARKDUP_OUT=$DIR_OUTPUT/out.markdup.bam

FASTQ_OUT=$DIR_OUTPUT/out.fastq
GROUP_OUT=$DIR_OUTPUT/groups.tsv

PICARD=/STARK/tools/picard/current/bin/picard.jar
BWA=/STARK/tools/bwa/current/bin/bwa
SAMTOOLS=/STARK/tools/samtools/current/bin/samtools
GENOME=/STARK/databases/genomes/hg19/hg19.fa
GZ=gzip

mkdir -p $DIR_OUTPUT


#$GZ -d -c $FASTQ | tr " " "|" | head -n 30
#$SAMTOOLS view -h $BAM_MARKDUP_OUT | head -n 30 | tr " " "|"
#exit 0

#java -jar $PICARD SamToFastq I=/STARK/data/SAMPLE/TEST/TEST.bam FASTQ=/STARK/data/SAMPLE/TEST/TEST.R1.fastq SECOND_END_FASTQ=/STARK/data/SAMPLE/TEST/TEST.R2.fastq

# FASTQ
echo $FASTQ
$GZ -d -c $FASTQ | head -n 10

# FASTQ R2
echo $FASTQ_R2
$GZ -d -c $FASTQ_R2 | head -n 10

# Extract
echo "## EXTRACT UMI"

time if ((1)); then

	# EXTRACT
	umi_tools extract --stdin=$FASTQ --bc-pattern=$BC_PATTERN --log=$LOG_UMI | awk '{
		if (NR%4==1) {
			n_head=split($0,head," ");
			n_read=split(head[1],read,"_");
			read_name=read[1];
			read_umi=read[2];
			printf read_name;
			if (n_head>=2) {
				for (i = 2; i <= n_head; ++i) {
					printf " "head[i];
				}
				printf "\t"
			} else {
				printf " "
			}
			#r2=((NR-1)/4)%10;
			printf "BX:Z:"read_umi;
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
	umi_tools extract --stdin=$FASTQ_R2 --bc-pattern=NNNNNNNN --log=$LOG_UMI | awk '{
		if (NR%4==1) {
			n_head=split($0,head," ");
			n_read=split(head[1],read,"_");
			read_name=read[1];
			read_umi=read[2];
			printf read_name;
			if (n_head>=2) {
				for (i = 2; i <= n_head; ++i) {
					printf " "head[i];
				}
				printf "\t"
			} else {
				printf " "
			}
			#r2=((NR-1)/4)%10;
			printf "BX:Z:"read_umi;
			print "";
		}
		else
		{
			print $0
		}
	}' | $GZ -c -1 > $FASTQ_R2_UMI

	echo $FASTQ_R2_UMI
	$GZ -d -c $FASTQ_R2_UMI | head -n 20

fi;



# Alignment
echo "## ALIGNMENT"

time if ((1)); then

	$BWA mem -C -M -t 3 -R "@RG\tID:1\tPL:ILLUMINA\tPU:PU\tLB:001\tSM:$SAMPLE_NAME" $GENOME $FASTQ_UMI $FASTQ_R2_UMI | \
	sed 's/[1-2]\:N\:0\:\([A-Z]*\)/BC\:Z\:\1/' | \
	$SAMTOOLS view -b -1 -S -T $GENOME - -@ 3 | \
	$SAMTOOLS sort - -l 1 -O BAM -o $BAM_OUT -T $BAM_OUT.SAMTOOLS_PREFIX -@ 3

	#echo $BAM_OUT
	#$SAMTOOLS view -h $BAM_OUT | head -n 30

fi;



# MARKDuplicates
echo "## MarkUMIDuplicates"

time if ((1)); then
	java -jar $PICARD MarkDuplicates \
	      I=$BAM_OUT \
	      O=$BAM_MARKDUP_OUT \
	      M=$BAM_MARKDUP_OUT.mx \
		  BARCODE_TAG=BX

	$SAMTOOLS view $BAM_MARKDUP_OUT | head -n 10
	grep "^#" -v $BAM_MARKDUP_OUT.mx | column -t

	$SAMTOOLS flagstat $BAM_MARKDUP_OUT

fi;








exit 0




if ((0)); then

	time $GZ -d -c $FASTQ | awk -v UMI=8 '{
		if (NR%4==1) {
			head=$0
		}
		if (NR%4==2) {
			seq=$0
			umi=substr($0,0,UMI)
			new_seq=substr($0,UMI+1)
		}
		if (NR%4==3) {
			three=$0
		}
		if (NR%4==0) {
			qual=$0
			new_qual=substr($0,UMI+1)
			print head"\tBX:Z:"umi;
			print new_seq;
			print three;
			print new_qual;
		}

	}' | $GZ -c -1 > $FASTQ_UMI

	$GZ -d -c $FASTQ_UMI | head -n 10


	time $GZ -d -c $FASTQ_R2 | awk -v UMI=8 '{
		if (NR%4==1) {
			head=$0
		}
		if (NR%4==2) {
			seq=$0
			umi=substr($0,0,UMI)
			new_seq=substr($0,UMI+1)
		}
		if (NR%4==3) {
			three=$0
		}
		if (NR%4==0) {
			qual=$0
			new_qual=substr($0,UMI+1)
			print head"\tBX:Z:"umi;
			print new_seq;
			print three;
			print new_qual;
		}

	}' | $GZ -c -1 > $FASTQ_R2_UMI

	$GZ -d -c $FASTQ_R2_UMI | head -n 10


fi;



if ((0)); then

	java -jar $PICARD MarkDuplicates \
		  I=$BAM_MARKDUP_OUT \
		  O=$BAM_MARKDUP_OUT.2.bam \
		  M=$BAM_MARKDUP_OUT.2.mx

	$SAMTOOLS view $BAM_MARKDUP_OUT.2.bam | head -n 10
	grep "^#" -v $BAM_MARKDUP_OUT.2.mx | column -t

	$SAMTOOLS flagstat $BAM_MARKDUP_OUT.2.bam

fi;

if ((0)); then
	java -jar $PICARD UmiAwareMarkDuplicatesWithMateCigar \
	      I=$BAM_OUT \
	      O=$BAM_MARKDUP_OUT.UMI.bam \
	      M=$BAM_MARKDUP_OUT.UMI.dup.mx \
	      UMI_METRICS=$BAM_MARKDUP_OUT.UMI.umi.mx \
		  UMI_TAG_NAME=BX

	$SAMTOOLS view $BAM_MARKDUP_OUT.UMI.bam | head -n 10
	grep "^#" -v $BAM_MARKDUP_OUT.UMI.dup.mx | column -t
	grep "^#" -v $BAM_MARKDUP_OUT.UMI.umi.mx | column -t

	$SAMTOOLS flagstat $BAM_MARKDUP_OUT.UMI.bam

fi;

if ((0)); then

	java -jar $PICARD MarkDuplicates \
		  I=$BAM_MARKDUP_OUT.UMI.bam \
		  O=$BAM_MARKDUP_OUT.UMI.bam.2.bam \
		  M=$BAM_MARKDUP_OUT.UMI.bam.2.mx

	$SAMTOOLS view $BAM_MARKDUP_OUT.UMI.bam.2.bam | head -n 10
	grep "^#" -v $BAM_MARKDUP_OUT.UMI.bam.2.mx | column -t

	$SAMTOOLS flagstat $BAM_MARKDUP_OUT.UMI.bam.2.bam


fi;


# FASTQ to BAM
#java -jar $PICARD FastqToSam F1=$FASTQ_UMI O=$BAM_OUT SM=$SAMPLE_NAME RG=rg0013

cat $FASTQ_UMI | awk '{
	if (NR%4==1) {
		n_head=split($0,head,"\t");
		n_read=split(head[1],read,"_");
		read_name=read[1]
		read_umi=read[2]
		printf read_name
		for (i = 2; i <= n_head; ++i) {
			print "\t"head[i]
		}
		printf "\tBX:Z:"read_umi
		print ""
	}
	else
	{
		print $0
	}
}' | head -n 10


exit 0

echo $BAM_OUT
samtools view -H $BAM_OUT | head


# GROUP
umi_tools group -I $BAM_OUT --output-bam -S $BAM_GROUP_OUT
# --paired
echo $BAM_GROUP_OUT
samtools view $BAM_GROUP_OUT | head




#
