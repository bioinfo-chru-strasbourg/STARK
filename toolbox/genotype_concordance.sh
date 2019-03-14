#! /bin/sh
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Authors: Amandine VELT
# Date: April 2016
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Define the function to print the usage of the script
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    sh genotype_concordance.sh -i initial-vcf -f final-vcf -n nb-workflow -p priority-list -v vcftools -b bcftools [-h]

		Description:
		    This script allows to add genotype concordance information on the final VCF file generated with the "PRIORITIZE" option of gatk 
		    by using the initial VCF file which has been generated with the "UNIQUIFY" option of gatk
			1. It verify if all the genotypes of all the pipelines in the initial VCF file are equal (with UNIQUIFY we have one column per pipeline)
			2. If yes, it print the GenotypeConcordance to TRUE in the final VCF file, else it print FALSE. If there is a unique pipeline, it print NONE.
			3. It count the number of pipelines which find the variant and print it in the FindByPipelines information.

		Options:
		 	-i, --initial-vcf the vcf file of the run which has been generated with the UNIQUIFY option of gatk combineVariants tool, 
		 	 it contains all the variants find by all the pipelines (union) with one column per pipeline.
		        	This option is required. 
		 	-f, --final-vcf the vcf file of the run which has been generated with the PRIORITIZE option of gatk combineVariants tool, 
		 	it contains all the variants find by all the pipelines (union) with an unique column, determined by the priority list.
		        	This option is required.
		 	-n, --nb-workflow The number of workflow(s) to consider
		      	This option is required.
		      	-p, --priority-list the priority list which contain all the pipeline to compare, the first writed pipeline being the most prioritize pipeline.
		      	This option is required.
		 	-v, --vcftools Path to vcftools
		 	This option is required.
		 	-b, --bcftools Path to bcftools
		        	This option is required.
		  	-h, --help
		        	Print this message and exit the program.
		__EOF__
}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ARGS=$(getopt -o "i:f:n:p:v:b:h" --long "initial-vcf:,final-vcf:,nb-workflow:,priority-list:,vcftools:,bcftools:,help" -- "$@" 2> /dev/null)
[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && \
	exit 1
eval set -- "$ARGS"
while true
do
	case "$1" in
		-i|--initial-vcf)
			INITIAL_VCF="$2"
			shift 2 
			;;
		-f|--final-vcf)
			FINAL_VCF="$2"
			shift 2 
			;;
		-n|--nb-workflow)
			NB_WORKFLOW="$2"
			shift 2 
			;;
		-p|--priority-list)
			PRIORITY_LIST="$2"
			shift 2 
			;;
		-v|--vcftools)
			VCFTOOLS="$2"
			shift 2 
			;;
		-b|--bcftools)
			BCFTOOLS="$2"
			shift 2 
			;;
		-h|--help)
			usage
			exit 0
			;;
		--) shift
			break 
			;;
		*) 	echo "Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Checking the input parameter
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
[ "$INITIAL_VCF" == "" ] || [ "$FINAL_VCF" == "" ] || [ "$NB_WORKFLOW" == "" ] || [ "$PRIORITY_LIST" == "" ] || [ "$VCFTOOLS" == "" ] || [ "$BCFTOOLS" == "" ] && \
	echo "Options --initial-vcf, --final-vcf --nb-workflow, --priority-list, --vcftools and --bcftools are required. " "Use -h or --help to display the help." && exit 1;
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# recover genotype values from all the pipelines
$VCFTOOLS/vcftools --vcf ${INITIAL_VCF} --extract-FORMAT-info GT --out ${INITIAL_VCF}.tmp1
sed -i '1d' ${INITIAL_VCF}.tmp1.GT.FORMAT
# recover the CHROM,POS,ID,REF and ALT from the prioritize pipeline, because sometimes, if there are several ALT, they may be not in the same order than in uniquify.
$BCFTOOLS query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ${FINAL_VCF} > ${INITIAL_VCF}.tmp2
# we verify that the two vcf are in the same order before paste them
cut -f1,2 ${INITIAL_VCF}.tmp1.GT.FORMAT > ${INITIAL_VCF}.a.tmp
cut -f1,2 ${INITIAL_VCF}.tmp2 > ${INITIAL_VCF}.b.tmp
DIFF=$(diff ${INITIAL_VCF}.a.tmp ${INITIAL_VCF}.b.tmp)
rm ${INITIAL_VCF}.a.tmp ${INITIAL_VCF}.b.tmp
if [ "$DIFF" != "" ]
then
	echo "${INITIAL_VCF}.tmp1.GT.FORMAT and ${INITIAL_VCF}.tmp2 are not in the same order !! Sort step running ..."
	sort -k1,1 ${INITIAL_VCF}.tmp1.GT.FORMAT > ${INITIAL_VCF}.tmp1.GT.FORMAT_ok; mv ${INITIAL_VCF}.tmp1.GT.FORMAT_ok ${INITIAL_VCF}.tmp1.GT.FORMAT
	sort -k1,1 ${INITIAL_VCF}.tmp2 > ${INITIAL_VCF}.tmp2_ok; mv ${INITIAL_VCF}.tmp2_ok ${INITIAL_VCF}.tmp2
else
	echo "${INITIAL_VCF}.tmp1.GT.FORMAT and ${INITIAL_VCF}.tmp2 are in the same order"
fi
# paste the CHROM,POS,ID,REF and ALT from prioritize and uniquify VCFs
# we used paste because the variants are in the same order in the two vcf
paste ${INITIAL_VCF}.tmp2 ${INITIAL_VCF}.tmp1.GT.FORMAT > ${INITIAL_VCF}.tmp.vcf
# keep the CHROM,POS,ID,REF and ALT from prioritize and the genotype values from all the pipelines
cut -f1-5,8- ${INITIAL_VCF}.tmp.vcf > ${INITIAL_VCF}.tmp.vcf_ok
mv ${INITIAL_VCF}.tmp.vcf_ok ${INITIAL_VCF}.tmp.vcf
# we recover the header (#) from the prioritize vcf and put it in a new vcf
grep "^#" ${FINAL_VCF} > ${FINAL_VCF}.tmp
# we add the information GenotypeConcordance in the INFO lines
tac ${FINAL_VCF}.tmp | awk '/##INFO/ && !seen{print "##INFO=<ID=GenotypeConcordance,Number=.,Type=String,Description=\"TRUE if genotypes are concordant between all the callers, else FALSE.\">"; seen++} 1' | tac > ${FINAL_VCF}.ok.tmp
mv ${FINAL_VCF}.ok.tmp ${FINAL_VCF}.tmp
# we add the information about the PRIORITY_LIST in the # lines
#tac ${FINAL_VCF}.tmp | awk -v r="$PRIORITY_LIST" '/##INFO/ && !seen{print "##Prioritize list is : " r; seen++} 1' | tac > ${FINAL_VCF}.ok.tmp
tac ${FINAL_VCF}.tmp | awk -v r="$PRIORITY_LIST" '/##INFO/ && !seen{print "##Prioritize_list_is=" r; seen++} 1' | tac > ${FINAL_VCF}.ok.tmp
mv ${FINAL_VCF}.ok.tmp ${FINAL_VCF}.tmp
# we add the information FindByPipelines in the INFO lines
tac ${FINAL_VCF}.tmp | awk '/##INFO/ && !seen{print "##INFO=<ID=FindByPipelines,Number=.,Type=String,Description=\"Number of pipelines which find the variant (if value = 0, that to say that the variant was filtered in by all pipelines).\">"; seen++} 1' | tac > ${FINAL_VCF}.ok.tmp
mv ${FINAL_VCF}.ok.tmp ${FINAL_VCF}.tmp
# here we check for each line of the ${INITIAL_VCF}.tmp.vcf if all the genotypes are equal. If yes, we write GenotypeConcordance=TRUE, else FALSE.
IFS='\t'; while read p
do
	min_to_select=$( expr `grep -v "#" "${INITIAL_VCF}.tmp.vcf" | awk -F'\t' '{print NF; exit}' ` - $NB_WORKFLOW + 1 )
	equal=""
	# we will go through all the genotype columns, by taking the first column ($min_to_select) as reference
	# we will look if the first genotype is equal to all others, if this is not the case, the genotype will be considered as discordant
	END=2
	END=$( grep -v "#" "${INITIAL_VCF}.tmp.vcf" | awk -F'\t' '{print NF; exit}' )	
	START=1
	START=$( expr `grep -v "#" "${INITIAL_VCF}.tmp.vcf" | awk -F'\t' '{print NF; exit}'` - $NB_WORKFLOW + 1 + 1 )
	i=$START
	# if there is only one pipeline, we don't define equal and genotype concordance will take NONE
	if [[ $i -gt $END ]]
	then
		equal="NONE"
	else
		while [[ $i -le $END ]]
		do
			# genotype1=$( echo $p | cut -f$min_to_select )
			# genotype2=$( echo $p | cut -f$i ) )
			# we select the two genotypes to compare
			# We increment the variable equal with TRUE if the genotypes are consistent and false otherwise
			if [ "$( echo $p | cut -f$min_to_select )" == "$( echo $p | cut -f$i )" ]
			then
				equal=$( echo "$equal TRUE" )
			else
				equal=$( echo "$equal FALSE" )
			fi
				let "i=i+1"
		done
	fi
	# we recovered the position/ref/alt informations of the variant we studied
	line=$( echo $p | cut -f1,2,3,4,5 )
	# if one of the genotypes is not equal to the other (if we have the presence of a "FALSE" in the equal variable); \
	if [[ "$( echo $equal )" =~ "FALSE" ]]
	then
		## we print the genotype as discordant in the info column
		awk -F"\t" -v pat="$line" '$0 ~ pat { $8 = $8 ";GenotypeConcordance=FALSE"; print }' ${FINAL_VCF} | tr '[[:space:]]' '\t' >> ${FINAL_VCF}.tmp
	elif [[ "$( echo $equal ) )" =~ "NONE" ]]
	then
		awk -F"\t" -v pat="$line" '$0 ~ pat { $8 = $8 ";GenotypeConcordance=NONE"; print }' ${FINAL_VCF} | tr '[[:space:]]' '\t' >> ${FINAL_VCF}.tmp
	else
		# else, we print the genotype as concordant in the info column
		awk -F"\t" -v pat="$line" '$0 ~ pat { $8 = $8 ";GenotypeConcordance=TRUE"; print }' ${FINAL_VCF} |  tr '[[:space:]]' '\t' >> ${FINAL_VCF}.tmp
	fi
done < ${INITIAL_VCF}.tmp.vcf

mv ${FINAL_VCF}.tmp ${FINAL_VCF}
sed -i 's/\tchr/\nchr/g' ${FINAL_VCF}
#echo '\n' >> ${FINAL_VCF}
echo "" >> ${FINAL_VCF}

## here we check for each line the number of pipelines which find the variant
IFS='\t'; while read p
do
	if [[ "$( echo $p )" =~ "#" ]]
	then
		echo "$p" >> ${FINAL_VCF}.tmp
	else
		# add of number of callers that find the variant
		nbset=$( echo "$p" | cut -f8 | sed 's/.*set/set/' | sed 's/;.*//' | sed 's/set=//' )
		### amandine code: missing 1 extra empty line at the end of $FINAL_VCF file  that why the laste ligne are not read by $p
		### try to check the script by adding 1 empty line
		line=$( echo $p | cut -f1,2,3,4,5 )
		if [[ "$( echo "$nbset" )" =~ "Intersection" ]]
		then
			# si c'est une intersection, tous les callers ont trouvé le variant
			awk -F"\t" -v pat="$line" '$0 ~ pat { $8 = $8 ";FindByPipelines='$NB_WORKFLOW'/'$NB_WORKFLOW'"; print }' $FINAL_VCF | tr '[[:space:]]' '\t' >> ${FINAL_VCF}.tmp
		elif [[ "$( echo "$nbset" )" =~ "FilteredInAll" ]]
		then
			awk -F"\t" -v pat="$line" '$0 ~ pat { $8 = $8 ";FindByPipelines=0/'$NB_WORKFLOW'"; print }' $FINAL_VCF | tr '[[:space:]]' '\t' >> ${FINAL_VCF}.tmp
		else
			nbcallers=$( expr $( expr `echo "$nbset" | grep -o "-" | wc -l` + 1 ) - `echo "$nbset" | grep -o "filterIn" | wc -l` )
			# else, we print the genotype as concordant in the info column
			awk -F"\t" -v pat="$line" '$0 ~ pat { $8 = $8 ";FindByPipelines='$nbcallers'/'$NB_WORKFLOW'"; print }' $FINAL_VCF | tr '[[:space:]]' '\t' >> ${FINAL_VCF}.tmp
		fi
	fi
done < ${FINAL_VCF}
mv ${FINAL_VCF}.tmp ${FINAL_VCF}
sed -i 's/\tchr/\nchr/g' ${FINAL_VCF}

rm -f ${INITIAL_VCF}.tmp*
rm -f ${INITIAL_VCF}.a.tmp ${INITIAL_VCF}.b.tmp
