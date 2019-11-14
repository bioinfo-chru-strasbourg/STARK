#!/bin/bash

#argument: path to the samplesheet in a run folder

#Path setup
#PICARD_PATH=/mnt/galaxy/data/galaxy/sw/picard-tools-2.0.1
PICARD_PATH=/home1/TOOLS/tools//picard/2.3.0/
CPU_COUNT=`lscpu -p | egrep -v '^#' | sort -u -t, -k 2,4 | wc -l`

#demultiplexer settings
MAX_MISMATCHES=2
MAX_NO_CALLS=0
MIN_MISMATCH_DELTA=2

#Java settings
FREE_MEMORY_KB=`grep -i memavailable /proc/meminfo | sed -e's/[^0-9]//g'`
if [[ -z "$FREE_MEMORY_KB" ]]; then #memavailable is only present in recent kernels, ~ free + cached
	FREE_MEMORY_KB=`grep -i memfree /proc/meminfo | sed -e's/[^0-9]//g'`;
	FREE_MEMORY_KB=$(( $FREE_MEMORY_KB + `grep -i '^cached:' /proc/meminfo | sed -e's/[^0-9]//g'` ))
fi
JAVA_OPTS="-Xmx${FREE_MEMORY_KB}k"
 

echo "Running on ${CPU_COUNT} cpus using ${FREE_MEMORY_KB}KB RAM"

if [ $# -eq 0 ]; then
	echo "specify the sample sheet"
	exit 1
fi

sample_sheet=`readlink -f "${1}"`
echo "Sample sheet: ${sample_sheet}"
run_path=`dirname "${sample_sheet}"`
echo "Run path: ${run_path}"

if [[ $( wc -l "${sample_sheet}" | cut -f 1 -d ' ' ) -eq 0 ]]; then
	echo "fixing line ending from \r to \n"
	sed -i 's/\r/\n/g' "${sample_sheet}"
fi

run_barcode=`echo "${run_path}" | rev | cut -f1 -d'/' | tr '_' '-' | rev | cut -f2 -d'-'`
echo "Run barcode: ${run_barcode}"

output_path=$run_path 
if [[ ! -d "${output_path}/fastq" ]] ; then
	mkdir -p -m 777 "${output_path}/fastq"
fi

if [[ ! -x `which xmllint` ]] ; then
	echo "xmllint program is not found on the path, cannot continue."
	exit 1
fi
num_reads=`echo 'xpath count(//Reads/Read)' | xmllint --shell "${run_path}/RunInfo.xml"  | grep ':' | awk -F ': ' '{print $2}'`
flowcell=`echo 'cat //Flowcell' | xmllint --shell "${run_path}/RunInfo.xml" | grep '<Flowcell' | sed -r 's/<[^>]+>//g'`
machine_name=`echo 'cat //Instrument' | xmllint --shell "${run_path}/RunInfo.xml" | grep '<Instrument' | sed -r 's/<[^>]+>//g'`
if  [[ $num_reads -lt 1 ]]; then
	echo "Failed to find any reads in ${run_path}/RunInfo.xml"
	exit 1
fi

bc1_cycles=0
bc2_cycles=0
read2_cycles=0
read1_cycles=0

if [[ $num_reads -ge 1 ]]; then
	read1_cycles=`echo 'cat //Read[@Number="1"]/@NumCycles' | xmllint -shell "${run_path}/RunInfo.xml"  | grep 'NumCycles=' | sed s/.*=// | sed s/\"//g`	
fi

if [[ $num_reads -ge 2 ]]; then
	bc1_cycles=`echo 'cat //Read[@Number="2"]/@NumCycles' | xmllint -shell "${run_path}/RunInfo.xml"  | grep 'NumCycles=' | sed s/.*=// | sed s/\"//g`
	if [[ $num_reads -eq 2 ]]; then
		if [[ $bc1_cycles -gt 12 ]]; then #read 2 is not a barcode, it's a second full read
                        read2_cycles=$bc1_cycles
                        bc1_cycles=0
		else
                        echo -n '' # no need to do anything here, already have the bc1 cycles
		fi
	elif [[ $num_reads -eq 3 ]]; then
		read2_cycles=`echo 'cat //Read[@Number="3"]/@NumCycles' | xmllint -shell "${run_path}/RunInfo.xml"  | grep 'NumCycles=' | sed s/.*=// | sed s/\"//g`	
		if [[ $read2_cycles -le 12 ]]; then #read3 is a barcode - maybe better to look at the isindexread attr...
		  bc2_cycles=$read2_cycles
		  read2_cycles=0
		fi
	elif [[ $num_reads -eq 4 ]]; then
		bc2_cycles=`echo 'cat //Read[@Number="3"]/@NumCycles' | xmllint -shell "${run_path}/RunInfo.xml"  | grep 'NumCycles=' | sed s/.*=// | sed s/\"//g`
		read2_cycles=`echo 'cat //Read[@Number="4"]/@NumCycles' | xmllint -shell "${run_path}/RunInfo.xml"  | grep 'NumCycles=' | sed s/.*=// | sed s/\"//g`	
	else
		echo "Unhandled number of reads ${num_reads}"
		exit 1
	fi
fi
read_structure="${read1_cycles}T"

if [[ $bc1_cycles  -gt 0 ]]; then
	read_structure="${read_structure}${bc1_cycles}B"
fi

if [[  $bc2_cycles  -gt 0 ]]; then
	#if read 2s contain only Ns don't try to demultiplex on this... instead just generate a parallel fastq file
	i5_indices_all_N=`cut -d',' -f 8 "${sample_sheet}" | grep -E '[GCATN]+' | uniq |  grep -Ei '^\s*[N]+\s*$' | wc -l`
	if [[ $i5_indices_all_N -eq 1 ]]; then
		read_structure="${read_structure}${bc2_cycles}T"
	else
		read_structure="${read_structure}${bc2_cycles}B"	
	fi
fi

if [[ $read2_cycles -gt 0 ]] ; then
	read_structure="${read_structure}${read2_cycles}T"
fi 

echo "read structure: ${read_structure}"
echo "writing fastq files to : ${output_path}/fastq"
pushd "${output_path}/fastq"

metrics_name="${MAX_NO_CALLS}nc_${MIN_MISMATCH_DELTA}mmd_${MAX_MISMATCHES}mis_bc_metrics.txt"

lanecount=`echo 'cat //FlowcellLayout/@LaneCount' | xmllint -shell "${run_path}/RunInfo.xml"  | grep 'LaneCount=' | sed s/.*=// | sed s/\"//g`
if [[ $lanecount -eq 1 ]]; then
	echo "Detected miseq format"
	is_miseq=true
elif [[ "$lanecount" -le 4 ]]; then
	is_miseq=true
else
	is_miseq=false
fi

for i in `seq 1 ${lanecount}`
do
	echo processing lane ${i} 

	multiplex_params="${run_path}/lane${i}_multiplex_params.txt" 
	barcode_params="${run_path}/lane${i}_barcode_params.txt"

	if [[ $bc2_cycles -gt 0 ]] ; then
		echo 'barcode_sequence_1	barcode_sequence_2	barcode_name	library_name' > "${barcode_params}"
	else
		echo 'barcode_sequence_1	barcode_name	library_name' > "${barcode_params}"
	fi
	regex=
	if $is_miseq ; then
		regex="/^([^,]+)(?:[^,]*,){3,4}([^,]+),([GCATN]+)(?:,([^,]*),(?:([GCATN]*),?))?.*$/"
		if [[ $bc2_cycles -gt 0 ]] ; then
			perl -nle "print \"\$3\t\$5\t\$2+\$4\t\$1\" if ${regex}" "${sample_sheet}" | tr ' ' '-' >> "${barcode_params}"
		else
			perl -nle "print \"\$3\t\$2\t\$1\" if ${regex}" "${sample_sheet}" | tr ' ' '-' >> "${barcode_params}"
		fi
	else
		regex="/^(?:[^,]+),${i},([^,]+),(?:[^,]*,)([GCATN]*),(?:[^,]*,){4}\w+\s*$/"
		perl -nle "print ((\$2 || 'N').\"\t\$1\t\$1\") if ${regex}" "${sample_sheet}" | tr ' ' '-' >> "${barcode_params}"
	fi

	barcode_count=$((`wc -l "${barcode_params}" | cut -d ' ' -f1` - 1))
	if [[ $barcode_count -eq 0 ]]; then
		echo "Warning: Failed to find any barcodes in ${sample_sheet}" 1>&2
	fi

	if [[ $bc2_cycles -gt 0 ]] ; then
		echo 'OUTPUT_PREFIX	BARCODE_1	BARCODE_2' > "${multiplex_params}"
	else
		echo 'OUTPUT_PREFIX	BARCODE_1' > "${multiplex_params}"
	fi

	#add an unassigned bin if there are any barcode cycles
	if [[ $bc2_cycles -gt 0 ]]; then #assumes if bc2 cycles is greater than 0 , there must also be bc1 cycles
		echo "L${i}_unassigned	N	N" >> "${multiplex_params}"
	else
		echo "L${i}_unassigned	N" >> "${multiplex_params}"
	fi
	
	if $is_miseq ; then
		if [[ $bc2_cycles -gt 0 ]] ; then		
			perl -nle "print \"L${i}_\$1\t\$3\t\$5\" if ${regex}" "${sample_sheet}" | tr ' ' '-' >> "${multiplex_params}"
		else 
			perl -nle "print \"L${i}_\$1\t\$3\" if ${regex}" "${sample_sheet}" | tr ' ' '-' >> "${multiplex_params}"
		fi
	else
		perl -nle "print \"L${i}_\$1\t\$2\" if ${regex} " "${sample_sheet}" | tr ' ' '-' >> "${multiplex_params}"
	fi
	
	if [[ $barcode_count -gt 0 ]]; then

  	    java  $JAVA_OPTS -jar $PICARD_PATH/picard.jar ExtractIlluminaBarcodes\
		MAX_NO_CALLS=$MAX_NO_CALLS MIN_MISMATCH_DELTA=$MIN_MISMATCH_DELTA \
		MAX_MISMATCHES=$MAX_MISMATCHES NUM_PROCESSORS=$CPU_COUNT \
		read_structure=$read_structure \
		LANE=${i} \
		BASECALLS_DIR="${run_path}/Data/Intensities/BaseCalls" \
		METRICS_FILE="L${i}_${metrics_name}" BARCODE_FILE="${barcode_params}"
	fi

	java  $JAVA_OPTS -jar $PICARD_PATH/picard.jar IlluminaBasecallsToFastq \
		NUM_PROCESSORS=$CPU_COUNT \
		read_structure=$read_structure \
		RUN_BARCODE=$run_barcode \
		LANE=${i} \
		MACHINE_NAME=$machine_name \
		FLOWCELL_BARCODE=$flowcell \
		BASECALLS_DIR="${run_path}/Data/Intensities/BaseCalls" \
		MULTIPLEX_PARAMS="${multiplex_params}" \
		MAX_READS_IN_RAM_PER_TILE=1500000 
done
popd
