#!/bin/bash
#################################
##
## STARK script
##
#################################

SCRIPT_NAME="STARKFixVCFHeader"
SCRIPT_DESCRIPTION="STARK Fix VCF Headner"
SCRIPT_RELEASE="0.9.0"
SCRIPT_DATE="13/02/2023"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="CHRU"
SCRIPT_LICENCE="GNU AGPLv3"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.0-13/02/2023: Script creation\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


# Header
function header () {
	echo "#######################################";
	echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
	echo "# $SCRIPT_DESCRIPTION ";
	echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © $SCRIPT_LICENCE";
	echo "#######################################";
}

# Release
function release () {
	echo -e $RELEASE_NOTES
}

# Usage
function usage {
	echo "# USAGE: $(basename $0) --vcf|input [options...]";
	echo "# --vcf|input=<FILE>              VCF input file to fix.";
	echo "#                                 Format: either '*.vcf' or '*.vcf.gz'.";
	echo "# --output=<FILE>                 VCF output file fixed.";
	echo "#                                 Compression depend on file extension.";
	echo "#                                 Default: VCF input file.";
	echo "# --reformat                      Reformat VCF header INFO tags";
	echo "#                                 Prevent '-' and first character as digit in INFO tags ID";
	echo "# --threads=<INTEGER>             Number of threads (for BCFTOOLS)";
	echo "# --bcftools=<STRING>             BCFTOOLS bin";
	echo "#                                 Default: 'bcftools'";

	echo "# --verbose                       VERBOSE option";
	echo "# --debug                         DEBUG option";
	echo "# --release                       RELEASE option";
	echo "# --help                          HELP option";
	echo "#";

}

# header
header;

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "i:o:rt:b:vdnh" --long "vcf:,input:,output:,reformat,threads:,bcftools:,verbose,debug,release,help" -- "$@" 2> /dev/null)


PARAM=$@

eval set -- "$ARGS"
while true
do
	#echo "$1=$2"
	#echo "Eval opts";
	case "$1" in
		-i|--vcf|--input)
			VCF="$2"
			shift 2
			;;
		-o|--output)
			VCF_OUTPUT="$2"
			shift 2
			;;
		-r|--reformat)
			REFORMAT=1
			shift 1
			;;
		-t|--threads)
			THREADS="$2"
			shift 2
			;;
		-b|--bcftools)
			BCFTOOLS="$2"
			shift 2
			;;
		-v|--verbose)
			VERBOSE=1
			shift 1
			;;
		-d|--debug)
			VERBOSE=1
			DEBUG=1
			shift 1
			;;
		-n|--release)
			release;
			exit 0
			;;
		-h|--help)
			usage
			exit 0
			;;
		--) shift
			break
			;;
		*) 	echo "# Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done



# OPTIONS
#########

# input
if [ -z "$VCF" ] || [ ! -e "$VCF" ]; then
    echo "#[ERROR] No VCF input"
    usage
    exit 1
fi;

# output
if [ -z "$VCF_OUTPUT" ]; then
    VCF_OUTPUT=$VCF
fi;

# threads
if [ -z "$THREADS" ]; then
    THREADS=1
fi;

# bcftools
if [ -z "$BCFTOOLS" ]; then
    BCFTOOLS=bcftools
fi;

# reformat
if [ -z "$REFORMAT" ]; then
    REFORMAT=0
fi;

# Extension
EXTENSION="${VCF_OUTPUT##*.}"


# tmp
OUTPUT_TMP=$VCF.$RANDOM.tmp


# Verbose
if (($VERBOSE)); then
    echo "#[INFO] Input VCF:    $VCF"
    echo "#[INFO] Output VCF:   $VCF_OUTPUT"
    echo "#[INFO] Reformat:     $REFORMAT"
    echo "#[INFO] Threads:      $THREADS"
    echo "#[INFO] bcftools:     $BCFTOOLS"
fi;



# Header TAGS fix
#################

(($VERBOSE)) && echo ""
(($VERBOSE)) && echo "#[INFO] Start header TAGS fixing"

$BCFTOOLS view --threads=$THREADS -h $VCF 2>$OUTPUT_TMP.input.header.err | awk '
# Reorder tags infos as Number, Type, Description
{
    # header INFO order tag infos
    if ( $0 ~ /^##INFO=/ || $0 ~ /^##FORMAT=/ || $0 ~ /^##FILTER=/ || $0 ~ /^##contig=/ ) {
        
        # Catch specific TAG infos
        match($0, /^##([^=]*)=.*[,<]ID=([^,]*).*[,>]$/, arr_ID);
        TAG = arr_ID[1]
        ID = arr_ID[2]
        match($0, /^##([^=]*)=.*[,<]Number=([^,]*).*[,>]$/, arr_Number);
        Number = arr_Number[2]
        match($0, /^##([^=]*)=.*[,<]Type=([^,]*).*[,>]$/, arr_Type);
        Type = arr_Type[2]
        match($0, /^##([^=]*)=.*[,<]Description=("[^"]*"|[^,]*).*[,>]$/, arr_Description);
        Description = arr_Description[2]
        match($0, /^##([^=]*)=.*[,<]Source=("[^"]*"|[^,]*).*[,>]$/, arr_Source);
        Source = arr_Source[2]
        match($0, /^##([^=]*)=.*[,<]Version=("[^"]*"|[^,]*).*[,>]$/, arr_Version);
        Version = arr_Version[2]
        # match($0, /^##([^=]*)=.*[,<]Source=(".*"|[^,]*).*[,>]$/, arr_Source);
        # Source = arr_Source[2]
        # match($0, /^##([^=]*)=.*[,<]Version=(".*"|[^,]*).*[,>]$/, arr_Version);
        #Version = arr_Version[2]
        match($0, /^##([^=]*)=.*[,<]length=([^,]*).*[,>]$/, arr_contig_length);
        contig_length = arr_contig_length[2]
        match($0, /^##([^=]*)=.*[,<]assembly=([^,]*).*[,>]$/, arr_assembly);
        assembly = arr_assembly[2]
        match($0, /^##([^=]*)=.*[,<]md5=([^,]*).*[,>]$/, arr_md5);
        md5 = arr_md5[2]
        match($0, /^##([^=]*)=.*[,<]species=([^,]*).*[,>]$/, arr_species);
        species = arr_species[2]
        match($0, /^##([^=]*)=.*[,<]taxonomy=([^,]*).*[,>]$/, arr_taxonomy);
        taxonomy = arr_taxonomy[2]
        match($0, /^##([^=]*)=.*[,<]URL=([^,]*).*[,>]$/, arr_URL);
        URL = arr_URL[2]
        
        # all TAGS infos
        match($0, /^##[^<]*<(.*)>$/, arr_TAGS);
        TAGS = arr_TAGS[1]
        
        # Check and fix
        if (Number == "") {
            Number = "."
        }
        if (Type == "") {
            Type = "String"
        }
        if (TAG == "INFO" && Type != "String" && Type != "Integer" && Type != "Float" && Type != "Flag") {
            Type = "String"
        }
        if (TAG == "FORMAT" && Type != "String" && Type != "Integer" && Type != "Float") {
            Type = "String"
        }
        if (Description == "") {
            Description = "\"\""
        }
        if (substr(Description, 1 , 1) != "\"" && substr(Description, length(Description) , 1) != "\"") {
            Description = "\"" Description "\""
        }
        if (Source == "") {
            Source = "\"unknown\""
        }
        if (substr(Source, 1 , 1) != "\"" && substr(Source, length(Source) , 1) != "\"") {
            Source = "\"" Source "\""
        }
        if (Version == "") {
            Version = "\"unknown\""
        }
        if (substr(Version, 1 , 1) != "\"" && substr(Version, length(Version) , 1) != "\"") {
            Version = "\"" Version "\""
        }
        # if (contig_length == "") {
        #     contig_length = "1"
        # }
        # if (assembly == "") {
        #     assembly = "unknown"
        # }
        # if (md5 == "") {
        #     md5 = "unknown"
        # }
        # if (species == "") {
        #     species = "\"unknown\""
        # }
        # if (substr(species, 1 , 1) != "\"" && substr(species, length(species) , 1) != "\"") {
        #     species = "\"" species "\""
        # }
        # if (taxonomy == "") {
        #     taxonomy = "unknown"
        # }

        # Create line
        if (TAG == "INFO") {
            line = "##" TAG "=<" "ID=" ID "," "Number=" Number "," "Type=" Type "," "Description=" Description "," "Source=" Source "," "Version=" Version ">"
        } else if (TAG == "FORMAT") {
            line = "##" TAG "=<" "ID=" ID "," "Number=" Number "," "Type=" Type "," "Description=" Description ">"
        } else if (TAG == "FILTER") {
            line = "##" TAG "=<" "ID=" ID "," "Description=" Description ">"
        } else if (TAG == "contig") {
            line = "##" TAG "=<" "ID=" ID 
            if (contig_length != "") {
                line = line  "," "length=" contig_length
            }
            if (assembly != "") {
                line = line  "," "assembly=" assembly
            }
            if (URL != "") {
                line = line  "," "URL=" URL
            }
            if (md5 != "") {
                line = line  "," "md5=" md5
            }
            if (species != "") {
                line = line  "," "species=" species
            }
            if (taxonomy != "") {
                line = line  "," "taxonomy=" taxonomy
            }
            split(TAGS, TAGS_array, ",")
            for (x in TAGS_array){
                TAGS_varval = TAGS_array[x]
                split(TAGS_varval, TAGS_varval_array, "=")
                var = TAGS_varval_array[1]
                val = TAGS_varval_array[2]
                if (var != "ID" && var != "length" && var != "assembly" && var != "URL" && var != "md5" && var != "species" && var != "taxonomy") {
                    line = line "," var "=" val
                }
            }
            line = line ">"
        }

        # print line with tag info ordered and fixed
        print line
        
    # other lines
    } else {
        print $0
    }
}' > $OUTPUT_TMP.header

(($VERBOSE)) && echo "#[INFO] Input VCF header errors (before fixes):" && cat $OUTPUT_TMP.input.header.err | awk -F"\t" '{print "#[INFO]    "$0}'

# Reheader with fixed header

$BCFTOOLS reheader --threads=$THREADS -h $OUTPUT_TMP.header $VCF > $OUTPUT_TMP 2>$OUTPUT_TMP.input.reheader.err;

$BCFTOOLS view -h $OUTPUT_TMP >/dev/null 2>$OUTPUT_TMP.input.header.fixed.err

(($VERBOSE)) && echo "#[INFO] Input VCF header errors (after fixes):" && cat $OUTPUT_TMP.input.header.fixed.err | awk -F"\t" '{print "#[INFO]    "$0}'


# Reformat VCF INFO TAGS ID
###########################

if (($REFORMAT)); then

    (($VERBOSE)) && echo ""
    (($VERBOSE)) && echo "#[INFO] Start reformat INFO TAGS ID"

    # Find error in bcftools
    $BCFTOOLS view -h $OUTPUT_TMP 1>/dev/null 2>$OUTPUT_TMP.invalid_tag_name.err

    # Find invalid tag name
    cat $OUTPUT_TMP.invalid_tag_name.err | grep -Po 'Invalid tag name: "\K.*?(?=")' > $OUTPUT_TMP.invalid_tag_name
    
    (($VERBOSE)) && echo "#[INFO] Invalid INFO TAGS ID (before reformat):" $(cat $OUTPUT_TMP.invalid_tag_name)


    # If invalid tag name
    if (( $(cat $OUTPUT_TMP.invalid_tag_name.err | wc -l) )); then

        # Create tab of invalid tag name and new tag
        cat $OUTPUT_TMP.invalid_tag_name | awk '
        {
            new = $0
            # Invalid characters "-"
            gsub("-","_",new)
            # Start with a digit
            if (substr(new,1,1) ~ /^[0-9]/ ) {
                new = "A" new
            }
            print "INFO/" $0 "\t" new
        }
        ' > $OUTPUT_TMP.invalid_tag_name.tab

        if (($VERBOSE)); then
            echo "#[INFO] Invalid INFO TAGS ID changed tab:"
            cat $OUTPUT_TMP.invalid_tag_name.tab | awk -F"\t" '{print "#[INFO]    From \"" $1 "\" To \"INFO/" $2 "\""}'
        fi;
        
        # Rename tag name and annotations
        $BCFTOOLS annotate --rename-annots $OUTPUT_TMP.invalid_tag_name.tab $OUTPUT_TMP > $OUTPUT_TMP.reformat 2>$OUTPUT_TMP.annotate.invalid_tag_name.err

        # Find error in bcftools
        $BCFTOOLS view -h $OUTPUT_TMP.reformat 1>/dev/null 2>$OUTPUT_TMP.invalid_tag_name.err

        # Find invalid tag name
        cat $OUTPUT_TMP.invalid_tag_name.err | grep -Po 'Invalid tag name: "\K.*?(?=")' > $OUTPUT_TMP.invalid_tag_name
    
        (($VERBOSE)) && echo "#[INFO] Invalid INFO TAGS ID (after reformat):" $(cat $OUTPUT_TMP.invalid_tag_name)

        # Move output VCF
        mv $OUTPUT_TMP.reformat $OUTPUT_TMP

    fi;

fi;

# OUTPUT VCF
############

if [ "$EXTENSION" == "gz" ]; then
    bgzip -c $OUTPUT_TMP > $OUTPUT_TMP.compressed
    mv $OUTPUT_TMP.compressed $OUTPUT_TMP
fi;

mv $OUTPUT_TMP $VCF_OUTPUT


# Clear
#######

rm -f $OUTPUT_TMP*

exit 0