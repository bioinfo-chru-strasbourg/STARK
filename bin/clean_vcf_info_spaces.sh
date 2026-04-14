#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARKCleanVCFInfoSpaces"
SCRIPT_DESCRIPTION="STARK clean VCF INFO field spaces"
SCRIPT_RELEASE="1.0.0"
SCRIPT_DATE="19/06/2025"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 01.0.0-19/06/2025: Script creation\n";

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
# Usage
function usage {
	echo "# USAGE: $(basename "$0") -i <input.vcf> [-o <output.vcf>] [-r <replace_char>] [-h]";
	echo "# -i, --input         Input VCF file (required)";
	echo "# -o, --output        Output cleaned VCF file (default: temp file)";
	echo "# -r, --replace-char  Character to replace spaces in INFO (default: '_')";

	echo "# -v/--verbose        VERBOSE option";
	echo "# -d/--debug          DEBUG option";
	echo "# -n/--release        RELEASE option";
	echo "# -h/--help           HELP option";
	echo "#";

}

# usage() {
#   cat <<EOF
# Usage: $(basename "$0") -i <input.vcf> [-o <output.vcf>] [-r <replace_char>] [-h]

# Options:
#   -i, --input         Input VCF file (required)
#   -o, --output        Output cleaned VCF file (default: temp file)
#   -r, --replace-char  Character to replace spaces in INFO (default: _)
  
#   -v/--verbose        VERBOSE option";
# 	-d/--debug          DEBUG option";
# 	-n/--release                   RELEASE option";
# 	-h/--help                      HELP option";
# EOF
#   exit 1
# }

# Parse args with getopt
ARGS=$(getopt -o i:o:r:h --long input:,output:,replace-char:,help -n "$(basename "$0")" -- "$@") || usage
eval set -- "$ARGS"

# Initialize variables
INPUT=""
OUTPUT=""
REPLACE="_"

# Parameters parsing
while true; do
  case "$1" in
    -i|--input)
      INPUT="$2";
      shift 2
      ;;
    -o|--output)
      OUTPUT="$2";
      shift 2
      ;;
    -r|--replace-char)
      REPLACE="$2";
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
    --) shift; break ;;
    *) usage ;;
  esac
done


# Check if input is provided and valid
[[ -z "$INPUT" || ! -f "$INPUT" ]] && echo "Error: valid input file required." && usage

# If output not set, create temp
if [[ -z "$OUTPUT" ]]; then
  OUTPUT=$INPUT
  echo "Output file not specified, using input file: $OUTPUT"
fi

# Check if input is a VCF file
if (( $(grep "^#" -v $INPUT | cut -f8 | grep " " | head -n1 | wc -l) )); then
  echo "Spaces found in INFO field, replacing with '$REPLACE'."
  FOUND_SPLACES=1
else
  echo "No spaces found in INFO field, nothing to replace."
  FOUND_SPLACES=0
fi

# If output = input, write to temp then overwrite
if [[ "$OUTPUT" == "$INPUT" ]]; then
  # Check if spaces were found
  if (( $FOUND_SPLACES )); then
    TMP=$(mktemp --suffix=".vcf")
    awk -F'\t' -v OFS='\t' -v r="$REPLACE" '/^#/ {print; next} {gsub(/ /, r, $8); print}' "$INPUT" > "$TMP" && mv "$TMP" "$INPUT" && echo "Input file overwritten with cleaned content."
  # If no spaces found, just copy input to output
  else
    echo "No spaces found, nothing to overwrite."
  fi
# If output is different from input
else
  # Check if spaces were found
  if (( $FOUND_SPLACES )); then
    awk -F'\t' -v OFS='\t' -v r="$REPLACE" '/^#/ {print; next} {gsub(/ /, r, $8); print}' "$INPUT" > "$OUTPUT" && echo "Cleaned VCF saved to: $OUTPUT"
  # If no spaces found, just copy input to output
  else
    cp "$INPUT" "$OUTPUT" && echo "No spaces found, copied input to output: $OUTPUT"
  fi
fi
