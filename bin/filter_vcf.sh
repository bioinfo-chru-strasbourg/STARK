
#!/usr/bin/env bash

set -Eeuo pipefail

# -----------------------------
# Default parameters
# -----------------------------
INPUT=""
OUTPUT=""
THREADS=1
EXCLUDE_GENOTYPES=""
EXCLUDE=""
INCLUDE=""
FILTER=""
NORM=""
COMPRESS=""
SORT=false
BCFTOOLS=$(command -v bcftools)

# -----------------------------
# Help
# -----------------------------
usage() {
    cat <<EOF
Usage:
  $0 --input input.vcf --output output.vcf [options]

Required:
  -i, --input FILE              Input VCF file
  -o, --output FILE             Output VCF file

Options:
  -t, --threads INT             Number of threads 
                                [default: 1]
  -g, --exclude-genotypes LIST  Genotypes to exclude, comma-separated
                                [default: '']
  -e, --exclude EXPR            bcftools exclusion expression
                                [default: '']
  -I, --include EXPR            bcftools inclusion expression
                                [default: '']
  -f, --filter LIST             bcftools FILTER selection (-f)
                                [default: '']
  -n, --norm EXPR               bcftools normalization expression
                                [default: '']
  -c, --compress INT            Compression level for output VCF
                                [default: '']
  -s, --sort                    Sort output with bcftools sort
                                [default: false]
  --bcftools PATH               Path to bcftools executable
                                [default: $(command -v bcftools)]

  -h, --help                    Show this help

Examples:
  $0 -i input.vcf -o output.vcf

  $0 -i input.vcf -o output.vcf \\
     --exclude-genotypes "0/0,./." \\
     --threads 8 --sort

  $0 -i input.vcf.gz -o output.vcf.gz \\
     --exclude 'QUAL<30' \\
     --include 'INFO/DP>10' \\
     --filter PASS \\
     --sort
EOF
}

# -----------------------------
# Argument parsing
# -----------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            INPUT="$2"; shift 2 ;;
        -o|--output)
            OUTPUT="$2"; shift 2 ;;
        -t|--threads)
            THREADS="$2"; shift 2 ;;
        -g|--exclude-genotypes)
            EXCLUDE_GENOTYPES="$2"; shift 2 ;;
        -e|--exclude)
            EXCLUDE="$2"; shift 2 ;;
        -I|--include)
            INCLUDE="$2"; shift 2 ;;
        -f|--filter)
            FILTER="$2"; shift 2 ;;
        -n|--norm)
            NORM="$2"; shift 2 ;;
        -c|--compress)
            COMPRESS="$2"; shift 2 ;;
        -s|--sort)
            SORT=true; shift ;;
        -h|--help)
            usage; exit 0 ;;
        -b|--bcftools)
            BCFTOOLS="$2"; shift 2 ;;
        *)
            echo "Unknown argument: $1" >&2
            usage
            exit 1 ;;
    esac
done

# -----------------------------
# Validation
# -----------------------------
if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo "Error: --input and --output are required." >&2
    usage
    exit 1
fi

if [[ ! -f "$INPUT" ]]; then
    echo "Error: input file does not exist: $INPUT" >&2
    exit 1
fi

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: threads must be a positive integer." >&2
    exit 1
fi

if ! [[ "$COMPRESS" =~ ^[1-9][0-9]*$ ]] && ! [[ -z "$COMPRESS" ]]; then
    echo "Error: compress must be a positive integer." >&2
    exit 1
fi

if ! command -v "$BCFTOOLS" >/dev/null 2>&1; then
    echo "Error: bcftools not found." >&2
    exit 1
fi

echo "compress: $COMPRESS"

# -----------------------------
# Build genotype exclusion expression
# -----------------------------
GENOTYPE_EXPR=""

if [[ -n "$EXCLUDE_GENOTYPES" ]]; then
    IFS=',' read -ra GENOTYPES <<< "$(echo "$EXCLUDE_GENOTYPES" | tr -d '"' | tr -d "'" )"

    for gt in "${GENOTYPES[@]}"; do
        [[ -z "$gt" ]] && continue

        if [[ -n "$GENOTYPE_EXPR" ]]; then
            GENOTYPE_EXPR+=" || "
        fi

        GENOTYPE_EXPR+="GT=\"$gt\""
    done
fi

# -----------------------------
# Build bcftools command
# -----------------------------
# -----------------------------
# Build bcftools pipeline
# -----------------------------

# Input format based on extension
if [[ "$INPUT" == *.vcf.gz ]]; then
    INPUT_TYPE="z"
elif [[ "$INPUT" == *.bcf ]]; then
    INPUT_TYPE="b"
elif [[ "$INPUT" == *.bcf.gz ]]; then
    INPUT_TYPE="b"
elif [[ "$INPUT" == *.vcf ]]; then
    INPUT_TYPE="v"
else
    echo "Error: unsupported input file extension: $INPUT" >&2
    exit 1
fi

# Output format based on extension
if [[ "$OUTPUT" == *.vcf.gz ]]; then
    OUTPUT_TYPE="z"
elif [[ "$OUTPUT" == *.bcf ]]; then
    OUTPUT_TYPE="b"
elif [[ "$OUTPUT" == *.bcf.gz ]]; then
    OUTPUT_TYPE="b"
elif [[ "$OUTPUT" == *.vcf ]]; then
    OUTPUT_TYPE="v"
else
    echo "Error: unsupported output file extension: $OUTPUT" >&2
    exit 1
fi

# Initial command
CMD="$BCFTOOLS view --threads $THREADS -Ou"

# -----------------------------
# Check if empty VCF - without variants
# -----------------------------


if [[ $($BCFTOOLS view -H "$INPUT" | grep -vc '^#') -eq 0 ]]; then

    # Input VCF contains no variants, just output as is (with compression if applicable)
    echo "Warning: input VCF contains no variants: $INPUT" >&2
    CMD+=" $INPUT -O $OUTPUT_TYPE -o \"$OUTPUT\""

else


    # Apply genotype exclusion if specified
    if [[ -n "$GENOTYPE_EXPR" ]]; then
        CMD+=" -e '$GENOTYPE_EXPR'"
    fi

    # Add input file to command
    CMD+=" $INPUT "

    # Apply custom exclusion if specified
    if [[ -n "$EXCLUDE" ]]; then
        CMD+=" | $BCFTOOLS view --threads $THREADS -e '$EXCLUDE' -Ou"
    fi

    # Apply inclusion if specified
    if [[ -n "$INCLUDE" ]]; then
        CMD+=" | $BCFTOOLS view --threads $THREADS -i '$INCLUDE' -Ou"
    fi

    # Apply filter if specified
    if [[ -n "$FILTER" ]]; then
        CMD+=" | $BCFTOOLS view --threads $THREADS -f '$FILTER' -Ou"
    fi

    # Apply filter if specified
    if [[ -n "$NORM" ]]; then
        CMD+=" | $BCFTOOLS norm --threads $THREADS $NORM -Ou"
    fi

    # Apply sorting if specified
    if [[ "$SORT" == true ]]; then
        CMD+=" | $BCFTOOLS sort -T $OUTPUT.$RANDOM.sort"
    fi

    # Determine final output command based on applied filters and sorting
    if [[ ! -n "$GENOTYPE_EXPR" ]] && [[ ! -n "$EXCLUDE" ]] && [[ ! -n "$INCLUDE" ]] && [[ ! -n "$FILTER" ]] && [[ ! -n "$NORM" ]] && [[ "$SORT" != true ]]; then
        if [[ "$INPUT_TYPE" == "$OUTPUT_TYPE" ]]; then
            CMD="cp \"$INPUT\" \"$OUTPUT\""
        else
            CMD+=" -O $OUTPUT_TYPE -o \"$OUTPUT\""
            if [[ -n "$COMPRESS" ]]; then
                CMD+=" -l $COMPRESS"
            fi
        fi
    else
        CMD+=" | $BCFTOOLS view --threads $THREADS -O $OUTPUT_TYPE -o \"$OUTPUT\""
        if [[ -n "$COMPRESS" ]]; then
            CMD+=" -l $COMPRESS"
        fi
    fi

fi


# -----------------------------
# Execute
# -----------------------------
echo "Input:  $INPUT"
echo "Output: $OUTPUT"
echo "Threads: $THREADS"
echo "Genotype exclusion: ${GENOTYPE_EXPR:-none}"
echo "Custom exclusion: ${EXCLUDE:-none}"
echo "Inclusion: ${INCLUDE:-none}"
echo "Filter: ${FILTER:-none}"
echo "Normalization: ${NORM:-none}"
echo "Compress: ${COMPRESS:-none}"
echo "Sort: $SORT"

echo "Executing command: $CMD"
if [[ -z "$CMD" ]]; then
    echo "Error: No command to execute" >&2
    exit 1
fi

if eval $CMD; then
    echo "Done: $OUTPUT"
    exit 0
else
    echo "Error: Command failed" >&2
    exit 1
fi

