VCF=$1
VCF_OUTPUT=$2
THREADS=$3
BCFTOOLS=$4
REFORMAT=$5

# input
if [ -z "$VCF" ]; then
    echo "#[ERROR] No VCF input"
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
OUTPUT_TMP=$VCF.tmp.$RANDOM

# Header with comma in Description
$BCFTOOLS view --threads=$THREADS -h $VCF | awk '
BEGIN {
    sep=""    # alternative separator 
}
{
    # header line with comma in Description
    if ( ($0 ~ /^##INFO/ || $0 ~ /^##FORMAT/ || $0 ~ /^##FILTER/) && $0 ~ /Description=\".*,.*\"/ ) {
        match($0, /^(.*)(Description=\".*\")(.*)$/, arr);
        gsub(",",sep,arr[2]);
        print arr[1] arr[2] arr[3]
    # other lines
    } else {
        print $0
    }
}' > $OUTPUT_TMP.header

# reheader
$BCFTOOLS reheader --threads=$THREADS -h $OUTPUT_TMP.header $VCF > $OUTPUT_TMP;

# Reformat
if (($REFORMAT)); then

    # Find error in bcftools
    $BCFTOOLS view -h $OUTPUT_TMP 1>/dev/null 2>$OUTPUT_TMP.invalid_tag_name.err

    # Find invalid tag name
    cat $OUTPUT_TMP.invalid_tag_name.err | grep -Po 'Invalid tag name: "\K.*?(?=")' > $OUTPUT_TMP.invalid_tag_name
    
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

        # Rename tag name and annotations
        $BCFTOOLS annotate --rename-annots $OUTPUT_TMP.invalid_tag_name.tab $OUTPUT_TMP > $OUTPUT_TMP.reformat

        # Move output VCF
        mv $OUTPUT_TMP.reformat $OUTPUT_TMP

    fi;

fi;

# Final VCF
if [ "$EXTENSION" == "gz" ]; then
    bgzip -c $OUTPUT_TMP > $OUTPUT_TMP.compressed
    mv $OUTPUT_TMP.compressed $OUTPUT_TMP
fi;

mv $OUTPUT_TMP $VCF_OUTPUT

# Clean
rm -f $OUTPUT_TMP*

exit 0

