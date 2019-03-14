

ENV=$1
BAMS=$2 # must be full path
BED=$3

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for BAM in $BAMS; do
	# Find FULL PATH
	BAM_FULL_PATH=$(echo "$(cd "$(dirname "$BAM")"; pwd)/$(basename "$BAM")")
	BED_FULL_PATH=$(echo "$(cd "$(dirname "$BED")"; pwd)/$(basename "$BED")")
	# PARAM
	PARAM="-f $BAM_FULL_PATH -s $(basename $BAM_FULL_PATH | cut -d. -f1) -r $(basename $(dirname $BAM_FULL_PATH))"
	if [ "$ENV" != "" ]; then PARAM=$PARAM" -e $ENV "; fi;
	if [ "$BED" != "" ]; then PARAM=$PARAM" -e $BED "; fi;
	# Launch STARK Command with COULSON
	$COULSON/add.sh '' "$SCRIPT_DIR/STARK $PARAM ";
done;
