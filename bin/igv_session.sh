#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="IGVSession"
SCRIPT_DESCRIPTION="Create igv_session.xml from files"
SCRIPT_RELEASE="0.9.0.0"
SCRIPT_DATE="17/11/2021"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-AGPL"

# Release note
#RELEASE_NOTES="#\n"
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.0.0-17/11/2021: Create script\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Configuration
ENV_CONFIG=$(find -L $SCRIPT_DIR/.. -name config.app)
source $ENV_CONFIG


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
	echo "# RELEASE NOTES:";
	echo -e $RELEASE_NOTES
}

# Usage
function usage {
	echo "# USAGE: $(basename $0) --files=<FILE1,FILE2,...> --output=<FILE> [options...]";
	echo "# -f|--files=<FILE1,FILE2,...>            List of files in igv session";
	echo "#                                         Formats: *bam|*cram|*vcf.gz|*bed";
	echo "# --igv_session|igv_session_xml=<FILE>    igv session xml file";
	echo "# --igv_session_json=<FILE>               igv session json file";
	echo "# --gmt=<FILE>                            gene list file";
	echo "# --folder=<FOLDER>                       Folder to find files for igv session";
	echo "# --folder_mindepth=<INTEGER>             Path min depth for folder for files for igv session";
	echo "#                                         Default: 0";
	echo "# --folder_maxdepth=<INTEGER>             Path min depth for folder for files for igv session";
	echo "#                                         Default: 1";
	echo "# --patterns=<STRING,STRING...>           Patterns for files for igv session";
	echo "# --ressources=<STRING>                   Ressources to add to igv session";
	echo "# --datapanel=<STRING>                    Data panel tracks to add to igv session";
	echo "# --featurepanel=<STRING>                 Feature panel tracks to add to igv session";
	echo "# --das_url=<URL>                         DAS URL for igv session json";
	echo "# --json_tracks_additional=<STRING>       JSON tracks additional for igv session json";
	
	echo "# -v|--verbose                            VERBOSE option";
	echo "# -d|--debug                              DEBUG option";
	echo "# -n|--release                            RELEASE option";
	echo "# -h|--help                               HELP option";
	echo "#";
}

# header
header;


# USAGE exemple


####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "e:f:o:vdnh" --long "env:,app:,application:,files:,folder:,folder_mindepth:,folder_maxdepth:,patterns:,ressources:,datapanel:,featurepanel:,output:,igv_session:,igv_session_xml:,igv_session_json:,gmt:,das_url:,json_tracks_additional:,verbose,debug,release,help" -- "$@" 2> /dev/null)
[ $? -ne 0 ] && \
	echo "#[ERROR] Error in the argument list." "Use -h or --help to display the help." >&2 && usage && \
	exit 1
eval set -- "$ARGS"

#echo $PARAM
PARAM=$@

while true
do
	#echo "$1=$2"
	#echo "Eval opts";
	case "$1" in
		-e|--env|--app|--application)
			APP="$2"
			shift 2
			;;
		-f|--files)
			FILES="$2"
			FILES=$(echo $FILES | tr "," " ")
			shift 2
			;;
		--folder)
			FOLDER="$2"
			shift 2
			;;
		--folder_maxdepth)
			FOLDER_MAXDEPTH="$2"
			shift 2
			;;
		--folder_mindepth)
			FOLDER_MINDEPTH="$2"
			shift 2
			;;
		--patterns)
			PATTERNS="$2"
			PATTERNS=$(echo $PATTERNS | tr "," " ")
			shift 2
			;;
		--ressources)
			RESSOURCES="$2"
			shift 2
			;;
		--datapanel)
			DATAPANEL="$2"
			shift 2
			;;
		--featurepanel)
			FEATUREPANEL="$2"
			shift 2
			;;
		-o|--igv_session|--igv_session_xml)
			IGV_SESSION="$2"
			shift 2
			;;
		--igv_session_json)
			IGV_SESSION_JSON="$2"
			shift 2
			;;
		--gmt)
			GMT="$2"
			shift 2
			;;
		--das_url)
			DAS_URL="$2"
			shift 2
			;;
		--json_tracks_additional)
			JSON_TRACKS_ADDITIONAL="$2"
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
		--no_header)
			NO_HEADER=1
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
		*) 	echo "#[WARNING] Option $1 is not recognized. Use -h or --help to display the help." && usage && exit 1
			;;
	esac
done
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Checking the input parameter
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if ((0))
then
	echo "Option --files is required. " "Use -h or --help to display the help." && usage && exit 1;
fi
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


# FUNCTIONS
#############

# function in_array
# input: $element $array
in_array () 
{ 
    param=$1;
    shift;
    for elem in "$@";
    do
        [[ "$param" = "$elem" ]] && return 0;
    done;
    return 1
}

VERBOSE=$(echo $VERBOSE | awk '{print $0+0}')
DEBUG=$(echo $DEBUG | awk '{print $0+0}')



### APPLICATION

ENV=$(find_app "$APP" "$STARK_FOLDER_APPS")
source_app "$APP" "$STARK_FOLDER_APPS" 1


### FOLDER DEPTH
if [ -z $FOLDER_MINDEPTH ] || ! [[ $FOLDER_MINDEPTH =~ ^[0-9]+$ ]]; then
	FOLDER_MINDEPTH=0
fi;

if [ -z $FOLDER_MAXDEPTH ] || ! [[ $FOLDER_MAXDEPTH =~ ^[0-9]+$ ]]; then
	FOLDER_MAXDEPTH=1
fi;


### FOLDER PATTERNS

if [ ! -z "$FOLDER" ] && [ ! -z "$PATTERNS" ]; then
	PATTERNS_PARAM='-name '$(echo $PATTERNS | sed "s/ / -or -name /gi")
	#echo $PATTERNS_PARAM
	#PATTERNS_PARAM=$(echo $PATTERNS_PARAM | sed 's/$SAMPLE/'$S'/gi')
	FILES=$(find $FOLDER -mindepth $FOLDER_MINDEPTH -maxdepth $FOLDER_MAXDEPTH $PATTERNS_PARAM | tac | sed "s#$FOLDER/\(.*\)#\1:$FOLDER/\1#gi" | tr "\n" " ")

fi;

### FILES


### IGV_SESSION
if [ -z "$IGV_SESSION" ]; then
	IGV_SESSION="igv_session.xml"
fi;



if (($DEBUG)); then
	echo "#[INFO] FILES=$FILES"
	echo "#[INFO] PATTERNS=$PATTERNS"
	echo "#[INFO] IGV_SESSION=$IGV_SESSION"
	#exit 0
fi;


Ressources_basename=""
Ressources=""
DataPanel=""
FeaturePanel=""
Regions=""

[ "$GMT" != "" ] && echo "#name=STARK" > $GMT

touch $IGV_SESSION.regions

JSON_TRACKS=""
JSON_TRACKS_ALIGNMENT=""
JSON_TRACKS_VARIANTS=""
JSON_TRACKS_BED=""


for file_infos in $FILES; do

	MyList=""

	file=$(echo $file_infos | awk -F: '{print $1}')
	file_path=$(echo $file_infos | awk -F: '{print $2}')

	(($VERBOSE)) && echo "#[INFO] file '$file'"

	if (($(echo "$file" | grep ".bam$\|.cram$" -c))); then
		FORMAT=${file##*.}
		INDEX="bai";
		if [ "$FORMAT" == "bam" ]; then
			INDEX="bai";
		elif [ "$FORMAT" == "cram" ]; then
			INDEX="crai";
		fi;
		(($VERBOSE)) && echo "#[INFO] file type: 'Alignment' [$FORMAT]"
		(($DEBUG)) && echo "#[INFO] $file_path.bai "$(basename $file)" $Ressources_basename"
		if ! in_array $(basename $file) $Ressources_basename; then
			if [ -e $file_path.$INDEX ]; then
				# XML
				file_index_attr='index="'$file'.'$INDEX'"'
				file_index_attr='' # bug
				Ressources=$Ressources' <Resource path="'$file'" '$file_index_attr' type="'$FORMAT'"/>'
				Ressources_basename=$Ressources_basename" "$(basename $file)
				AlignmentPanel=$AlignmentPanel'
				<Panel name="Panel'$RANDOM'">
					<Track attributeKey="'$file' Coverage" autoScale="true" clazz="org.broad.igv.sam.CoverageTrack" color="175,175,175" colorScale="ContinuousColorScale;0.0;1494.0;255,255,255;175,175,175" fontSize="10" id="'$file'_coverage" name="'$file' Coverage" snpThreshold="1" visible="true">
						<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" type="LINEAR"/>
					</Track>
					<Track attributeKey="'$file' Junctions" clazz="org.broad.igv.sam.SpliceJunctionTrack" fontSize="10" groupByStrand="false" height="60" id="'$file'_junctions" name="'$file' Junctions" visible="false"/>
					<Track attributeKey="'$file'" clazz="org.broad.igv.sam.AlignmentTrack" displayMode="SQUISHED" experimentType="OTHER" fontSize="10" id="'$file'" name="'$file'" visible="true">
						<RenderOptions showAllBases="false"/>
					</Track>
				</Panel>'
				(($VERBOSE)) && echo "#[INFO] file XML '$FORMAT' processed"
				# JSON
				if [ "$IGV_SESSION_JSON" != "" ] && [ "$DAS_URL" != "" ]; then
					# TODO
					[ "$JSON_TRACKS_ALIGNMENT" != "" ] && JSON_TRACKS_SEP="," || JSON_TRACKS_SEP=""
					JSON_TRACKS_ALIGNMENT=$JSON_TRACKS_ALIGNMENT$JSON_TRACKS_SEP'{
						"url": "'$DAS_URL'/'$file'",
						"indexURL": "'$DAS_URL'/'$file'.'$INDEX'",
						"name": "'$file'",
						"sourceType": "file",
						"format": "'$INDEX'",
						"displayMode": "SQUISHED",
						"type": "alignment"
					}'
				else
					(($VERBOSE)) && echo "#[WARN] file JSON '$FORMAT' not processed!"
				fi;
			else
				(($VERBOSE)) && echo "#[WARN] file '$FORMAT' index not found!"
			fi;
		else
			(($VERBOSE)) && echo "#[WARN] file '$FORMAT' already processed!"
		fi;

	elif (($(echo "$file" | grep ".vcf.gz$" -c))); then
		FORMAT="vcf";
		(($VERBOSE)) && echo "#[INFO] file type: 'Variant'"
		if ! in_array $(basename $file) $Ressources_basename; then
			INDEX="tbi"
			file_index_attr=''
			if [ -e $file_path.$INDEX ]; then
				file_index_attr='index="'$file'.'$INDEX'"'
			fi;
			file_index_attr='' # bug
			if (("$($BCFTOOLS view $file_path | grep '^#' -vc)")); then
				# XML
				Ressources=$Ressources' <Resource '$file_index_attr' path="'$file'" type="'$FORMAT'"/>'
				Ressources_basename=$Ressources_basename" "$(basename $file)
				DataPanel=$DataPanel' <Track attributeKey="'$file'" clazz="org.broad.igv.variant.VariantTrack" displayMode="COLLAPSED" fontSize="10" groupByStrand="false" id="'$file'" name="'$file'" siteColorMode="ALLELE_FREQUENCY" squishedHeight="1" visible="true"/>'
				(($VERBOSE)) && echo "#[INFO] file '$FORMAT' processed"
				# JSON
				if [ "$IGV_SESSION_JSON" != "" ] && [ "$DAS_URL" != "" ]; then
					# TODO
					[ "$JSON_TRACKS_VARIANTS" != "" ] && JSON_TRACKS_SEP="," || JSON_TRACKS_SEP=""
					JSON_TRACKS_VARIANTS=$JSON_TRACKS_VARIANTS$JSON_TRACKS_SEP'{
						"url": "'$DAS_URL'/'$file'",
						"indexURL": "'$DAS_URL'/'$file'.'$INDEX'",
						"name": "'$file'",
						"sourceType": "file",
						"format": "vcf",
						"displayMode": "COLLAPSED",
						"height": "40",
						"type": "variant"
					}'
				else
					(($VERBOSE)) && echo "#[WARN] file JSON '$FORMAT' not processed!"
				fi;
			else
				#(($VERBOSE)) && echo "#[WARN] file index not found!"
				(($VERBOSE)) && echo "#[WARN] file '$FORMAT' do not contain variants!"
			fi;
		else
			(($VERBOSE)) && echo "#[WARN] file '$FORMAT' already processed!"
		fi;
	elif (($(echo "$file" | grep ".bed$" -c))); then
		FORMAT="bed";
		(($VERBOSE)) && echo "#[INFO] file type: 'Region'"
		if ! in_array $(basename $file) $Ressources_basename; then
			# XML
			Ressources=$Ressources' <Resource path="'$file'" type="'$FORMAT'"/>'
			Ressources_basename=$Ressources_basename" "$(basename $file)
			FeaturePanel=$FeaturePanel' <Track attributeKey="'$file'" clazz="org.broad.igv.track.FeatureTrack" colorScale="ContinuousColorScale;0.0;58.0;255,255,255;0,0,178" fontSize="10" displayMode="COLLAPSED" groupByStrand="false" id="'$file'" name="'$file'" visible="true"/>'
			(($VERBOSE)) && echo "#[INFO] file '$FORMAT' processed"
			# regions
			cat $file_path | awk '{print "<Region chromosome=\""$1"\" description=\""$4"\" start=\""$2"\" end=\""$3"\"/>"}' >> $IGV_SESSION.regions
			(($VERBOSE)) && echo "#[INFO] file '$FORMAT' regions processed"
			# GMT
			if [ "$GMT" != "" ]; then
				#Regions=$Regions$"\n"$(cat $file_path | awk '{print "<Region chromosome=\""$1"\" description=\""$4"\" start=\""$2"\" end=\""$3"\"/>'$"\n"'"}')
				echo -e "$file		"$(cat $file_path | cut -f4 | sort -u | tr "\n" "\t") >> $GMT
				(($VERBOSE)) && echo "#[INFO] file '$FORMAT' GMT processed"
			else
				(($VERBOSE)) && echo "#[INFO] file '$FORMAT' GMT not processed"
			fi;
			# JSON
			if [ "$IGV_SESSION_JSON" != "" ] && [ "$DAS_URL" != "" ]; then
				# TODO
				[ "$JSON_TRACKS_BED" != "" ] && JSON_TRACKS_SEP="," || JSON_TRACKS_SEP=""
				JSON_TRACKS_BED=$JSON_TRACKS_BED$JSON_TRACKS_SEP'{
					"url": "'$DAS_URL'/'$file'",
					"indexURL": null,
					"name": "'$file'",
					"sourceType": "file",
					"format": "",
					"displayMode": "COLLAPSED",
					"height": "40",
					"type": "annotation"
				}'
			else
				(($VERBOSE)) && echo "#[WARN] file JSON '$FORMAT' not processed!"
			fi;
		else
			(($VERBOSE)) && echo "#[WARN] file '$FORMAT' already processed!"
		fi;
	else
		FORMAT="unknown";
		(($VERBOSE)) && echo "#[WARN] file '$FORMAT' not processed!"
	fi;

done;


# Regions
Regions=$(cat $IGV_SESSION.regions | sort -u)
rm -f $IGV_SESSION.regions


# JSON
JSON_TRACKS=$JSON_TRACKS_VARIANTS
if [ "$JSON_TRACKS_ALIGNMENT" != "" ]; then
	[ "$JSON_TRACKS" != "" ] && JSON_TRACKS=$JSON_TRACKS","$JSON_TRACKS_ALIGNMENT || JSON_TRACKS=$JSON_TRACKS$JSON_TRACKS_ALIGNMENT
fi;
if [ "$JSON_TRACKS_BED" != "" ]; then
	[ "$JSON_TRACKS" != "" ] && JSON_TRACKS=$JSON_TRACKS","$JSON_TRACKS_BED || JSON_TRACKS=$JSON_TRACKS$JSON_TRACKS_BED
fi;
if [ "$JSON_TRACKS_ADDITIONAL" != "" ]; then
	[ "$JSON_TRACKS" != "" ] && JSON_TRACKS=$JSON_TRACKS","$JSON_TRACKS_ADDITIONAL || JSON_TRACKS=$JSON_TRACKS$JSON_TRACKS_ADDITIONAL
fi;


# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz


if (($DEBUG)); then
	echo $Ressources
	echo $DataPanel
	echo $AlignmentPanel
	echo $FeaturePanel
	echo $Regions
	cat $GMT
fi;


echo '<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="hg19" hasGeneTrack="false" hasSequenceTrack="true" locus="All" version="8">
    <Resources>
        '$Ressources'
		'$RESSOURCES'
    </Resources>
    <Panel name="DataPanel">
        '$DATAPANEL'
		'$DataPanel'
    </Panel>
    '$AlignmentPanel'
    <Panel name="FeaturePanel">
		<Track attributeKey="Reference sequence" clazz="org.broad.igv.track.SequenceTrack" fontSize="10" id="Reference sequence" name="Reference sequence" sequenceTranslationStrandValue="POSITIVE" shouldShowTranslation="false" visible="true"/>
            '$FEATUREPANEL'
			'$FeaturePanel'
    </Panel>
    <PanelLayout dividerFractions="0.20,0.80"/>
	<Regions>
        '$Regions'
    </Regions>
    <HiddenAttributes>
        <Attribute name="DATA FILE"/>
        <Attribute name="DATA TYPE"/>
        <Attribute name="NAME"/>
    </HiddenAttributes>
</Session>' > $IGV_SESSION

# <Resource index="https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz.tbi" name="Refseq Genes" path="https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz" type="refgene"/>
# <Track attributeKey="Refseq Genes" clazz="org.broad.igv.track.FeatureTrack" fontSize="10" groupByStrand="false" id="https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz" name="Refseq Genes" visible="true"/>

if [ "$IGV_SESSION_JSON" != "" ] && [ "$DAS_URL" != "" ]; then

	echo '
	{
		"genome": "'$ASSEMBLY'",
		"tracks": [
					'$JSON_TRACKS'
				]
	}
	' > $IGV_SESSION_JSON

fi;


exit 0
