#!/bin/bash
#################################
## STARK environment
#################################

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

# APPLICATION INFOS
############

# Name of the APP. usefull to include specific rules if exists (in $STARK/$APP_NAME.rules.mk/*rules.mk)
# AUTO detect: $(basename ${BASH_SOURCE[0]} | sed 's/^env//' | sed 's/\.sh$//gi' | cut -d. -f2-)
if [ "$APP_NAME" == "" ]; then
	APP_NAME="DEFAULT"
fi;
export APP_NAME

# Release of the APP
# AUTO detect: $(basename ${BASH_SOURCE[0]} | sed 's/^env//' | sed 's/\.sh$//gi' | cut -d. -f2-)
if [ "$APP_RELEASE" == "" ]; then
	APP_RELEASE=""
fi;
export APP_RELEASE

# GROUP and PROJECT Associated with the APP
# Use to structure data in the repository folder
if [ "$APP_GROUP" == "" ]; then
	APP_GROUP="UNKNOWN"
fi;
export GROUP
if [ "$APP_PROJECT" == "" ]; then
	APP_PROJECT="UNKNOWN"
fi;
export APP_PROJECT



# SCRIPT DIR
##############

# Folder of STARK ENV
#export STARK_FOLDER_CONFIG="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )"



# FOLDERS
#################
# Checking and adaptation of variables for the workflow


# MAIN FOLDERS
[ -z $STARK_FOLDER_CONFIG ] && export STARK_FOLDER_CONFIG=$STARK_FOLDER_ROOT/config
export STARK_FOLDER_ROOT=$(dirname $STARK_FOLDER_CONFIG)
#export STARK_FOLDER_CONFIG=$STARK_FOLDER_ROOT/config
export STARK_FOLDER_BIN=$STARK_FOLDER_ROOT/bin
export STARK_FOLDER_DOCS=$STARK_FOLDER_ROOT/docs
export STARK_FOLDER_TOOLBOX=$STARK_FOLDER_ROOT/toolbox
export STARK_FOLDER_APPS=$STARK_FOLDER_CONFIG/apps
export STARK_FOLDER_RULES=$STARK_FOLDER_CONFIG/rules


# FOLDER INFRASTRUCTURE

# MAIN STARK FOLDER
STARK_FOLDER_MAIN="/STARK"

# INPUT FOLDER
[ "$INPUT" != "" ] && FOLDER_INPUT=$INPUT && FOLDER_RUN="$FOLDER_INPUT/runs" && FOLDER_MANIFEST="$FOLDER_INPUT/manifests" && FOLDER_MANIFEST="$FOLDER_INPUT/pedigree"
if [ "$FOLDER_INPUT" == "" ]; then
	FOLDER_INPUT="$STARK_FOLDER_MAIN/input"
fi;

# RUN FOLDER
if [ "$FOLDER_RUN" == "" ]; then
	FOLDER_RUN="$FOLDER_INPUT/runs"
fi;
if [ "$FOLDER_MANIFEST" == "" ]; then
	FOLDER_MANIFEST="$FOLDER_INPUT/manifests"
fi;
if [ "$FOLDER_PEDIGREE" == "" ]; then
	FOLDER_PEDIGREE="$FOLDER_INPUT/pedigree"
fi;

# OUTPUT FOLDER
[ "$OUTPUT" != "" ] && FOLDER_OUTPUT=$OUTPUT
[ "$OUTPUT" != "" ] && FOLDER_OUTPUT=$OUTPUT && FOLDER_RESULTS="$FOLDER_OUTPUT/results" && FOLDER_DEMULTIPLEXING="$FOLDER_OUTPUT/demultiplexing"  && FOLDER_LOG="$FOLDER_OUTPUT/log"  && FOLDER_TMP="$FOLDER_OUTPUT/tmp"

if [ "$FOLDER_OUTPUT" == "" ]; then
	FOLDER_OUTPUT="$STARK_FOLDER_MAIN/output"
fi;

# RESULTS FOLDER
[ "$RESULTS" != "" ] && FOLDER_RESULTS=$RESULTS
if [ "$FOLDER_RESULTS" == "" ]; then
	FOLDER_RESULTS="$FOLDER_OUTPUT/results"
fi;

# DEMULTIPLEXING FOLDER
[ "$DEMULTIPLEXING" != "" ] && FOLDER_DEMULTIPLEXING=$DEMULTIPLEXING
if [ "$FOLDER_DEMULTIPLEXING" == "" ]; then
	FOLDER_DEMULTIPLEXING="$FOLDER_OUTPUT/demultiplexing"
fi;

# ANALYSIS/LOG FOLDER
[ "$LOG" != "" ] && FOLDER_LOG=$LOG
if [ "$FOLDER_LOG" == "" ]; then
	FOLDER_LOG="$FOLDER_OUTPUT/log"
fi;

# TMP FOLDER
[ "$TMP" != "" ] && FOLDER_TMP=$TMP
if [ "$FOLDER_TMP" == "" ]; then
	FOLDER_TMP="$FOLDER_OUTPUT/tmp"
fi;

# EXPORT FOLDER REPOSITORY
export FOLDER_REPOSITORY

# EXPORT FOLDER ARCHIVES
export FOLDER_ARCHIVES

# EXPORT FOLDER FAVORITES
export FOLDER_FAVORITES


# CONFIG FOLDERS

if [ "$FOLDER_TOOLS" == "" ]; then
	FOLDER_TOOLS="$STARK_FOLDER_MAIN/tools"
fi;

[ "$DATABASES" != "" ] && [ -d $DATABASES ] && FOLDER_DATABASES=$DATABASES
if [ "$FOLDER_DATABASES" == "" ]; then
	FOLDER_DATABASES="$STARK_FOLDER_MAIN/databases"
fi;


# RUNS FOLDER
export MISEQ_FOLDER=$FOLDER_RUN				# MISEQ RAW folder
MSR_SUBFOLDER=Data/Intensities/BaseCalls		# SUBFOLDER Illumina Sequencer Data subfolder

# MANIFEST FOLDER
export MANIFEST_FOLDER=$FOLDER_MANIFEST			# MANIFEST folder

# PEDIGREE FOLDER
export PEDIGREE_FOLDER=$FOLDER_PEDIGREE			# PEDIGREE folder


# RESULTS FOLDER
#export DEMULTIPLEXING_FOLDER=$FOLDER_RESULTS/DEM	# DEMULTIPLEXING folder
export DEMULTIPLEXING_FOLDER=$FOLDER_DEMULTIPLEXING	# DEMULTIPLEXING folder
export RESULTS_FOLDER=$FOLDER_RESULTS		# RESULTS folder ALL
#export ANALYSIS_FOLDER=$FOLDER_RESULTS/ANA		# ANALYSIS folder
export ANALYSIS_FOLDER=$FOLDER_LOG		# ANALYSIS folder
#export TMP_FOLDER_TMP=$FOLDER_RESULTS/TMP		# TEMPORARY folder tmp
export TMP_FOLDER_TMP=$FOLDER_TMP		# TEMPORARY folder tmp
#export VALIDATION_FOLDER=$FOLDER_RESULTS/VAL		# VALIDATION folder
#export VALIDATION_FOLDER=$FOLDER_VALIDATION		# VALIDATION folder
export TMP_FOLDER=$FOLDER_OUTPUT			# TEMPORARY folder = MAIN V2 folder
export MAIN_FOLDER=$FOLDER_OUTPUT			# MAIN V2 Folder (RES, DEM...)

# Create directories if does not exist
mkdir -p $DEMULTIPLEXING_FOLDER $RESULTS_FOLDER $ANALYSIS_FOLDER $TMP_FOLDER_TMP $TMP_FOLDER $MAIN_FOLDER 2>/dev/null
# TOOLS FOLDER
export NGS_TOOLS=$FOLDER_TOOLS				# TOOLS
export DBFOLDER=$FOLDER_DATABASES			# DB
export NGS_GENOMES=$DBFOLDER/genomes			# GENOMES folder
export APPS=$STARK_FOLDER_APPS				# APPS
#export GENOMES=$NGS_GENOMES				# Genomes folder in NGS folder
export TMP_SYS_FOLDER=/tmp				# TEMPORARY SYSTEM folder (tmpfs)
export NGS_FOLDER=$FOLDER_TOOLS
export NGS_SCRIPTS=$STARK_FOLDER_BIN			#"$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )"


# REPOSITORY AND ARCHIVES FOLDER
export RESULTS_FOLDER_BY_GROUP_PROJECT_COPY=$FOLDER_REPOSITORY	# Copy data into group and project (if defined)
export RESULTS_SUBFOLDER_DATA="STARK";				# Copy sample results in a SUBFOLDER

# CORE Files to copy to REPOSITORY
export REPOSITORY_FILE_PATTERNS_CORE='$SAMPLE.reports/$SAMPLE.final.vcf.gz $SAMPLE.reports/$SAMPLE.final.vcf.gz.tbi $SAMPLE.reports/$SAMPLE.final.Panel*.vcf.gz $SAMPLE.reports/$SAMPLE.final.Panel*.vcf.gz.tbi $SAMPLE.reports/*.*.config $SAMPLE.reports/*.report*.html $SAMPLE.reports/*.report.html.folder:FOLDER $SAMPLE.*.bam.metrics/$SAMPLE.*.validation.flags.Panel*.bed $SAMPLE*.tag $SAMPLE.reports/$SAMPLE.final.Panel*.tsv ';
# CORE Files to exclude to REPOSITORY results subfolder
export REPOSITORY_FILE_SUBFOLDER_PATTERNS_CORE='';
# CORE Files to copy to ARCHIVES
export ARCHIVES_FILE_PATTERNS_CORE='$SAMPLE.reports/$SAMPLE.final.vcf.gz $SAMPLE.reports/$SAMPLE.final.vcf.gz.tbi $SAMPLE.reports/$SAMPLE.final.Panel*.vcf.gz $SAMPLE.reports/$SAMPLE.final.Panel*.vcf.gz.tbi $SAMPLE.reports/*.*.config $SAMPLE.reports/*.report*.html $SAMPLE.reports/*.report.html.folder:FOLDER $SAMPLE.*.bam.metrics/$SAMPLE.*.validation.flags.Panel*.bed $SAMPLE*.tag $SAMPLE.bed $SAMPLE.manifest $SAMPLE*.genes $SAMPLE*.transcripts  $SAMPLE.archive.cram $SAMPLE.archive.cram.crai $SAMPLE.analysis.json';
# CORE Files to copy to FAVORITES
export FAVORITES_FILE_PATTERNS_CORE='';

# Copy some sample results files in the root sample folder, if any SUBFOLDER defined
export REPOSITORY_FILE_PATTERNS=$(echo $REPOSITORY_FILE_PATTERNS' '$REPOSITORY_FILE_PATTERNS_CORE | tr "," " " | tr " " "\n" | sort -u | tr "\n" " ");
export REPOSITORY_FILE_SUBFOLDER_PATTERNS=$(echo $REPOSITORY_FILE_SUBFOLDER_PATTERNS' '$REPOSITORY_FILE_SUBFOLDER_PATTERNS_CORE | tr "," " " | tr " " "\n" | sort -u | tr "\n" " ");
export ARCHIVES_FILE_PATTERNS=$(echo $ARCHIVES_FILE_PATTERNS' '$ARCHIVES_FILE_PATTERNS_CORE | tr "," " " | tr " " "\n" | sort -u | tr "\n" " ");
export FAVORITES_FILE_PATTERNS=$(echo $FAVORITES_FILE_PATTERNS' '$FAVORITES_FILE_PATTERNS_CORE | tr "," " " | tr " " "\n" | sort -u | tr "\n" " ");

# Useful file patterns to add in repository (variable REPOSITORY_FILE_PATTERNS)
# $SAMPLE.reports/$SAMPLE.full.Design.vcf.gz $SAMPLE.reports/$SAMPLE.full.Design.tsv $SAMPLE.*.validation.bam $SAMPLE.*.validation.bam.bai
# Useful file patterns to add in archives (variable ARCHIVES_FILE_PATTERNS)
# $SAMPLE.reports/$SAMPLE.full.vcf.gz $SAMPLE.reports/$SAMPLE.final.tsv $SAMPLE.*.bam.metrics/$SAMPLE.*.validation.flags.*.bed


# IGV SESSION
# IGV session xml file (*igv_session.xml) is generated with Samples and Runs files (BAM, VCF, BED)
# Parameters select files through patterns and folder depth
# Be aware to coordonate these parameters with availabled files within repositories folders

# SAMPLE
IGV_SESSION_PATTERNS_SAMPLE_DEFAULT=' $SAMPLE.final.vcf.gz $SAMPLE*Design*bed $SAMPLE*Panel*bed *Design*vcf.gz *Panel*vcf.gz *bam *cram '
IGV_SESSION_MINDEPTH_SAMPLE_DEFAULT=0
IGV_SESSION_MAXDEPTH_SAMPLE_DEFAULT=1

[ "$REPOSITORY_SAMPLE_IGV_SESSION_PATTERNS" == "" ] && REPOSITORY_SAMPLE_IGV_SESSION_PATTERNS=$IGV_SESSION_PATTERNS_SAMPLE_DEFAULT
export REPOSITORY_SAMPLE_IGV_SESSION_PATTERNS
[ "$REPOSITORY_SAMPLE_IGV_SESSION_MINDEPTH" == "" ] && REPOSITORY_SAMPLE_IGV_SESSION_MINDEPTH=$IGV_SESSION_MINDEPTH_SAMPLE_DEFAULT
export REPOSITORY_SAMPLE_IGV_SESSION_MINDEPTH
[ "$REPOSITORY_SAMPLE_IGV_SESSION_MAXDEPTH" == "" ] && REPOSITORY_SAMPLE_IGV_SESSION_MAXDEPTH=$IGV_SESSION_MAXDEPTH_SAMPLE_DEFAULT
export REPOSITORY_SAMPLE_IGV_SESSION_MAXDEPTH

[ "$ARCHIVES_SAMPLE_IGV_SESSION_PATTERNS" == "" ] && ARCHIVES_SAMPLE_IGV_SESSION_PATTERNS=$IGV_SESSION_PATTERNS_SAMPLE_DEFAULT
export ARCHIVES_SAMPLE_IGV_SESSION_PATTERNS
[ "$ARCHIVES_SAMPLE_IGV_SESSION_MINDEPTH" == "" ] && ARCHIVES_SAMPLE_IGV_SESSION_MINDEPTH=$IGV_SESSION_MINDEPTH_SAMPLE_DEFAULT
export ARCHIVES_SAMPLE_IGV_SESSION_MINDEPTH
[ "$ARCHIVES_SAMPLE_IGV_SESSION_MAXDEPTH" == "" ] && ARCHIVES_SAMPLE_IGV_SESSION_MAXDEPTH=$IGV_SESSION_MAXDEPTH_SAMPLE_DEFAULT
export ARCHIVES_SAMPLE_IGV_SESSION_MAXDEPTH

[ "$FAVORITES_SAMPLE_IGV_SESSION_PATTERNS" == "" ] && FAVORITES_SAMPLE_IGV_SESSION_PATTERNS=$IGV_SESSION_PATTERNS_SAMPLE_DEFAULT
export FAVORITES_SAMPLE_IGV_SESSION_PATTERNS
[ "$FAVORITES_SAMPLE_IGV_SESSION_MINDEPTH" == "" ] && FAVORITES_SAMPLE_IGV_SESSION_MINDEPTH=$IGV_SESSION_MINDEPTH_SAMPLE_DEFAULT
export FAVORITES_SAMPLE_IGV_SESSION_MINDEPTH
[ "$FAVORITES_SAMPLE_IGV_SESSION_MAXDEPTH" == "" ] && FAVORITES_SAMPLE_IGV_SESSION_MAXDEPTH=$IGV_SESSION_MAXDEPTH_SAMPLE_DEFAULT
export FAVORITES_SAMPLE_IGV_SESSION_MAXDEPTH

# RUN
IGV_SESSION_PATTERNS_RUN_DEFAULT=' *.final.vcf.gz *Design*bed *Panel*bed *Design*vcf.gz *Panel*vcf.gz *bam *cram '
IGV_SESSION_MINDEPTH_RUN_DEFAULT=0
IGV_SESSION_MAXDEPTH_RUN_DEFAULT=2

[ "$REPOSITORY_RUN_IGV_SESSION_PATTERNS" == "" ] && REPOSITORY_RUN_IGV_SESSION_PATTERNS=$IGV_SESSION_PATTERNS_RUN_DEFAULT
export REPOSITORY_RUN_IGV_SESSION_PATTERNS
[ "$REPOSITORY_RUN_IGV_SESSION_MINDEPTH" == "" ] && REPOSITORY_RUN_IGV_SESSION_MINDEPTH=$IGV_SESSION_MINDEPTH_RUN_DEFAULT
export REPOSITORY_RUN_IGV_SESSION_MINDEPTH
[ "$REPOSITORY_RUN_IGV_SESSION_MAXDEPTH" == "" ] && REPOSITORY_RUN_IGV_SESSION_MAXDEPTH=$IGV_SESSION_MAXDEPTH_RUN_DEFAULT
export REPOSITORY_RUN_IGV_SESSION_MAXDEPTH

[ "$ARCHIVES_RUN_IGV_SESSION_PATTERNS" == "" ] && ARCHIVES_RUN_IGV_SESSION_PATTERNS=$IGV_SESSION_PATTERNS_RUN_DEFAULT
export ARCHIVES_RUN_IGV_SESSION_PATTERNS
[ "$ARCHIVES_RUN_IGV_SESSION_MINDEPTH" == "" ] && ARCHIVES_RUN_IGV_SESSION_MINDEPTH=$IGV_SESSION_MINDEPTH_RUN_DEFAULT
export ARCHIVES_RUN_IGV_SESSION_MINDEPTH
[ "$ARCHIVES_RUN_IGV_SESSION_MAXDEPTH" == "" ] && ARCHIVES_RUN_IGV_SESSION_MAXDEPTH=$IGV_SESSION_MAXDEPTH_RUN_DEFAULT
export ARCHIVES_RUN_IGV_SESSION_MAXDEPTH

[ "$FAVORITES_RUN_IGV_SESSION_PATTERNS" == "" ] && FAVORITES_RUN_IGV_SESSION_PATTERNS=$IGV_SESSION_PATTERNS_RUN_DEFAULT
export FAVORITES_RUN_IGV_SESSION_PATTERNS
[ "$FAVORITES_RUN_IGV_SESSION_MINDEPTH" == "" ] && FAVORITES_RUN_IGV_SESSION_MINDEPTH=$IGV_SESSION_MINDEPTH_RUN_DEFAULT
export FAVORITES_RUN_IGV_SESSION_MINDEPTH
[ "$FAVORITES_RUN_IGV_SESSION_MAXDEPTH" == "" ] && FAVORITES_RUN_IGV_SESSION_MAXDEPTH=$IGV_SESSION_MAXDEPTH_RUN_DEFAULT
export FAVORITES_RUN_IGV_SESSION_MAXDEPTH


# IGV Display Mode
# IGV Display mode for files found in patterns
# Depend on file format: either BEM, VCF or BED formats
# There are 3 different options for viewing the feature track.  These allow you to display overlapping features, such as different transcripts of a gene, on one line or multiple lines
# Dislpay mode options available in IGV: SQUISHED, COLLAPSED and EXPANDED
DISPLAYMODE_ARRAY="SQUISHED COLLAPSED EXPANDED"
if [ -z "$DISPLAYMODE_BAM" ] || ! in_array $DISPLAYMODE_BAM $DISPLAYMODE_ARRAY; then
	DISPLAYMODE_BAM="SQUISHED"
fi;
export DISPLAYMODE_BAM
if [ -z "$DISPLAYMODE_VCF" ] || ! in_array $DISPLAYMODE_VCF $DISPLAYMODE_ARRAY; then
	DISPLAYMODE_VCF="COLLAPSED"
fi;
export DISPLAYMODE_VCF
if [ -z "$DISPLAYMODE_BED" ] || ! in_array $DISPLAYMODE_BED $DISPLAYMODE_ARRAY; then
	DISPLAYMODE_BED="COLLAPSED"
fi;
export DISPLAYMODE_BED


# IGV SESSION DB
# Additionnal databases for IGV session (see IGV doc)
# Example:
# IGV_SESSION_RESSOURCES='
# 	<Resource index="https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz.tbi" path="https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz" type="vcf"/>
# 	<Resource index="http://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi" path="http://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz" type="vcf"/>
# 	<Resource name="GC Percentage" path="http://www.broadinstitute.org/igvdata/annotations/hg19/hg19.gc5base.tdf" type="tdf"/>
# 	<Resource index="https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz.tbi" name="Refseq Genes" path="https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz" type="refgene"/>
# 	<Resource hyperlink="http://www.gencodegenes.org/" index="https://s3.amazonaws.com/igv.broadinstitute.org/annotations/hg19/genes/gencode.v18.annotation.sorted.gtf.gz.tbi" name="Gencode Genes (v18)" path="https://s3.amazonaws.com/igv.broadinstitute.org/annotations/hg19/genes/gencode.v18.annotation.sorted.gtf.gz" trackLine="visibilitywindow:10000000" type="gtf"/>
# '
# IGV_SESSION_DATAPANEL='
# 	<Track attributeKey="00-All.vcf.gz" clazz="org.broad.igv.variant.VariantTrack" featureVisibilityWindow="100100" fontSize="10" groupByStrand="false" id="https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz" name="dbSNP" siteColorMode="ALLELE_FREQUENCY" displayMode="COLLAPSE" visible="true"/>
# 	<Track attributeKey="clinvar.vcf.gz" clazz="org.broad.igv.variant.VariantTrack" featureVisibilityWindow="100100" fontSize="10" groupByStrand="false" id="http://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz" name="ClinVar" siteColorMode="ALLELE_FREQUENCY" displayMode="COLLAPSE" visible="true"/>
# '
# IGV_SESSION_FEATUREPANEL='
# 	<Track attributeKey="GC Percentage" autoScale="false" clazz="org.broad.igv.track.DataSourceTrack" colorScale="ContinuousColorScale;0.0;80.0;255,255,255;0,0,178" fontSize="10" height="20" id="http://www.broadinstitute.org/igvdata/annotations/hg19/hg19.gc5base.tdf" name="GC Percentage" renderer="HEATMAP" visible="true" windowFunction="mean">
# 		<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="80.0" minimum="0.0" type="LINEAR"/>
# 	</Track>
# 	<Track attributeKey="Refseq Genes" clazz="org.broad.igv.track.FeatureTrack" featureVisibilityWindow="63761233" fontSize="10" groupByStrand="false" id="https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz" name="Refseq Genes" visible="true"/>
# 	<Track attributeKey="Gencode Genes (v18)" clazz="org.broad.igv.track.FeatureTrack" featureVisibilityWindow="1000000" fontSize="10" groupByStrand="false" id="https://s3.amazonaws.com/igv.broadinstitute.org/annotations/hg19/genes/gencode.v18.annotation.sorted.gtf.gz" name="Gencode Genes (v18)" visible="true"/>
# '
# Default Refseq Gene
if [ "$IGV_SESSION_RESSOURCES" == "" ] && [ "$IGV_SESSION_DATAPANEL" == "" ] && [ "$IGV_SESSION_FEATUREPANEL" == "" ]; then

	IGV_SESSION_RESSOURCES='<Resource index="https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz.tbi" name="Refseq Genes" path="https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz" type="refgene"/>'
	IGV_SESSION_DATAPANEL=''
	IGV_SESSION_FEATUREPANEL='<Track attributeKey="Refseq Genes" clazz="org.broad.igv.track.FeatureTrack" featureVisibilityWindow="63761233" fontSize="10" groupByStrand="false" id="https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz" name="Refseq Genes" visible="true"/>'

fi;

export IGV_SESSION_RESSOURCES
export IGV_SESSION_DATAPANEL
export IGV_SESSION_FEATUREPANEL


# IGV SESSION JSON
# Example:
# IGV_SESSION_JSON_TRACKS_ADDITIONAL='
# {
# 	"type": "bed",
# 	"url": "https://s3.amazonaws.com/igv.org.test/data/gencode.v18.collapsed.bed",
# 	"indexURL": "https://s3.amazonaws.com/igv.org.test/data/gencode.v18.collapsed.bed.idx",
# 	"name": "Gencode V18",
#   "type": "annotation",
#   "displayMode": "COLLAPSED"
# }'
if [ "$IGV_SESSION_DAS" == "" ]; then
	IGV_SESSION_DAS="http://localhost:4201/static/data/public/repositories"
fi;
if [ "$IGV_SESSION_DAS_REPOSITORIES" == "" ]; then
	IGV_SESSION_DAS_REPOSITORIES="$IGV_SESSION_DAS/Repository $IGV_SESSION_DAS/Archives $IGV_SESSION_DAS/Favorites"
fi;
if [ "$IGV_SESSION_JSON_TRACKS_ADDITIONAL" == "" ]; then
	IGV_SESSION_JSON_TRACKS_ADDITIONAL=""
fi;
export IGV_SESSION_DAS
export IGV_SESSION_DAS_REPOSITORIES
export IGV_SESSION_JSON_TRACKS_ADDITIONAL



#Rules defined in the app
for RA in $(echo "$RULES_APP" | tr "," " " | tr " " "\n"); do
	RULES_APP=$RULES_APP" "$(find $STARK_FOLDER_APPS/$RA -maxdepth 0 -type f 2>/dev/null)
done;

# RULES for the application
[ "$APP_NAME" != "" ] && RULES_APP="$RULES_APP $STARK_FOLDER_RULES/$APP_NAME.app/*rules.mk";

# RULES for the group
[ "$APP_GROUP" != "" ] && RULES_APP="$RULES_APP $STARK_FOLDER_RULES/$APP_GROUP/*rules.mk";

# RULES for the project
[ "$APP_PROJECT" != "" ] && RULES_APP="$RULES_APP $STARK_FOLDER_RULES/$APP_PROJECT/*rules.mk";

#SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )"
#echo $SCRIPT_DIR
#echo $APP_NAME
#echo $APP_GROUP
#echo $APP_PROJECT

# TOOLS
#########

#source $STARK_FOLDER_CONFIG/tools.app
source $CONFIG_TOOLS

# DATABASES
#############

#source $STARK_FOLDER_CONFIG/databases.app
source $CONFIG_DATABASES

# databases release
DATABASES_RELEASE_FILE=$FOLDER_DATABASES/RELEASE
[ -s $DATABASES_RELEASE_FILE ] && DATABASES_RELEASE=$(grep "^RELEASE" $DATABASES_RELEASE_FILE | cut -f2) || DATABASES_RELEASE="UNKNOWN";
export DATABASES_RELEASE


# PIPELINES
##############


# PIPELINES (default "bwamem.gatkHC.howard")
# pipelines to use for the analysis
# Format: "ALIGNER1.CALLER1.ANNOTATOR ALIGNER1.CALLER2.ANNOTATOR1 ALIGNER2.CALLER1.ANNOTATOR1"
# If variables ALIGNERS, CALLERS, and ANNOTATORS are defined (see below), PIPELINES variable will be automatically generated
# This variable can be a additionnal pipeline that those defined by the combinasion of ALIGNERS, CALLERS, and ANNOTATORS
# If no pipeline is finally defined, the default pipeline will be applied
#PIPELINES=""

# ALIGNERS (default "")
# Aligners to use for the analysis
# Example of available aligners: bwamem  bwasw bwaaln
#ALIGNERS=""
export ALIGNERS

# CALLERS (default "")
# Callers to use for the analysis
# Example of available callers: gatkHC gatkUG VarScan samtools
#CALLERS=""
export CALLERS

# ANNOTATORS (default "")
# Annotator to use for the analysis
# Example of available annotators: howard snpeff
#ANNOTATORS=""
export ANNOTATORS

# CREATE PIPELINES from ALIGNERS/CALLERS/ANNOTATORS
for ALIGNER in $ALIGNERS; do
	for CALLER in $CALLERS; do
		for ANNOTATOR in $ANNOTATORS; do
			PIPELINES="$PIPELINES $ALIGNER.$CALLER.$ANNOTATOR"
		done;
	done;
done;
#PIPELINES=$(echo $PIPELINES | tr "," "\n" | tr " " "\n" | sort | uniq | tr "\n" " " | sed -e 's/^[[:space:]]*//'  | sed -e 's/[[:space:]]*$//')
PIPELINES=$(echo $PIPELINES | tr "," "\n" | tr " " "\n" | awk '!x[$0]++' | tr "\n" " " | sed -e 's/^[[:space:]]*//'  | sed -e 's/[[:space:]]*$//')
if [ -z "$PIPELINES" ]; then
	PIPELINES="bwamem.gatkHC.howard"
fi;
export PIPELINES

# MANIFEST
export MANIFEST

# BLANK samples
export BLANK

# BARCODE_MISMATCHES (default 1)
# Used to demultiplex
if [ -z $BARCODE_MISMATCHES ]; then
	BARCODE_MISMATCHES=1
fi;
export BARCODE_MISMATCHES

# BAM METRICS (default 1/TRUE/YES/Y)
# Performs BAM METRICS (1/TRUE/YES/Y or 0/FALSE/NO/N). Time and space consuming. Switch off for exome/genome for better perfomances
if [ "${BAM_METRICS^^}" == "FALSE" ] || [ "${BAM_METRICS^^}" == "NO" ] || [ "${BAM_METRICS^^}" == "N" ] || [ "$BAM_METRICS" == "0" ]; then
	BAM_METRICS=0
elif [ -z $BAM_METRICS ] || [ "${BAM_METRICS^^}" == "TRUE" ] || [ "${BAM_METRICS^^}" == "YES" ] || [ "${BAM_METRICS^^}" == "Y" ] || [ "$BAM_METRICS" == "1" ]; then
	BAM_METRICS=1
else
	BAM_METRICS=1
fi;
export BAM_METRICS


# GLOBAL METRICS VARIABLES
# Values for BAM metrics
# Minimum mapping quality to consider in the metrics BAM
if [ -z "$METRICS_MINIMUM_MAPPING_QUALITY" ] || ! [[ $METRICS_MINIMUM_MAPPING_QUALITY =~ ^[0-9]+$ ]]; then
	METRICS_MINIMUM_MAPPING_QUALITY=20
fi;
export METRICS_MINIMUM_MAPPING_QUALITY

# Minimum bases quality to consider in the metrics BAM
if [ -z "$METRICS_MINIMUM_BASE_QUALITY" ] || ! [[ $METRICS_MINIMUM_BASE_QUALITY =~ ^[0-9]+$ ]]; then
	METRICS_MINIMUM_BASE_QUALITY=20
fi;
export METRICS_MINIMUM_BASE_QUALITY

# Clipping overlapping reads in the metrics BAM
if [ -z "$CLIP_OVERLAPPING_READS" ]; then
	CLIP_OVERLAPPING_READS=1
fi;
export CLIP_OVERLAPPING_READS

# Flagged reads in the metrics BAM (mpileup format)
# default UNMAP,SECONDARY,QCFAIL,DUP
if [ -z "$METRICS_FLAGS" ]; then
	METRICS_FLAGS="UNMAP,SECONDARY,QCFAIL,DUP"
fi;
METRICS_FLAGS=$(echo $METRICS_FLAGS | tr -d " ")
export METRICS_FLAGS

# Flagged reads in the metrics BAM (samtools format)
# Generated from mpileup format (if empty)
if [ -z "$SAMTOOLS_METRICS_FLAG_PARAM" ]; then
	#SAMTOOLS_METRICS_FLAG_PARAM=" -F 0x4 -F 0x100 -F 0x200 -F 0x400"
	SAMTOOLS_METRICS_FLAG_PARAM=""
	for F in $(echo $METRICS_FLAGS | tr "," " "); do
		[ "$F" == "UNMAP" ] && SAMTOOLS_METRICS_FLAG_PARAM=$SAMTOOLS_METRICS_FLAG_PARAM" -F 0x4 "
		[ "$F" == "SECONDARY" ] && SAMTOOLS_METRICS_FLAG_PARAM=$SAMTOOLS_METRICS_FLAG_PARAM" -F 0x100 "
		[ "$F" == "QCFAIL" ] && SAMTOOLS_METRICS_FLAG_PARAM=$SAMTOOLS_METRICS_FLAG_PARAM" -F 0x200 "
		[ "$F" == "DUP" ] && SAMTOOLS_METRICS_FLAG_PARAM=$SAMTOOLS_METRICS_FLAG_PARAM" -F 0x400 "
	done;
fi;
export SAMTOOLS_METRICS_FLAG_PARAM


# VARIANT RECALIBRATION (default 0/FALSE/NO/N)
# Performs Variant recalibration after Calling (1/TRUE/YES/Y or 0/FALSE/NO/N).
# If the recalibration fail (usually due to lack of data for statistic calculation), nothing will be done
# if [ -z "$VARIANT_RECALIBRATION" ] || [ "${VARIANT_RECALIBRATION^^}" == "FALSE" ] || [ "${VARIANT_RECALIBRATION^^}" == "NO" ] || [ "${VARIANT_RECALIBRATION^^}" == "N" ]  || [ "$VARIANT_RECALIBRATION" == "0" ]; then
# 	VARIANT_RECALIBRATION=0
# elif [ "${VARIANT_RECALIBRATION^^}" == "TRUE" ] || [ "${VARIANT_RECALIBRATION^^}" == "YES" ] || [ "${VARIANT_RECALIBRATION^^}" == "Y" ]  || [ "$VARIANT_RECALIBRATION" == "1" ]; then
# 	VARIANT_RECALIBRATION=1

# else
# 	VARIANT_RECALIBRATION=0
# fi;
# export VARIANT_RECALIBRATION

# INTERVAL_PADDING (default 0)
# Add some “padding” to the intervals used (manifest) in order to include the flanking regions (typically ~100 bp)
if [ -z "$INTERVAL_PADDING" ] || ! [[ $INTERVAL_PADDING =~ ^[0-9]+$ ]]; then
	INTERVAL_PADDING=0
fi;
export INTERVAL_PADDING

# COVERAGE CRITERIA (default "1,30")
# For gene coverage metrics
# the criteria to calculate the percent of bases over $COVERAGE_CRITERIA X (eg 30 for 30X)
if [ -z "$COVERAGE_CRITERIA" ]; then
	COVERAGE_CRITERIA="1,5,10,20,30,50,100,200,300"
fi;
export COVERAGE_CRITERIA


# COVERAGE DP THRESHOLD (default "30" "100" "1")
# For gene coverage metrics
# the criteria test if genes failed (or just warning) the coverage threshold
# Sequencing depth threshold
if [ -z "$SEQUENCING_DEPTH" ]; then
	SEQUENCING_DEPTH="1"
fi;
export SEQUENCING_DEPTH
# Sequencing coverage threshold
if [ -z "$SEQUENCING_COVERAGE_THRESHOLD" ]; then
	SEQUENCING_COVERAGE_THRESHOLD="1"
fi;
export SEQUENCING_COVERAGE_THRESHOLD
# fail DP threshold (default 30X)
if [ -z "$MINIMUM_DEPTH" ]; then
	MINIMUM_DEPTH="30"
fi;
export MINIMUM_DEPTH
# warn DP threshold (default 100X)
if [ -z "$EXPECTED_DEPTH" ]; then
	EXPECTED_DEPTH="100"
fi;
export EXPECTED_DEPTH
# threshold percentage of bases over the DP threshold
if [ -z "$DEPTH_COVERAGE_THRESHOLD" ]; then
	DEPTH_COVERAGE_THRESHOLD="0.95"
fi;
export DEPTH_COVERAGE_THRESHOLD


# COVERAGE HARMONIZATION
COVERAGE_CRITERIA=$(echo "$SEQUENCING_DEPTH,$MINIMUM_DEPTH,$EXPECTED_DEPTH,$COVERAGE_CRITERIA" | tr "," "\n" | tr " " "\n" | grep -v "^$" | sort -u -n | tr "\n" ","  | sed "s/,$//")


# CLIP_OVERLAPPING_READS (default 1)
# From PICARD: For paired reads, soft clip the 3' end of each read if necessary so that it does not extend past the 5' end of its mate
if [ -z "$CLIP_OVERLAPPING_READS" ]; then
	CLIP_OVERLAPPING_READS=1
fi;
export CLIP_OVERLAPPING_READS


# NB_BASES_AROUND (default 0)
# For gene coverage metrics
# the number of bases to look around the exons from the given bed file
if [ -z $NB_BASES_AROUND ] || ! [[ $NB_BASES_AROUND =~ ^[0-9]+$ ]]; then
	NB_BASES_AROUND=0
fi;
export NB_BASES_AROUND


# GENESCOVERAGE_PRECISION (default 2)
# Genes Coverage calculation precision
if [ -z $GENESCOVERAGE_PRECISION ] || ! [[ $GENESCOVERAGE_PRECISION =~ ^[0-9]+$ ]]; then
	GENESCOVERAGE_PRECISION=2
fi;
export GENESCOVERAGE_PRECISION


# BEDFILE_GENES (default "")
# For gene coverage metrics
# the bed file containing the 5'UTR, 3'UTR and genomic coding coordinates.
#BEDFILE_GENES=""
#if [ "$BEDFILE_GENES" != "" ] && [ ! -s $BEDFILE_GENES ]; then
	export BEDFILE_GENES
#fi;

# VARANK ANALYSIS (default 0/FALSE/NO/N)
# Performs VARANK ANALYSIS with Alamut (1/TRUE/YES/Y or 0/FALSE/NO/N).
if [ -z $VARANK_ANALYSIS ] || [ "${VARANK_ANALYSIS^^}" == "FALSE" ] || [ "${VARANK_ANALYSIS^^}" == "NO" ] || [ "${VARANK_ANALYSIS^^}" == "N" ]  || [ "$VARANK_ANALYSIS" == "0" ]; then
	VARANK_ANALYSIS=0
elif [ "${VARANK_ANALYSIS^^}" == "TRUE" ] || [ "${VARANK_ANALYSIS^^}" == "YES" ] || [ "${VARANK_ANALYSIS^^}" == "Y" ]  || [ "$VARANK_ANALYSIS" == "1" ]; then
	VARANK_ANALYSIS=1
else
	VARANK_ANALYSIS=0
fi;
export VARANK_ANALYSIS
# if yes, define a folder to store the results
#VARANK_FOLDER=$FOLDER_RESULTS/VARANK
export VARANK_FOLDER

# BAM CHECK (default 0)
# Check BAM for each manipulation step (clipping, realignment...)
# Time consuming, but will stop the analysis in case of missing reads.
# A Metrics on each BAM check the BAM decrependy in any case
if [ -z $BAM_CHECK_STEPS ] || [ "${BAM_CHECK_STEPS^^}" == "FALSE" ] || [ "${BAM_CHECK_STEPS^^}" == "NO" ] || [ "${BAM_CHECK_STEPS^^}" == "N" ]  || [ "$BAM_CHECK_STEPS" == "0" ]; then
	BAM_CHECK_STEPS=0
elif [ "${BAM_CHECK_STEPS^^}" == "TRUE" ] || [ "${BAM_CHECK_STEPS^^}" == "YES" ] || [ "${BAM_CHECK_STEPS^^}" == "Y" ]  || [ "$BAM_CHECK_STEPS" == "1" ]; then
	BAM_CHECK_STEPS=1
else
	BAM_CHECK_STEPS=0
fi;
export BAM_CHECK_STEPS



# METRICS SNPEFF (default 0)
# Generate snpEff variant metrics from VCF
# Only for Report final VCF
if [ -z $METRICS_SNPEFF ] || [ "${METRICS_SNPEFF^^}" == "FALSE" ] || [ "${METRICS_SNPEFF^^}" == "NO" ] || [ "${METRICS_SNPEFF^^}" == "N" ]  || [ "$METRICS_SNPEFF" == "0" ]; then
	METRICS_SNPEFF=0
elif [ "${METRICS_SNPEFF^^}" == "TRUE" ] || [ "${METRICS_SNPEFF^^}" == "YES" ] || [ "${METRICS_SNPEFF^^}" == "Y" ]  || [ "$METRICS_SNPEFF" == "1" ]; then
	METRICS_SNPEFF=1
else
	METRICS_SNPEFF=0
fi;
export METRICS_SNPEFF



# PIPELINES PRIORITIZATION
# List of pipelines to prioritize for the report (final.vcf)
#PRIORITIZE_PIPELINES_LIST=""
export PRIORITIZE_PIPELINES_LIST



### FASTQ processing

# Set mask for demultiplexing
# e.g. "" (auto from SampleSheet), "Y150,I10,Y10,Y150", "Y150,I8,Y10,Y150" (UMI index2)
# default ""
export STARK_DEMULTIPLEXING_BASES_MASK

# Set short read size for demultiplexing
# If demultiplexing UMI within a read must be set to 0 for Agilent XTHS kits
# default "" (auto)
export STARK_DEMULTIPLEXING_MASK_SHORT_ADAPTATER_READ

# Set read mapping
# Redefine FASTQ files in order to identify R1, R2, I1 and I2
# Order: R1 I1 I2 R2
# e.g. "" (default, corresponding to "R1 I1 I2 R2"), "R1 I1 R2 R3" (UMI index2)
# default "R1 I1 I2 R2"
export STARK_DEMULTIPLEXING_READS_MAPPING

# Demultiplexing adaptated stringency
# For BCL2FASTQ demultiplexing (see doc)
[ "$ADAPTER_STRINGENCY" == "" ] && ADAPTER_STRINGENCY=0.9
export ADAPTER_STRINGENCY

# Demultiplexing options
# For BCL2FASTQ demultiplexing (see doc)
# Usually: "--no-lane-splitting --create-fastq-for-index-reads"
# Default: ""
export STARK_DEMULTIPLEXING_BCL2FASTQ_OPTIONS

# FASTQ compression level for demultiplexing FASTQ files
# zlib compression level (1-9) used for FASTQ files during demultiplexing
# Used by BCL2FASTQ
# If FASTQ_DEMULTIPLEXING_KEEP=1, we suggest a high level of compression (at least 5)
[ "$FASTQ_DEMULTIPLEXING_COMPRESSION_LEVEL" == "" ] && FASTQ_DEMULTIPLEXING_COMPRESSION_LEVEL=1
export FASTQ_DEMULTIPLEXING_COMPRESSION_LEVEL

# FASTQ compression level for main FASTQ files
# zlib compression level (1-9) used for FASTQ files
# Used by FASTP
[ "$FASTQ_COMPRESSION_LEVEL" == "" ] && FASTQ_COMPRESSION_LEVEL=5
export FASTQ_COMPRESSION_LEVEL

# ENABLE_ADAPTER_TRIMMING
# Trim adapter and autodetect adapter for paired end
# Either 0 or 1
# Default: 0 (i.e. adapter trimming is disable)
[ "$ENABLE_ADAPTER_TRIMMING" == "" ] && ENABLE_ADAPTER_TRIMMING=0
export ENABLE_ADAPTER_TRIMMING

# FASTQ Read quality filtering
# Read Quality threshold. Read quality below will be removed
# Default: 0 (i.e. disable)
[ "$FASTQ_QUALITY_FILTERING" == "" ] && FASTQ_QUALITY_FILTERING=0
export FASTQ_QUALITY_FILTERING

# polyG tail trimming
# force polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
# the minimum length to detect polyG in the read tail.
# Default: 0 (i.e. disable)
[ "$POLY_G_MIN_LEN" == "" ] && POLY_G_MIN_LEN=0
export POLY_G_MIN_LEN

# Read length filtering
# reads shorter than length_required will be discarded
# Default: 0 (i.e. disable)
[ "$POLY_G_MIN_LEN" == "" ] && POLY_G_MIN_LEN=0
export READ_LENGTH_REQUIRED

# UMI extract location
# Set the UMI location
# If not null, NO UMI extraction and analysis
# Available locations: 
#    index1: the first index is used as UMI. If the data is PE, this UMI will be used for both read1/read2.
#    index2 the second index is used as UMI. PE data only, this UMI will be used for both read1/read2.
#    read1 the head of read1 is used as UMI. If the data is PE, this UMI will be used for both read1/read2.
#    read2 the head of read2 is used as UMI. PE data only, this UMI will be used for both read1/read2.
#    per_index read1 will use UMI extracted from index1, read2 will use UMI extracted from index2.
#    per_read read1 will use UMI extracted from the head of read1, read2 will use UMI extracted from the head of read2.
# UMI_TYPE will be calculated from UMI_LOC ("simplex" for index1, index2, read1, read2, and "duplex" for per_index and per_read)
# e.g.: UMI_LOC="index2"
# See FASTP/UMI TOOLS documentation for more information
UMI_LOC=${UMI_LOC,,}
if [ "$UMI_LOC" == "" ] || [ "$UMI_LOC" == "index1" ] || [ "$UMI_LOC" == "index2" ] || [ "$UMI_LOC" == "per_index" ] || [ "$UMI_LOC" == "read1" ] || [ "$UMI_LOC" == "read2" ] || [ "$UMI_LOC" == "per_read" ]; then
	if [ "$UMI_LOC" == "" ]; then
		UMI_TYPE=""
	elif [[ "$UMI_LOC" =~ ^per_* ]]; then
		UMI_TYPE="duplex"
	else
		UMI_TYPE="simplex"
	fi;
else
	echo "#[ERROR] UMI LOC '$UMI_LOC' not available"
	return 1
fi;
export UMI_LOC
export UMI_TYPE

# UMI extract tag
# Set the UMI Barcode pattern
# If not null, STARK will prepare fastq containg UMIs +/- cell barcodes for alignment
# e.g.: UMI_BARCODE_PATTERN="NNNNNNNNNN" for simplex
# e.g.: UMI_BARCODE_PATTERN="NNNNN-NNNNN" for duplex
# if UMI_LOC is "per_index" or "per_read", and UMI_BARCODE_PATTERN is defined as "NNNNN", it with be redefined as "NNNNN-NNNNN"
# Only length of the first part of the duplex barcode will be used with FASTP
# See FASTP/UMI TOOLS documentation for more information
if [ "$UMI_LOC" == "per_index" ] || [ "$UMI_LOC" == "per_read" ]; then
	if ! [[ "$UMI_BARCODE_PATTERN" =~ .*"-".* ]]; then
		UMI_BARCODE_PATTERN="$UMI_BARCODE_PATTERN-$UMI_BARCODE_PATTERN"
	fi;
fi; 
export UMI_BARCODE_PATTERN

# Barcode tag
# Barcode to use for Mark Duplicates
# If not null, Mark Duplicates will consider this tag (default null)
# e.g.: BARCODE_TAG="BC" for 10X Genomics, BARCODE_TAG="BX" for UMI
# See PICARD documentation for more information
export BARCODE_TAG

# set variable to "READ_NAME_REGEX=null" to disable optical deduplication
# PICARD_MARKDUP_OPTICAL_DEDUP="-READ_NAME_REGEX null"
export PICARD_MARKDUP_OPTICAL_DEDUP

# Keep demultiplexing FASTQ
# Keep fastq demultiplexed or from input reads/reads2
[ "$FASTQ_DEMULTIPLEXING_KEEP" == "" ] && FASTQ_DEMULTIPLEXING_KEEP=0
export FASTQ_DEMULTIPLEXING_KEEP

# Folder for demultiplexing FASTQ
# within $SAMPLE.sequencing folder
[ "$SEQUENCING_DEMULTIPLEXING_FOLDER" == "" ] && SEQUENCING_DEMULTIPLEXING_FOLDER=demultiplexing
export SEQUENCING_DEMULTIPLEXING_FOLDER

# Additional FASTP options
# see FASTP doc
[ "$FASTP_ADDITIONAL_OPTIONS" == "" ] && FASTP_ADDITIONAL_OPTIONS=""
export FASTP_ADDITIONAL_OPTIONS


# FASTQ_PROCESSING_STEPS
# All steps to process input FASTQ files, after sequencing and demultiplexing (if any)
# Format: "step1 step2 step3"
# Example (default): fastq_reheader sort fastp fastq_clean_header compress
# Example (UMItools): fastq_reheader sort umi_tools fastp fastq_clean_header compress
# Available steps:
#    fastq_reheader: FASTQ reheader to integreate index within FASTQ comment Illumina tag (e.g. 1:N:0:xxx). Nothing done if already integrated (same header or tag BC or RX exists)
#    fastq_clean_header: FASTQ read head formatting, especially SAMTOOLS tags. Nothing done if no needs
#    compress: FASTQ files compression (see FASTQ_COMPRESSION_LEVEL)
#    sort: sort FASTQ using read name
#    fastp: process FASTP algorithm and report, UMI extraction (if any, see UMI_LOC and UMI_BARCODE_PATTERN), quality filtration...
#    umi_tools: process UMITools algorithm for UMI extraction (if any, see UMI_LOC and UMI_BARCODE_PATTERN)
# Usually:
#    "fastq_reheader sort fastp fastq_clean_header compress" for UMI technology
# dafault:
#    "sort compress" for sorting and compression

if [ -z "$FASTQ_PROCESSING_STEPS" ]; then
	FASTQ_PROCESSING_STEPS="sort compress"
fi;

# Create FASTQ_PROCESSED_STEPS variable
[ "$FASTQ_PROCESSING_STEPS" != "" ] && FASTQ_PROCESSED_STEPS=$(echo -n "." && echo $FASTQ_PROCESSING_STEPS | tr -d "." | tr " " ".") || FASTQ_PROCESSED_STEPS=""
export FASTQ_PROCESSED_STEPS



# POST SEQUENCING STEPS (default '')
# All steps and before alignment
# This sequence correspond to the FASTQ file processing before the alignemnt (trimming...)
# Format: "step1 step2 step3"
# The steps are defined as makefiles rules
# Check available steps by using the command: STARK --pipelines_infos
# Available steps (not up-to-date):
# Usually:
#    "" nothing to do

if [ -z "$POST_SEQUENCING_STEPS" ]; then
	POST_SEQUENCING_STEPS=""
fi;

# Create POST_ALIGNEMENT variable
[ "$(echo $POST_SEQUENCING_STEPS | sed 's/[[:blank:]]*$//' | sed 's/^[[:blank:]]*//')" != "" ] && POST_SEQUENCING="."$(echo $POST_SEQUENCING_STEPS | tr "," " " | tr "." " " | tr " " "\n" | sed '/^$/d' | tac  | tr "\n" " " | sed s/\.$//g | tr " " ".") || POST_SEQUENCING=""
export POST_SEQUENCING



# POST ALIGNEMENT STEPS (default "sorting realignment clipping compress")
# All steps after alignement and before calling
# This sequence correspond to the BAM file generated jsut after the alignemnt
# Format: "step1 step2 step3"
# Example: "sorting realignment clipping compress"
#    This sequence will generate the file $ALIGNER.compress.clipping.realigned.sorting.bam whose will be processed
#    Then, this BAM file will be 1/ sorted, 2/ realigned, 3/ clipped and 4/ compressed
# The steps are defined as makefiles rules
# Check available steps by using the command: STARK --pipelines_infos
# Available steps (not up-to-date):
#    sorting: BAM sorting
#    compress: BAM compression (see $BAM_COMPRESSION variable)
#    realignment: local realignment
#    recalibration: reads recalibration
#    gencore: gencore is a tool for fast and powerful deduplication for paired-end next-generation sequencing
#    markduplicates: Mark duplicated reads in BAM with PICARD MarkDuplicates. Use BARCODE_TAG to specify tag
#    clipping: BAM Clipping according to primer definition in manifest file, if any
# Usually:
#    "sorting realignment clipping compress" for Amplicon technology
#    "sorting markduplicates realignment compress" for Capture technology
#    "sorting gencore markduplicates realignment compress" for UMI technology
#POST_ALIGNMENT_STEPS="sorting realignment recalibration clipping compress"

if [ -z "$POST_ALIGNMENT_STEPS" ]; then
	POST_ALIGNMENT_STEPS="sorting realignment clipping compress"
fi;

# Create POST_ALIGNEMENT variable
[ "$(echo $POST_ALIGNMENT_STEPS | sed 's/[[:blank:]]*$//' | sed 's/^[[:blank:]]*//')" != "" ] && POST_ALIGNMENT="."$(echo $POST_ALIGNMENT_STEPS | tr "," " " | tr "." " " | tr " " "\n" | sed '/^$/d' | tac  | tr "\n" " " | sed s/\.$//g | tr " " ".") || POST_ALIGNMENT=""
export POST_ALIGNMENT



# POST CALLING STEPS (default " ")
# All steps after calling
# This sequence correspond to the VCF file generated just after the calling
# Format: "step1 step2 step3"
# Example: "recalibration filtration"
#    This sequence will generate the file $CALLER.filtration.recalibration.vcf whose will be processed
#    Then, this VCF file will be 1/ recalibrated, 2/ filtrered
# The steps are defined as makefiles rules
# Check available steps by using the command: STARK --pipelines_infos
# Available steps (not up-to-date):
#    sorting: VCF sort
#    normalization: VCF normalization
#    variantrecalibration: VCF recalibration (using GATK4). Include variantfiltration if no recalibration possible
#    variantfiltration: VCF filtration (using GATK4)
# Usually:
#    " " to avoid at this pipeline step (see POST_CALLING_MERGING_STEPS)
#    "normalization variantfiltration" for gene panel
#    "normalization variantrecalibration" for exome or genome
if [ -z "$POST_CALLING_STEPS" ]; then
	POST_CALLING_STEPS=" "
fi;

# Create POST_CALLING variable
#[ "$POST_CALLING_STEPS" != "" ] && POST_CALLING="."$(echo $POST_CALLING_STEPS | tr "," " " | tr "." " " | tr " " "\n" | sed '/^$/d' | tac  | tr "\n" " " | sed s/\.$//g | tr " " ".") || POST_CALLING=""
[ "$(echo $POST_CALLING_STEPS | sed 's/[[:blank:]]*$//' | sed 's/^[[:blank:]]*//')" != "" ] && POST_CALLING="."$(echo $POST_CALLING_STEPS | tr "," " " | tr "." " " | tr " " "\n" | sed '/^$/d' | tac  | tr "\n" " " | sed s/\.$//g | tr " " ".") || POST_CALLING=""
export POST_CALLING



# POST CALLING MERGING STEPS (default "sorting normalization variantrecalibration")
# All steps after merging calling
# This sequence correspond to the VCF file generated after the merge of VCF calling
# Format: "step1 step2 step3"
# Example: "recalibration filtration"
#    This sequence will generate the file $CALLER.filtration.recalibration.vcf whose will be processed
#    Then, this VCF file will be 1/ recalibrated, 2/ filtrered
# The steps are defined as makefiles rules
# Check available steps by using the command: STARK --pipelines_infos
# Available steps (not up-to-date):
#    sorting: VCF sort
#    normalization: VCF normalization
#    variantrecalibration: VCF recalibration (using GATK4). Include variantfiltration if no recalibration possible
#    variantfiltration: VCF filtration (using GATK4)
# Usually:
#    "sorting normalization variantrecalibration" for exome or genome
#    "sorting normalization variantfiltration" for gene panel
if [ -z "$POST_CALLING_MERGING_STEPS" ]; then
	POST_CALLING_MERGING_STEPS="sorting normalization variantrecalibration"
fi;

# Create POST_CALLING variable
[ "$(echo $POST_CALLING_MERGING_STEPS | sed 's/[[:blank:]]*$//' | sed 's/^[[:blank:]]*//')" != "" ] && POST_CALLING_MERGING="."$(echo $POST_CALLING_MERGING_STEPS | tr "," " " | tr "." " " | tr " " "\n" | sed '/^$/d' | tac  | tr "\n" " " | sed s/\.$//g | tr " " ".") || POST_CALLING_MERGING=""
export POST_CALLING_MERGING



# POST ANNOTATION STEPS (default " ")
# All steps after annotation
# This sequence correspond to the VCF file generated just after the annotation
# Format: "step1 step2 step3"
# Example: "sorting normalization"
#    This sequence will generate the file $ANNOTATION.sorting.vcf whose will be processed
#    Then, this VCF file will be 1/ sorted, 2/ normalized
# The steps are defined as makefiles rules
# Check available steps by using the command: STARK --pipelines_infos
# Available steps (not up-to-date):
#    sorting: VCF sorting
# Usually:
#    " " to avoid at this pipeline step (see POST_CALLING_MERGING_STEPS)
#    "sorting"
if [ -z "$POST_ANNOTATION_STEPS" ]; then
	POST_ANNOTATION_STEPS=" "
fi;

# Create POST_ANNOTATION variable
[ "$(echo $POST_ANNOTATION_STEPS | sed 's/[[:blank:]]*$//' | sed 's/^[[:blank:]]*//')" != "" ] && POST_ANNOTATION="."$(echo $POST_ANNOTATION_STEPS | tr "," " " | tr "." " " | tr " " "\n" | sed '/^$/d' | tac  | tr "\n" " " | sed s/\.$//g | tr " " ".") || POST_ANNOTATION=""
export POST_ANNOTATION



# BAM COMPRESSION
# Final BAM compression level (ALIGNER.bam)
if [ -z $BAM_COMPRESSION ] || ! [[ $BAM_COMPRESSION =~ ^[0-9]$ ]]; then
	BAM_COMPRESSION=9
fi;
export BAM_COMPRESSION


# BAM VALIDATION COMPRESSION
# Validation BAM compression level (ALIGNER.validation.bam)
if [ -z $BAM_VALIDATION_COMPRESSION ] || ! [[ $BAM_VALIDATION_COMPRESSION =~ ^[0-9]$ ]]; then
	BAM_VALIDATION_COMPRESSION=5
fi;
export BAM_VALIDATION_COMPRESSION



### GENCORE

# supporting_reads ; set to 2 to keep the clusters with 2 or more supporting reads / ultrasensible filter
# set to 1 to replace Picard Markduplicates dedup
# default "" (corresponding to "--supporting_reads 1", see gencore doc)
export GENCORE_SUP_READS

# --score_threshold
# set to 8 recommanded for dup-rate < 50% if you want to keep all the DNA fragments, and for each output read you want to discard all the low quality unoverlapped mutations to obtain a relative clean data
# default "" (corresponfing to "--score_threshold 6", see gencore doc)
export GENCORE_SCORE_THREESHOLD

# --ratio_threshold
# if the ratio of the major base in a cluster is less than <ratio_threshold>, it will be further compared to the reference.
# The valud should be 0.5~1.0, and the default value is 0.8
# default "" (corresponding to "--ratio_threshold 0.8", see gencore doc)
export GENCORE_RATIO_THREESHOLD

# --umi_diff_threshold
# if two reads with identical mapping position have UMI difference <= <umi_diff_threshold>, then they will be merged to generate a consensus read
# Default "" (corresponding to "--umi_diff_threshold 2", see gencore doc)
export GENCORE_DIFF_THREESHOLD

# Quality Threeshold
# --high_qual (Q30) ; --moderate_qual (Q20) ; --low_qual (Q15)
# --high_qual : the threshold for a quality score to be considered as high quality. Default 30 means Q30. (int [=30])
# --moderate_qual : the threshold for a quality score to be considered as moderate quality. Default 20 means Q20. (int [=20])
# --low_qual : the threshold for a quality score to be considered as low quality. Default 15 means Q15. (int [=15])
# e.g. "--moderate_qual 20", "--high_qual 20 --moderate_qual 15 --low_qual 10"
# default "" (corresponding to default quality, see gencore doc)
export GENCORE_QUAL_THREESHOLD

# --coverage_sampling
# the sampling rate for genome scale coverage statistics. Default 10000 means 1/10000
# for statistics purpose, can be reduce to gain performance
# not included in the mk file
# default "" (corresponding to "--coverage_sampling=10000", see gencore doc)
export GENCORE_COVERAGE_SAMPLING


### CRAM

# CRAM OPTIONS
# Final CRAM options for compression (archive.cram)
# example: CRAM_OPTIONS="version=3.0,level=9,no_ref,use_lzma,seqs_per_slice=100000"
export CRAM_OPTIONS


# CRAM REMOVE TAGS
# Final CRAM options for tags (archive.cram)
# example: CRAM_REMOVE_TAGS="BD,BI,OQ"
export CRAM_REMOVE_TAGS


# ANNOTATION
#############

# DEFAULT

# Default annotation with HOWARD for HOWARD for intermediate VCF (for each caller) used by default with annotation rule "howard"
if [ -z "$HOWARD_ANNOTATION" ]; then
	#HOWARD_ANNOTATION="core,symbol,location,outcome,hgvs,snpeff,snpeff_hgvs"
	#HOWARD_ANNOTATION="core,snpeff_hgvs"
	HOWARD_ANNOTATION="symbol,location,outcome,hgvs"
	# ANNOTATINO FULL
	#HOWARD_ANNOTATION="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs"
fi;
export HOWARD_ANNOTATION


# DEJAVU
# if [ -s $ANNOVAR_DATABASES/$ASSEMBLY"_dejavu."$APP_GROUP.$APP_PROJECT.txt ]; then
# 	HOWARD_ANNOTATION="$HOWARD_ANNOTATION,dejavu.$APP_GROUP.$APP_PROJECT"
# fi
# # DEJAVU for all the GROUP
# for DEJAVU_DATABASE_ONE in $(find $ANNOVAR_DATABASES -name $ASSEMBLY"_dejavu."$APP_GROUP".*txt" 2>/dev/null); do
# 	HOWARD_ANNOTATION=$HOWARD_ANNOTATION","$(basename $DEJAVU_DATABASE_ONE | sed s/^$ASSEMBLY"_"//g | sed s/.txt$//g );
# done;
# export HOWARD_ANNOTATION

# DEJAVU
if [ -s $DEJAVU_ANNOVAR_DATABASES/$ASSEMBLY"_dejavu."$APP_GROUP.$APP_PROJECT.txt ]; then
	HOWARD_DEJAVU_ANNOTATION="dejavu.$APP_GROUP.$APP_PROJECT"
fi
# DEJAVU for all the GROUP
for DEJAVU_DATABASE_ONE in $(find $DEJAVU_ANNOVAR_DATABASES -name $ASSEMBLY"_dejavu."$APP_GROUP".*txt" 2>/dev/null); do
	HOWARD_DEJAVU_ANNOTATION=$HOWARD_DEJAVU_ANNOTATION","$(basename $DEJAVU_DATABASE_ONE | sed s/^$ASSEMBLY"_"//g | sed s/.txt$//g );
done;
HOWARD_DEJAVU_ANNOTATION=$(echo $HOWARD_DEJAVU_ANNOTATION | tr "," "\n" | sort -u | tr "\n" "," | sed s/,$//)
export HOWARD_DEJAVU_ANNOTATION


# MINIMAL

# Minimal annotation with HOWARD for minimal VCF annotation  (rule howard_minimal)
if [ -z "$HOWARD_ANNOTATION_MINIMAL" ]; then
	HOWARD_ANNOTATION_MINIMAL=$HOWARD_ANNOTATION
fi;
export HOWARD_ANNOTATION_MINIMAL

# REPORT
# Default annotation with HOWARD for Report (/fullfinal VCF)
if [ -z "$HOWARD_ANNOTATION_REPORT" ]; then
	HOWARD_ANNOTATION_REPORT=$HOWARD_ANNOTATION
fi;
export HOWARD_ANNOTATION_REPORT

# ANALYSIS
# Default annotation with HOWARD for whole analysis (annotation NOT forced)
if [ -z "$HOWARD_ANNOTATION_ANALYSIS" ]; then
	HOWARD_ANNOTATION_ANALYSIS=$HOWARD_ANNOTATION
fi;
export HOWARD_ANNOTATION_ANALYSIS



# CALCULATION
###############

# DEFAULT
# Default calculation with HOWARD
if [ -z "$HOWARD_CALCULATION" ]; then
	HOWARD_CALCULATION="VARTYPE,NOMEN"
fi;
export HOWARD_CALCULATION

# MINIMAL
# Minimal calculation with HOWARD
if [ -z "$HOWARD_CALCULATION_MINIMAL" ]; then
	HOWARD_CALCULATION_MINIMAL=$HOWARD_CALCULATION
fi;
export HOWARD_CALCULATION_MINIMAL

# REPORT
# Default calculation with HOWARD for Report (full/final VCF)
if [ -z "$HOWARD_CALCULATION_REPORT" ]; then
	HOWARD_CALCULATION_REPORT=$HOWARD_CALCULATION
fi;
export HOWARD_CALCULATION_REPORT

# ANALYSIS
# Default calculation with HOWARD for whole analysis (prioritization forced)
if [ -z "$HOWARD_CALCULATION_ANALYSIS" ]; then
	HOWARD_CALCULATION_ANALYSIS=$HOWARD_CALCULATION
fi;
export HOWARD_CALCULATION_ANALYSIS


# List of annotation fields to extract NOMEN annotation (default 'hgvs', see HOWARD docs)
if [ -z "$HOWARD_NOMEN_FIELDS" ]; then
	HOWARD_NOMEN_FIELDS="hgvs"
fi;
export HOWARD_NOMEN_FIELDS



# PRIORITIZATION
###################
# This option create ranking scores in VCF and comment in TXT (after translation).
# Scores can be used to sort variant in the TXT

# DEFAULT
if [ -z "$HOWARD_PRIORITIZATION_DEFAULT" ]; then
	HOWARD_PRIORITIZATION_DEFAULT="null"
fi;
export HOWARD_PRIORITIZATION_DEFAULT


# # DEFAULT (deprecated)
# # Default filter to prioritize/rank variant.
# if [ -z "$HOWARD_PRIORITIZATION" ]; then
# 	HOWARD_PRIORITIZATION=$HOWARD_PRIORITIZATION_DEFAULT
# fi;
# if [ ! -z "$APP_NAME" ] && [ ${APP_NAME^^} != "DEFAULT" ]; then
# 	HOWARD_PRIORITIZATION="$APP_NAME,$HOWARD_PRIORITIZATION"
# fi;
# # Keep fist as first and sort the rest
# HOWARD_PRIORITIZATION=$(echo $(echo $HOWARD_PRIORITIZATION | tr "," " " | cut -d" " -f1 | tr " " "," && echo $HOWARD_PRIORITIZATION | tr "," " " | tr " " "\n" | sort -u | grep -v "^$(echo $HOWARD_PRIORITIZATION | tr "," " " | cut -d" " -f1)$" | tr "\n" "," | sed s/,$//) | tr " " "," )
# export HOWARD_PRIORITIZATION


# # MINIMAL (deprecated)
# # Minimal filter to prioritize/rank variant.
# if [ -z "$HOWARD_PRIORITIZATION_MINIMAL" ]; then
# 	HOWARD_PRIORITIZATION_MINIMAL=$HOWARD_PRIORITIZATION
# fi;
# if [ ! -z "$APP_NAME" ] && [ ${APP_NAME^^} != "DEFAULT" ]; then
# 	HOWARD_PRIORITIZATION_MINIMAL="$APP_NAME,$HOWARD_PRIORITIZATION_MINIMAL"
# fi;
# # Keep fist as first and sort the rest
# HOWARD_PRIORITIZATION_MINIMAL=$(echo $(echo $HOWARD_PRIORITIZATION_MINIMAL | tr "," " " | cut -d" " -f1 | tr " " "," && echo $HOWARD_PRIORITIZATION_MINIMAL | tr "," " " | tr " " "\n" | sort -u | grep -v "^$(echo $HOWARD_PRIORITIZATION_MINIMAL | tr "," " " | cut -d" " -f1)$" | tr "\n" "," | sed s/,$//) | tr " " "," )
# export HOWARD_PRIORITIZATION_MINIMAL


# Report
# Default filter to prioritize/rank variant for Report
if [ -z "$HOWARD_PRIORITIZATION_REPORT" ]; then
	HOWARD_PRIORITIZATION_REPORT=$HOWARD_PRIORITIZATION_DEFAULT
fi;
if [ ! -z "$APP_NAME" ] && [ ${APP_NAME^^} != "DEFAULT" ]; then
	HOWARD_PRIORITIZATION_REPORT="$APP_NAME,$HOWARD_PRIORITIZATION_REPORT"
fi;
# Keep fist as first and sort the rest
HOWARD_PRIORITIZATION_REPORT=$(echo $(echo $HOWARD_PRIORITIZATION_REPORT | tr "," " " | cut -d" " -f1 | tr " " "," && echo $HOWARD_PRIORITIZATION_REPORT | tr "," " " | tr " " "\n" | sort -u | grep -v "^$(echo $HOWARD_PRIORITIZATION_REPORT | tr "," " " | cut -d" " -f1)$" | tr "\n" "," | sed s/,$//) | tr " " "," )
export HOWARD_PRIORITIZATION_REPORT


# ANALYSIS
# Default filter to prioritize/rank variant for whole analysis (calculation forced)
if [ -z "$HOWARD_PRIORITIZATION_ANALYSIS" ]; then
	HOWARD_PRIORITIZATION_ANALYSIS=$HOWARD_PRIORITIZATION_DEFAULT
fi;
if [ ! -z "$APP_NAME" ] && [ ${APP_NAME^^} != "DEFAULT" ]; then
	HOWARD_PRIORITIZATION_ANALYSIS="$APP_NAME,$HOWARD_PRIORITIZATION_ANALYSIS"
fi;
# Keep fist as first and sort the rest
HOWARD_PRIORITIZATION_ANALYSIS=$(echo $(echo $HOWARD_PRIORITIZATION_ANALYSIS | tr "," " " | cut -d" " -f1 | tr " " "," && echo $HOWARD_PRIORITIZATION_ANALYSIS | tr "," " " | tr " " "\n" | sort -u | grep -v "^$(echo $HOWARD_PRIORITIZATION_ANALYSIS | tr "," " " | cut -d" " -f1)$" | tr "\n" "," | sed s/,$//) | tr " " "," )
export HOWARD_PRIORITIZATION_ANALYSIS


# VARANK
# Default prioritization with HOWARD for VaRank score mode
if [ -z "$HOWARD_PRIORITIZATION_VARANK" ]; then
	HOWARD_PRIORITIZATION_VARANK=$HOWARD_PRIORITIZATION_DEFAULT
fi;
if [ ! -z "$APP_NAME" ] && [ ${APP_NAME^^} != "DEFAULT" ]; then
	HOWARD_PRIORITIZATION_VARANK="$APP_NAME,$HOWARD_PRIORITIZATION_VARANK"
fi;
# Keep fist as first and sort the rest
HOWARD_PRIORITIZATION_VARANK=$(echo $(echo $HOWARD_PRIORITIZATION_VARANK | tr "," " " | cut -d" " -f1 | tr " " "," && echo $HOWARD_PRIORITIZATION_VARANK | tr "," " " | tr " " "\n" | sort -u | grep -v "^$(echo $HOWARD_PRIORITIZATION_VARANK | tr "," " " | cut -d" " -f1)$" | tr "\n" "," | sed s/,$//) | tr " " "," )
export HOWARD_PRIORITIZATION_VARANK


# TRANSLATION
################
# List of fields to show in the TXT file
# use ALL to show ALL "other" annotations

# DEFAULT
# Fields to show after translation
if [ -z $HOWARD_FIELDS ]; then
	HOWARD_FIELDS="NOMEN,PZFlag,PZScore,PZComment,CNOMEN,PNOMEN,location,outcome,snpeff_impact,VAF_average,dbSNP,dbSNPNonFlagged,popfreq,ALL"
fi;
export HOWARD_FIELDS
# Sort variant in the TXT using 2 fields
if [ -z $HOWARD_SORT_BY ]; then
	HOWARD_SORT_BY="PZFlag,PZScore"
fi;
export HOWARD_SORT_BY
# Order fields in variant ranking
if [ -z $HOWARD_ORDER_BY ]; then
	HOWARD_ORDER_BY="DESC,DESC"
fi;
export HOWARD_ORDER_BY

# MINIMAL
# Fields to show after minimal translation
if [ -z $HOWARD_FIELDS_MINIMAL ]; then
	#HOWARD_FIELDS_MINIMAL=$HOWARD_FIELDS
	HOWARD_FIELDS_MINIMAL=$(echo $HOWARD_FIELDS | sed "s/ALL//")
fi;
export HOWARD_FIELDS_MINIMAL
# Sort variant in the TXT using 2 fields
if [ -z $HOWARD_SORT_BY_MINIMAL ]; then
	HOWARD_SORT_BY_MINIMAL=$HOWARD_SORT_BY
fi;
export HOWARD_SORT_BY_MINIMAL
# Order fields in variant ranking
if [ -z $HOWARD_ORDER_BY_MINIMAL ]; then
	HOWARD_ORDER_BY_MINIMAL=$HOWARD_ORDER_BY
fi;
export HOWARD_ORDER_BY_MINIMAL

# REPORT
# Fields to show after minimal translation
if [ -z $HOWARD_FIELDS_REPORT ]; then
	#HOWARD_FIELDS_REPORT=$HOWARD_FIELDS
	HOWARD_FIELDS_REPORT=$(echo $HOWARD_FIELDS | sed "s/ALL//")
fi;
export HOWARD_FIELDS_REPORT
# Sort variant in the TXT using 2 fields
if [ -z $HOWARD_SORT_BY_REPORT ]; then
	HOWARD_SORT_BY_REPORT=$HOWARD_SORT_BY
fi;
export HOWARD_SORT_BY_REPORT
# Order fields in variant ranking
if [ -z $HOWARD_ORDER_BY_REPORT ]; then
	HOWARD_ORDER_BY_REPORT=$HOWARD_ORDER_BY
fi;
export HOWARD_ORDER_BY_REPORT



# INFO to FORMAT
# Transfers INFO annotation to FORMAT annotation
# Useful for annotations on full VCF to final VCF on each sample
if [ -z $INFO_TO_FORMAT_ANNOTATIONS ]; then
	INFO_TO_FORMAT_ANNOTATIONS=""
fi;
export INFO_TO_FORMAT_ANNOTATIONS


# Recalibration and Filtration
################################

# Variant Filtration
# Filter variant calls based on INFO and/or FORMAT annotations

# Variant Filtration main option (see documentation guide for more info)
# default: VARIANTFILTRATION_OPTIONS=
#VARIANTFILTRATION_OPTIONS=
export VARIANTFILTRATION_OPTIONS

# One or more expression used with INFO fields to filter SNP (see documentation guide for more info)
# default: VARIANTFILTRATION_SNP_FILTER_OPTION=''
# example: VARIANTFILTRATION_SNP_FILTER_OPTION='--filter-name "SNP_filter_QD" --filter-expression "QD < 2.0" --filter-name "SNP_filter_FS" --filter-expression "FS > 60.0" --filter-name "SNP_filter_MQ" --filter-expression "MQ < 40.0" --filter-name "SNP_filter_MQRankSum" --filter-expression "MQRankSum < -12.5" --filter-name "SNP_filter_ReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0"'
# example: VARIANTFILTRATION_SNP_FILTER_OPTION='--filter-name "HARD_TO_VALIDATE" --filter-expression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)"    --filter-name "VeryVeryLowDepth" --filter-expression "DP == 0"    --filter-name "VeryLowDepth" --filter-expression "DP > 0 && DP < 10"    --filter-name "LowDepth" --filter-expression "DP >= 10 && DP < 30"    --filter-name "VeryVeryLowQual" --filter-expression "QUAL == 0"    --filter-name "VeryLowQual" --filter-expression "QUAL > 0 && QUAL < 30.0"    --filter-name "LowQual" --filter-expression "QUAL >= 30.0 && QUAL < 50.0"    --filter-name "LowQD" --filter-expression "QD >= 0.0 && QD < 1.5"'
if [ -z "$VARIANTFILTRATION_SNP_FILTER_OPTION" ]; then
	VARIANTFILTRATION_SNP_FILTER_OPTION=''
fi;
export VARIANTFILTRATION_SNP_FILTER_OPTION

# One or more expression used with FORMAT (sample/genotype-level) fields to filter SNP (see documentation guide for more info)
# default: VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION=''
# example: VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION='--genotype-filter-expression "GQ == 0" --genotype-filter-name "genotype_GQ_filter_VeryVeryLow" --genotype-filter-expression "GQ > 0 && GQ < 50.0" --genotype-filter-name "genotype_GQ_filter_VeryLow"  --genotype-filter-expression "GQ >= 50.0 && GQ < 90.0" --genotype-filter-name "genotype_GQ_filter_Low" --genotype-filter-expression "DP == 0" --genotype-filter-name "genotype_DP_filter_VeryVeryLow" --genotype-filter-expression "DP >= 0 && DP < 10" --genotype-filter-name "genotype_DP_filter_VeryLow"  --genotype-filter-expression "DP >= 10 && DP < 30" --genotype-filter-name "genotype_GQ_filter_Low"'
# example: VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION='--genotype-filter-name "VeryVeryLowGQ" --genotype-filter-expression "GQ == 0"    --genotype-filter-name "VeryLowGQ" --genotype-filter-expression "GQ > 0 && GQ < 50.0"    --genotype-filter-name "LowGQ" --genotype-filter-expression "GQ >= 50.0 && GQ < 90.0"    --genotype-filter-name "VeryVeryLowDP" --genotype-filter-expression "DP == 0"    --genotype-filter-name "VeryLowDP" --genotype-filter-expression "DP >= 0 && DP < 10"    --genotype-filter-name "LowDP" --genotype-filter-expression "DP >= 10 && DP < 30" '
if [ -z "$VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION" ]; then
	VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION=''
fi;
export VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION

# One or more expression used with INFO fields to filter INDEL (see documentation guide for more info)
# default: VARIANTFILTRATION_INDEL_FILTER_OPTION=''
# example: VARIANTFILTRATION_INDEL_FILTER_OPTION='--filter-name "INDEL_filter_QD" --filter-expression "QD < 2.0" --filter-name "INDEL_filter_FS" --filter-expression "FS > 200.0" --filter-name "INDEL_filter_ReadPosRankSum" --filter-expression "ReadPosRankSum < -20.0"'
# example: VARIANTFILTRATION_INDEL_FILTER_OPTION='--filter-name "HARD_TO_VALIDATE" --filter-expression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)"    --filter-name "VeryVeryLowDepth" --filter-expression "DP == 0"    --filter-name "VeryLowDepth" --filter-expression "DP > 0 && DP < 10"    --filter-name "LowDepth" --filter-expression "DP >= 10 && DP < 30"    --filter-name "VeryVeryLowQual" --filter-expression "QUAL == 0"    --filter-name "VeryLowQual" --filter-expression "QUAL > 0 && QUAL < 30.0"    --filter-name "LowQual" --filter-expression "QUAL >= 30.0 && QUAL < 50.0"    --filter-name "LowQD" --filter-expression "QD >= 0.0 && QD < 1.5"'
if [ -z "$VARIANTFILTRATION_INDEL_FILTER_OPTION" ]; then
	VARIANTFILTRATION_INDEL_FILTER_OPTION=''
fi;
export VARIANTFILTRATION_INDEL_FILTER_OPTION

# One or more expression used with FORMAT (sample/genotype-level) fields to filter INDEL (see documentation guide for more info)
# default: VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION=''
# example: VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION='--genotype-filter-expression "GQ == 0" --genotype-filter-name "genotype_GQ_filter_VeryVeryLow" --genotype-filter-expression "GQ > 0 && GQ < 50.0" --genotype-filter-name "genotype_GQ_filter_VeryLow"  --genotype-filter-expression "GQ >= 50.0 && GQ < 90.0" --genotype-filter-name "genotype_GQ_filter_Low" --genotype-filter-expression "DP == 0" --genotype-filter-name "genotype_DP_filter_VeryVeryLow" --genotype-filter-expression "DP >= 0 && DP < 10" --genotype-filter-name "genotype_DP_filter_VeryLow"  --genotype-filter-expression "DP >= 10 && DP < 30" --genotype-filter-name "genotype_GQ_filter_Low"'
# example: VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION='--genotype-filter-name "VeryVeryLowGQ" --genotype-filter-expression "GQ == 0"    --genotype-filter-name "VeryLowGQ" --genotype-filter-expression "GQ > 0 && GQ < 50.0"    --genotype-filter-name "LowGQ" --genotype-filter-expression "GQ >= 50.0 && GQ < 90.0"    --genotype-filter-name "VeryVeryLowDP" --genotype-filter-expression "DP == 0"    --genotype-filter-name "VeryLowDP" --genotype-filter-expression "DP >= 0 && DP < 10"    --genotype-filter-name "LowDP" --genotype-filter-expression "DP >= 10 && DP < 30" '
if [ -z "$VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION" ]; then
	VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION=''
fi;
export VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION

# Remove previous filters applied to the VCF
# within makefile rule: VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS_OPTION?=$(shell if (( $(VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS) )); then echo " --invalidate-previous-filters "; fi )
# default: VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS=0
if [ -z "$VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS" ]; then
	VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS=0
fi;
export VARIANTFILTRATION_INVALIDATE_PREVIOUS_FILTERS


# Variant Recalibrator

# Variant Recalibrator main option (see documentation guide for more info)
# default: VARIANTRECALIBRATOR_OPTIONS=
#VARIANTRECALIBRATOR_OPTIONS=
export VARIANTRECALIBRATOR_OPTIONS

# Variant Recalibrator SNP resources option (see documentation guide for more info)
# These resources need to be available on STARK Databases folder for GATK
# default:
if [ "$ASSEMBLY" == 'hg19' ]; then
VARIANTRECALIBRATION_SNP_RESOURCES="
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.vcf.gz
    -resource:omni,known=false,training=true,truth=true,prior=12.0 1000G_omni2.5.b37.vcf.gz
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.b37.vcf.gz 
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.b37.vcf.gz
"
fi;
if [ "$ASSEMBLY" == 'hg38' ]; then
VARIANTRECALIBRATION_SNP_RESOURCES="
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz
    -resource:omni,known=false,training=true,truth=true,prior=12.0 1000G_omni2.5.hg38.vcf.gz
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf.gz
"
fi;
export VARIANTRECALIBRATION_SNP_RESOURCES

# Variant Recalibrator INDEL resources option (see documentation guide for more info)
# These resources need to be available on STARK Databases folder for GATK
# default:
if [ "$ASSEMBLY" == 'hg19' ]; then
VARIANTRECALIBRATION_INDEL_RESOURCES="
    -resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.b37.vcf.gz
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.b37.vcf.gz
"
fi;
if [ "$ASSEMBLY" == 'hg38' ]; then
VARIANTRECALIBRATION_INDEL_RESOURCES="
    -resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf.gz
"
fi;
export VARIANTRECALIBRATION_INDEL_RESOURCES

# Variant Recalibrator SNP annotations option (see documentation guide for more info)
# default: VARIANTRECALIBRATION_SNP_ANNOTATIONS="-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP"
VARIANTRECALIBRATION_SNP_ANNOTATIONS="-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP"

# Variant Recalibrator INDEL annotations option (see documentation guide for more info)
# default: VARIANTRECALIBRATION_INDEL_ANNOTATIONS="-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum"
if [ -z "$VARIANTRECALIBRATION_INDEL_ANNOTATIONS" ]; then
	VARIANTRECALIBRATION_INDEL_ANNOTATIONS="-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum"
fi;
export VARIANTRECALIBRATION_INDEL_ANNOTATIONS

# Variant Recalibrator SNP tranches option (see documentation guide for more info)
# default: VARIANTRECALIBRATION_SNP_TRANCHES="-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.1 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0"
if [ -z "$VARIANTRECALIBRATION_SNP_TRANCHES" ]; then
	VARIANTRECALIBRATION_SNP_TRANCHES="-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.1 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0"
fi;
export VARIANTRECALIBRATION_SNP_TRANCHES

# Variant Recalibrator INDEL tranches option (see documentation guide for more info)
# default: VARIANTRECALIBRATION_INDEL_TRANCHES="-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.1 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0"
if [ -z "$VARIANTRECALIBRATION_INDEL_TRANCHES" ]; then
	VARIANTRECALIBRATION_INDEL_TRANCHES="-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.1 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0"
fi;
export VARIANTRECALIBRATION_INDEL_TRANCHES


# Variant Recalibrator Optional Variant Filtration 
# If Variant Recalibrator failed, usually due to lack of variant in the input callset, a Variant Filtration is optional
# default: no filtration
# example: Use empty value to switch off this option (no filtration)
# VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_OPTION=
# VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION=
# VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_OPTION=
# VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION=
# example: Use same options than Variant Filtration in stand alone (see above)
# VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_OPTION=$VARIANTFILTRATION_SNP_FILTER_OPTION
# VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION=$VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION
# VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_OPTION=$VARIANTFILTRATION_INDEL_FILTER_OPTION
# VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION=$VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION
# Filter variant calls based on INFO and/or FORMAT annotations for SNP and INDEL 
if [ -z "$VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_OPTION" ]; then
	VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_OPTION=""
fi;
export VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_OPTION
if [ -z "$VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION" ]; then
	VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION=""
fi;
export VARIANTRECALIBRATOR_VARIANTFILTRATION_SNP_FILTER_EXPRESSION_OPTION
if [ -z "$VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_OPTION" ]; then
	VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_OPTION=""
fi;
export VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_OPTION
if [ -z "$VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION" ]; then
	VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION=""
fi;
export VARIANTRECALIBRATOR_VARIANTFILTRATION_INDEL_FILTER_EXPRESSION_OPTION



# DAEMON CONFIG
##########

# Config
export TMP_REMOVE_TIME=10
export LOG_FILE=analysis.log
export RUN_ANALYSIS_CONDITIONS="RTAComplete.txt SampleSheet.csv"


# DATABASE PEPPER CONFIG (depreciated)
########################################

CONFIG=$PEPPER_CONFIG

# RELEASE
###########

README_FILE=$STARK_FOLDER_CONFIG/README
RELEASE_FILE=$STARK_FOLDER_CONFIG/RELEASE

# INFRASTRUCTURE
##################

# THREADS is CORES - 1

if [ -z $THREADS ] || [ "${THREADS^^}" == "AUTO" ] || ! [[ $THREADS =~ ^[0-9]+$ ]]; then
	#export CORES=$(ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w)	# NB of cores in the server
	export CORES=$(nproc)	# NB of cores in the server
	export CORES_FREE=1							# Number of threads free for other command
	export CORES_TO_USE=$(($CORES-$CORES_FREE))				# Nb of cores to use
	THREADS=$CORES_TO_USE
	#THREADS=2
	#if ! [[ $THREADS =~ "^\d+" ]]; then # || [ "$THREADS" == "" ]; then
	#if [[ -z $THREADS ]]; then # || [ "$THREADS" == "" ]; then
	#	THREADS=1
	#fi;
	#echo "|$THREADS|"
	if [ $THREADS -gt $CORES ]; then
		THREADS=$CORES
	fi;

fi;
export THREADS
export CORES_TO_USE=$THREADS

#export THREADS_BY_SAMPLE=1				# Number of threads for a sample's command
#export THREADS_BY_PIPELINE=1				# Number of threads for a pipeline's command
#export THREADS_BY_ALIGNER=1				# Number of threads for a aligner's command
#export THREADS_BY_CALLER=1				# Number of threads for a caller's command
#export THREADS=$(($CORES_TO_USE/$THREADS_BY_SAMPLE))	# Total nb of threads to use for Makefile
#export THREADS_BWA=$THREADS_BY_SAMPLE			# NB of threads for BWA command
#export THREADS_RTC=$THREADS_BY_SAMPLE			# NB of threads for GATK realigment RealignerTargetCreator command

export THREADS_AUTO=1					# Automatize nb of threads depending on the input and number of cores to use (TRUE if = 1)


# THREADS_LOADING (default THREADS)
# Number of threads used for loading demultiplexed data
# AUTO will considere THREADS threads to use
# The number of threads need to be between 1 and the total number of cores available (autoadjusting if bad value)
if ! [[ $THREADS_LOADING =~ $re ]] || [ -z "$THREADS_LOADING" ] || [ "$THREADS_LOADING" == "" ] || [ $THREADS_LOADING -gt $CORES ] ; then
	THREADS_LOADING=$THREADS
fi;
export THREADS_LOADING

# THREADS_WRITING (default THREADS)
# Number of threads used for writing demultiplexed data
# AUTO will considere THREADS threads to use
# The number of threads need to be between 1 and the total number of cores available (autoadjusting if bad value)
if ! [[ $THREADS_WRITING =~ $re ]] || [ -z "$THREADS_WRITING" ] || [ "$THREADS_WRITING" == "" ] || [ $THREADS_WRITING -gt $CORES ] ; then
	THREADS_WRITING=$THREADS
fi;
export THREADS_WRITING


# THREADS COPY
# Threads to use for results copy in repositories (repository, archives, favorites)
if ! [[ $THREADS_COPY =~ $re ]] || [ -z "$THREADS_COPY" ] || [ "$THREADS_COPY" == "" ] || [ $THREADS_COPY -gt $CORES ] ; then
	THREADS_COPY=1
fi;
export THREADS_COPY


# MEMORY
MEMTOTAL=$(cat /proc/meminfo 2>/dev/null | grep MemTotal | awk '{print $2}')	# MEMORY in octet
#export MEMORY=$(($MEMTOTAL/$THREADS/1024/1024))			# MEMORY in Go

if [ "$MEMORY" == "" ] || [ $MEMORY -lt 1 ]; then MEMORY=$(($MEMTOTAL/$CORES_TO_USE/1024/1024)); fi;
if [ "$MEMORY" == "" ] || [ $MEMORY -lt 1 ]; then MEMORY=1; fi;
export MEMORY				# MEMORY in Go

# for naming...
if [ "$JAVA_MEMORY" == "" ] || [ $JAVA_MEMORY -lt 1 ]; then JAVA_MEMORY=$MEMORY; fi;
export JAVA_MEMORY		


# JAVA_MEMORY by sample and by caller
###############

if ((1)); then

	[ "$NB_SAMPLE" == "" ] && NB_SAMPLE=1
	[ "$NB_PIPELINES" == "" ] && NB_PIPELINES=1
	[ "$NB_ALIGNERS" == "" ] && NB_ALIGNERS=1
	[ "$NB_CALLERS" == "" ] && NB_CALLERS=1
	# NB_SAMPLE=1
	# NB_PIPELINES=1
	# NB_ALIGNERS=1
	# NB_CALLERS=1

	# Default
	JAVA_MEMORY_BY_SAMPLE=$JAVA_MEMORY
	JAVA_MEMORY_BY_CALLER=$JAVA_MEMORY

	# Calcul
	JAVA_MEMORY_BY_SAMPLE=$(($MEMTOTAL/$NB_SAMPLE/1024/1024))				# Number of threads by sample
	JAVA_MEMORY_BY_PIPELINE=$(($JAVA_MEMORY_BY_SAMPLE/$NB_PIPELINES))	# Number of threads for a pipeline's command
	JAVA_MEMORY_BY_ALIGNER=$(($JAVA_MEMORY_BY_SAMPLE/$NB_ALIGNERS))	# Number of threads for a aligner's command
	JAVA_MEMORY_BY_CALLER=$(($JAVA_MEMORY_BY_SAMPLE/$NB_CALLERS))		# Number of threads for a caller's command

	# Test
	if [ "$JAVA_MEMORY_BY_SAMPLE" == "" ] || [ $JAVA_MEMORY_BY_SAMPLE -lt 1 ]; then JAVA_MEMORY_BY_SAMPLE=1; fi;
	if [ "$JAVA_MEMORY_BY_PIPELINE" == "" ] || [ $JAVA_MEMORY_BY_PIPELINE -lt 1 ]; then JAVA_MEMORY_BY_PIPELINE=1; fi;
	if [ "$JAVA_MEMORY_BY_ALIGNER" == "" ] || [ $JAVA_MEMORY_BY_ALIGNER -lt 1 ]; then JAVA_MEMORY_BY_ALIGNER=1; fi;
	if [ "$JAVA_MEMORY_BY_CALLER" == "" ] || [ $JAVA_MEMORY_BY_CALLER -lt 1 ]; then JAVA_MEMORY_BY_CALLER=1; fi;

fi;


# JAVA FLAGS
##############
export JAVA_MEMORY_BY_SAMPLE
export JAVA_MEMORY_BY_PIPELINE
export JAVA_MEMORY_BY_ALIGNER
export JAVA_MEMORY_BY_CALLER
export JAVA_FLAGS_TMP_FOLDER=" -Dorg.xerial.snappy.tempdir=$TMP_FOLDER_TMP -Djava.io.tmpdir=$TMP_FOLDER_TMP";
export JAVA_FLAGS_OTHER_PARAM=" -Dsnappy.disable=true -Dsamjdk.try_use_intel_deflater=false ";
export JAVA_FLAGS_DEFAULT=" -Xmx"$JAVA_MEMORY"g $JAVA_FLAGS_OTHER_PARAM $JAVA_FLAGS_TMP_FOLDER";

# JAVA_FLAGS
if [ "$JAVA_FLAGS" == "" ]; then
	JAVA_FLAGS=$JAVA_FLAGS_DEFAULT;
fi;
#JAVA_FLAGS=$JAVA_FLAGS_DEFAULT;
export JAVA_FLAGS;
#if [ "$JAVA_FLAGS_BY_SAMPLE" == "" ]; then
#	JAVA_FLAGS_BY_SAMPLE=" -Xmx"$JAVA_MEMORY_BY_SAMPLE"g $JAVA_FLAGS_OTHER_PARAM $JAVA_FLAGS_TMP_FOLDER";
#fi;
#export JAVA_FLAGS_BY_SAMPLE;


# RESOURCES MANAGMENT
#######################

### Picard CollectHsMetrics
# If validation bam is small enough (lower than MAX_VALIDATION_BAM_SIZE Kb), launch CollectHsMetrics the classic way
# Otherwise increase RAM to MAX_CONCURRENT_HSMETRICS_RAM and limit command to MAX_CONCURRENT_HSMETRICS concurrent launches
# default: MAX_VALIDATION_BAM_SIZE=1000000000 (exageratly big to switch off)
# example:
# - MAX_VALIDATION_BAM_SIZE=2097152		# 2Go
# - MAX_CONCURRENT_HSMETRICS=1 			# No concurrence, only 1 HsMetrics at once
# - MAX_CONCURRENT_HSMETRICS_RAM=16g 	# 16Go of RAM to Collect HsMetrics
if [ "$MAX_VALIDATION_BAM_SIZE" == "" ]; then
	MAX_VALIDATION_BAM_SIZE=1000000000;
fi;
export MAX_VALIDATION_BAM_SIZE
if [ "$MAX_CONCURRENT_HSMETRICS" == "" ]; then
	MAX_CONCURRENT_HSMETRICS=1;
fi;
export MAX_CONCURRENT_HSMETRICS
if [ "$MAX_CONCURRENT_HSMETRICS_RAM" == "" ]; then
	MAX_CONCURRENT_HSMETRICS_RAM=16g;
fi;
export MAX_CONCURRENT_HSMETRICS_RAM




# REPORT
# Report variables

# Report Sections
# List of sections to show in the report (default "ALL")
# Sections :
#   Report Sections: results_summary sequencing_mapping depth coverage variant_calling variant_stats
#   Report Annex Sections: annex_coverage annex_depth annex_genes_coverage annex_variants annex_annotations
if [ -z "$REPORT_SECTIONS" ]; then
	REPORT_SECTIONS="ALL"
fi;
export REPORT_SECTIONS



### REPORT for run Files
# Generate variants files from run with full VCF (include all calling information)
export REPORT_VARIANTS_FULL



# RULES export has to be after PIPELINE is defined and no longer changed, otherwise some rules might not be loaded

# Main rules
RULES=$(ls $STARK_FOLDER_RULES/*.rules.mk 2>/dev/null | tr '\n' ' ')

# APP rules
RULES+=" "$(ls $RULES_APP/*.rules.mk 2>/dev/null | tr '\n' ' ')

# Pipelines rules
# For developpers:
#    - beware of similar rules and file names within subfolders pipelines rules
#    - if a pipeline rule is a common rule (e.g. for vcf filtration which is a rule used in multiple steps), just do not create file, or let it empty
#    - use "_" to gather rules 
for rule in $(echo $PIPELINES" "$POST_SEQUENCING_STEPS" "$POST_ALIGNMENT_STEPS" "$POST_CALLING_STEPS" "$POST_ANNOTATION_STEPS" "$POST_CALLING_MERGING_STEPS | tr '.' ' ' | sort -u); do 
	rule_root=$(echo $rule | cut -d_ -f1)
	# Add root rule file
	RULES+=" "$(find $STARK_FOLDER_RULES/*/ -name $rule_root.rules.mk | tr '\n' ' ')
	# Add rule file
	RULES+=" "$(find $STARK_FOLDER_RULES/*/ -name $rule.rules.mk | tr '\n' ' ')
done;

# Uniqify
RULES=$(echo $RULES | tr " " "\n" | sort -u | tr "\n" " ")

# Export
export RULES
