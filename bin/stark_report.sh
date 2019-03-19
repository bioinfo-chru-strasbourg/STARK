#! /bin/sh
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Authors: Amandine VELT
# Date: May 2016
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# in run_analysis.sh, this report is launch with the command : stark_report.sh -f $RUN -p $SAMPLE_PROJECT -g $SAMPLE_GROUP -u $SAMPLE_USER -s $SAMPLE -e $ENV

#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARKReport"
SCRIPT_DESCRIPTION="STARK Report"
SCRIPT_RELEASE="0.9.1b"
SCRIPT_DATE="16/03/2019"
SCRIPT_AUTHOR="Amandine VELT, Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU AGPL V3"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-01/05/2016: Script creation\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.1b-16/03/2019: Add features, tables, new VCF tables\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Configuration
#ENV_CONFIG=$(find $SCRIPT_DIR/.. -name config.app)
#source $ENV_CONFIG



####################################################################################################################################
# Define the function to print the usage of the script
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    sh stark_report.sh -f flowcell -p project -g group -u user -s sample -e env -d date [-h]

		Description:
		    This script allows to generates a readable pdf to see the quality and metrics of each sample.

		Options:
		 	-f, --flowcell The flowcell ID
		 	This option is required.
		 	-p, --project The project ID
		 	This option is NOT required.
		 	-g --group The group name
		 	This option is NOT required.
		 	-u --user The user name
		 	This option is NOT required.
		 	-s, --sample The sample ID
		 	This option is required.
		 	-e --env The environment file
		 	This option is required.
		 	-i --pipelines Pipelines if different from Pipelines option in ENV
		 	This option is not required.
		 	-d --date The date of the current run
		 	This option is required.
		  	-o --output The output PDF
		 	This option is not required.
		  	-h, --help
		 	Print this message and exit the program.
		__EOF__
}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ARGS=$(getopt -o "r:f:p:s:e:i:g:u:d:o:h" --long "results:flowcell:,project:,sample:,env:,pipelines:,group:,user:,date:,output:,help" -- "$@" 2> /dev/null)
[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && \
	exit 1
eval set -- "$ARGS"
while true
do
	case "$1" in
		-r|--results)
			RESULTS_FOLDER_INPUT="$2"
			shift 2
			;;
		-f|--flowcell)
			FLOWCELL="$2"
			shift 2
			;;
		-p|--project)
			PROJECT_INPUT="$2"
			shift 2
			;;
		-g|--group)
			GROUP_INPUT="$2"
			shift 2
			;;
		-u|--user)
			USER_INPUT="$2"
			shift 2
			;;
		-s|--sample)
			SAMPLE="$2"
			shift 2
			;;
		-e|--env)
			ENV="$2"
			shift 2
			;;
		-i|--pipelines)
			PIPELINES_INPUT="$2"
			PIPELINES_INPUT=$(echo $PIPELINES_INPUT | tr "," " ")
			shift 2
			;;
		-d|--date)
			DATE="$2"
			shift 2
			;;
		-o|--output)
			OUTPUT_INPUT="$2"
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
#[ "$FLOWCELL" == "" ] || [ "$PROJECT" == "" ] || [ "$SAMPLE" == "" ] || [ "$ENV" == "" ] || [ "$GROUP" == "" ] || [ "$USER" == "" ] || [ "$DATE" == "" ] &&
[ "$FLOWCELL" == "" ] || [ "$SAMPLE" == "" ] || [ "$ENV" == "" ] || [ "$DATE" == "" ] && \
	echo "Options --flowcell, --project, --sample, --group, --user, --env and --date are required. " "Use -h or --help to display the help." && exit 1;
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# source of environement file or exit
source $ENV || (echo "$ENV can't be sourced ! Exit." && exit 1);

#for E in $ENV; do
#	source $ENV || (echo "$ENV can't be sourced ! Exit." && exit 1);
#done;

#ENV=$(find_app "$ENV")
#echo "ENV=$ENV"
#source_app "$ENV"

# PIPELINES
if [ "$PIPELINES_INPUT" != "" ]; then
	PIPELINES=$PIPELINES_INPUT
fi;

if [ -z "$LATEX" ]; then
	LATEX="latex";
fi;

if [ ! -z $GROUP_INPUT ]; then
	GROUP=$GROUP_INPUT
fi;
if [ -z $GROUP ]; then
	GROUP="UNKNOWN"
fi;
if [ ! -z $PROJECT_INPUT ]; then
	PROJECT=$PROJECT_INPUT
fi;
if [ -z $PROJECT ]; then
	PROJECT="UNKNOWN"
fi;
if [ ! -z $USER_INPUT ]; then
	USER=$USER_INPUT
fi;
if [ -z $USER ]; then
	USER="UNKNOWN"
fi;
if [ ! -z $TMP_FOLDER_TMP ]; then
	TMP_REPORT=$TMP_FOLDER_TMP;
else
	TMP_REPORT=/tmp;
fi;
if [ ! -z $NB_VARIANT_LIMIT ]; then
	NB_VARIANT_LIMIT=$NB_VARIANT_LIMIT;
else
	NB_VARIANT_LIMIT=20;
fi;



# logo must be sourced from env, thn it can be specific to a project etc ...
#LOGO="/home1/DIAG/DATA/NGS/veltaman/test_report/structures_photo_logo_1777.jpg"
LOGO="$STARK/logo_report.jpg"
if [ ! -e $LOGO ]; then LOGO=""; fi;

if [ ! -z $RESULTS_FOLDER_INPUT ] && [ -d $RESULTS_FOLDER_INPUT ]; then
	RESULTS_FOLDER=$RESULTS_FOLDER_INPUT;
fi;

# $MISEQ_FOLDER contains all the raw data (in the $FLOWCELL folder)
# DEMULTIPLEXING_FOLDER contains all the demultiplexing results (in the $FLOWCELL folder)
# RESULTS_FOLDER contains all the other results (in the $FLOWCELL folder)

[ "$MISEQ_FOLDER" == "" ] || [ "$DEMULTIPLEXING_FOLDER" == "" ] || [ "$RESULTS_FOLDER" == "" ] && \
	echo "One of the necessary folder is not defined ! (MISEQ_FOLDER : $MISEQ_FOLDER, DEMULTIPLEXING_FOLDER :$DEMULTIPLEXING_FOLDER and RESULTS_FOLDER : $RESULTS_FOLDER)" && exit 1;

# creation and verification of the directory containing the sample report
#REPORTDIR="$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.reports/latex_report.$DATE"
REPORTDIR="$TMP_REPORT/STARK_REPORT_LATEX_$RANDOM"
mkdir -p $REPORTDIR
if [ ! -e $REPORTDIR ]; then
	echo "$REPORTDIR can't be created ! Exit."
	exit 1;
fi

# creation and verification of the log file
#LOGFILE="$RESULTS_FOLDER/$FLOWCELL/latex_report.V$DATE.log"
#LOGFILE="$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.reports/latex_report.$DATE.log"
LOGFILE="$REPORTDIR/latex_report.$DATE.log"
if [ ! -e "$LOGFILE" ] ; then
    touch "$LOGFILE" || (echo "$LOGFILE can't be created ! Exit." && exit 1);
fi

# creation and verification of the tex file report
#TEXFILE="$REPORTDIR/report.$SAMPLE.V$DATE.tex"
#PDFFILE="$REPORTDIR/report.$SAMPLE.V$DATE.pdf"
TEXFILE="$REPORTDIR/$SAMPLE.$DATE.stark.report.tex"
PDFFILE="$REPORTDIR/$SAMPLE.$DATE.stark.report.pdf"
if [ -z $OUTPUT_INPUT ]; then
	OUTPUT_PDFFILE=$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.reports/$SAMPLE.$DATE.stark.report.pdf
else
	OUTPUT_PDFFILE=$OUTPUT_INPUT
fi;
#STATSFILE="$RESULTS_FOLDER/$FLOWCELL/global_stats_$DATE.txt"
#STATSFILE="$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.reports/global_stats_$DATE.txt"
STATSFILE="$REPORTDIR/global_stats_$DATE.txt"
touch $TEXFILE || (echo "$TEXFILE can't be created ! Exit." && exit 1);
touch $STATSFILE || (echo "$STATSFILE can't be created ! Exit." && exit 1);

# ALIGNERS

ALIGNERS=$(echo $PIPELINES | tr " " "\n"  | cut -d. -f1 | sort | uniq | tr "\n" " ")
CALLERS=$(echo $PIPELINES | tr " " "\n"  | cut -d. -f2 | sort | uniq | tr "\n" " ")
ANNOTATORS=$(echo $PIPELINES | tr " " "\n"  | cut -d. -f3 | sort | uniq | tr "\n" " ")

# SAMPLEID
SAMPLEID=$(echo $SAMPLE | sed "s/\_/\\\_/gi")
#echo "SAMPLEID=$SAMPLEID"; exit 0;

# INFOS
APP=$(echo $(basename $ENV) | sed "s/^env.//gi" | sed "s/.sh$//gi" | sed "s/sh$//gi")
#SAMPLE_GROUP=$(echo $APP | awk -F- '{print $1}'); if [ "$SAMPLE_GROUP" == "" ]; then SAMPLE_GROUP="UNKNOWN"; fi;
#SAMPLE_PROJECT=$(echo $APP | awk -F- '{print $2}'); if [ "$SAMPLE_PROJECT" == "" ]; then SAMPLE_PROJECT="UNKNOWN"; fi;
#SAMPLE_USER=$(echo $APP| awk -F- '{print $3}'); if [ "$SAMPLE_USER" == "" ]; then SAMPLE_USER="UNKNOWN"; fi;
SAMPLE_GROUP=$GROUP; if [ "$SAMPLE_GROUP" == "" ]; then SAMPLE_GROUP="UNKNOWN"; fi;
SAMPLE_PROJECT=$PROJECT; if [ "$SAMPLE_PROJECT" == "" ]; then SAMPLE_PROJECT="UNKNOWN"; fi;
SAMPLE_USER=$USER; if [ "$SAMPLE_USER" == "" ]; then SAMPLE_USER="UNKNOWN"; fi;
if [ -z "$APP" ]; then APP="DEFAULT"; fi;


# PANEL name
MANIFEST_NAME_FILE=$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.manifest_name
BED_NAME_FILE=$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.bed_name
MANIFEST_NAME=""
MANIFEST_SOURCE=""
BED_NAME=""
BED_SOURCE=""
if [ -s "$MANIFEST_NAME_FILE" ]; then
	MANIFEST_NAME=$(basename $(awk -F"\t" '{print $1}' $MANIFEST_NAME_FILE)) # awk -F\t '{print $2}'
	MANIFEST_SOURCE=$(awk -F"\t" '{print $2}' $MANIFEST_NAME_FILE)
fi;
if [ -s "$BED_NAME_FILE" ]; then
	BED_NAME=$(basename $(awk -F"\t" '{print $1}' $BED_NAME_FILE))
	BED_SOURCE=$(awk -F"\t" '{print $2}' $BED_NAME_FILE)
fi;
if [ -z "$MANIFEST_NAME" ]; then
	MANIFEST_NAME="unknown";
fi;
if [ -z "$BED_NAME" ]; then
	BED_NAME="unknown";
fi;
if [ ! -z "$MANIFEST_SOURCE" ]; then
	MANIFEST_SOURCE=" ($MANIFEST_SOURCE)";
fi;
if [ ! -z "$BED_SOURCE" ]; then
	BED_SOURCE=" ($BED_SOURCE)";
fi;



echo "\documentclass[a4paper]{report}
\usepackage[a4paper]{geometry}
\usepackage{longtable}
\usepackage{supertabular}
\usepackage{lipsum}
\usepackage[group-separator={,}]{siunitx}
\usepackage[dvips]{graphicx}
\selectlanguage{english}
\usepackage{fancybox}
\usepackage[english]{babel}
\usepackage{listings}
\usepackage{hyperref}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage{amsmath}
\usepackage{color}
\usepackage{xcolor}
\usepackage[table]{xcolor}% http://ctan.org/pkg/xcolor
\usepackage{tikz,array,collcell}
\usepackage{listings}
\usepackage{tcolorbox}
\usepackage{lastpage}
\usepackage{tabularx}
\usepackage{longtable}
\usepackage{array}% http://ctan.org/pkg/array
\makeatletter
\g@addto@macro{\endtabular}{\rowfont{}}% Clear row font
\makeatother
\newcommand{\rowfonttype}{}% Current row font
\newcommand{\rowfont}[1]{% Set current row font
   \gdef\rowfonttype{#1}#1%
}
\newcolumntype{L}{>{\rowfonttype}l}
\usepackage{multirow}
\usepackage[justification=justified,singlelinecheck=off]{caption}
\setlength{\abovecaptionskip}{2pt}
\usepackage{fancyhdr}
\addtolength{\headheight}{40pt}
\addtolength{\topmargin}{-10pt}
\addtolength{\textheight}{-10pt}
\pagestyle{fancy}
\fancyhf{}% clear all headers and footer
%\renewcommand{\headrulewidth}{0pt} % no line for the header
\renewcommand{\footrulewidth}{1pt}

\chead{
	\begin{tabularx}{\textwidth}{cXc}
	\begin{minipage}[c]{0.15\textwidth}
		\vspace{3pt}
		\centering \\includegraphics[width=0.8\\textwidth]{$LOGO}
		\vspace{3pt}
	\end{minipage} &
	\centering $SAMPLEID Report &
	\begin{minipage}[c]{0.35\textwidth}
		\centering [$APP]
	\end{minipage}
	\end{tabularx}
}

\lfoot{\vspace{-9pt}\scriptsize This report was generated automatically}
\cfoot{\thepage\ / \pageref{LastPage}}
\rfoot{\today}
\fancypagestyle{plain}{} % prevent chapter page style to differ from the rest, works if nothing is inside (strange :s)


\usepackage{graphicx}
\usepackage{subfig}
\usepackage{float}
\usepackage[autolanguage]{numprint}
\usepackage{xcolor, colortbl}
\usepackage{tocloft}

\setcounter{tocdepth}{6}
\setcounter{secnumdepth}{4}
\renewcommand\thesection{\arabic{section}.}
\renewcommand\thesubsection{\hspace{3mm}\thesection\arabic{subsection}.}
\renewcommand\thesubsubsection{\hspace{6mm}\thesubsection\arabic{subsubsection}.}
\setlength{\cftbeforetoctitleskip}{2mm}  % vertical space after and before the summary title
\setlength{\cftaftertoctitleskip}{1mm}
\setlength{\cftsecindent}{10pt}
\setlength{\cftsecnumwidth}{10pt}
\setlength{\cftsubsecindent}{20pt}
\setlength{\cftsubsecnumwidth}{27pt}
\setlength{\cftsubsubsecindent}{48pt}
\setlength{\cftsubsubsecnumwidth}{52pt}
\renewcommand{\cfttoctitlefont}{\large\bfseries}

\usepackage{setspace}
\usepackage[nomessages]{fp}% http://ctan.org/pkg/fp
\newcommand{\maxnum}{100.00}
\newlength{\maxlen}
\newcommand{\databar}[2][green]{%
  \settowidth{\maxlen}{\maxnum}%
  \addtolength{\maxlen}{\tabcolsep}%
  \FPeval\result{round(#2/\maxnum:4)}%
  \rlap{\color{green}\hspace*{-.5\tabcolsep}\rule[-.05\ht\strutbox]{\result\maxlen}{.95\ht\strutbox}}%
  \makebox[\dimexpr\maxlen-\tabcolsep][r]{#2}%
}
\begin{document}
\tableofcontents

\section{General information}
\begin{tabular}{l l}
Sample: & ` echo $SAMPLEID`\\\\
Run: & ` echo $FLOWCELL | sed 's/\\_/\\\_/g'`\\\\
Group: & ` echo $GROUP | sed 's/\\_/\\\_/g'`\\\\
Project: & ` echo $PROJECT | sed 's/\\_/\\\_/g'`\\\\
Assembly: & ` echo $ASSEMBLY | sed 's/\\_/\\\_/g'`\\\\
STARK: & ` echo $STARK_VERSION | sed 's/\\_/\\\_/g'`\\\\
Application: & ` echo $APP_NAME | sed 's/\\_/\\\_/g'` (` echo $APP_RELEASE | sed 's/\\_/\\\_/g'`)\\\\
Manifest: & ` echo $MANIFEST_NAME | sed 's/\\_/\\\_/g'` ` echo $MANIFEST_SOURCE | sed 's/\\_/\\\_/g'`\\\\
Bed: & ` echo $BED_NAME | sed 's/\\_/\\\_/g'` ` echo $BED_SOURCE | sed 's/\\_/\\\_/g'`\\\\
\end{tabular}



\section{FastQC results}
FastQC performs analyses to assess the quality of the data. Informative controls are shown in Figure~\ref{quality} to~\ref{adapter}. Table~\ref{fastqc_recap} shows general information" >> $TEXFILE

#Application: & $APP\\\\
# & (`basename $ENV`)\\\\

#FASTQCDATA="$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.fastqc/$SAMPLE.unaligned_fastqc/fastqc_data.txt"
FASTQCDATA="$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.fastqc/metrics.fastqc.txt"
FASTQCIMAGESDIR="$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.fastqc/$SAMPLE.unaligned_fastqc/Images"

Filename=$(grep "Filename" $FASTQCDATA | cut -f2  | sed "s/\_/\\\_/gi")
TotalSequences=`grep "Total Sequences" $FASTQCDATA | cut -f2`
SequenceLength=`grep "Sequence length" $FASTQCDATA | cut -f2`
if [ -s $RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.fastqc/metrics.counts.txt ]; then
	UniqueNumberReads=$(grep "^unique" $RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.fastqc/metrics.counts.txt | cut -f2);
	PercentUniqueNumberReads=$(grep "^%unique" $RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.fastqc/metrics.counts.txt | cut -f2);
elif [ -s $RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.unaligned.bam ]; then
	UniqueNumberReads=$( $SAMTOOLS view $RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.unaligned.bam | cut -f10 | sort -T "$TMP_FOLDER_TMP" -u | wc -l )
	PercentUniqueNumberReads=$( echo "scale=2; $UniqueNumberReads*100/$TotalSequences" | bc )
else
	UniqueNumberReads=$(zcat $RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE*.fastq.gz 2>/dev/null | awk 'NR % 4 == 2' | sort -T "$TMP_FOLDER_TMP" -u | wc -l);
	PercentUniqueNumberReads=$( echo "scale=2; $UniqueNumberReads*100/$TotalSequences" | bc )
#else
#	UniqueNumberReads=0
fi;
#UniqueNumberReads=$(grep "Total Sequences" $FASTQCDATA | cut -f2)
#HORIZON_bam_env_CPSGEN_sh.fastqc\HORIZON_bam_env_CPSGEN_sh.unaligned_fastqc

# Q30
Q30_ALL=$(cat $FASTQCDATA | awk '/>>Per sequence quality scores/,/>>END_MODULE/' | head -n -1 | tail -n+3 | awk '{s+=$2}END{print s}');
Q30_UNTIL=$(cat $FASTQCDATA | awk '/>>Per sequence quality scores/,/>>END_MODULE/' | awk '/>>Per sequence quality scores/,/30\t/' | head -n -1 | tail -n+3 | awk '{s+=$2}END{print s}');
Q30=$(bc <<< "scale = 4; (($Q30_ALL - $Q30_UNTIL) / $Q30_ALL * 100)");
Q30_print=$(echo $Q30 | sed s/00$//)"\%"


echo "\begin{table}[!h]
\begin{center}
\begin{tabular}{| c | c | c | c | c | c |}
\hline
Filename & $Filename \\\\
\hline
Total number of reads & \num[group-separator={,}]{$TotalSequences} \\\\
\hline
Unique number of reads & \num[group-separator={,}]{$UniqueNumberReads} \\\\
\hline
\% of unique number reads & $PercentUniqueNumberReads \\\\
\hline
Sequence length & $SequenceLength \\\\
\hline
Q30 & $Q30_print \\\\
\hline
\end{tabular}
\end{center}
\caption{\label{fastqc_recap}Basic information}
\end{table}

Note : the reads length varying because we trimm the adapters from the reads with bcl2fastq2 of Illumina.

\begin{figure}[H]
\centering
\includegraphics[width=0.8\textwidth]{$FASTQCIMAGESDIR/per_base_quality.png}
\caption{\label{quality}Per base quality.}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=0.8\textwidth]{$FASTQCIMAGESDIR/per_base_n_content.png}
\caption{\label{n_content}Proportion of undetermined nucleotide (N) along the read.}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=0.8\textwidth]{$FASTQCIMAGESDIR/per_base_sequence_content.png}
\caption{Proportion of each nucleotide (A, C, G, T) along the read.}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=0.8\textwidth]{$FASTQCIMAGESDIR/duplication_levels.png}
\caption{\label{duplication}Duplicate sequences.}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=0.8\textwidth]{$FASTQCIMAGESDIR/adapter_content.png}
\caption{\label{adapter}Adapter content.}
\end{figure}

\section{Alignments}
Sequenced reads were mapped to the $ASSEMBLY assembly version using the \"$ALIGNERS\" aligner(s). Alignment statistics are provided in following table(s) (one per aligner)." >> $TEXFILE

# one table of alignments statistics per aligner.
for aligner in `echo $ALIGNERS`
do
	STATFILE="$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.$aligner.bam.metrics/$SAMPLE.$aligner.flagstat"
	Total_qc_passed_reads=`grep "in total" $STATFILE | cut -d"+" -f1 | sed 's/ //'`
	#Percent_mapped_reads=`grep "mapped (" $STATFILE | cut -d"(" -f2 | sed 's/%:.*//'`
	Percent_mapped_reads=`grep "mapped (" $STATFILE | cut -d"(" -f2 | sed 's/[%:].*//'`
	Total_mapped_reads=`grep "mapped (" $STATFILE | cut -d"+" -f1 | sed 's/ //'`
	echo "\begin{table}[!h]
	\begin{center}
	\begin{tabular}{| c | c | c | c |}
	\hline
	Sample & QC-passed Reads & Mapped Reads & \\% Mapped Reads \\\\
	\hline
	$SAMPLEID & \num[group-separator={,}]{$Total_qc_passed_reads} & \num[group-separator={,}]{$Total_mapped_reads} & $Percent_mapped_reads \\\\
	\hline
	\end{tabular}
	\end{center}
	\caption{\label{stats_$aligner}Alignment statistics with $aligner on $SAMPLEID. Statistics were determined with Samtools flagstat on the $SAMPLEID.$aligner.bam file.}
	\end{table}" >> $TEXFILE
done

if [[ -f "$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.$aligner.markDuplicates.metrics.txt" ]]
then
for aligner in `echo $ALIGNERS`
do
	MARKDUPLICATESFILE="$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.$aligner.markDuplicates.metrics.txt"
	READ_PAIRS_EXAMINED=`cat $MARKDUPLICATESFILE | grep -A 2 "## METRICS CLASS" | tail -n 1 | cut -f3 | sed 's/,/./'`
	PERCENT_DUPLICATION=`cat $MARKDUPLICATESFILE | grep -A 2 "## METRICS CLASS" | tail -n 1 | cut -f8 | sed 's/,/./'`
	echo "\section{MarkDuplicates}
	\begin{table}[!h]
	\begin{center}
	\begin{tabular}{| c | c | c |}
	\hline
	Sample & Read pairs examined & Percent duplication \\\\
	\hline
	$SAMPLEID & \num[group-separator={,}]{$READ_PAIRS_EXAMINED} & `echo "scale=2;($PERCENT_DUPLICATION * 100)/1" | bc` \\\\
	\hline
	\end{tabular}
	\end{center}
	\caption{\label{stats_markduplicates_$aligner}Markduplicate statistics on $aligner on $SAMPLEID. Statistics were determined with GATK Markduplicates on the $SAMPLEID.$aligner.bam file.}
	\end{table}" >> $TEXFILE
done
fi


	echo "\section{Coverage}" >> $TEXFILE

for aligner in `echo $ALIGNERS`
do

	# METRICS File
	COVERAGEFILESTATS_SAMTOOLS="$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.$aligner.bam.metrics/$SAMPLE.$aligner.depthbed"
	COVERAGEFILESTATSOFF_SAMTOOLS="$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.$aligner.bam.metrics/$SAMPLE.$aligner.off.depthbed"
	COVERAGEFILESTATS_ONREADS_SAMTOOLS="$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.$aligner.bam.metrics/$SAMPLE.$aligner.on.nbreads"
	COVERAGEFILESTATS_OFFREADS_SAMTOOLS="$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.$aligner.bam.metrics/$SAMPLE.$aligner.off.nbreads"
	COVERAGEFILESTATS_HSMETRICS="$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.$aligner.bam.metrics/$SAMPLE.$aligner.HsMetrics"
	COVERAGEFILESTATS_GENES_COVERAGE_MSG="$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.$aligner.bam.metrics/$SAMPLE.msg"

	# COVERAGE SUMMARY AVERAGE on all sequencing with SAMTOOLS
	if [ -s $COVERAGEFILESTATS_SAMTOOLS ]; then
		COV_CALC=$( cat $COVERAGEFILESTATS_SAMTOOLS | awk '{SUM+=$3; LINES+=1} END {print SUM/LINES}' )
		echo "Sequencing coverage estimation on target ( total number of sequenced bases / number bases of target) : `echo "scale=2;($COV_CALC)/1" | bc` X" >> $TEXFILE
	fi


	# COVERAGE SUMMARY ON/OFF target with HsMetrics
	if [ -s $COVERAGEFILESTATS_HSMETRICS ] && ((1)); then
		PCT_SELECTED_BASES=$( C=1; for i in $(grep -v "#" $COVERAGEFILESTATS_HSMETRICS | sed '/^\s*$/d' | head -n 2) ; do if [ $i == "PCT_SELECTED_BASES" ] ; then break ; else C=$(( $C + 1 )) ; fi ; done ; grep -v "#" $COVERAGEFILESTATS_HSMETRICS | sed '/^\s*$/d'  | head -n2 | tail -n1 | cut -f$C | sed 's/,/./' )
		PCT_OFF_BAIT=$( C=1; for i in $(grep -v "#" $COVERAGEFILESTATS_HSMETRICS | sed '/^\s*$/d' | head -n 2) ; do if [ $i == "PCT_OFF_BAIT" ] ; then break ; else C=$(( $C + 1 )) ; fi ; done ; grep -v "#" $COVERAGEFILESTATS_HSMETRICS | sed '/^\s*$/d'  | head -n2 | tail -n1 | cut -f$C | sed 's/,/./' )
		#echo "PCT_SELECTED_BASES=$PCT_SELECTED_BASES PCT_OFF_BAIT=$PCT_OFF_BAIT"

		if [ "$PCT_SELECTED_BASES" != "" ] && [ "$PCT_OFF_BAIT" != "" ]; then

			echo "\begin{table}[H]
			\begin{center}
			\begin{tabular}{| c | c | c | }
			\hline
			Sample & \% ON target & \% OFF target \\\\
			\hline
			$SAMPLEID & `echo "scale=2;($PCT_SELECTED_BASES*100)/1" | bc` & `echo "scale=2;($PCT_OFF_BAIT*100)/1" | bc` \\\\
			\hline
			\end{tabular}
			\end{center}
			\caption{\label{stats_on_off_$aligner}Coverage statistics between on and off targets reads on $aligner on $SAMPLEID. Statistics were determined with Picard HSMetrics on the $SAMPLEID.$aligner.bam file.}
			\end{table}" >> $TEXFILE

		fi;
	fi;

	if [ -e $COVERAGEFILESTATS_ONREADS_SAMTOOLS ] && [ -e $COVERAGEFILESTATS_OFFREADS_SAMTOOLS ]; then
		READS_ON_TARGET=$(echo $(cat $COVERAGEFILESTATS_ONREADS_SAMTOOLS)" +0 " | bc);
		READS_OFF_TARGET=$(echo $(cat $COVERAGEFILESTATS_OFFREADS_SAMTOOLS)" +0 " | bc);
		READS_ALL_TARGET=$(echo "$READS_ON_TARGET+$READS_OFF_TARGET" | bc)
		if [ "$READS_ON_TARGET" != "" ] && [ "$READS_OFF_TARGET" != "" ]; then

			echo "\begin{table}[H]
			\begin{center}
			\begin{tabular}{| c | c | c | }
			\hline
			Sample & ON target & OFF target  \\\\
			\hline
			$SAMPLEID & $READS_ON_TARGET (`echo "scale=2;(($READS_ON_TARGET/$READS_ALL_TARGET)*100)/1" | bc`\%) & $READS_OFF_TARGET (`echo "scale=2;(($READS_OFF_TARGET/$READS_ALL_TARGET)*100)/1" | bc`\%) \\\\
			\hline
			\end{tabular}
			\end{center}
			\caption{\label{stats_on_off_$aligner}Coverage statistics between on and off targets reads on $aligner on $SAMPLEID. Statistics were determined with SAMTOOLS/BEDTOOLS (duplicates and bad quality reads removed) on the $SAMPLEID.$aligner.bam file.}
			\end{table}" >> $TEXFILE

		fi;



	fi;

	# COVERAGE SUMMARY with SAMTOOLS
	if [ -s $COVERAGEFILESTATS_SAMTOOLS ]
	then

		if [ -z $COVS ]; then
			COVS="5,10,20,30,50,100,200,300";
		fi;
		MDP=$(echo $COVS | tr "," "\n" | sort -n | tail -n 1);

		cat $COVERAGEFILESTATS_SAMTOOLS | awk -v MDP=$MDP  '{SUM++} { if ($3>MDP) {DP[MDP]++} else {DP[$3]++} } END { for (i=MDP; i>=0; i-=1) {print i" "DP[i]" SUM"SUM}}' | sort -g -r | awk -v COVS=$COVS '{SUM+=$2} {CUM[$1]=SUM} {split(COVS,C,",")}  END  { for (j in C) {print C[j]" "CUM[C[j]]" SUM "(CUM[C[j]]/SUM)} }' | sort -g > $COVERAGEFILESTATS_SAMTOOLS.summary.tmp

		#cat $COVERAGEFILESTATS_SAMTOOLS.summary.tmp

		cat $COVERAGEFILESTATS_SAMTOOLS.summary.tmp | awk '{print $1"X & "($4*100)"\\% \\\\ \\hline"}' > $COVERAGEFILESTATS_SAMTOOLS.summary.tmp.latex


		echo "\begin{longtable}{|r|l|}

			\hline \multicolumn{1}{|c|}{\textbf{Coverage}} & \multicolumn{1}{c|}{\textbf{\\%targeted bases}}  \\\\ \hline
			\endfirsthead

			\multicolumn{2}{c}%
			{{\bfseries \tablename\ \thetable{} -- continued from previous page}} \\\\
			\hline \multicolumn{1}{|c|}{\textbf{Coverage}} & \multicolumn{1}{c|}{\textbf{\\%targeted bases}} \\\\ \hline
			\endhead

			\hline \multicolumn{2}{|r|}{{Continued on next page}} \\\\ \hline
			\endfoot

			\hline
			\endlastfoot
			" >> $TEXFILE
		cat $COVERAGEFILESTATS_SAMTOOLS.summary.tmp.latex >> $TEXFILE
		echo "
			\caption{Coverage statistics on targeted regions on $aligner on $SAMPLEID. Statistics were determined with SAMTOOLS (duplicates and bad quality reads removed) on the $SAMPLEID.$aligner.bam file.} \label{tab:long} \\\\
			\end{longtable}

			" >> $TEXFILE


		rm $COVERAGEFILESTATS_SAMTOOLS.summary.tmp $COVERAGEFILESTATS_SAMTOOLS.summary.tmp.latex;

	fi;
	#fi;

	# COVERAGE SUMMARY with SAMTOOLS
	if [ -s $COVERAGEFILESTATS_GENES_COVERAGE_MSG ]
	then
		echo "\begin{itemize}" >> $TEXFILE
		echo "\item "$(awk -F"\t" '$2=="PASS" {print $1}' $COVERAGEFILESTATS_GENES_COVERAGE_MSG | wc -w)" genes passing the coverage threshold (more than "$(echo "$DP_THRESHOLD * 100" | bc)" percent bases with DP greater than $DP_WARN)" >> $TEXFILE
		echo "\item "$(awk -F"\t" '$2=="WARN" {print $1}' $COVERAGEFILESTATS_GENES_COVERAGE_MSG | wc -w)" genes with a warning coverage (less than "$(echo "$DP_THRESHOLD * 100" | bc)" percent bases with DP greater than $DP_WARN) " >> $TEXFILE
		echo $(awk -F"\t" '$2=="WARN" {print $1}' $COVERAGEFILESTATS_GENES_COVERAGE_MSG) >> $TEXFILE
		echo "\item "$(awk -F"\t" '$2=="FAIL" {print $1}' $COVERAGEFILESTATS_GENES_COVERAGE_MSG | wc -w)" genes with a failed coverage (less than "$(echo "$DP_THRESHOLD * 100" | bc)" percent bases with DP greater than $DP_FAIL) " >> $TEXFILE
		echo $(awk -F"\t" '$2=="FAIL" {print $1}' $COVERAGEFILESTATS_GENES_COVERAGE_MSG) >> $TEXFILE
		echo "\end{itemize}" >> $TEXFILE
	fi;


done


#rm  ${STATSFILE}
if grep -q "Sample" ${STATSFILE}; then
	:
else
	echo -ne "Sample ID\tTotal number of reads\tSequences length\tUnique number of reads\tDuplication percentage" >> ${STATSFILE}
	for aligner in `echo $ALIGNERS`
	do
		echo -ne "\tTotal QC passed reads ($aligner)\tPercent mapped reads ($aligner)\tTotal mapped reads ($aligner)" >> ${STATSFILE}
		for caller in `echo $CALLERS`
		do
			for annotator in `echo $ANNOTATORS`
			do
				if [[ "$caller" != "canoe"  ]]
				then
					echo -ne "\tNumber find variants (${aligner}-${caller})" >> ${STATSFILE}
				fi
			done
		done
	done

	echo -ne "\tSequencing coverage estimation\tPercent target covered > 10 X\tPercent target covered > 20 X\tPercent target covered > 30 X\tPercent target covered > 40 X\tPercent target covered > 50 X\tPercent target covered > 100 X"  >> ${STATSFILE}


fi
#echo $PERCENT_DUPLICATION; exit 0;
if [ ! -z $PERCENT_DUPLICATION ]; then
	PERCENT_DUPLICATION_SHOW=$(echo "scale=2;($PERCENT_DUPLICATION * 100)/1" | bc)
else
	PERCENT_DUPLICATION_SHOW="?"
fi;
#exit 0;

echo -ne "\n${SAMPLE}\t${TotalSequences}\t${SequenceLength}\t${UniqueNumberReads}\t$PERCENT_DUPLICATION_SHOW\t" >> ${STATSFILE}

#exit 0;





echo "\section{Variant calling}" >> $TEXFILE
echo "\begin{itemize}" >> $TEXFILE

# FINAL VCF
FINALVCFFILE=$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.reports/$SAMPLE.final.vcf
FINALTXTFILE=$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.reports/$SAMPLE.final.tsv
number_variants_final=$(grep -cv ^# $FINALVCFFILE)
echo "\item Total number of variants: $number_variants_final" >> $TEXFILE

for aligner in `echo $ALIGNERS`
do
	STATFILE="$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.$aligner.bam.metrics/$SAMPLE.$aligner.flagstat"
	Total_qc_passed_reads=`grep "in total" $STATFILE | cut -d"+" -f1 | sed 's/ //'`
	Percent_mapped_reads=`grep "mapped (" $STATFILE | cut -d"(" -f2 | sed 's/%:.*//'`
	Total_mapped_reads=`grep "mapped (" $STATFILE | cut -d"+" -f1 | sed 's/ //'`
	echo -ne "${Total_qc_passed_reads}\t${Percent_mapped_reads}\t${Total_mapped_reads}" >> ${STATSFILE}
	for caller in `echo $CALLERS`
	do
		for annotator in `echo $ANNOTATORS`
		do
			if [[ "$caller" != "canoe"  ]]
			then
				#VARIANTSFILESTATS="$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.$aligner.$caller.$annotator.vcf.metrics/$SAMPLE.$aligner.$caller.$annotator.vcf.snpeff.metrics.html"
				#if [ -e $VARIANTSFILESTATS ]; then
				#number_variants=$( grep -A 1 "Number of variants (before filter)" $VARIANTSFILESTATS | tail -n1 | sed 's/<td> //' | sed 's/ <\/td>//' | tr -d '\t' )
				VCFFILE=$RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.$aligner.$caller.$annotator.vcf
				if [ -e $VCFFILE ]; then
					number_variants=$(grep -cv ^# $VCFFILE)

				else
					number_variants="?"
				fi;
				echo -ne "\t${number_variants}" >> ${STATSFILE}
				echo "\item Number of variants find by pipeline ` echo $aligner.$caller.$annotator | sed 's/\\_/\\\_/g'`: $number_variants" >> $TEXFILE
				#fi;
			fi
		done
	done
done

echo "\end{itemize}" >> $TEXFILE



if ((1)); then



	echo "\section{Variants prioritization}" >> $TEXFILE

	#NB_FILTERED=$(echo $(awk -F "\t" '$9 != "PASS" {print $0}' $FINALTXTFILE | wc -l)" -1 " | bc)
	#NB_PASS=$(awk -F "\t" '$9 == "PASS" {print $0}' $FINALTXTFILE | wc -l)

	echo "\begin{itemize}" >> $TEXFILE
	echo "\item Prioritization filters (only first is listed below): "$(echo $HOWARD_PRIORITIZATION_REPORT | tr "," " " | sed 's/\_/\\_/g') >> $TEXFILE
	#echo "\item Number of variants filtered by the prioritization filter '"$(echo $HOWARD_PRIORITIZATION_REPORT | cut -d, -f1 | sed 's/\\_/\\\_/g')"': $NB_FILTERED" >> $TEXFILE
	#echo "\item Number of variants passing the prioritization filter '"$(echo $HOWARD_PRIORITIZATION_REPORT | cut -d, -f1 | sed 's/\\_/\\\_/g')"': $NB_PASS" >> $TEXFILE
	echo "\end{itemize}" >> $TEXFILE

	#if (($NB_PASS)); then

		#for $HOWARD_FILTER_one in $(echo $HOWARD_FILTER | tr "," " "); do echo $H;

		#echo "LIST: "$(echo "DEFAULT "$(echo $HOWARD_PRIORITIZATION_REPORT | cut -d, -f2- | tr "," " "))

		for HOWARD_PRIORITIZATION_REPORT_ONE in $(echo "DEFAULT "$(echo $HOWARD_PRIORITIZATION_REPORT | cut -d, -f2- | tr "," " ")); do


			if [ $HOWARD_PRIORITIZATION_REPORT_ONE == "DEFAULT" ]; then
				PZFLAG_FIELD="PZFlag"
				HOWARD_PRIORITIZATION_REPORT_NAME=$(echo $HOWARD_PRIORITIZATION_REPORT | cut -d, -f1 | sed 's/\_/\\_/g')
			else
				PZFLAG_FIELD="PZFlag-$HOWARD_PRIORITIZATION_REPORT_ONE"
				HOWARD_PRIORITIZATION_REPORT_NAME=$(echo $HOWARD_PRIORITIZATION_REPORT_ONE | sed 's/\_/\\_/g')
			fi;

			NB_FILTERED=$(cat $FINALTXTFILE | awk -f $SCRIPT_DIR/tsv_extract.awk -F'\t' -v COLS="$PZFLAG_FIELD" | grep ^# -v | awk -F'\t' '$1 != "PASS" {print $0}' | wc -l)
			NB_PASS=$(cat $FINALTXTFILE | awk -f $SCRIPT_DIR/tsv_extract.awk -F'\t' -v COLS="$PZFLAG_FIELD"  | grep ^# -v | awk -F'\t' '$1 == "PASS" {print $0}' | wc -l)

			#echo "$HOWARD_PRIORITIZATION_REPORT_ONE $PZFLAG_FIELD $HOWARD_PRIORITIZATION_REPORT_NAME $NB_FILTERED $NB_PASS"

			echo "\begin{itemize}" >> $TEXFILE
			#echo "\item Prioritization filters (only first is listed below): "$(echo $HOWARD_PRIORITIZATION_REPORT | tr "," " " | sed 's/\\_/\\\_/g') >> $TEXFILE
			echo "\item Number of variants filtered by the prioritization filter '$HOWARD_PRIORITIZATION_REPORT_NAME': $NB_FILTERED" >> $TEXFILE
			echo "\item Number of variants passing the prioritization filter '$HOWARD_PRIORITIZATION_REPORT_NAME': $NB_PASS" >> $TEXFILE
			echo "\end{itemize}" >> $TEXFILE

			if (($NB_PASS)); then

				echo "\begin{longtable}{|p{8cm}|p{5cm}|p{1.2cm}|}
				\caption{Variants passing prioritization filter '$HOWARD_PRIORITIZATION_REPORT_NAME' (first $NB_VARIANT_LIMIT variants).} \label{tab:long} \\\\

				\hline \multicolumn{1}{|p{7cm}|}{\textbf{HGVS}} & \multicolumn{1}{p{5.5cm}|}{\textbf{Location Outcome Impact}} & \multicolumn{1}{p{1cm}|}{\textbf{VAF}} \\\\ \hline
				\endfirsthead

				\multicolumn{3}{c}%
				{{\bfseries \tablename\ \thetable{} -- continued from previous page}} \\\\
				\hline \multicolumn{1}{|p{7cm}|}{\textbf{HGVS}} & \multicolumn{1}{p{5.5cm}|}{\textbf{Location Outcome Impact}} & \multicolumn{1}{p{1cm}|}{\textbf{VAF}} \\\\ \hline
				\endhead

				\hline \multicolumn{3}{|r|}{{Continued on next page}} \\\\ \hline
				\endfoot

				\hline \hline
				\endlastfoot

					" >> $TEXFILE
					#awk -F "\t" '$9 == "PASS" {print $11" & "$13" & "$14" & \begin{left} "$12" \end{left} \\\\ \\hline"}' $FINALTXTFILE | sed "s/,/ /gi" | sed 's/_/\\_/gi' >> $TEXFILE
					#awk -F "\t" '$9 == "PASS" {print $11" & "$13" & "$14" & \\hline "$12" \\\\ \\hline"}' $FINALTXTFILE | cut -d, -f1 | sed 's/_/\\_/gi' >> $TEXFILE
					#awk -F "\t" '{split($8,H,",")} {V=$16*100} length(H[1])>50 {N=substr(H[1],1,47)"..." } length(H[1])<=50 {N=H[1]} $9 == "PASS" {print "N" & "$14"/"$15"} & "V"\\% \\\\ \\hline"}' $FINALTXTFILE | sed 's/_/\\_/gi' | sed 's/>/\\textgreater /gi' >> $TEXFILE
					#awk -F "\t" '{split($8,H,",")} {split($15,Ol,",")} {V=$16*100} length(H[1])>50 {N=substr(H[1],1,47)"..." } length(H[1])<=50 {N=H[1]} {O=Ol[1]} $9 == "PASS" {print "{\\small "N"} & {\\small "$14"/"$O"} & {\\small "V"\\%} \\\\ \\hline"}' $FINALTXTFILE | sed 's/_/\\_/gi' | sed 's/>/\\textgreater /gi' >> $TEXFILE

					#cat $FINALTXTFILE | sed 's/&/ /gi' | awk -F "\t" '{split($8,H,",")} {split($15,Ol,",")} {V=$16*100} length(H[1])>50 {N=substr(H[1],1,47)"..." } length(H[1])<=50 {N=H[1]} {O=Ol[1]} {L=$14} {I=$20} $9 == "PASS" {print "{\\small "N"} & {\\small "L" "O" "I"} & {\\small "V"\\%} \\\\ \\hline"}' | sed 's/_/\\_/gi' | sed 's/>/\\textgreater /gi' >> $TEXFILE

					#echo "cat $FINALTXTFILE | awk -f $SCRIPT_DIR/tsv_extract.awk -F'\t' -v COLS='NOMEN,location,outcome,snpeff_impact,VAF_median,$PZFLAG_FIELD,#CHROM,POS,REF,ALT'"
					#cat $FINALTXTFILE | awk -f $SCRIPT_DIR/tsv_extract.awk -F"\t" -v COLS="NOMEN,location,outcome,snpeff_impact,VAF_median,$PZFLAG_FIELD,#CHROM,POS,REF,ALT" -v DEFAULT_NA=""
					#cat $FINALTXTFILE | awk -f $SCRIPT_DIR/tsv_extract.awk -F"\t" -v COLS="NOMEN,location,outcome,snpeff_impact,VAF_median,$PZFLAG_FIELD,#CHROM,POS,REF,ALT" -v DEFAULT_NA="" | sed 's/&/ /gi' | awk -F "\t" '{split($1,NOMEN,",")} {CHROMPOSREFALT=$7":"$8$9">"$10} length(NOMEN[1])==0 {NOMEN[1]=CHROMPOSREFALT} length(NOMEN[1])<=50 {HGVS=NOMEN[1]} length(NOMEN[1])>50 {HGVS=substr(NOMEN[1],1,47)"..." } {OUTPUT=$2} {LOCATION=$3} {IMPACT=$4} {VAF=$5*100} {PZFlag=$6} PZFlag == "PASS" {print "{\\small "HGVS"} & {\\small "LOCATION" "OUTPUT" "IMPACT"} & {\\small "VAF"\\%} \\\\ \\hline"}' | sed 's/_/\\_/gi' | sed 's/>/\\textgreater /gi'
					#cat $FINALTXTFILE | awk -f $SCRIPT_DIR/tsv_extract.awk -F"\t" -v COLS="NOMEN,location,outcome,snpeff_impact,VAF_median,$PZFLAG_FIELD"

					#cat $FINALTXTFILE | awk -f $SCRIPT_DIR/tsv_extract.awk -F"\t" -v COLS="NOMEN,location,outcome,snpeff_impact,VAF_median,$PZFLAG_FIELD" | sed 's/&/ /gi' | awk -F "\t" '{split($1,NOMEN,",")} length(NOMEN[1])<=50 {HGVS=NOMEN[1]} length(NOMEN[1])>50 {HGVS=substr(NOMEN[1],1,47)"..." } {OUTPUT=gsub(","," ",$2)} {LOCATION=gsub(","," ",$3)} {IMPACT=gsub(","," ",$4)} {VAF=$5*100} {PZFlag=$6} $PZFlag == "PASS" {print "{\\small "HGVS"} & {\\small "LOCATION" "OUTPUT" "IMPACT"} & {\\small "VAF"\\%} \\\\ \\hline"}' | sed 's/_/\\_/gi' | sed 's/>/\\textgreater /gi' >> $TEXFILE
					head -n$(echo "$NB_VARIANT_LIMIT +1" | bc) $FINALTXTFILE | awk -f $SCRIPT_DIR/tsv_extract.awk -F"\t" -v COLS="NOMEN,location,outcome,snpeff_impact,VAF_median,$PZFLAG_FIELD,#CHROM,POS,REF,ALT" -v DEFAULT_NA="" | sed 's/&/ /gi' | awk -F "\t" '{split($1,NOMEN,",")} {CHROMPOSREFALT=$7":"$8$9">"$10} length(NOMEN[1])==0 {NOMEN[1]=CHROMPOSREFALT} length(NOMEN[1])<=50 {HGVS=NOMEN[1]} {HGVS=NOMEN[1]} length(NOMEN[1])>50 {HGVS=substr(NOMEN[1],1,47)"..." } {OUTPUT=$2} {LOCATION=$3} {IMPACT=$4} {VAF=$5*100} {PZFlag=$6} PZFlag == "PASS" {print "{\\small "HGVS"} & {\\small "LOCATION" "OUTPUT" "IMPACT"} & {\\small "VAF"\\%} \\\\ \\hline"}' | sed 's/_/\\_/gi' | sed 's/>/\\textgreater /gi' >> $TEXFILE

				#substr($0,7,20)

					#echo " \\\\ \hline";
				echo "
						\end{longtable}" >> $TEXFILE

			fi;

		done;

	#fi;
	#done

fi;


#######
### Gene coverage section
#########

echo "\section{Gene Coverage}" >> $TEXFILE


for aligner in `echo $ALIGNERS`
do
	for metrics in $(ls $RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.$aligner.bam.metrics/$SAMPLE*latex.txt 2>/dev/null)
	do
		if [ -s "${metrics}" ]
		then

				NB_ERROR_COVERED=$(grep "orange\|red\|yellow" ${metrics} | cut -d" " -f1 | sort -u | wc -l)
				NB_GENES=$(grep "^#" ${metrics} -vc)
				NB_FULL_COVERED=$(echo $(grep -v "orange\|red\|yellow" ${metrics} | cut -d" " -f1 | sort -u | wc -l)" -1 " | bc)

				number_col=`awk -F'&' '{print NF; exit}' $metrics`
				first=1
				while read line
				do
					if [[ "$first" == "1" ]]
					then
						# table for genes coverage withUTR
						echo "`echo $line | sed 's/#/\\\\tablefirsthead{\\\\hline \\\\multicolumn{1}{|c}{/' | sed 's/&/} \& \\\\multicolumn{1}{|p{1cm}|}{/g' | sed 's/$/} \\\\\\\\ \\\\hline}/'`" >> $TEXFILE
						echo "`echo $line | sed 's/#/\\\\tablehead{\\\\hline \\\\multicolumn{'"${number_col}"'}{|l|}{\\\\small ... next} \\\\\\\\ \\\\hline \\\\multicolumn{1}{|c}{/' | sed 's/&/} \& \\\\multicolumn{1}{|p{1cm}|}{/g' | sed 's/$/} \\\\\\\\ \\\\hline}/'`" >> $TEXFILE
						echo "\par
\bottomcaption{Genes coverage statistics on reads aligned with $aligner from $SAMPLEID from `basename "${metrics}" | sed 's/\\_/\\\_/g'` file. Genes coverage was calculated with a padding of ${NB_BASES_AROUND} bases around the coordinates.}
\begin{center}
\small
\begin{supertabular}{*{${number_col}}{|l}|}" >> $TEXFILE
					else
						echo "`echo $line | sed 's/$/ \\\\\\\\ \\\\hline/'`" >> $TEXFILE
					fi
				let "first=first+1"
				done < ${metrics}
				echo "\end{supertabular}
				\end{center}"  >> $TEXFILE

				echo "\begin{itemize}" >> $TEXFILE
				echo "\item$NB_FULL_COVERED (out of $NB_GENES) genes fully covered" >> $TEXFILE
				echo "\item$NB_ERROR_COVERED (out of $NB_GENES) genes not fully covered" >> $TEXFILE
				echo "\end{itemize}" >> $TEXFILE

		else
			echo "${metrics} doesn't exists ! No coverage genes table write in the report."
		fi
	done
done


#######
### COVERAGE PLUS section
#########

echo "\section{Coverage Annexe}" >> $TEXFILE

# Targeted Region Coverage
############################

if ((1)); then

for aligner in `echo $ALIGNERS`
do
	for metrics in $(ls $RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.$aligner.bam.metrics/$SAMPLE.$aligner.HsMetrics.per_target_coverage 2>/dev/null)
	do
		if [ -s "${metrics}" ]
		then

			echo "\begin{longtable}{|c|c|c|c|}

				\hline \multicolumn{1}{|c|}{\textbf{Region}} & \multicolumn{1}{c|}{\textbf{Min}} & \multicolumn{1}{c|}{\textbf{Mean}} & \multicolumn{1}{c|}{\textbf{Max}}  \\\\ \hline
				\endfirsthead

				\multicolumn{4}{c}%
				{{\bfseries \tablename\ \thetable{} -- continued from previous page}} \\\\
				\hline \multicolumn{1}{|c|}{\textbf{Region}} & \multicolumn{1}{c|}{\textbf{Min}} & \multicolumn{1}{c|}{\textbf{Mean}} & \multicolumn{1}{c|}{\textbf{Max}} \\\\ \hline
				\endhead

				\hline \multicolumn{4}{|r|}{{Continued on next page}} \\\\ \hline
				\endfoot

				\hline
				\endlastfoot
				" >> $TEXFILE
			#cat $COVERAGEFILESTATS_SAMTOOLS.summary.tmp.latex >> $TEXFILE

			tail -n +2 ${metrics} |  sort -k5 | awk -F "\t" '{split($7,M1,".")} {Mean=M1[1]} length($5)>50 {N=substr($5,1,47)"..." } length($5)<=50 {N=$5} {O=Ol[1]} {print "{\\small "N"} & {\\small "$11"} & {\\small "Mean"} & {\\small "$12"} \\\\ \\hline"}' | sed 's/_/\\_/gi' | sed 's/>/\\textgreater /gi' >> $TEXFILE;

			echo "
				\caption{Coverage statistics on each targeted region (primers excluded if any) on $aligner on $SAMPLEID. Statistics were determined with PICARD on the $SAMPLEID.$aligner.bam file.} \label{tab:long} \\\\
				\end{longtable}



				" >> $TEXFILE


		fi;
	done;
done;

fi;

echo "\\\\" >> $TEXFILE



#AMPLICON Coverage
###################

if ((1)); then

for aligner in `echo $ALIGNERS`
do
	for metrics in $(ls $RESULTS_FOLDER/$FLOWCELL/$SAMPLE/$SAMPLE.$aligner.bam.metrics/$SAMPLE.$aligner.amplicon_coverage 2>/dev/null)
	do
		if [ -s "${metrics}" ]
		then

			echo "\begin{longtable}{|c|c|c|c|}

				\hline \multicolumn{1}{|c|}{\textbf{Region}} & \multicolumn{1}{c|}{\textbf{Min}} & \multicolumn{1}{c|}{\textbf{Mean}} & \multicolumn{1}{c|}{\textbf{Max}}  \\\\ \hline
				\endfirsthead

				\multicolumn{4}{c}%
				{{\bfseries \tablename\ \thetable{} -- continued from previous page}} \\\\
				\hline \multicolumn{1}{|c|}{\textbf{Region}} & \multicolumn{1}{c|}{\textbf{Min}} & \multicolumn{1}{c|}{\textbf{Mean}} & \multicolumn{1}{c|}{\textbf{Max}} \\\\ \hline
				\endhead

				\hline \multicolumn{4}{|r|}{{Continued on next page}} \\\\ \hline
				\endfoot

				\hline
				\endlastfoot
				" >> $TEXFILE
			#cat $COVERAGEFILESTATS_SAMTOOLS.summary.tmp.latex >> $TEXFILE

			tail -n +2 ${metrics} |  sort -k5 | awk -F "\t" '{split($7,M1,".")} {Mean=M1[1]} length($5)>50 {N=substr($5,1,47)"..." } length($5)<=50 {N=$5} {O=Ol[1]} {print "{\\small "N"} & {\\small "$11"} & {\\small "Mean"} & {\\small "$12"} \\\\ \\hline"}' | sed 's/_/\\_/gi' | sed 's/>/\\textgreater /gi' >> $TEXFILE;

			echo "
				\caption{Coverage statistics on each amplicon (primers excluded) on $aligner on $SAMPLEID. Statistics were determined with FATBAM/PICARD on the $SAMPLEID.$aligner.bam file.} \label{tab:long} \\\\
				\end{longtable}



				" >> $TEXFILE


		fi;
	done;
done;

fi;

echo "\\\\" >> $TEXFILE
echo "\\\\" >> $TEXFILE
echo "\\vskip 0.5cm" >> $TEXFILE






#grep -Po '.*hgvs=\K.*?(?=;.*)' $FINALVCFFILE
#cat $FINALVCFFILE | awk -F "\t" '{print $8}' | tr ";" "\n"  | awk -F "=" '/hgvs/ {print $2} | awk -F ","  '{print $1}'
#grep -v ^# T_HORIZON/T_HORIZON.reports/T_HORIZON.final.vcf | awk -F "\t" '{print $1" "$2" "$4" "$5" "}' ; cat T_HORIZON/T_HORIZON.reports/T_HORIZON.final.vcf | awk -F "\t" '{print $8}' | tr ";" "\n"  | awk -F "=" '/hgvs/ {print $2}'

echo -ne "\t`echo "scale=0;($COV_CALC)/1" | bc` X\t`echo "scale=2;($PCT_TARGET_BASES_10X*100)/1" | bc`\t`echo "scale=2;($PCT_TARGET_BASES_20X*100)/1" | bc`\t`echo "scale=2;($PCT_TARGET_BASES_30X*100)/1" | bc`\t`echo "scale=2;($PCT_TARGET_BASES_40X*100)/1" | bc`\t`echo "scale=2;($PCT_TARGET_BASES_50X*100)/1" | bc`\t`echo "scale=2;($PCT_TARGET_BASES_100X*100)/1" | bc`" >> ${STATSFILE}
echo "\end{itemize}" >> $TEXFILE

echo "\end{document}" >> $TEXFILE
# compilation and creation of the report in pdf format
# because pdf latex write A LOT of messages during compilation, the option "-file-line-error" allows to do a grep command on the log to see only the ERRORS, with regexp ".*:[0-9]*:.*"
$LATEX --output-format=pdf -interaction=nonstopmode -output-directory=$REPORTDIR -file-line-error $TEXFILE >> $LOGFILE
$LATEX --output-format=pdf -interaction=nonstopmode -output-directory=$REPORTDIR -file-line-error $TEXFILE >> $LOGFILE

if [ ! -z $OUTPUT_PDFFILE ] && [ -e $PDFFILE ]; then
	cp $PDFFILE $OUTPUT_PDFFILE;
fi;


# cleaning
#rm -f $STATSFILE
#echo $REPORTDIR
rm -rf $REPORTDIR
