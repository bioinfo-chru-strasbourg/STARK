#!/bin/bash
#################################
## STARK environment
#################################

# TOOLS
########

TOOLS_LIST=""


# GZIP
export GZ="gzip"			# BIN
export UNGZ="gzip -d"		# BIN
export GZIP=""				# PARAM GZIP


# JAVA
export JAVA=$NGS_TOOLS/java/current/bin/java		# BIN
export JAVA_PATH=$NGS_TOOLS/java/current/bin		# BIN
export JAVA_VERSION=current							# VER
export JAVA_DESCRIPTION="A high-level programming language developed by Sun Microsystems"
export JAVA_REF="http://java.com"
TOOLS_LIST=$TOOLS_LIST" JAVA"


# JAVA7
export JAVA7=$NGS_TOOLS/java/1.7.0/bin/java			# BIN
export JAVA7_PATH=$NGS_TOOLS/java/1.7.0/bin			# BIN
export JAVA7_VERSION=1.7.0							# VER
export JAVA7_DESCRIPTION="A high-level programming language developed by Sun Microsystems"
export JAVA7_REF="http://java.com"
TOOLS_LIST=$TOOLS_LIST" JAVA7"


# PYTHON
export PYTHON=$NGS_TOOLS/python/current/bin/python	# BIN
export PYTHON_PATH==$NGS_TOOLS/python/current/bin	# FOLDER
export PYTHON_VERSION=current						# VER
export PYTHON_DESCRIPTION="Python is a programming language that lets you work quickly and integrate systems more efficiently"
export PYTHON_REF="http://python.com"
TOOLS_LIST=$TOOLS_LIST" PYTHON"


# PYTHON2
export PYTHON2=$NGS_TOOLS/python/2/bin/python2		# BIN
export PYTHON2_PATH==$NGS_TOOLS/python/2/bin		# FOLDER
export PYTHON2_VERSION=current						# VER
export PYTHON2_DESCRIPTION="Python is a programming language that lets you work quickly and integrate systems more efficiently"
export PYTHON2_REF="http://python.com"
TOOLS_LIST=$TOOLS_LIST" PYTHON2"


# PYTHON3
export PYTHON3=$NGS_TOOLS/python/3/bin/python3		# BIN
export PYTHON3_PATH==$NGS_TOOLS/python/3/bin		# FOLDER
export PYTHON3_VERSION=current						# VER
export PYTHON3_DESCRIPTION="Python is a programming language that lets you work quickly and integrate systems more efficiently"
export PYTHON3_REF="http://python.com"
TOOLS_LIST=$TOOLS_LIST" PYTHON3"


# BCL2FASTQ
export BCL2FASTQ=$NGS_TOOLS/bcl2fastq/current/bin/bcl2fastq		# BIN
export BCL2FASTQ_VERSION=2.20.0									# VER
export BCL2FASTQ_DESCRIPTION="BCL to FASTQ conversion"
export BCL2FASTQ_REF="http://support.illumina.com/sequencing/sequencing_software/casava.html"
TOOLS_LIST=$TOOLS_LIST" BCL2FASTQ"


# SAMTOOLS
export SAMTOOLS=$NGS_TOOLS/samtools/current/bin/samtools		# BIN
export SAMTOOLS_VERSION=1.8										# VER
export SAMTOOLS_DESCRIPTION="Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format"
export SAMTOOLS_REF="Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]. Li H A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93. Epub 2011 Sep 8. [PMID: 21903627]"
TOOLS_LIST=$TOOLS_LIST" SAMTOOLS"


# VCFUTILS
export VCFUTILS=$NGS_TOOLS/bcftools/current/bin/vcfutils.pl		# BIN-SCRIPT
export VCFUTILS_VERSION=1.8										# VER
export VCFUTILS_DESCRIPTION="fix a compatibility issue with the new bcftools"
export VCFUTILS_REF="unknown"
TOOLS_LIST=$TOOLS_LIST" VCFUTILS"


# HTSlib
export HTSLIB_DESCRIPTION="A C library for reading/writing high-throughput sequencing data"
export HTSLIB_REF="http://www.htslib.org/"


# TABIX
export TABIX=$NGS_TOOLS/htslib/current/bin/tabix			# BIN
export TABIX_PATH=$(dirname $TABIX)							# BIN
export TABIX_VERSION=1.8									# VER
export TABIX_DESCRIPTION="Indexing VCF files"
export TABIX_REF=$HTSLIB_REF
TOOLS_LIST=$TOOLS_LIST" TABIX"


# BGZIP
export BGZIP=$NGS_TOOLS/htslib/current/bin/bgzip			# BIN
export BGZIP_VERSION=1.8									# VER
export BGZIP_DESCRIPTION="Compressing VCF files"
export BGZIP_REF=$HTSLIB_REF
TOOLS_LIST=$TOOLS_LIST" BGZIP"


# BCFTOOLS
export BCFTOOLS=$NGS_TOOLS/bcftools/current/bin/bcftools		# BIN
export BCFTOOLS_VERSION=1.8										# VER
export BCFTOOLS_DESCRIPTION="Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP and short indel sequence variants"
export BCFTOOLS_REF=$HTSLIB_REF
TOOLS_LIST=$TOOLS_LIST" BCFTOOLS"


# PICARD
export PICARD=$NGS_TOOLS/picard/current/bin/picard.jar		# BIN
export PICARD_VERSION=2.18.5								# VER
export PICARD_DESCRIPTION="Java command line tools for manipulating high-throughput sequencing data (HTS) data and formats"
export PICARDLIB=$NGS_TOOLS/picard/2.18.5/bin				# DIR
export PICARD_REF="http://broadinstitute.github.io/picard/"
TOOLS_LIST=$TOOLS_LIST" PICARD"


# IGV
export IGV=$NGS_TOOLS/igv/current				# BIN-JAR
export IGV_VERSION=2.4.10						# VER
export IGV_DESCRIPTION="high-performance visualization tool for interactive exploration of large, integrated genomic datasets"
export IGV_REF="James T. Robinson, Helga Thorvaldsdóttir, Wendy Winckler, Mitchell Guttman, Eric S. Lander, Gad Getz, Jill P. Mesirov. Integrative Genomics Viewer. Nature Biotechnology 29, 2426 (2011). Helga Thorvaldsdóttir, James T. Robinson, Jill P. Mesirov. Integrative Genomics Viewer (IGV): high-performance genomics data visualization and exploration.  Briefings in Bioinformatics 14, 178-192 (2013)."
TOOLS_LIST=$TOOLS_LIST" IGV"


# IGV TOOLS
export IGVTOOLS=$NGS_TOOLS/igvtools/current/bin/lib/igvtools.jar	# BIN-JAR
export IGVTOOLS_VERSION=2.4.16										# VER
export IGVTOOLS_DESCRIPTION="provides a set of tools for pre-processing data files"
export IGVTOOLS_REF="https://www.broadinstitute.org/igv/igvtools"
TOOLS_LIST=$TOOLS_LIST" IGVTOOLS"


# GATK
export GATK=$NGS_TOOLS/gatk/current/bin/GenomeAnalysisTK.jar	# BIN-JAR
export GATK_VERSION=3.8-1-0										# VER
export GATK_DESCRIPTION="The toolkit offers a wide variety of tools, with a primary focus on variant discovery and genotyping as well as strong emphasis on data quality assurance."
export GATK_REF="The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA, 2010 GENOME RESEARCH 20:1297-303 [Article] [Pubmed]"
TOOLS_LIST=$TOOLS_LIST" GATK"


# GATK4
export GATK4=$NGS_TOOLS/gatk/4.1.4.0/bin/gatk-package-4.1.4.0-local.jar	# BIN-JAR
export GATK4_VERSION=4.1.4.0											# VER
export GATK4_DESCRIPTION="The toolkit offers a wide variety of tools, with a primary focus on variant discovery and genotyping as well as strong emphasis on data quality assurance."
export GATK4_REF="The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA, 2010 GENOME RESEARCH 20:1297-303 [Article] [Pubmed]"
TOOLS_LIST=$TOOLS_LIST" GATK4"


# R
export R=$NGS_TOOLS/anaconda/miniconda2/bin/R	# BIN-R
export R_VERSION=3.2.2							# VER
export R_DESCRIPTION="R: A Language and Environment for Statistical Computing."
export R_REF="R Development Core Team (2008). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. ISBN 3-900051-07-0, URL http://www.R-project.org."
TOOLS_LIST=$TOOLS_LIST" R"


# MUTECT
export MUTECT=$NGS_TOOLS/mutect/current/bin/mutect.jar		# BIN-JAR
export MUTECT_VERSION=1.1.7									# VER
export MUTECT_DESCRIPTION="MuTect is a method developed at the Broad Institute for the reliable and accurate identification of somatic point mutations in next generation sequencing data of cancer genomes."
export MUTECT_REF="Cibulskis, K. et al. Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples. Nat Biotechnology (2013).doi:10.1038/nbt.2514"
TOOLS_LIST=$TOOLS_LIST" MUTECT"


# OUTLYZER
export OUTLYZER=$NGS_TOOLS/outlyzer/current/bin/outLyzer.py		# BIN
export OUTLYZER_VERSION=2										# VER
export OUTLYZER_DESCRIPTION="outLyzer is a computer program whose purpose is to detect variations, specifically low allele frequency variation, in next generation sequencing data (tumor samples, mosaïc mutation)."
export OUTLYZER_REF="https://github.com/EtieM/outLyzer"
TOOLS_LIST=$TOOLS_LIST" OUTLYZER"


# FASTQC
export FASTQC=$NGS_TOOLS/fastqc/current/bin/fastqc		# BIN
export FASTQC_VERSION=0.11.8							# VER
export FASTQC_DESCRIPTION="A quality control tool for high throughput sequence data."
export FASTQC_REF="http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc"
TOOLS_LIST=$TOOLS_LIST" FASTQC"


# FASTP
export FASTP=$NGS_TOOLS/fastp/current/bin/fastp		# BIN
export FASTP_VERSION=0.20.0							# VER
export FASTP_DESCRIPTION="A tool designed to provide fast all-in-one preprocessing for FastQ files."
export FASTP_REF="https://github.com/OpenGene/fastp"
TOOLS_LIST=$TOOLS_LIST" FASTP"


# UMI TOOLS
export UMITOOLS=$NGS_TOOLS/umi_tools/current/bin/umi_tools		# BIN
export FASTQC_VERSION=1.0.0										# VER
export FASTQC_DESCRIPTION="UMI-tools contains tools for dealing with Unique Molecular Identifiers (UMIs)/Random Molecular Tags (RMTs) and single cell RNA-Seq cell barcodes."
export FASTQC_REF="https://github.com/CGATOxford/UMI-tools"
TOOLS_LIST=$TOOLS_LIST" UMITOOLS"


# BWA
export BWA=$NGS_TOOLS/bwa/current/bin/bwa			# BIN
export BWA_VERSION=0.7.17							# VER
export BWA_DESCRIPTION="package for mapping low-divergent sequences against a large reference genome"
export BWA_REF="Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]"
TOOLS_LIST=$TOOLS_LIST" BWA"


# BOWTIE_0_7_15
export BOWTIE=$NGS_TOOLS/bowtie2/current/bin/bowtie2			# BIN
export BOWTIE_VERSION=2.3.4.3									# VER
export BOWTIE_DESCRIPTION="Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences."
export BOWTIE_REF="Langmead B1, Trapnell C, Pop M, Salzberg SL. (2009) Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol. 2009;10(3):R25. doi: 10.1186/gb-2009-10-3-r25. Epub 2009 Mar 4. [PMID: 19261174]"
TOOLS_LIST=$TOOLS_LIST" BOWTIE"


# BEDTOOLS
export BEDTOOLS=$NGS_TOOLS/bedtools/current/bin/bedtools	# BIN
export BEDTOOLS_DIR=$NGS_TOOLS/bedtools/current/bin			# DIR
export BEDTOOLS_VERSION=2.29.0								# VER
export BEDTOOLS_DESCRIPTION="a powerful toolset for genome arithmetic"
export BEDTOOLS_REF="http://bedtools.readthedocs.org/"
TOOLS_LIST=$TOOLS_LIST" BEDTOOLS"


# ANNOVAR
export ANNOVAR=$NGS_TOOLS/annovar/current/bin			# DIR
export ANNOVAR_VERSION=2019Oct24						# VER
export ANNOVAR_DESCRIPTION="an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes"
export ANNOVAR_REF="Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research, 38:e164, 2010"
TOOLS_LIST=$TOOLS_LIST" ANNOVAR"


# VARSCAN
export VARSCAN=$NGS_TOOLS/varscan/current/bin/VarScan.jar		# BIN-JAR
export VARSCAN_VERSION=2.4.4									# VER
export VARSCAN_DESCRIPTION="variant detection in massively parallel sequencing data"
export VARSCAN_REF="VarScan 2: Koboldt, D., Zhang, Q., Larson, D., Shen, D., McLellan, M., Lin, L., Miller, C., Mardis, E., Ding, L., & Wilson, R. (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing Genome Research DOI: 10.1101/gr.129684.111 "
TOOLS_LIST=$TOOLS_LIST" VARSCAN"


# ITDSEEK
export ITDSEEK=$NGS_TOOLS/itdseek/current/bin/itdseek.sh	# BIN-JAR
export ITDSEEK_VERSION=1.2-2								# VER
export ITDSEEK_DESCRIPTION="FLT3 ITD detection algorithm"
export ITDSEEK_REF="Chun Hang Au, Anna Wa, Dona N. Ho, Tsun Leung Chan and Edmond S. K. Ma. Clinical evaluation of panel testing by next-generation sequencing (NGS) for gene mutations in myeloid neoplasms. Diagn Pathol. 2016 Jan 22;11:11. doi: 10.1186/s13000-016-0456-8."
TOOLS_LIST=$TOOLS_LIST" ITDSEEK"


# TRIMMOMATIC
export TRIMMOMATIC=$NGS_TOOLS/trimmomatic/current/bin/trimmomatic.jar	# BIN-JAR
export TRIMMOMATIC_VERSION=0.35											# VER
export TRIMMOMATIC_DESCRIPTION="A flexible read trimming tool for Illumina NGS data"
export TRIMMOMATIC_REF="Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.     A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.,  (Austin). 2012 Apr-Jun;6(2):80-92."
TOOLS_LIST=$TOOLS_LIST" TRIMMOMATIC"


# SNPEFF
export SNPEFF_FOLDER=$NGS_TOOLS/snpeff/current/bin		# FOLDER
export SNPEFF=$SNPEFF_FOLDER/snpEff.jar					# BIN-JAR
export SNPEFF_VERSION=4.3t								# VER
export SNPEFF_DESCRIPTION="Genetic variant annotation and effect prediction toolbox. It annotates and predicts the effects of variants on genes (such as amino acid changes)"
export SNPEFF_REF="A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3., Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (Austin). 2012 Apr-Jun;6(2):80-92 "
TOOLS_LIST=$TOOLS_LIST" SNPEFF"


# STARK
export STARK=$NGS_TOOLS/stark/$ENV_RELEASE/bin			# DIR
if [ ! -d $STARK ]; then
	export STARK=$STARK_FOLDER_BIN;
fi;
if [ ! -d $STARK ]; then
	export STARK=$STARK_FOLDER;
fi;
export STARK_VERSION=$ENV_RELEASE						# VER
export STARK_QUEUED=STARKQueued.txt
export STARK_RUNNING=STARKRunning.txt
export STARK_COMPLETE=STARKComplete.txt
export STARK_DESCRIPTION="Stellar Tools for varaints Analysis and RanKing"
export STARK_REF="inhouse"
TOOLS_LIST=$TOOLS_LIST" STARK"


# CAP
export CAP_FOLDER=$NGS_TOOLS/cap/current/bin				# DIR
export CAP=$CAP_FOLDER/CAP									# BIN
export CAP_VERSION=0.9.12									# VER
export CAP_SOFTCLIPTOQ0=$CAP_FOLDER/CAP.SoftClipToQ0.pl		# BIN-SCRIPT
export CAP_SOFTCLIPTOQ0_VERSION=$CAP_VERSION				# VER
export CAP_ManifestToBED=$CAP_FOLDER/CAP.ManifestToBED.pl	# BIN-SCRIPT
export CAP_ManifestToBED_VERSION=$CAP_VERSION				# VER
export CAP_DESCRIPTION="Clipping Amplicons Primers"
export CAP_REF="inhouse"
TOOLS_LIST=$TOOLS_LIST" CAP"


# HOWARD
export HOWARD_FOLDER=$NGS_TOOLS/howard/current			# DIR
export HOWARD_FOLDER_BIN=$HOWARD_FOLDER/bin				# DIR
#export HOWARD_FOLDER_CONFIG=$HOWARD_FOLDER/config		# DIR
export HOWARD_FOLDER_DOCS=$HOWARD_FOLDER/docs			# DIR
export HOWARD_VERSION=0.9.15.1b							# VER
export HOWARD=$HOWARD_FOLDER_BIN/HOWARD					# BIN-SCRIPT
export HOWARD_RELEASE=$HOWARD_VERSION
export HOWARDDIR=$HOWARD_FOLDER_BIN
export HOWARD_DESCRIPTION="Highly Open and Valuable tool for Variant Annotation & Ranking"
export HOWARD_REF="inhouse"
TOOLS_LIST=$TOOLS_LIST" HOWARD"


# SCRIPTS
export STARK_BED_NORMALIZATION=$STARK_FOLDER_BIN/bed_normalization.awk
export STARK_RUN_METRICS=$STARK_FOLDER_BIN/runmetrics.py

# PERL5LIB
ENV_PERLLIB=$NGS_FOLDER/tools/perl/lib


# PATH
########


# ADD JAVA
export PATH=$PATH:$TABIX_PATH:$JAVA_PATH

# ADD PERL
if (($(grep -c " 6." /etc/centos-release 2>/dev/null))); then
	export PERL5LIB=$BCL2FASTQ_PERLLIB:$ENV_PERLLIB
else
	export PERL5LIB=
fi;
