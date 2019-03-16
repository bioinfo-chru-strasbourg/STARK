#!/bin/bash
#################################
## STARK environment
#################################

# TOOLS
########

TOOLS_LIST=""


# GZIP
export GZ="gzip"			# BIN
export UNGZ="gzip -d"			# BIN
export GZIP=""				# PARAM GZIP


# JAVA
export JAVA=$NGS_TOOLS/java/current/bin/java						# BIN
export JAVA_PATH=$NGS_TOOLS/java/current/bin						# BIN
export JAVA_VERSION=1.8.0_65					# VER
#export JAVA6=$NGS_TOOLS/java/1.6.0_45/bin/java			# BIN
#export JAVA6_VERSION=1.6.0_45					# VER
export JAVA_DESCRIPTION="A high-level programming language developed by Sun Microsystems"
export JAVA_REF="http://java.com"
TOOLS_LIST=$TOOLS_LIST" JAVA"


# CASAVA
#export CASAVA=$NGS_TOOLS/casava/1.8.2/bin			# DIR
#export CASAVA_PERLLIB=$NGS_TOOLS/casava/1.8.2/perl/lib			# LIB
#export CASAVA_BCLTOFASTQ=$CASAVA/configureBclToFastq.pl		# BIN
#export CASAVA_VERSION=1.8.2					# VER
#export CASAVA_DESCRIPTION="Processes sequencing reads provided by RTA or OLB, especially demultiplexing"
#export CASAVA_REF="http://support.illumina.com/sequencing/sequencing_software/casava.html"


# BCL2FASTQ
#export BCL2FASTQ=$CASAVA_BCLTOFASTQ				# BIN
#export BCL2FASTQ=$NGS_TOOLS/bcl2fastq/1.8.4compfromint/bin	# BIN
#export BCL2FASTQ_PERLLIB=$NGS_TOOLS/bcl2fastq/1.8.4compfromint/lib/bcl2fastq-1.8.4/perl			# LIB
export BCL2FASTQ=$NGS_TOOLS/bcl2fastq/current/bin		# BIN
export BCL2FASTQ_VERSION=2.20.0				# VER
#export BCL2FASTQ_PERLLIB=$NGS_TOOLS/bcl2fastq/1.8.4/lib/bcl2fastq-1.8.4/perl			# LIB
export BCL2FASTQ_BCLTOFASTQ=$BCL2FASTQ/bcl2fastq		# BIN
#export BCL2FASTQ_VERSION=1.8.2					# VER
export BCL2FASTQ_DESCRIPTION="BCL to FASTQ conversion"
export BCL2FASTQ_REF="http://support.illumina.com/sequencing/sequencing_software/casava.html"
TOOLS_LIST=$TOOLS_LIST" BCL2FASTQ"


# SAMTOOLS
#export SAMTOOLS=$NGS_TOOLS/samtools/1.3.1/bin/samtools		# BIN
#export SAMTOOLS_VERSION=1.3.1					# VER
export SAMTOOLS=$NGS_TOOLS/samtools/current/bin/samtools		# BIN
export SAMTOOLS_VERSION=1.8					# VER
export SAMTOOLS_DESCRIPTION="Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format"
export SAMTOOLS_REF="Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]. Li H A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93. Epub 2011 Sep 8. [PMID: 21903627]"
TOOLS_LIST=$TOOLS_LIST" SAMTOOLS"


# VCFUTILS
export VCFUTILS=$NGS_TOOLS/bcftools/current/bin/vcfutils.pl		# BIN-SCRIPT
export VCFUTILS_VERSION=1.8					# VER
export VCFUTILS_DESCRIPTION="fix a compatibility issue with the new bcftools"
export VCFUTILS_REF="unknown"
TOOLS_LIST=$TOOLS_LIST" VCFUTILS"


# HTSlib - A C library for reading/writing high-throughput sequencing data
export HTSLIB_DESCRIPTION="A C library for reading/writing high-throughput sequencing data"
export HTSLIB_REF="http://www.htslib.org/"


# TABIX
#export TABIX=$NGS_TOOLS/htslib/1.3.2/bin/tabix			# BIN
#export TABIX_PATH=$(dirname $TABIX)				# BIN
export TABIX=$NGS_TOOLS/htslib/current/bin/tabix			# BIN
export TABIX_PATH=$(dirname $TABIX)				# BIN
export TABIX_VERSION=1.8					# VER
export TABIX_DESCRIPTION="Indexing VCF files"
export TABIX_REF=$HTSLIB_REF
TOOLS_LIST=$TOOLS_LIST" TABIX"

# BGZIP
#export BGZIP=$NGS_TOOLS/htslib/1.3.2/bin/bgzip			# BIN
#export BGZIP_VERSION=1.3.2					# VER
export BGZIP=$NGS_TOOLS/htslib/current/bin/bgzip			# BIN
export BGZIP_VERSION=1.8					# VER
export BGZIP_DESCRIPTION="Compressing VCF files"
export BGZIP_REF=$HTSLIB_REF
TOOLS_LIST=$TOOLS_LIST" BGZIP"


# BCFTOOLS
#export BCFTOOLS=$NGS_TOOLS/bcftools/1.3.1/bin/bcftools		# BIN
#export BCFTOOLS_VERSION=1.3.1					# VER
export BCFTOOLS=$NGS_TOOLS/bcftools/current/bin/bcftools		# BIN
export BCFTOOLS_VERSION=1.8					# VER
export BCFTOOLS_DESCRIPTION="Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP and short indel sequence variants"
export BCFTOOLS_REF=$HTSLIB_REF
TOOLS_LIST=$TOOLS_LIST" BCFTOOLS"


# VCFTOOLS
#export VCFTOOLS=$NGS_TOOLS/vcftools/0.1.12b/bin/		# DIR
#export VCFTOOLSLIB=$NGS_TOOLS/vcftools/0.1.12b/perl/		# DIR
export VCFTOOLS=$NGS_TOOLS/vcftools/current/bin			# DIR
export VCFTOOLSLIB=$NGS_TOOLS/vcftools/current/perl		# DIR
export VCFTOOLS_VERSION=0.1.14
export VCFTOOLS_DESCRIPTION="provide easily accessible methods for working with complex genetic variation data in the form of VCF files"
export VCFTOOLS_REF="http://vcftools.sourceforge.net/"
TOOLS_LIST=$TOOLS_LIST" VCFTOOLS"


# PICARD
export PICARD=$NGS_TOOLS/picard/current/bin/picard.jar		# BIN
export PICARD_VERSION=2.18.5					# VER
export PICARD_DESCRIPTION="Java command line tools for manipulating high-throughput sequencing data (HTS) data and formats"
export PICARDLIB=$NGS_TOOLS/picard/2.18.5/bin			# DIR
#export PICARDDEV=$NGS_TOOLS/picard/dev/bin			# DIR
#export PICARDLIB_VERSION=1.95					# VER
#export PICARDLIB_DESCRIPTION="Java command line tools for manipulating high-throughput sequencing data (HTS) data and formats"
export PICARD_REF="http://broadinstitute.github.io/picard/"
TOOLS_LIST=$TOOLS_LIST" PICARD"


# IGV
export IGV=$NGS_TOOLS/igv/current				# BIN-JAR
export IGV_VERSION=2.4.10					# VER
export IGV_DESCRIPTION="high-performance visualization tool for interactive exploration of large, integrated genomic datasets"
export IGV_REF="James T. Robinson, Helga Thorvaldsdóttir, Wendy Winckler, Mitchell Guttman, Eric S. Lander, Gad Getz, Jill P. Mesirov. Integrative Genomics Viewer. Nature Biotechnology 29, 2426 (2011). Helga Thorvaldsdóttir, James T. Robinson, Jill P. Mesirov. Integrative Genomics Viewer (IGV): high-performance genomics data visualization and exploration.  Briefings in Bioinformatics 14, 178-192 (2013)."
TOOLS_LIST=$TOOLS_LIST" IGV"


# IGV TOOLS
export IGVTOOLS=$NGS_TOOLS/igvtools/current/bin/lib/igvtools.jar	# BIN-JAR
export IGVTOOLS_VERSION=2.4.16					# VER
export IGVTOOLS_DESCRIPTION="provides a set of tools for pre-processing data files"
export IGVTOOLS_REF="https://www.broadinstitute.org/igv/igvtools"
TOOLS_LIST=$TOOLS_LIST" IGVTOOLS"


# GATK
#export GATK=$NGS_TOOLS/gatk/3.5-0/bin/GenomeAnalysisTK.jar	# BIN-JAR
#export GATK_VERSION=3.5-0					# VER
#export GATK_DESCRIPTION="The toolkit offers a wide variety of tools, with a primary focus on variant discovery and genotyping as well as strong emphasis on data quality assurance."
#export GATK_REF="The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA, 2010 GENOME RESEARCH 20:1297-303 [Article] [Pubmed]"
#export GATK=$NGS_TOOLS/gatk/3.8-1_2017-12-04-1/bin/GenomeAnalysisTK.jar	# BIN-JAR
#export GATK_VERSION=3.8-1_2017-12-04-1					# VER
export GATK=$NGS_TOOLS/gatk/current/bin/GenomeAnalysisTK.jar	# BIN-JAR
export GATK_VERSION=3.8-1-0					# VER
export GATK_DESCRIPTION="The toolkit offers a wide variety of tools, with a primary focus on variant discovery and genotyping as well as strong emphasis on data quality assurance."
export GATK_REF="The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA, 2010 GENOME RESEARCH 20:1297-303 [Article] [Pubmed]"
TOOLS_LIST=$TOOLS_LIST" GATK"


# CANOE
export CANOE=$NGS_TOOLS/canoe/current # BIN-R
export CANOE_VERSION=1.0.0					# VER
export CANOE_DESCRIPTION="An algorithm for the detection of rare copy number variants from exome sequencing data. CANOES models read counts using a negative binomial distribution and estimates variance of the read counts using a regression-based approach based on selected reference samples in a given dataset."
export CANOE_REF="Backenroth D, Homsy J, Murillo LR, Glessner J, Lin E, Brueckner M, Lifton R, Goldmuntz E, Chung WK, Shen Y, (2014) CANOES: Detecting rare copy number variants from whole exome sequencing data, Nucleic Acids Research,  doi: 10.1093/nar/gku345 PMID: 24771342"
TOOLS_LIST=$TOOLS_LIST" CANOE"


# R
export R=$NGS_TOOLS/anaconda/miniconda2/bin/R	# BIN-R
export R_VERSION=3.2.2					# VER
export R_DESCRIPTION="R: A Language and Environment for Statistical Computing."
export R_REF="R Development Core Team (2008). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. ISBN 3-900051-07-0, URL http://www.R-project.org."
TOOLS_LIST=$TOOLS_LIST" R"


# MUTECT
export MUTECT=$NGS_TOOLS/mutect/current/bin/muTect-1.1.4.jar	# BIN-JAR
export MUTECT_VERSION=1.1.4					# VER
export MUTECT_DESCRIPTION="MuTect is a method developed at the Broad Institute for the reliable and accurate identification of somatic point mutations in next generation sequencing data of cancer genomes."
export MUTECT_REF="Cibulskis, K. et al. Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples. Nat Biotechnology (2013).doi:10.1038/nbt.2514"
TOOLS_LIST=$TOOLS_LIST" MUTECT"


# FASTQC
export FASTQC=$NGS_TOOLS/fastqc/current/bin/fastqc		# BIN
export FASTQC_VERSION=0.11.8					# VER
export FASTQC_DESCRIPTION="A quality control tool for high throughput sequence data."
export FASTQC_REF="http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc"
TOOLS_LIST=$TOOLS_LIST" FASTQC"


# # BWA
#export BWA_OLD=$NGS_TOOLS/bwa/0.7.12/bin/bwa			# BIN
#export BWA_VERSION=0.7.12					# VER
#export BWA_DESCRIPTION="package for mapping low-divergent sequences against a large reference genome"
#export BWA_REF="Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]"


# BWA_0_7_15
export BWA=$NGS_TOOLS/bwa/current/bin/bwa			# BIN
export BWA_VERSION=0.7.17					# VER
export BWA_DESCRIPTION="package for mapping low-divergent sequences against a large reference genome"
export BWA_REF="Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]"
TOOLS_LIST=$TOOLS_LIST" BWA"


# BOWTIE_0_7_15
export BOWTIE=$NGS_TOOLS/bowtie2/current/bin/bowtie2			# BIN
export BOWTIE_VERSION=2.3.4.3					# VER
export BOWTIE_DESCRIPTION="Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences."
export BOWTIE_REF="Langmead B1, Trapnell C, Pop M, Salzberg SL. (2009) Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol. 2009;10(3):R25. doi: 10.1186/gb-2009-10-3-r25. Epub 2009 Mar 4. [PMID: 19261174]"
TOOLS_LIST=$TOOLS_LIST" BOWTIE"


# BEDTOOLS
export BEDTOOLS=$NGS_TOOLS/bedtools/current/bin			# DIR
export BEDTOOLS2=$NGS_TOOLS/bedtools/2.27.1/bin			# DIR
export BEDTOOLS_VERSION=2.27.1					# VER
export BEDTOOLS2_VERSION=2.27.1					# VER
export BEDTOOLS_DESCRIPTION="a powerful toolset for genome arithmetic"
export BEDTOOLS2_DESCRIPTION="a powerful toolset for genome arithmetic"
export BEDTOOLS_REF="http://bedtools.readthedocs.org/"
export BEDTOOLS2_REF="http://bedtools.readthedocs.org/"
#export BEDTOOLS2="/home1/DIAG/TOOLS/tools/bedtools/2.25.0/bin/"
TOOLS_LIST=$TOOLS_LIST" BEDTOOLS BEDTOOLS2"

# ANNOVAR
export ANNOVAR=$NGS_TOOLS/annovar/current/bin			# DIR
export ANNOVAR_VERSION=2018Apr16				# VER
export ANNOVAR_DESCRIPTION="an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes"
export ANNOVAR_REF="Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research, 38:e164, 2010"
TOOLS_LIST=$TOOLS_LIST" ANNOVAR"


# VARSCAN
export VARSCAN=$NGS_TOOLS/varscan/current/bin/VarScan.jar		# BIN-JAR
export VARSCAN_VERSION=2.4.3					# VER
export VARSCAN_DESCRIPTION="variant detection in massively parallel sequencing data"
export VARSCAN_REF="VarScan 2: Koboldt, D., Zhang, Q., Larson, D., Shen, D., McLellan, M., Lin, L., Miller, C., Mardis, E., Ding, L., & Wilson, R. (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing Genome Research DOI: 10.1101/gr.129684.111 "
TOOLS_LIST=$TOOLS_LIST" VARSCAN"


# ITDSEEK
export ITDSEEK=$NGS_TOOLS/itdseek/current/bin/itdseek.sh	# BIN-JAR
export ITDSEEK_VERSION=1.2-2				# VER
export ITDSEEK_DESCRIPTION="FLT3 ITD detection algorithm"
export ITDSEEK_REF="Chun Hang Au, Anna Wa, Dona N. Ho, Tsun Leung Chan and Edmond S. K. Ma. Clinical evaluation of panel testing by next-generation sequencing (NGS) for gene mutations in myeloid neoplasms. Diagn Pathol. 2016 Jan 22;11:11. doi: 10.1186/s13000-016-0456-8."
TOOLS_LIST=$TOOLS_LIST" ITDSEEK"


# TRIMMOMATIC
export TRIMMOMATIC=$NGS_TOOLS/trimmomatic/current/bin/trimmomatic.jar	# BIN-JAR
export TRIMMOMATIC_VERSION=0.35						# VER
export TRIMMOMATIC_DESCRIPTION="A flexible read trimming tool for Illumina NGS data"
export TRIMMOMATIC_REF="Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.     A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.,  (Austin). 2012 Apr-Jun;6(2):80-92."
TOOLS_LIST=$TOOLS_LIST" TRIMMOMATIC"


# SNPEFF
export SNPEFF_FOLDER=$NGS_TOOLS/snpeff/current/bin		# FOLDER
export SNPEFF=$SNPEFF_FOLDER/snpEff.jar			# BIN-JAR
export SNPEFF_VERSION=4.3t				# VER
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
export STARK_VERSION=$ENV_RELEASE				# VER
export STARK_QUEUED=STARKQueued.txt
export STARK_RUNNING=STARKRunning.txt
export STARK_COMPLETE=STARKComplete.txt
export STARK_DESCRIPTION="Stellar Tools for varaints Analysis and RanKing"
export STARK_REF="inhouse"
TOOLS_LIST=$TOOLS_LIST" STARK"


# FATBAM
#export FATBAM=$NGS_TOOLS/fatbam/0.9.6d/bin			# DIR
export FATBAM=$NGS_TOOLS/fatbam/current/bin			# DIR
#export FATBAM_VERSION=0.9.6d					# VER
export FATBAM_VERSION=0.9.9b					# VER
export FATBAM_CLIPPING=$FATBAM/FATBAM.clipping.sh		# BIN-SCRIPT
export FATBAM_COVERAGE=$FATBAM/FATBAM.coverage.sh		# BIN-SCRIPT
export FATBAM_ManifestToBED=$FATBAM/FATBAM.ManifestToBED.pl
export CLIPPING=1						# Process clipping by default
export FATBAM_DESCRIPTION="Clipping Amplicons' Primers"
export FATBAM_REF="inhouse"
TOOLS_LIST=$TOOLS_LIST" FATBAM"


# HOWARD
export HOWARD_FOLDER_BIN=$NGS_TOOLS/howard/current/bin			# DIR
export HOWARD_FOLDER_CONFIG=$NGS_TOOLS/howard/current/config	# DIR
export HOWARD_VERSION=0.9.14b				# VER
export HOWARD=$HOWARD_FOLDER_BIN/HOWARD.sh			# BIN-SCRIPT
export HOWARD_RELEASE=$HOWARD_VERSION
export HOWARDDIR=HOWARD_FOLDER_BIN
export HOWARD_DESCRIPTION="Highly Open and Valuable tool for Variant Annotation & Ranking"
export HOWARD_REF="inhouse"
# Configuration
TOOLS_LIST=$TOOLS_LIST" HOWARD"


# VaRank
export VARANK=$NGS_TOOLS/varank/VaRank_1.4.2			# DIR
export VARANK_VERSION=1.4.2			# VER
export VARANK_DESCRIPTION="VaRank is a program for genetic Variant Ranking from NGS data"
export VARANK_REF="Geoffroy V.*, Pizot C.*, Redin C., Piton A., Vasli N., Stoetzel C., Blavier A., Laporte J. and Muller J. VaRank: a simple and powerful tool for ranking genetic variants. PeerJ. 2015."
TOOLS_LIST=$TOOLS_LIST" VARANK"


# PEPPER # PrEtty PiPEline Repository
export PEPPER=$NGS_TOOLS/pepper/current/bin			# DIR
export PEPPER_VERSION=0.9d					# VER
export PEPPER_INTEGRATION=$PEPPER/DBintegration.sh		# BIN-SCRIPT
export PEPPER_EXPORT=$PEPPER/DBexport.pl			# BIN-SCRIPT
export PEPPER_CONFIG=$PEPPER/config.ini				# INI
export PEPPER_INTEGRATION_RELEASE=$PEPPER_VERSION
export PEPPER_EXPORT_RELEASE=$PEPPER_VERSION
export PEPPER_DESCRIPTION="PrEtty PiPEline Repository - dedicated to store VCF/Variants and metadata in a database"
export PEPPER_REF="inhouse"
TOOLS_LIST=$TOOLS_LIST" PEPPER"


# ALAMUT
export ALAMUT=$NGS_TOOLS/alamut_batch/alamut-batch-standalone-1.9
export ALAMUT_VERSION=1.9
export ALAMUT_DESCRIPTION="unknown"
export ALAMUT_REF="unknown"
TOOLS_LIST=$TOOLS_LIST" ALAMUT"


# AnnotSV
export ANNOTSV=$NGS_TOOLS/AnnotSV/AnnotSV_1.0
export ANNOTSV_VERSION=1.0
export ANNOTSV_DESCRIPTION="unknown"
export ANNOTSV_REF="unknown"
TOOLS_LIST=$TOOLS_LIST" ANNOTSV"


# AnnotCNV
export ANNOTCNV=$NGS_TOOLS/AnnotSV/AnnotCNV_1.0
#export ANNOTCNV=$NGS_TOOLS/AnnotSV/AnnotSV_1.0
export ANNOTCNV_VERSION=1.0
export ANNOTCNV_DESCRIPTION="unknown"
export ANNOTCNV_REF="unknown"
TOOLS_LIST=$TOOLS_LIST" ANNOTCNV"


# LATEX
#export LATEX=$NGS_TOOLS/texlive/current/bin/x86_64-linux/latex	# BIN
export LATEX=$NGS_TOOLS/texlive/current/bin/latex	# BIN
export LATEX_VERSION=2018						# VER
export LATEX_DESCRIPTION="unknown"
export LATEX_REF="unknown"
export LATEX_DIR=$NGS_TOOLS/texlive3/20161110/bin/x86_64-linux		# DIR
export PATH=$LATEX_DIR:$PATH
TOOLS_LIST=$TOOLS_LIST" LATEX"


# COULSON
export COULSON=$NGS_TOOLS/coulson/current/bin		# DIR
export COULSON_VERSION=0.9.1b
export COULSON_DESCRIPTION="Command Organiser Using LiSt Of LiNes"
export COULSON_REF="in-house"
TOOLS_LIST=$TOOLS_LIST" COULSON"



# DAEMON
export DAEMON=$NGS_TOOLS/daemon/current/bin			# DIR
export DAEMON_VERSION=0.9b
export DAEMON_DESCRIPTION="daemon tool"
export DAEMON_REF="in-house"
TOOLS_LIST=$TOOLS_LIST" DAEMON"


# PERL5LIB
ENV_PERLLIB=$NGS_FOLDER/tools/perl/lib


# PATH
########


# ADD JAVA
export PATH=$PATH:$TABIX_PATH:$JAVA_PATH

# ADD PERL
if (($(grep -c " 6." /etc/centos-release))); then
	export PERL5LIB=$VCFTOOLSLIB:$BCL2FASTQ_PERLLIB:$ENV_PERLLIB
else
	export PERL5LIB=$VCFTOOLSLIB
fi;
