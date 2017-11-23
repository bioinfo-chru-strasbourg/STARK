#################################
## STARK environment
#################################

# TOOLS
########

# GZIP
export GZIP=gzip						# BIN
export GUNZIP=gunzip						# BIN

# JAVA
#export JAVA=$NGS_TOOLS/java/1.8.0_65/bin/java						# BIN
export JAVA=java						# BIN
export JAVA_PATH=$NGS_TOOLS/java/1.8.0_65/bin						# BIN
export JAVA_VERSION=1.8.0_65					# VER
export JAVA6=$NGS_TOOLS/java/1.6.0_45/bin/java			# BIN
export JAVA6_VERSION=1.6.0_45					# VER
export JAVA_DESCRIPTION="A high-level programming language developed by Sun Microsystems"
export JAVA_REF="http://java.com"

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
export BCL2FASTQ=$NGS_TOOLS/bcl2fastq/2.17.1.14/bin		# BIN
export BCL2FASTQ_VERSION=2.17.1.14				# VER
#export BCL2FASTQ_PERLLIB=$NGS_TOOLS/bcl2fastq/1.8.4/lib/bcl2fastq-1.8.4/perl			# LIB
export BCL2FASTQ_BCLTOFASTQ=$BCL2FASTQ/bcl2fastq		# BIN
#export BCL2FASTQ_VERSION=1.8.2					# VER
export BCL2FASTQ_DESCRIPTION="BCL to FASTQ conversion"
export BCL2FASTQ_REF="http://support.illumina.com/sequencing/sequencing_software/casava.html"

# SAMTOOLS
export SAMTOOLS=$NGS_TOOLS/samtools/1.3.1/bin/samtools		# BIN
export SAMTOOLS_VERSION=1.3.1					# VER
export SAMTOOLS_DESCRIPTION="Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format"
export SAMTOOLS_REF="Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]. Li H A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93. Epub 2011 Sep 8. [PMID: 21903627]"
# VCFUTILS
export VCFUTILS=$NGS_TOOLS/vcfutils/1.2/bin/vcfutils.pl		# BIN-SCRIPT
export VCFUTILS_VERSION=1.2					# VER
export VCFUTILS_DESCRIPTION="fix a compatibility issue with the new bcftools"
# HTSlib - A C library for reading/writing high-throughput sequencing data
export HTSLIB_DESCRIPTION="A C library for reading/writing high-throughput sequencing data"
export HTSLIB_REF="http://www.htslib.org/"
# TABIX
export TABIX=$NGS_TOOLS/htslib/1.3.2/bin/tabix			# BIN
export TABIX_PATH=$(dirname $TABIX)				# BIN
export TABIX_VERSION=1.3.2					# VER
export TABIX_DESCRIPTION="Indexing VCF files"
# BGZIP
export BGZIP=$NGS_TOOLS/htslib/1.3.2/bin/bgzip			# BIN
export BGZIP_VERSION=1.3.2					# VER
export BGZIP_DESCRIPTION="Compressing VCF files"
# BCFTOOLS
export BCFTOOLS=$NGS_TOOLS/bcftools/1.3.1/bin/bcftools		# BIN
export BCFTOOLS_VERSION=1.3.1					# VER
export BCFTOOLS_DESCRIPTION="Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP and short indel sequence variants"
export BCFTOOLS_REF="http://www.htslib.org/"
# VCFTOOLS
#export VCFTOOLS=$NGS_TOOLS/vcftools/0.1.12b/bin/		# DIR
#export VCFTOOLSLIB=$NGS_TOOLS/vcftools/0.1.12b/perl/		# DIR
export VCFTOOLS=$NGS_TOOLS/vcftools/0.1.13/bin			# DIR
export VCFTOOLSLIB=$NGS_TOOLS/vcftools/0.1.13/perl		# DIR
export VCFTOOLS_VERSION=0.1.13
export VCFTOOLS_DESCRIPTION="provide easily accessible methods for working with complex genetic variation data in the form of VCF files"
export VCFTOOLS_REF="http://vcftools.sourceforge.net/"

# PICARD
export PICARD=$NGS_TOOLS/picard/2.3.0/bin/picard.jar		# BIN
export PICARD_VERSION=2.3.0					# VER
export PICARD_DESCRIPTION="Java command line tools for manipulating high-throughput sequencing data (HTS) data and formats"
export PICARDLIB=$NGS_TOOLS/picard/1.95/bin			# DIR
export PICARDDEV=$NGS_TOOLS/picard/dev/bin			# DIR
export PICARDLIB_VERSION=1.95					# VER
export PICARDLIB_DESCRIPTION="Java command line tools for manipulating high-throughput sequencing data (HTS) data and formats"
export PICARD_REF="http://broadinstitute.github.io/picard/"

# IGV
export IGV=$NGS_TOOLS/igv/2.3.52				# BIN-JAR
export IGV_VERSION=2.3.52					# VER
export IGV_DESCRIPTION="high-performance visualization tool for interactive exploration of large, integrated genomic datasets"
export IGV_REF="James T. Robinson, Helga Thorvaldsdóttir, Wendy Winckler, Mitchell Guttman, Eric S. Lander, Gad Getz, Jill P. Mesirov. Integrative Genomics Viewer. Nature Biotechnology 29, 2426 (2011). Helga Thorvaldsdóttir, James T. Robinson, Jill P. Mesirov. Integrative Genomics Viewer (IGV): high-performance genomics data visualization and exploration.  Briefings in Bioinformatics 14, 178-192 (2013)."
export IGVTOOLS=$NGS_TOOLS/igvtools/2.3.52/bin/igvtools.jar	# BIN-JAR
export IGVTOOLS_VERSION=2.3.52					# VER
export IGVTOOLS_DESCRIPTION="provides a set of tools for pre-processing data files"
export IGVTOOLS_REF="https://www.broadinstitute.org/igv/igvtools"

# GATK
export GATK=$NGS_TOOLS/gatk/3.5-0/bin/GenomeAnalysisTK.jar	# BIN-JAR
export GATK_VERSION=3.5-0					# VER
export GATK_DESCRIPTION="The toolkit offers a wide variety of tools, with a primary focus on variant discovery and genotyping as well as strong emphasis on data quality assurance."
export GATK_REF="The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA, 2010 GENOME RESEARCH 20:1297-303 [Article] [Pubmed]"

# CANOE
#export CANOE=$NGS_TOOLS/canoe/1.0.0/CANOES_without_SEX.R # BIN-R
export CANOE_DIR=$NGS_TOOLS/canoe/1.0.0				# DIR

export CANOE_VERSION=1.0.0					# VER
export CANOE_DESCRIPTION="An algorithm for the detection of rare copy number variants from exome sequencing data. CANOES models read counts using a negative binomial distribution and estimates variance of the read counts using a regression-based approach based on selected reference samples in a given dataset."
export CANOE_REF="Backenroth D, Homsy J, Murillo LR, Glessner J, Lin E, Brueckner M, Lifton R, Goldmuntz E, Chung WK, Shen Y, (2014) CANOES: Detecting rare copy number variants from whole exome sequencing data, Nucleic Acids Research,  doi: 10.1093/nar/gku345 PMID: 24771342"

# R
export R=$NGS_TOOLS/anaconda/miniconda2/bin/R	# BIN-R
export R_VERSION=3.2.2					# VER
export R_DESCRIPTION="R: A Language and Environment for Statistical Computing."
export R_REF="R Development Core Team (2008). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. ISBN 3-900051-07-0, URL http://www.R-project.org."

# MUTECT
export MUTECT=$NGS_TOOLS/mutect/1.1.4/bin/muTect-1.1.4.jar	# BIN-JAR
export MUTECT_VERSION=1.1.4					# VER
export MUTECT_DESCRIPTION="MuTect is a method developed at the Broad Institute for the reliable and accurate identification of somatic point mutations in next generation sequencing data of cancer genomes."
export MUTECT_REF="Cibulskis, K. et al. Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples. Nat Biotechnology (2013).doi:10.1038/nbt.2514"

# FASTQC
export FASTQC=$NGS_TOOLS/fastqc/0.11.3/bin/fastqc		# BIN
export FASTQC_VERSION=0.11.3					# VER
export FASTQC_DESCRIPTION="A quality control tool for high throughput sequence data."
export FASTQC_REF="http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc"

# # BWA
#export BWA_OLD=$NGS_TOOLS/bwa/0.7.12/bin/bwa			# BIN
#export BWA_VERSION=0.7.12					# VER
#export BWA_DESCRIPTION="package for mapping low-divergent sequences against a large reference genome"
#export BWA_REF="Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]"

# BWA_0_7_15
export BWA=$NGS_TOOLS/bwa/0.7.15/bwa			# BIN
export BWA_VERSION=0.7.15					# VER
export BWA_DESCRIPTION="package for mapping low-divergent sequences against a large reference genome"
export BWA_REF="Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]"

# BOWTIE_0_7_15
export BOWTIE=$NGS_TOOLS/bowtie2/2.2.8/bowtie2			# BIN
export BWA_VERSION=2.2.8					# VER
export BWA_DESCRIPTION="Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences."
export BWA_REF="Langmead B1, Trapnell C, Pop M, Salzberg SL. (2009) Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol. 2009;10(3):R25. doi: 10.1186/gb-2009-10-3-r25. Epub 2009 Mar 4. [PMID: 19261174]"

# BEDTOOLS
export BEDTOOLS=$NGS_TOOLS/bedtools/2.17.0/bin			# DIR
export BEDTOOLS2=$NGS_TOOLS/bedtools/2.25.0/bin			# DIR
export BEDTOOLS_VERSION=2.17.0					# VER
export BEDTOOLS2_VERSION=2.25.0					# VER
export BEDTOOLS_DESCRIPTION="a powerful toolset for genome arithmetic"
export BEDTOOLS_REF="http://bedtools.readthedocs.org/"
#export BEDTOOLS2="/home1/DIAG/TOOLS/tools/bedtools/2.25.0/bin/"

# ANNOVAR
export ANNOVAR=$NGS_TOOLS/annovar/2015Mar22/bin		# DIR
export ANNOVAR_VERSION=2015Mar22				# VER
export ANNOVAR_DATABASES=$NGS_FOLDER/annovar_sources		# DATA
export ANNOVAR_DESCRIPTION="an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes"
export ANNOVAR_REF="Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research, 38:e164, 2010"

# VARSCAN
export VARSCAN=$NGS_TOOLS/varscan/v2.3.7/bin/VarScan.jar	# BIN-JAR
export VARSCAN_VERSION=v2.3.7					# VER
export VARSCAN_DESCRIPTION="variant detection in massively parallel sequencing data"
export VARSCAN_REF="VarScan 2: Koboldt, D., Zhang, Q., Larson, D., Shen, D., McLellan, M., Lin, L., Miller, C., Mardis, E., Ding, L., & Wilson, R. (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing Genome Research DOI: 10.1101/gr.129684.111 "

# TRIMMOMATIC
export TRIMMOMATIC=$NGS_TOOLS/trimmomatic/0.35/bin/trimmomatic.jar	# BIN-JAR
export TRIMMOMATIC_VERSION=0.35						# VER
export TRIMMOMATIC_DESCRIPTION="A flexible read trimming tool for Illumina NGS data"
export TRIMMOMATIC_REF="Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.     A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.,  (Austin). 2012 Apr-Jun;6(2):80-92."

# SNPEFF
export SNPEFF_FOLDER=$NGS_TOOLS/snpeff/4.2/bin		# FOLDER
export SNPEFF=$SNPEFF_FOLDER/snpEff.jar			# BIN-JAR
export SNPEFF_VERSION=4.2				# VER
export SNPEFF_CONFIG=$SNPEFF_FOLDER/snpeff.config	# CONFIG # NOT USED !!! # CHANGE CONFIG FILE in SNPEFF TOOL if necessary
export SNPEFF_DATABASES=$NGS_FOLDER/snpeff_sources	# DATA # NOT USED !!! # CHANGE DATABASE location in CONFIG FILE in SNPEFF TOOL if necessary
export SNPEFF_DESCRIPTION="Genetic variant annotation and effect prediction toolbox. It annotates and predicts the effects of variants on genes (such as amino acid changes)"
export SNPEFF_REF="A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3., Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (Austin). 2012 Apr-Jun;6(2):80-92 "

# STARK
export STARK=$NGS_TOOLS/stark/$ENV_RELEASE/bin			# DIR
export STARK_VERSION=$ENV_RELEASE				# VER
export STARK_QUEUED=STARKQueued.txt
export STARK_RUNNING=STARKRunning.txt
export STARK_COMPLETE=STARKComplete.txt
export STARK_DESCRIPTION="Stellar Tools for varaints Analysis and RanKing"
export STARK_REF="inhouse"

# FATBAM
#export FATBAM=$NGS_TOOLS/fatbam/0.9.6d/bin			# DIR
export FATBAM=$NGS_TOOLS/fatbam/0.9.7b/bin			# DIR
#export FATBAM_VERSION=0.9.6d					# VER
export FATBAM_VERSION=0.9.7b					# VER
export FATBAM_CLIPPING=$FATBAM/FATBAM.clipping.sh		# BIN-SCRIPT
export FATBAM_COVERAGE=$FATBAM/FATBAM.coverage.sh		# BIN-SCRIPT
export FATBAM_ManifestToBED=$FATBAM/FATBAM.ManifestToBED.pl
export CLIPPING=1						# Process clipping by default
export FATBAM_DESCRIPTION="Clipping Amplicons' Primers"
export FATBAM_REF="inhouse"

# HOWARD
#export HOWARDDIR=$NGS_TOOLS/howard/0.9.5d/bin			# DIR
export HOWARDDIR=$NGS_TOOLS/howard/0.9.9b/bin			# DIR
#export HOWARD_VERSION=0.9.5d					# VER
export HOWARD_VERSION=0.9.9b					# VER
export HOWARD=$HOWARDDIR/HOWARD.sh				# BIN-SCRIPT
#export HOWARD_ANNOTATION=$HOWARDDIR/VCFannotation.pl		# BIN-SCRIPT
#export HOWARD_PRIORITIZATION=$HOWARDDIR/VCFprioritization.pl	# BIN-SCRIPT
#export HOWARD_TRANSLATION=$HOWARDDIR/VCFtranslation.pl		# BIN-SCRIPT
export HOWARD_CONFIG=$HOWARDDIR/config.ini			# INI
export HOWARD_RELEASE=$HOWARD_VERSION
#export HOWARD_PRIORITIZATION_RELEASE=$HOWARD_VERSION
#export HOWARD_TRANSLATION_RELEASE=$HOWARD_VERSION
export HOWARD_DESCRIPTION="Highly Open and Valuable tool for Variant Annotation & Ranking"
export HOWARD_REF="inhouse"
# Configuration
export HOWARD_CONFIG_FILTER=$HOWARDDIR/config.filter.ini

# VaRank
export VARANK=$NGS_TOOLS/varank/current			# DIR
export VARANK_VERSION=1.3.4					# VER
export VARANK_DESCRIPTION="VaRank is a program for genetic Variant Ranking from NGS data"
export VARANK_REF="Geoffroy V.*, Pizot C.*, Redin C., Piton A., Vasli N., Stoetzel C., Blavier A., Laporte J. and Muller J. VaRank: a simple and powerful tool for ranking genetic variants. PeerJ. 2015."

# PEPPER # PrEtty PiPEline Repository
export PEPPER=$NGS_TOOLS/pepper/0.9d/bin			# DIR
export PEPPER_VERSION=0.9d					# VER
export PEPPER_INTEGRATION=$PEPPER/DBintegration.sh		# BIN-SCRIPT
export PEPPER_EXPORT=$PEPPER/DBexport.pl			# BIN-SCRIPT
export PEPPER_CONFIG=$PEPPER/config.ini				# INI
export PEPPER_INTEGRATION_RELEASE=$PEPPER_VERSION
export PEPPER_EXPORT_RELEASE=$PEPPER_VERSION
export PEPPER_DESCRIPTION="PrEtty PiPEline Repository - dedicated to store VCF/Variants and metadata in a database"
export PEPPER_REF="inhouse"

# ALAMUT
export ALAMUT=$NGS_TOOLS/alamut_batch/current

# AnnotSV
export ANNOTSV=$NGS_TOOLS/AnnotSV/AnnotSV_1.0
export ANNOTCNV=$NGS_TOOLS/AnnotSV/AnnotCNV_1.0

# LATEX
export LATEX=$NGS_TOOLS/texlive/20161110/bin/x86_64-linux/latex	# BIN
#export LATEX=latex							# BIN
#export LATEX_VERSION=20161110						# VER
export LATEX_VERSION=3.5.6						# VER
export LATEX_DESCRIPTION=""
export LATEX_REF=""
#export LATEX_DIR=$NGS_TOOLS/texlive/20161110/bin/x86_64-linux		# DIR
export PATH=$LATEX_DIR:$PATH

# COULSON
export COULSON=$NGS_TOOLS/coulson/0.9.1b/bin		# DIR
export COULSON_VERSION=0.9.1b

# DAEMON
export DAEMON=$NGS_TOOLS/daemon/0.9b/bin			# DIR
export DAEMON_VERSION=0.9b

# PERL5LIB
ENV_PERLLIB=$NGS_FOLDER/tools/perl/lib




# PATH
########


# ADD JAVA
export PATH=$PATH:$TABIX_PATH:$JAVA_PATH

# ADD PERL
export PERL5LIB=$VCFTOOLSLIB:$BCL2FASTQ_PERLLIB:$ENV_PERLLIB




