#!/bin/bash
#################################
## STARK environment
#################################

# TOOLS
########

TOOLS_LIST=""
DOCKER_LIST=""


# GZIP
# export GZ="gzip"			# BIN
# export UNGZ="gzip -d"		# BIN
# export GZIP=""				# PARAM GZIP
export GZ="pigz"			# BIN
export UNGZ="pigz -d"		# BIN
export GZIP=""				# PARAM GZIP


# JAVA
export JAVA_VERSION=21 									# VER
export JAVA=$NGS_TOOLS/java/$JAVA_VERSION/bin/java		# BIN
export JAVA_PATH=$NGS_TOOLS/java/$JAVA_VERSION/bin		# BIN
export JAVA_VERSION=current								# VER
export JAVA_DESCRIPTION="A high-level programming language developed by Sun Microsystems"
export JAVA_REF="http://java.com"
TOOLS_LIST=$TOOLS_LIST" JAVA"


# JAVA8 (for GATK3)
export JAVA8_VERSION=1.8.0 										# VER
export JAVA8=$NGS_TOOLS/java/$JAVA8_VERSION/bin/java			# BIN
export JAVA8_PATH=$NGS_TOOLS/java/$JAVA8_VERSION/bin			# BIN
export JAVA8_DESCRIPTION="A high-level programming language developed by Sun Microsystems"
export JAVA8_REF="http://java.com"
TOOLS_LIST=$TOOLS_LIST" JAVA8"

# JAVA7 (for MuTect)
export JAVA7_VERSION=1.7.0 										# VER
export JAVA7=$NGS_TOOLS/java/$JAVA7_VERSION/bin/java			# BIN
export JAVA7_PATH=$NGS_TOOLS/java/$JAVA7_VERSION/bin			# BIN
export JAVA7_DESCRIPTION="A high-level programming language developed by Sun Microsystems"
export JAVA7_REF="http://java.com"
TOOLS_LIST=$TOOLS_LIST" JAVA7"


# # JAVA17
# export JAVA17_VERSION=17 							# VER
# export JAVA17=$NGS_TOOLS/java/$JAVA17_VERSION/bin/java			# BIN
# export JAVA17_PATH=$NGS_TOOLS/java/$JAVA17_VERSION/bin			# BIN
# export JAVA17_DESCRIPTION="A high-level programming language developed by Sun Microsystems"
# export JAVA17_REF="http://java.com"
# TOOLS_LIST=$TOOLS_LIST" JAVA17"


# PYTHON
export PYTHON_VERSION=3.10									# VER
export PYTHON=$NGS_TOOLS/python/$PYTHON_VERSION/bin/python	# BIN
export PYTHON_DESCRIPTION="Python is a programming language that lets you work quickly and integrate systems more efficiently"
export PYTHON_REF="http://python.com"
TOOLS_LIST=$TOOLS_LIST" PYTHON"


# PYTHON3
export PYTHON3_VERSION=3.10											# VER
export PYTHON3=$NGS_TOOLS/python/$PYTHON3_VERSION/bin/python		# BIN
export PYTHON3_DESCRIPTION="Python is a programming language that lets you work quickly and integrate systems more efficiently"
export PYTHON3_REF="http://python.com"
TOOLS_LIST=$TOOLS_LIST" PYTHON3"


# BCL2FASTQ
export BCL2FASTQ_VERSION=2.20.0									# VER
export BCL2FASTQ=$NGS_TOOLS/bcl2fastq/current/bin/bcl2fastq		# BIN
export BCL2FASTQ_DESCRIPTION="BCL to FASTQ conversion"
export BCL2FASTQ_REF="http://support.illumina.com/sequencing/sequencing_software/casava.html"
TOOLS_LIST=$TOOLS_LIST" BCL2FASTQ"


# SAMTOOLS
export SAMTOOLS_VERSION=1.23									# VER
export SAMTOOLS=$NGS_TOOLS/samtools/current/bin/samtools		# BIN 
export SAMTOOLS_DESCRIPTION="Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format"
export SAMTOOLS_REF="Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]. Li H A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93. Epub 2011 Sep 8. [PMID: 21903627]"
TOOLS_LIST=$TOOLS_LIST" SAMTOOLS"


# VCFUTILS
export VCFUTILS_VERSION=1.23									# VER
export VCFUTILS=$NGS_TOOLS/bcftools/current/bin/vcfutils.pl		# BIN-SCRIPT
export VCFUTILS_DESCRIPTION="fix a compatibility issue with the new bcftools"
export VCFUTILS_REF="unknown"
TOOLS_LIST=$TOOLS_LIST" VCFUTILS"


# HTSlib
export HTSLIB_DESCRIPTION="A C library for reading/writing high-throughput sequencing data"
export HTSLIB_REF="http://www.htslib.org/"


# TABIX
export TABIX_VERSION=1.23									# VER
export TABIX=$NGS_TOOLS/htslib/current/bin/tabix			# BIN 
export TABIX_PATH=$(dirname $TABIX)							# BIN
export TABIX_DESCRIPTION="Indexing VCF files"
export TABIX_REF=$HTSLIB_REF
TOOLS_LIST=$TOOLS_LIST" TABIX"


# BGZIP
export BGZIP_VERSION=1.23									# VER
export BGZIP=$NGS_TOOLS/htslib/current/bin/bgzip			# BIN   
export BGZIP_DESCRIPTION="Compressing VCF files"
export BGZIP_REF=$HTSLIB_REF
TOOLS_LIST=$TOOLS_LIST" BGZIP"


# BCFTOOLS
export BCFTOOLS_VERSION=1.23								# VER
export BCFTOOLS=$NGS_TOOLS/bcftools/current/bin/bcftools	# BIN $NGS_TOOLS/bcftools/current/bin/ 
export BCFTOOLS_DOCKER=dceoy/bcftools:latest				# DOCKER
export BCFTOOLS_DESCRIPTION="Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP and short indel sequence variants"
export BCFTOOLS_REF=$HTSLIB_REF
TOOLS_LIST=$TOOLS_LIST" BCFTOOLS"
DOCKER_LIST=$DOCKER_LIST" $BCFTOOLS_DOCKER"


# PICARD
export PICARD_VERSION=3.4.0									# VER
export PICARD=$NGS_TOOLS/picard/current/bin/picard.jar		# BIN
export PICARD_DESCRIPTION="Java command line tools for manipulating high-throughput sequencing data (HTS) data and formats"
export PICARDLIB=$NGS_TOOLS/picard/$PICARD_VERSION/bin				# DIR
export PICARD_REF="http://broadinstitute.github.io/picard/"
TOOLS_LIST=$TOOLS_LIST" PICARD"


# IGV TOOLS
export IGVTOOLS_VERSION=2.17.3-0								# VER
export IGVTOOLS=$NGS_TOOLS/igvtools/current/bin/igvtools		# BIN-JAR
export IGVTOOLS_DESCRIPTION="high-performance visualization tool for interactive exploration of large, integrated genomic datasets"
export IGVTOOLS_REF="James T. Robinson, Helga Thorvaldsdóttir, Wendy Winckler, Mitchell Guttman, Eric S. Lander, Gad Getz, Jill P. Mesirov. Integrative Genomics Viewer. Nature Biotechnology 29, 2426 (2011). Helga Thorvaldsdóttir, James T. Robinson, Jill P. Mesirov. Integrative Genomics Viewer (IGV): high-performance genomics data visualization and exploration.  Briefings in Bioinformatics 14, 178-192 (2013). https://www.broadinstitute.org/igv"
TOOLS_LIST=$TOOLS_LIST" IGVTOOLS"


# GATK3
export GATK3_VERSION=3.8											# VER
export GATK3=$NGS_TOOLS/gatk/current/bin/GenomeAnalysisTK.jar		# BIN-JAR
export GATK3_DESCRIPTION="The toolkit offers a wide variety of tools, with a primary focus on variant discovery and genotyping as well as strong emphasis on data quality assurance."
export GATK3_REF="The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA, 2010 GENOME RESEARCH 20:1297-303"
TOOLS_LIST=$TOOLS_LIST" GATK3"


# GATK4
export GATK4_VERSION=4.6.1.0										# VER
export GATK4=$NGS_TOOLS/gatk4/current/bin/GenomeAnalysisTK4.jar		# JAR
export GATK4_BIN=$NGS_TOOLS/gatk4/$GATK4_VERSION/bin/gatk			# BIN
export GATK4_DESCRIPTION="The toolkit offers a wide variety of tools, with a primary focus on variant discovery and genotyping as well as strong emphasis on data quality assurance."
export GATK4_REF="The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA, 2010 GENOME RESEARCH 20:1297-303"
TOOLS_LIST=$TOOLS_LIST" GATK4"


# GATK (default)
export GATK=$GATK4								# BIN-JAR
export GATK_VERSION=$GATK4_VERSION				# VER
export GATK_DESCRIPTION=$GATK4_DESCRIPTION
export GATK_REF=$GATK4_REF
TOOLS_LIST=$TOOLS_LIST" GATK"


# GENCORE
export GENCORE_VERSION=0.17.2							# VER
export GENCORE=$NGS_TOOLS/gencore/current/bin/gencore	# BIN-JAR $NGS_TOOLS/gencore/current/bin/
export GENCORE_DESCRIPTION="An efficient tool to remove sequencing duplications and eliminate sequencing errors by generating consensus reads."
export GENCORE_REF="Chen, S., Zhou, Y., Chen, Y. et al. Gencore: an efficient tool to generate consensus reads for error suppressing and duplicate removing of NGS data. BMC Bioinformatics 20, 606 (2019) doi:10.1186/s12859-019-3280-9"
TOOLS_LIST=$TOOLS_LIST" GENCORE"

# MUTECT
export MUTECT_VERSION=1.1.6									# VER
export MUTECT=$NGS_TOOLS/mutect/current/bin/muTect.jar		# BIN-JAR
export MUTECT_DESCRIPTION="MuTect is a method developed at the Broad Institute for the reliable and accurate identification of somatic point mutations in next generation sequencing data of cancer genomes."
export MUTECT_REF="Cibulskis, K. et al. Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples. Nat Biotechnology (2013).doi:10.1038/nbt.2514"
TOOLS_LIST=$TOOLS_LIST" MUTECT"


# OUTLYZER
export OUTLYZER_VERSION=3.2															# VER
export OUTLYZER=$NGS_TOOLS/outlyzer/current/bin/outLyzer_V$OUTLYZER_VERSION.py		# BIN
export OUTLYZER_DESCRIPTION="outLyzer is a computer program whose purpose is to detect variations, specifically low allele frequency variation, in next generation sequencing data (tumor samples, mosaïc mutation)."
export OUTLYZER_REF="https://github.com/EtieM/outLyzer"
TOOLS_LIST=$TOOLS_LIST" OUTLYZER"


# FASTQC
# export FASTQC=$NGS_TOOLS/fastqc/current/bin/fastqc		# BIN
# export FASTQC_VERSION=0.11.8								# VER
# export FASTQC_DESCRIPTION="A quality control tool for high throughput sequence data."
# export FASTQC_REF="http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc"
# TOOLS_LIST=$TOOLS_LIST" FASTQC"


# FASTP
export FASTP_VERSION=1.1.0							# VER
export FASTP=$NGS_TOOLS/fastp/current/bin/fastp		# BIN $NGS_TOOLS/fastp/current/bin/
export FASTP_DESCRIPTION="A tool designed to provide fast all-in-one preprocessing for FastQ files."
export FASTP_REF="https://github.com/OpenGene/fastp"
TOOLS_LIST=$TOOLS_LIST" FASTP"


# UMI TOOLS
export UMITOOLS_VERSION=1.1.6									# VER
export UMITOOLS=$NGS_TOOLS/umi_tools/current/bin/umi_tools		# BIN
export UMITOOLS_DESCRIPTION="UMI-tools contains tools for dealing with Unique Molecular Identifiers (UMIs)/Random Molecular Tags (RMTs) and single cell RNA-Seq cell barcodes."
export UMITOOLS_REF="https://github.com/CGATOxford/UMI-tools"
TOOLS_LIST=$TOOLS_LIST" UMITOOLS"


# BWA
export BWA1_VERSION=0.7.19							# VER
export BWA1=$NGS_TOOLS/bwa/current/bin/bwa			# BIN
export BWA1_DESCRIPTION="package for mapping low-divergent sequences against a large reference genome"
export BWA1_REF="Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]"
TOOLS_LIST=$TOOLS_LIST" BWA1"


# BWA2
export BWA2_VERSION=2.3										# VER
export BWA2=$NGS_TOOLS/bwa-mem2/current/bin/bwa-mem2		# BIN
export BWA2_DESCRIPTION="package for mapping low-divergent sequences against a large reference genome"
export BWA2_REF="Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru. Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019. https://github.com/bwa-mem2/bwa-mem2"
TOOLS_LIST=$TOOLS_LIST" BWA2"


# BWA
export BWA=$BWA1							# BIN
export BWA_VERSION=$BWA1_VERSION			# VER
export BWA_DESCRIPTION=$BWA1_DESCRIPTION
export BWA_REF=$BWA1_REF
TOOLS_LIST=$TOOLS_LIST" BWA"


# STAR
export STAR_VERSION=2.7.11b										# VER
export STAR=$NGS_TOOLS/star-fusion/current/bin/STAR-plain		# BIN
export STAR_DESCRIPTION="Spliced Transcripts Alignment to a Reference"
export STAR_REF="https://github.com/alexdobin/STAR"
TOOLS_LIST=$TOOLS_LIST" STAR"


# BOWTIE2
export BOWTIE_VERSION=2.5.4										# VER
export BOWTIE=$NGS_TOOLS/bowtie2/current/bin/bowtie2			# BIN
export BOWTIE_DESCRIPTION="Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences."
export BOWTIE_REF="Langmead B1, Trapnell C, Pop M, Salzberg SL. (2009) Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol. 2009;10(3):R25. doi: 10.1186/gb-2009-10-3-r25. Epub 2009 Mar 4. [PMID: 19261174]"
TOOLS_LIST=$TOOLS_LIST" BOWTIE"


# BEDTOOLS
export BEDTOOLS_VERSION=2.31.1								# VER
export BEDTOOLS=$NGS_TOOLS/bedtools/current/bin/bedtools	# BIN
export BEDTOOLS_DIR=$NGS_TOOLS/bedtools/current/bin			# DIR
export BEDTOOLS_DESCRIPTION="a powerful toolset for genome arithmetic"
export BEDTOOLS_REF="http://bedtools.readthedocs.org/"
TOOLS_LIST=$TOOLS_LIST" BEDTOOLS"


# ANNOVAR
export ANNOVAR_VERSION=2025May02						# VER
export ANNOVAR=$NGS_TOOLS/annovar/current/bin			# DIR
export ANNOVAR_DESCRIPTION="an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes"
export ANNOVAR_REF="Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research, 38:e164, 2010"
TOOLS_LIST=$TOOLS_LIST" ANNOVAR"


# VARSCAN
export VARSCAN_VERSION=2.4.6-6									# VER
export VARSCAN=$NGS_TOOLS/varscan/current/bin/VarScan.jar		# BIN-JAR
export VARSCAN_DESCRIPTION="variant detection in massively parallel sequencing data"
export VARSCAN_REF="VarScan 2: Koboldt, D., Zhang, Q., Larson, D., Shen, D., McLellan, M., Lin, L., Miller, C., Mardis, E., Ding, L., & Wilson, R. (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing Genome Research DOI: 10.1101/gr.129684.111 "
TOOLS_LIST=$TOOLS_LIST" VARSCAN"


# ITDSEEK
export ITDSEEK_VERSION=1.2-2								# VER
export ITDSEEK=$NGS_TOOLS/itdseek/current/bin/itdseek.sh	# BIN-JAR
export ITDSEEK_DESCRIPTION="FLT3 ITD detection algorithm"
export ITDSEEK_REF="Chun Hang Au, Anna Wa, Dona N. Ho, Tsun Leung Chan and Edmond S. K. Ma. Clinical evaluation of panel testing by next-generation sequencing (NGS) for gene mutations in myeloid neoplasms. Diagn Pathol. 2016 Jan 22;11:11. doi: 10.1186/s13000-016-0456-8."
TOOLS_LIST=$TOOLS_LIST" ITDSEEK"


# RNASeQC
export RNASEQC_VERSION=2.4.2								# VER
export RNASEQC=$NGS_TOOLS/rnaseqc/current/bin/rnaseqc		# BIN-JAR
export RNASEQC_DESCRIPTION="RNA-SeQC is a tool for evaluating the quality of RNA-Seq data."
export RNASEQC_REF="DeLuca DS, Levin JZ, Sivachenko A, Fennell T, Nazaire M-D, Williams C, Reich M, Winckler W, Getz G. RNA-SeQC: RNA-seq metrics for quality control and process optimization. Bioinformatics. 2012;28(11):1530-1532. doi:10.1093/bioinformatics/bts196 PMID: 22539670"
TOOLS_LIST=$TOOLS_LIST" RNASEQC"


# RSCRIPT
export RSCRIPT_VERSION=4.5.2										# VER
export RSCRIPT=$NGS_TOOLS/r-base/current/bin/Rscript				# BIN-JAR
export RSCRIPT_DESCRIPTION="R is a free software environment for statistical computing and graphics. "
export RSCRIPT_REF="https://www.r-project.org/ "
TOOLS_LIST=$TOOLS_LIST" RSCRIPT"


# SNPEFF
export SNPEFF_VERSION=5.3.0a-0							# VER - 5.3.0a
export SNPEFF_FOLDER=$NGS_TOOLS/snpeff/current			# FOLDER
export SNPEFF=$SNPEFF_FOLDER/bin/snpEff.jar				# BIN-JAR
export SNPEFF_DESCRIPTION="Genetic variant annotation and effect prediction toolbox. It annotates and predicts the effects of variants on genes (such as amino acid changes)"
export SNPEFF_REF="A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3., Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (Austin). 2012 Apr-Jun;6(2):80-92 "
TOOLS_LIST=$TOOLS_LIST" SNPEFF"


# ARRIBA
export ARRIBA_VERSION=2.5.1							# VER
export ARRIBA=$NGS_TOOLS/arriba/current/bin/arriba	# BIN-JAR
export ARRIBA_DESCRIPTION="Arriba is a command-line tool for the detection of gene fusions from RNA-Seq data."
export ARRIBA_REF="https://github.com/suhrig/arriba"
TOOLS_LIST=$TOOLS_LIST" ARRIBA"


# STARFUSION
export STARFUSION_VERSION=1.15.1												# VER
export STARFUSION=$NGS_TOOLS/star-fusion/$STARFUSION_VERSION/bin/STAR-Fusion	# BIN
export STARFUSION_ENV=$NGS_TOOLS/star-fusion/$STARFUSION_VERSION				# ENV
export STARFUSION_DESCRIPTION="STAR-Fusion uses the STAR aligner to identify candidate fusion transcripts supported by Illumina reads. STAR-Fusion further processes the output generated by the STAR aligner to map junction reads and spanning reads to a reference annotation set."
export STARFUSION_REF="https://github.com/STAR-Fusion/STAR-Fusion"
TOOLS_LIST=$TOOLS_LIST" STARFUSION"


# VARIANTCONVERT
export VARIANTCONVERT_VERSION=2.0.1															# VER
export VARIANTCONVERT=$NGS_TOOLS/variantconvert/$VARIANTCONVERT_VERSION/bin/variantconvert	# BIN
export VARIANTCONVERT_ENV=$NGS_TOOLS/variantconvert/$VARIANTCONVERT_VERSION					# ENV
#export VARIANTCONVERT_CONFIGS=$NGS_TOOLS/variantconvert/$VARIANTCONVERT_VERSION/src/variantconvert/configs	# CONFIG
#export VARIANTCONVERT_CONFIGS=$NGS_TOOLS/variantconvert/$VARIANTCONVERT_VERSION/configs		# CONFIG
export VARIANTCONVERT_DESCRIPTION="VariantConvert is a tool for converting variant call formats."
export VARIANTCONVERT_REF="https://github.com/SamuelNicaise/variantconvert"
TOOLS_LIST=$TOOLS_LIST" VARIANTCONVERT"


# DEEPVARIANT
export DEEPVARIANT_VERSION=1.10.0										# VER
export DEEPVARIANT_DOCKER=google/deepvariant:$DEEPVARIANT_VERSION		# DOCKER
export DEEPVARIANT_DESCRIPTION="DeepVariant is a deep learning-based variant caller developed by Google."
export DEEPVARIANT_REF="https://github.com/google/deepvariant"
TOOLS_LIST=$TOOLS_LIST" DEEPVARIANT"
DOCKER_LIST=$DOCKER_LIST" $DEEPVARIANT_DOCKER"


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
export CAP_VERSION=0.9.13									# VER
export CAP_FOLDER=$NGS_TOOLS/cap/current/bin				# DIR
export CAP=$CAP_FOLDER/CAP									# BIN
export CAP_SOFTCLIPTOQ0=$CAP_FOLDER/CAP.SoftClipToQ0.pl		# BIN-SCRIPT
export CAP_SOFTCLIPTOQ0_VERSION=$CAP_VERSION				# VER
export CAP_ManifestToBED=$CAP_FOLDER/CAP.ManifestToBED.pl	# BIN-SCRIPT
export CAP_ManifestToBED_VERSION=$CAP_VERSION				# VER
export CAP_DESCRIPTION="Clipping Amplicons Primers"
export CAP_REF="inhouse"
TOOLS_LIST=$TOOLS_LIST" CAP"


# # HOWARD
# export HOWARD_VERSION=0.9.15.6							# VER
# export HOWARD_FOLDER=$NGS_TOOLS/howard/$HOWARD_VERSION			# DIR 
# export HOWARD_FOLDER_BIN=$HOWARD_FOLDER/bin				# DIR
# export HOWARD_FOLDER_DOCS=$HOWARD_FOLDER/docs			# DIR
# export HOWARD=$HOWARD_FOLDER_BIN/HOWARD					# BIN-SCRIPT
# export HOWARD_RELEASE=$HOWARD_VERSION
# export HOWARDDIR=$HOWARD_FOLDER_BIN
# export HOWARD_DESCRIPTION="Highly Open and Valuable tool for Variant Annotation & Ranking"
# export HOWARD_REF="inhouse"
# TOOLS_LIST=$TOOLS_LIST" HOWARD"

# HOWARD (devel)
export HOWARD_VERSION=devel								# VER
export HOWARD=$NGS_TOOLS/howard/current/bin/howard		# BIN
export HOWARD_DESCRIPTION="Highly Open and Valuable tool for Variant Annotation & Ranking"
export HOWARD_REF="inhouse"
TOOLS_LIST=$TOOLS_LIST" HOWARD"


# SCRIPTS
export STARK_BED_NORMALIZATION=$STARK_FOLDER_BIN/bed_normalization.awk
export FASTQ_CLEAN_HEADER=$STARK_FOLDER_BIN/fastq_clean_header.awk
export FASTQ_REHEADER=$STARK_FOLDER_BIN/fastq_reheader.awk
export RELOCATE_UMI=$STARK_FOLDER_BIN/relocate_umi.awk
export STARK_RUN_METRICS=$STARK_FOLDER_BIN/runmetrics.py


# PERL5LIB
ENV_PERLLIB=$NGS_FOLDER/tools/perl/lib


# DOCKER
export DOCKER_VERSION=29.2.1												# VER
export DOCKER=$NGS_TOOLS/docker/$DOCKER_VERSION/bin/docker					# BIN
export DOCKER_DESCRIPTION="Docker is a platform for developing, shipping, and running applications."
export DOCKER_REF="https://www.docker.com/"
TOOLS_LIST=$TOOLS_LIST" DOCKER"

# Docker Config
if [ -z "$DOCKER_MOUNTS" ]; then
	export DOCKER_MOUNTS=$($PYTHON $STARK/extract_mounts.py --docker=$DOCKER --volumes_from_enable)	# MOUNTS
fi;
export DOCKER_MOUNTS
export DOCKER_RUN="$DOCKER run $DOCKER_MOUNTS"								# RUN


# PATH
########

# ADD JAVA
export PATH=$PATH:$TABIX_PATH:$JAVA_PATH

# ADD PERL
if (($(grep -c " 6." /etc/centos-release 2>/dev/null))); then
	export PERL5LIB=$BCL2FASTQ_PERLLIB:$ENV_PERLLIB
else
	export PERL5LIB
fi;