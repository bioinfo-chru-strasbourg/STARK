# Help

```bash
# USAGE: STARK --analysis=<FILE>|--reads=<FASTQ>|--run=<RUN> [options...]

### Launch an analysis through (by order of priority):
# --analysis=<FILE1,FILE2...>              List of configuration file in JSON format defining options (see below).
#                                          Format: { "option1":"value1", "option2":"value2"...}
# --reads=<FILE1,FILE2...>                 List of Reads files
#                                          For each samples, allowed formats FASTQ|BAM|SAM|CRAM:
#                                          *fastq.gz|*fq.gz|*bam|*ubam|*cram|*ucram|*sam|*usam
# --run=<FOLDER1,FOLDER2...>               List of RUNs to analyse
#                                          From Illumina sequencers (BCL), demultiplexing (FASTQ), or any folder containing FASTQ|BAM|SAM|CRAM.
#                                          For folder containing FASTQ, corresponding Paired-End Read2 will be autodetected: /.*[._-])R1([._-].*/.*[._-])R2([._-].*/
#                                          RUN can be a folder path, or the run folder in defined runs folders (see Application configuration)
#                                          Analysis name of each RUN can be defined using ':' after each RUN (default: run folder name)
#                                          Format: RUN1:RUN1_ANALYSIS_NAME,RUN2:RUN2_ANALYSIS_NAME...
#                                          Example: /path_of_my_runs/my_run:my_run1,run_name_folder_in_configured_runs_folder:my_run2

### SAMPLE Analysis options
# --sample=<STRING1,STRING2...>            List of corresponding SAMPLE Name
#                                          Automatically detected from input files.
# --sample_tag=<STRING1,STRING2...>        List of corresponding SAMPLE Tags.
#                                          Format: [TYPE]#TAG[#TAG]{!}[[TYPE]#TAG[#TAG]]
# --reads2=<FILE1,FILE2...>                List of Reads files with corresponding Paired-End Read2 (beware of correspondance).
#                                          For each samples, allowed formats FASTQ:
#                                          *fastq.gz|*fq.gz
# --index1=<FILE1,FILE2...>                List of corresponding Index1 files (beware of correspondance).
#                                          For each samples, allowed formats FASTQ:
#                                          *fastq.gz|*fq.gz
# --index2=<FILE1,FILE2...>                List of corresponding Index2 files (beware of correspondance).
#                                          For each samples, allowed formats FASTQ:
#                                          *fastq.gz|*fq.gz
# --other_files=<FILE1,FILE2...>           List of corresponding other files (beware of correspondance).
#                                          Files will be copied into each SAMPLE folder

### RUN Analysis options
# --samplesheet=<FILE1,FILES2...>          List of corresponding Illumina SampleSheet.csv file to use with RUNS
#                                          Default: found in RUN folder
# --sample_filter=<STRING1,STRING2...>     List of SAMPLE name in Illumina SampleSheet.csv to analyse
#                                          Default: all samples in the SampleSheet
# --demultiplexing_only                    Perform only demultiplexing (do NOT analysis FASTQ)

### Folder options
# --input=<FOLDER>                         Input folder (Default: Defined in Application).
#                                          Contains folders 'runs' and 'manifests' to analyse Illumina RUN (see --run option)
# --output=<FOLDER>                        Output folder (Default: Defined in Application).
#                                          Contains folders 'demultiplexing', 'results', 'log' and 'tmp'
# --demultiplexing=<FOLDER>                Output/demultiplexing folder
# --results=<FOLDER>                       Output/results directory to generate RESULTS|RUN|SAMPLE|* files
# --log=<FOLDER>                           Output/log directory to generate RESULTS|RUN|SAMPLE|* files
# --tmp=<FOLDER>                           Output/tmp directory to generate RESULTS|RUN|SAMPLE|* files
#                                          Default: Defined in Application, or first --fastq file folder
# --repository=<FOLDER>                    Repository directory to generate GROUP|PROJECT|RUN|SAMPLE|* specific files
#                                          Default: no copy in a repository
# --archives=<FOLDER>                      Archives directory to generate GROUP|PROJECT|RUN|SAMPLE|* specific files
#                                          Default: no copy in a archives
# --favorites=<FOLDER>                     Favorites directory to generate GROUP|PROJECT|RUN|* specific files from repository folder
#                                          Default: no copy in a favorites
# --databases=<FOLDER>                     Databases folder (requires STARK databases folder structure)

### Other options
# --analysis_name=<STRING>                 Analysis name.
#                                          Default: 'date' (format YYYYMMDD) for a SAMPLE list analysis .
# --analysis_tag=<STRING>                  List of ANALYSIS Tags.
#                                          Format: TYPE#TAG[#TAG]{!}[TYPE#TAG[#TAG]]
# --application=<STRING|FILE>              APP name or APP file configuration of the APPLICATION.
#                                          Must be in the STARK APPS folder if relative path
#                                          Default: defined in the RUN SampleSheet, or default.app if not defined
# --design=<FILE1,FILE2...>                List of corresponding design for SAMPLE analysis
#                                          OR force design for RUN analysis (use only one design).
#                                          If not *.bed file (BED format), considered as Illumina manifest.
#                                          Default: first BED or empty
# --panels=<FILE1+FILE2,FILE3+FILE4...>    List of corresponding GENES files (Panels).
#                                          OR force GENES file for RUN analysis (use only one GENES file).
#                                          Format: <FILE1+FILE2,FILE3+FILE4...>, multiple panels for each sample with '+' separator.
#                                          File format: BED (chr<TAB>start<TAB>stop<TAB>strand<TAB>gene_name).
#                                          Default first GENES file or empty
# --transcripts=<FILE1,FILE2...>           List of corresponding TRANSCRIPTS files.
#                                          OR force TRANSCRIPTS file for RUN analysis (use only one TRANSCRIPTS file).
#                                          Format: TSV (transcript<TAB>geneID).
#                                          Default first TRANSCRIPTS file or empty
# --pedigree=<FILE1,FILE2...>              List of corresponding PEDIGREE files.
#                                          OR force PEDIGREE file for RUN analysis.
#                                          Format: PED (see GATK doc).
#                                          Default first PEDIGREE file or empty
# --threads=<INTEGER>                      Number of thread to use
#                                          Default: all cores in your system minus one
# --by_sample                              Split analysis by SAMPLE, all threads on each sample, one by one.
# --keep_alignment                         Keep alignment from input read file, for format BAM|SAM|CRAM only (alignment autodetected).
#                                          CRAM Archives integrity will not be checked

### Launch STARK with Docker
# --docker-compose-file=<FILE>             Docker compose file.
#                                          Default: 'docker-compose.yml' if any, or ''
#                                          Warning: Docker env file '.env' in the Docker compose file folder will be used to build images and populate databases
# --docker-env-file=<FILE>                 Docker environment file.
#                                          Default: found in Docker compose file folder if any, or ''
# --docker-parameters=<STRING>             Docker parameters, added if Docker compose file and Docker env file exists.
#                                          Format: see Docker doc
#                                          Default: ''
# --docker-stark-image=<STRING>            Docker STARK image.
#                                          Default: found in Docker env file if any, or 'stark:latest'
# --docker-stark-container=<STRING>        Docker STARK container name. Start container as deamon/detached (if any) and execute command
#                                          Default: no container started

### Information options
# --applications_infos                     Applications informations.
#                                          Use --runs=<MyRun> to detect RAW FOLDER and APPLICATION of <MyRun>.
#                                          Use --application=<MyApp> to show only <MyApp> informations.
# --applications_infos_all                 Applications informations with all variables (--runs and --application option available).
# --pipelines_infos                        Pipelines informations.
# --release_infos                          Pipelines with tools and databases information.
# --tools_infos                            Tools release information.
# --databases_infos                        Databases information.
# --verbose                                VERBOSE
# --debug                                  DEBUG
# --release                                RELEASE
# --help                                   HELP
```