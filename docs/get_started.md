# STARK Get started

See [STARK User guide](user_guides/USER-GUIDE.en.19.0.0.md) for more information.

## Quick installation

Use curl from GitHub bioinfo-chru-strasbourg to setup STARK environment by default. This setup will build STARK docker image and setup folders (in "${HOME}/STARK") and needed databases.

```bash
mkdir -p ${HOME}/STARK && cd ${HOME}/STARK && curl https://raw.githubusercontent.com/bioinfo-chru-strasbourg/STARK/master/setup.sh | bash
```

See [STARK Installation](documentation/installation.md) for more information.

## Quick access

STARK CLI (Command Line Interface) is started as a container to execute custom analyses with data and runs, available in STARK main folder (default `${HOME}/STARK`, with `${HOME}/STARK/data` corresponds to `/STARK/data`, with `${HOME}/STARK/input/runs` corresponds to `/STARK/input/runs`, etc.).
Connect to STARK CLI with Docker command:
```bash
docker exec -ti stark-module-stark-submodule-stark-service-cli bash
STARK --help
# USAGE: STARK --analysis=<FILE>|--reads=<FASTQ>|--run=<RUN> [options...]
...
```

In addition, STARK XTerm is a web interface providing a terminal to launch STARK command within STARK CLI service.
STARK XTerm is available through URI <http://localhost:4299> (by default, account 'stark/password').

![panel](documentation/images/XTerm.png)

See [STARK services](documentation/services.md) for more information

## Quick docs

STARK Docs provides documentation such as a user guide and many other information to configure STARK analysis. STARK Docs is available through URI <http://localhost:4298> (by default).

![panel](documentation/images/mkdocs.png)

For more information about STARK Command, use `help` option.

```bash
docker exec stark-module-stark-submodule-stark-service-cli STARK --help
```

For informations about available applications (e.g. EXOME, GENOME, GERMLINE, SOMATIC), pipelines, tools, databases, use options:

```bash
### Information options
# --applications_infos                     Applications informations
#                                          Use --runs=<MyRun> to detect RAW FOLDER and APPLICATION of <MyRun>
#                                          Use --application=<MyApp> to show only <MyApp> informations
# --applications_infos_all                 Applications informations with all variables (--runs and --application option available)
# --pipelines_infos                        Pipelines informations
# --release_infos                          Pipelines with tools and databases information
# --tools_infos                            Tools release information
# --databases_infos                        Databases information
```

See [STARK Help](documentation/help.md) for more information.

## Quick command

Either through STARK CLI or STARK XTerm, STARK command can be executed with your own data and runs (run names will be automatically found in input folder). STARL CLI can be used for a one shot command, or as a interactive mode.

STARK CLI with a one shot command

```bash
# Analyze FASTQ with a design file and generated results on an output folder, with a one shot command
docker exec stark-module-stark-submodule-stark-service-cli STARK --reads="/STARK/data/my_sample.R1.fastq.gz" --reads2="/STARK/data/my_sample.R2.fastq.gz" --design="/STARK/data/my_sample.bed" --application="default.app" --output="/STARK/data/my_output"
```

STARK CLI in interactive mode to launch multiple command in terminal

```bash
# Connect to STARK CLI
docker exec -ti stark-module-stark-submodule-stark-service-cli bash
# Launch a STARK analysis with FASTQ and a design
STARK --reads="/STARK/data/my_sample.R1.fastq.gz" --reads2="/STARK/data/my_sample.R2.fastq.gz" --design="/STARK/data/my_sample.bed" --application="default.app" --output="/STARK/data/my_output"
# Launch a STARK analysis for an entire Illumina run available in ${HOME}/STARK/input/runs, for results in repository folder ${HOME}/STARK/repository
STARK --run="MY_RUN" 
# Launch another tool such as BCFTools with your data
bcftools view /STARK/data/my_sample.vcf.gz
```

## Quick container

A specific STARK container can be created with Docker for custom commands. Mount any volumes you need, such as `/STARK/data` and `/STARK/databases` (needed for STARK analysis), to let them available in the container.

```bash
# Launch an docker terminal interface with a STARK container 
docker run --rm --name "STARK-terminal" -v ${HOME}/STARK/data:/STARK/data -v ${HOME}/STARK/databases:/STARK/databases --entrypoint="bash" -ti stark/stark:19.0.0
```

```bash
# Launch a STARK command with an docker STARK container 
docker run --rm --name "STARK-command" -v ${HOME}/STARK/data:/STARK/data -v ${HOME}/STARK/databases:/STARK/databases stark/stark:19.0.0 --reads="/STARK/data/my_sample.R1.fastq.gz" --reads2="/STARK/data/my_sample.R2.fastq.gz" --design="/STARK/data/my_sample.bed" --application="default.app" --output="/STARK/data/my_output"
```

```bash
# Launch a tool command with an docker STARK container 
docker run --rm --name "STARK-bcftools" -v ${HOME}/STARK/data:/STARK/data --entrypoint="bcftools" stark/stark:19.0.0 view /STARK/data/my_sample.vcf.gz
```
