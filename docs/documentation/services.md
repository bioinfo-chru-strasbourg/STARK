# Services

In addition to STARK CLI, STARK provides services to automatically launch runs newly sequenced and available in `/STARK/input/runs`, through a listener and an API. STARK also provides a DAS to share data publically.

## STARK Command Line Interface (CLI)

A STARK Command Line Interface (CLI) is started as a container to execute custom analyses with data and runs, available in STARK main folder (default `${HOME}/STARK`, with `${HOME}/STARK/data` corresponds to `/STARK/data`, with `${HOME}/STARK/input/runs` corresponds to `/STARK/input/runs`, etc.).

Use STARK Command Line Interface with command `docker exec stark-module-stark-submodule-stark-service-cli STARK`, to execute a STARK command with your own data and runs (run names will be automatically found in input folder).

```bash
# Analyze FASTQ with a design file and generated results on an output folder
docker exec stark-module-stark-submodule-stark-service-cli STARK --reads="/STARK/data/my_sample.R1.fastq.gz" --reads2="/STARK/data/my_sample.R2.fastq.gz" --design="/STARK/data/my_sample.bed" --application="default.app" --output="/STARK/data/my_output"
```

```bash
# Analyze a run as a folder (of multiple FASTQ files or an Illumina folder)
docker exec stark-module-stark-submodule-stark-service-cli STARK --run="/STARK/data/my_data/"
```

```bash
# Analyze a run as a Illumina folder by name in input/run folder (e.g. /STARK/input/runs/MY_RUN)
docker exec stark-module-stark-submodule-stark-service-cli STARK --run="MY_RUN"
```

STARK Command Line Interface can be used in interactive mode ('-ti' option).

```bash
# Connect to STARK CLI
docker exec -ti stark-module-stark-submodule-stark-service-cli bash
# Launch an analysis
STARK --run="MY_RUN"
```

```bash
# Connect to STARK CLI and use available tools
docker exec -ti stark-module-stark-submodule-stark-service-cli bash
# Explore your data
cd /STARK/data
ls -lah
# Launch a tool command
samtools view my_sample.bam
bcftools view my_sample.vcf.gz
```

All tools used by STARK can be executed as they are in the PATH environment variable (e.g. samtools, bcftools). Available tools can be found in 'STARK/tools' folder.

```bash
# Launch tools 
docker exec stark-module-stark-submodule-stark-service-cli samtools view /STARK/data/my_sample.bam
docker exec stark-module-stark-submodule-stark-service-cli bcftools view /STARK/data/my_sample.vcf.gz
# List available tools 
docker exec stark-module-stark-submodule-stark-service-cli bash -c "find /STARK/tools -mindepth 2 -maxdepth 2 -type d"
```

## STARK XTerm

STARK XTerm is a web interface providing a terminal to launch STARK command within STARK CLI service. STARK XTerm is available through URI http://localhost:4299 (by default, account 'stark/password').

![panel](images/XTerm.png)

## STARK Docs

STARK Docs provides documentation such as a user guide and many other information to configure STARK analysis. STARK Docs is available through URI <http://localhost:4298> (by default).

![panel](images/mkdocs.png)

## STARK Application Program Interface (API)

A STARK Application Program Interface (API) is available through URI `http://<ip>:<port>` (default <http://localhost:4200>, help with an internet browser). This service prodives an interface to run STARK analysis with parameters in JSON format through URI (`http://<ip>:<port>/analysis`), and to manage analyses queue (`http://<ip>:<port>/queue`)

```bash
# STARK analysis with curl in POST method
curl -X POST -H 'Content-Type: application/json' -d '{"run":"MY_RUN"}' http://<ip>:<port>/analysis 
# List of analysis running, queued and finished
curl http://<ip>:<port>/queue?list 
```

## STARK listener

A STARK listener service is started as a daemon, listening for new sequenced NGS run (new folder in input/runs) and well configured (RTAComplete.txt and SampleSheet.csv), and send a request to STARK API. A STARK listener clear service is checking (once at services start) STARK listener and STARK API log files to reload requests if needed (useful after a server stop/crash).

## STARK DAta Sharing (DAS)

A STARK DAta Sharing (DAS) web server provides data publically through URI `http://<ip>:<port>` (default <http://localhost:4201/>`<path>`). This server may be used with application able to open file through URI (such as IGV), or to share data and files between other STARK modules and services.

By default, data available are (`<path>`):

- runs: inputs/Input/runs
- repository: repositories/Repository
- archives: repositories/Archives
- data: data
- databases: databases
