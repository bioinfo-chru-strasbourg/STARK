STARK
============
STARK is a Next-Generation Sequencing data analysis pipeline for clinical diagnosis
* Stellar Tools for variants Analysis and RanKing
* Author: Antony Le Béchec
* Copyright: HUS/CPS
* License: GNU GPLA V3
* Release : 0.9.18.2
* Date : 20200820



Getting Started
---------------


Use curl from BioInfoDiag GitLab to setup STARK environment by default.

```
$ mkdir -p ~/STARK && cd ~/STARK && curl https://gitlab.bioinfo-diag.fr/Strasbourg/STARK/raw/master/setup.sh | bash
```


Use STARK Command Line Interface (CLI) to execute custom analyses with data in \~/STARK/data (/STARK/data whtin the container).

```
$ docker exec STARK-CLI STARK --help
```


Complete setup
---------------


---
**1. Download**

Download STARK script from BioInfoDiag GitLab.

```
$ git clone https://gitlab.bioinfo-diag.fr/Strasbourg/STARK.git .
```


---
**2. Configuration**

Edit ".env" file to configure STARK environment with ".env". Basically, change the STARK main host folder with the variable "DOCKER_STARK_MAIN_FOLDER". All sub-folders (input, output, databases location... and create/configure them by yourself) can be configured, such as STARK variables (see ".env" file comments). The default configuration is adequate for a standard environment, but all variables in ".env" file and services in "docker-compose.yml" file can be modified to fit infrastructure specificity.


---
**3. Build**

Build all docker images needed by STARK environment.

```
$ docker-compose build
```


---
**4. Setup**

The setup step will create folders (if not exist), populate databases folder if needed, and incrementally archive tools setup sources and binaries. Use --project-name if STARK scripts are not in a folder named "STARK". Variable DOCKER_STARK_MAIN_FOLDER corresponds to variable in ".env" configuration file.

```
$ DOCKER_STARK_MAIN_FOLDER=<STARK_main_folder>
$ mkdir -p $DOCKER_STARK_MAIN_FOLDER
$ docker-compose --project-name STARK up stark-folders
$ docker-compose --project-name STARK up stark-databases
$ docker-compose --project-name STARK up stark-sources-archive
```

---
**5. Services**

Services are located in the folder 'services', and are organized in separated modules (folders), containing 'STARK.docker-compose.yml' file describing services, 'STARK.env' file including all parameters, and 'STARK.module' file describing the module and all services, especially to share information and access to other modules.

To automatically start all services modules:

```
$ services/services.sh --module=* --command=up
```

Main STARK services in the folder 'services/STARK' contains a CLI (Command Line Interface), an API (Application Program Interface), a Listener and its cleaner, and a DAS service (DAta Sharing).

```
$ services/services.sh --module=STARK --command=up
```


Analysis
--------


---
**1. STARK Command Line Interface (CLI)**


A STARK Command Line Interface (CLI) is started as a container to execute custom analyses with data and runs, available in inner main folder (default /STARK/data and /STARK/data, resp.).

Use STARK Command Line Interface with command 'docker exec STARK-CLI STARK', to execute a STARK command with data and runs (run names will be automatically found in input folder). For more information, use HELP option.

```
$ docker exec STARK-CLI STARK --help
$ docker exec STARK-CLI STARK --run=<my_run>
$ docker exec STARK-CLI STARK --reads=/STARK/data/<my_data>/<my_fastq> --design=/STARK/data/<my_data>/<my_design> --application=<my_application>
```

STARK Command Line Interface can be used in interactive mode ('-ti' option). All tools used by STARK can be executed as they are in the PATH environment variable (e.g. samtools, bcftools). Available tools can be found in 'STARK/tools' folder.

```
$ docker exec -ti STARK-CLI bash
$ docker exec STARK-CLI samtools
$ docker exec STARK-CLI bcftools
$ docker exec STARK-CLI bash -c "find /STARK/tools -mindepth 2 -maxdepth 2 -type d"
```

---
**2. STARK Application Program Interface (API)**


A STARK Application Program Interface (API) is available through URI http://\<ip\>:\<port\> (default http://localhost:4200, help with an internet browser). This service prodives an interface to run STARK analysis with parameters in JSON format through URI (http://\<ip\>:\<port\>/analysis), and to manage analyses queue (http://\<ip\>:\<port\>/queue)

```
$ curl -X POST -H 'Content-Type: application/json' -d '{"run":"MY_RUN"}' http://<ip>:<port>/analysis # STARK analysis with curl in POST method:
$ curl http://<ip>:<port>/queue?list # List of analysis running, queued and finished
```

---
**3. STARK DAta Sharing (DAS)**


A STARK DAta Sharing (DAS) web server provides data publically through URI http://\<ip\>:\<port\> (default http://localhost:4201/static/data/public/<path\>). This server may be used with application able to open file through URI (such as IGV), or to share data and files between other STARK modules and services.


By default, data available are (\<path\>):
- runs: inputs/Input/runs
- repository: repositories/Repository
- archive: repositories/Archive
- data: data
- databases: databases



---
**4. STARK listener**


A STARK listener service is started as a daemon, listening for new sequenced NGS run (new folder in input/runs) and well configured (RTAComplete.txt and SampleSheet.csv), and send a request to STARK API. A STARK listener clear service is checking (once at services start) STARK listener and STARK API log files to reload requests if needed (useful after a server stop/crash).
