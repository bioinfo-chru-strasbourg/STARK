STARK
============
STARK is a Next-Generation Sequencing data analysis pipeline for clinical diagnosis
* Stellar Tools for variants Analysis and RanKing
* Author: Antony Le Béchec
* Copyright: HUS/CPS
* License: GNU GPLA V3
* Release : 0.9.18
* Date : 20191216



Getting Started
---------------


---
**1. Download**

Download STARK script from BioInfoDiag GitLab.

```
$ git clone https://gitlab.bioinfo-diag.fr/Strasbourg/STARK.git
$ cd STARK
```


---
**2. Configuration**

Edit ".env" file to configure your STARK environment with ".env". Basically, change the STARK main host folder with the variable "DOCKER_STARK_MAIN_FOLDER".  You can also setup all sub-folders (input, output, databases location...), STARK services port pattern, and all specific variables (see ".env" file comments). The default configuration is adequate for a standard environment, but all variables in ".env" file and services in "docker-compose.yml" file can be modified to fit your infrastructure specificity.


---
**3. Build**

Build all docker images needed by STARK environment.

```
$ docker-compose build
```


---
**4. Setup**

The setup step will create needed folders and populate databases folder

```
$ source .env
$ mkdir -p $DOCKER_STARK_MAIN_FOLDER
$ docker-compose up stark-folders
$ docker-compose up stark-databases
```


---
**5. Start**

Launch docker compose to start all services.

```
$ docker-compose up -d
```


---
**6. Dashboard**

Finally, use a browser to open the dashboard (default "http://localhost:4200")


---
**7. Analysis**

In order to analysis a NGS run, copy target design and raw data folder to "DOCKER_STARK_MAIN_FOLDER/input" folders ("manifests" and "runs", resp.).

For command line execution, just launch STARK command and follow the instructions to analyze your data:

```
$ bin/STARK --help
```
