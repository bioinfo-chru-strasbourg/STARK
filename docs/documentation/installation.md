# Installation

## Quick installation

Use curl from GitHub bioinfo-chru-strasbourg to setup STARK environment by default. This setup will build STARK docker image and setup folders and needed databases.

```bash
mkdir -p ${HOME}/STARK && cd ${HOME}/STARK && curl https://raw.githubusercontent.com/bioinfo-chru-strasbourg/STARK/master/setup.sh | bash
```

Use STARK Command Line Interface (CLI) to execute custom analyses with data in ${HOME}/STARK/data (/STARK/data whtin the container).

```bash
docker exec stark-module-stark-submodule-stark-service-cli STARK --help
```

## Complete installation

### Download

Download STARK script from BioInfoDiag GitLab.

```bash
git clone https://github.com/bioinfo-chru-strasbourg/STARK.git .
```

### Configuration

Edit ".env" file to configure STARK environment with ".env". Basically, change the STARK main host folder with the variable "DOCKER_STARK_MAIN_FOLDER" (default ${HOME}/STARK). All sub-folders (input, output, databases location... and create/configure them by yourself) can be configured, such as STARK variables (see ".env" file comments). The default configuration is adequate for a standard environment, but all variables in ".env" file and services in "docker-compose.yml" file can be modified to fit infrastructure specificity.

### Build

Build all docker images needed by STARK environment.

```bash
docker-compose build
```

### Setup

The setup step will create folders (if not exist), populate databases folder if needed, and incrementally archives tools setup sources and binaries. Use `--project-name` if STARK scripts are not in a folder named "STARK". Variable DOCKER_STARK_MAIN_FOLDER corresponds to variable in ".env" configuration file (default "$HOME/STARK").

```bash
DOCKER_STARK_MAIN_FOLDER=<STARK_main_folder>
mkdir -p $DOCKER_STARK_MAIN_FOLDER
docker-compose --project-name STARK up stark-setup
docker-compose --project-name STARK up stark-databases
docker-compose --project-name STARK up stark-sources-archives
```

### Services

Services are located in the folder 'services', and are organized in separated modules (folders), containing 'STARK.docker-compose.yml' file describing services, 'STARK.env' file including all parameters, and 'STARK.module' file describing the module and all services, especially to share information and access to other modules.

To automatically start all services modules (detached):

```bash
services/services.sh --modules=* --command=up
```

Main STARK services in the folder 'services/STARK' contains a CLI (Command Line Interface), an API (Application Program Interface), a Listener and its cleaner, and a DAS service (DAta Sharing).

```bash
services/services.sh --modules=stark --command=up
```
