STARK
============

* Stellar Tools for variants Analysis and RanKing
* Author: Antony Le Béchec
* Copyright: IRC
* Licence: GNU GPLA V3
* Release : 0.9.18
* Date : 20191216

STARK is a Next-Generation Sequencing data analysis pipeline for clinical diagnosis


Getting Started
---------------

Welcome to STARK environment!!!

First, configure ".env" file to setup input and output folders, databases location and analyses log folders.

Then, launch docker compose to build images and launch all services:

```
$ docker-compose up -d
```

Finally, use a browser to open the dashboard (default "http://<your_ip>:4200")


For command line execution, just launch STARK command and follow the instructions to analyze your data:

```
$ bin/STARK --help
```
