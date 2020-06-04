STARK module
============
STARK module prodives services to run analyses
* Stellar Tools for variants Analysis and RanKing module
* Author: Antony Le Béchec
* Copyright: HUS/CPS
* License: GNU GPLA V3
* Release : 0.9.18.1
* Date : 20200602



Description
-----------


Main STARK services in the folder 'services/STARK' contains a CLI (Command Line Interface), an API (Application Program Interface), a Listener and its cleaner.



Services
---------


---
**STARK Command Line Interface (CLI)**


A STARK Command Line Interface (CLI) is started as a container to execute custom analyses with data and runs, available in inner main folder (default /STARK/data and /STARK/input/runs, resp.).

Use STARK Command Line Interface with command 'docker exec STARK-module-STARK-service-CLI STARK', to execute a STARK command with data and runs (run names will be automatically found in input folder). For more information, use HELP option.

```
$ docker exec STARK-module-STARK-service-CLI STARK --help
$ docker exec STARK-module-STARK-service-CLI STARK --run=<my_run>
$ docker exec STARK-module-STARK-service-CLI STARK --reads=/STARK/data/<my_data>/<my_fastq> --design=/STARK/data/<my_data>/<my_design> --application=<my_application>
```

STARK Command Line Interface can be used in interactive mode ('-ti' option). All tools used by STARK can be executed as they are in the PATH environment variable (e.g. samtools, bcftools). Available tools can be found in 'STARK/tools' folder.

```
$ docker exec -ti STARK-module-STARK-service-CLI bash
$ docker exec STARK-module-STARK-service-CLI samtools
$ docker exec STARK-module-STARK-service-CLI bcftools
$ docker exec STARK-module-STARK-service-CLI bash -c "find /STARK/tools -mindepth 2 -maxdepth 2 -type d"
```

---
**STARK Application Program Interface (API)**


A STARK Application Program Interface (API) is available through URI http://\<ip\>:\<port\> (default http://localhost:4200, help with an internet browser). This service prodives an interface to run STARK analysis with parameters in JSON format through URI (http://\<ip\>:\<port\>/analysis), and to manage analyses queue (http://\<ip\>:\<port\>/queue)

```
$ curl -X POST -H 'Content-Type: application/json' -d '{"run":"MY_RUN"}' http://<ip>:<port>/analysis # STARK analysis with curl in POST method: 
$ curl http://<ip>:<port>/queue?list # List of analysis running, queued and finished
```

---
**STARK DAta Sharing (DAS)**


A STARK DAta Sharing (DAS) web server provides data publically through URI http://\<ip\>:\<port\> (default http://localhost:4201/static/data/public/<path\>). This server may be used with application able to open file through URI (such as IGV), or to share data and files between other STARK modules and services.


By default, data available are (\<path\>):
- runs: inputs/Input/runs
- repository: repositories/Repository 
- archive: repositories/Archive 
- data: data
- databases: databases


---
**STARK listener**


A STARK listener service is started as a daemon, listening for new sequenced NGS run (new folder in input/runs) and well configured (RTAComplete.txt and SampleSheet.csv), and send a request to STARK API. A STARK listener clear service is checking (once at services start) STARK listener and STARK API log files to reload requests if needed (useful after a server stop/crash).

