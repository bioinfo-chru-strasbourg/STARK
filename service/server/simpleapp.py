#!/usr/bin/env python

from flask import Flask, request, jsonify
import sys
import optparse
import time
import os
import subprocess
import urllib
import json
import re
import string
import random


# DEF

def randomStringDigits(stringLength=6):
    """Generate a random string of letters and digits """
    lettersAndDigits = string.ascii_letters + string.digits
    return ''.join(random.choice(lettersAndDigits) for i in range(stringLength))


# Task Spooler env param

#TS_SOCKET="STARK"
#TS_SAVELIST="/tmp/STARKTSLIST";
#TS_ENV=" export TS_SOCKET=" + TS_SOCKET + " && export TS_SAVELIST=" + TS_SAVELIST + " ";
#TS_ENV=" export TS_SAVELIST=" + TS_SAVELIST + " ";
#TS_BIN="ts";
#TS=TS_ENV + " && " + TS_BIN + " ";
#TS=""


app = Flask(__name__)

#app.url_map.strict_slashes = True

start = int(round(time.time()))

@app.route("/analysis", methods = ['POST']) #methods = ['GET', 'POST', 'DELETE'])
#def hello_world(json=None):
def stark_launch():


    # Input JSON data
    json_input = request.get_json(force=True)
    json_dump=json.dumps(json_input)
    #print json_dump

    # analysis ID & NAME default
    analysesID=randomStringDigits(12)
    analysesNAME=analysesID


    # Find analysisID from RUN info in JSON
    runID=""
    if 'run' in json_input.keys():
        #print json_input['run']
        for i, name in enumerate(json_input['run'].split(":")):
            #print i, name
            #runID=name #run_split[i]
            #runID=re.sub("/", '_', name)
            runID=os.path.basename(name)
    if runID != "":
        analysesNAME=runID

    # Find analysisID from ANALYSIS_NAME info in JSON
    analysis_name=""
    if 'analysis_name' in json_input.keys():
        #print json_input['run']
        analysis_name=json_input['analysis_name']
    if analysis_name != "":
        analysesNAME=analysis_name

    #print "analysesID:"+analysesID

    # DOCKER MOUNT PARAMETERS
    docker_mount=""
    docker_stark_folder_pattern="DOCKER_STARK_FOLDER_"
    for docker_stark_folder in os.environ:
        if docker_stark_folder.startswith(docker_stark_folder_pattern):
            #print "Found",docker_stark_folder
            docker_stark_folder_basename=re.sub('_', '/', re.sub(docker_stark_folder_pattern, '', docker_stark_folder)).lower()
            #print docker_stark_folder_basename
            docker_mount+=" -v "+os.environ[docker_stark_folder]+":"+"/STARK/"+docker_stark_folder_basename
        #print docker_stark_folder
    #print docker_mount

    docker_env=""
    if 'DOCKER_STARK_ENV' in os.environ:
        docker_env=os.environ['DOCKER_STARK_ENV']

    # Docker name
    docker_name=" --name " + analysesID + '.' + analysesNAME + " "

    # docker parameters
    docker_parameters=' --rm ' + docker_mount + ' ' + docker_env + ' ' + docker_name + ' '



    analysisFILE="/STARK/analyses/analysis."+analysesID+".json"
    analysisLOG="/STARK/analyses/analysis."+analysesID+".log"
    analysisERR="/STARK/analyses/analysis."+analysesID+".err"
    analysisLOGERR_PARAM=' 1>' + analysisLOG + ' 2>' + analysisERR
    f = open(analysisFILE, "w")
    f.write(json_dump)
    f.close()

    # Prepare and Launch TS Docker run command
    myCmd = 'ts -L ' + analysesID + '.' + analysesNAME + ' docker run ' + docker_parameters + ' stark --analysis_name=' + analysesNAME + ' --analysis=' + analysisFILE #+ analysisLOGERR_PARAM
    print myCmd
    getCmd = subprocess.check_output(myCmd, shell=True);
    #return "Hello world from Distelli & Docker!"+json+getCmd
    if getCmd:
        return analysesID + '.' + analysesNAME, 200
    else:
        return "KO"



@app.route("/queue", methods = ['GET','POST']) #methods = ['GET', 'POST', 'DELETE'])
def queue():

    # queue action
    if 'log' in request.args:
        id=request.args['log']
        myCmd = 'cat $(ts -o ' + str(id) + ')'
    elif 'kill' in request.args:
        id=request.args['kill']
        myCmd = 'ts -k ' + str(id) + ' && echo "#[INFO] task #' + id + ' killed" ' + ' &&  ts -l'
    elif 'info' in request.args:
        id=request.args['info']
        myCmd = 'ts -i ' + str(id)
    elif 'prioritize' in request.args:
        id=request.args['prioritize']
        myCmd = 'ts -u ' + str(id) + ' && echo "#[INFO] task #' + id + ' prioritized" ' + ' &&  ts -l'
    elif 'remove' in request.args:
        id=request.args['remove']
        myCmd = 'ts -r ' + str(id) + ' && echo "#[INFO] task #' + id + ' removed"'  + ' &&  ts -l'
    elif 'swap' in request.args:
        id=request.args['swap']
        myCmd = 'ts -U ' + str(id) + ' && echo "#[INFO] tasks #' + id + ' swapped" ' + ' &&  ts -l'
    else:
        myCmd = 'ts -l'

    #print myCmd

    # format
    if 'format' in request.args:
        format=request.args['format']
    else:
        format=""

    # command
    try:
        getCmd = subprocess.check_output(myCmd, shell=True)
        ret_val=200
    except:
        getCmd = "#[ERROR] Command error"
        ret_val=500


    # output
    if format == "html":
        ret="<pre>"+getCmd+"</pre>"
    else:
        ret=getCmd

    return ret, ret_val




if __name__ == '__main__':
    parser = optparse.OptionParser(usage="python simpleapp.py -p ")
    parser.add_option('-p', '--port', action='store', dest='port', help='The port to listen on.')
    (args, _) = parser.parse_args()
    if args.port == None:
        print "Missing required argument: -p/--port"
        sys.exit(1)
    app.run(host='0.0.0.0', port=int(args.port), debug=False)
