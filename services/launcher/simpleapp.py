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



# VARIABLES
if 'DOCKER_STARK_IMAGE' in os.environ:
    docker_stark=os.environ['DOCKER_STARK_IMAGE']
else:
    docker_stark="stark"

if 'TS' in os.environ:
    ts=os.environ['TS']
else:
    ts=""

if 'TS_SOCKET' in os.environ:
    ts_socket=os.environ['TS_SOCKET']
else:
    ts_socket="/ts-tmp/STARK"

if 'TS_SAVELIST' in os.environ:
    ts_savelist=os.environ['TS_SAVELIST']
else:
    ts_savelist="/ts-tmp"

#ts_env=" TS_SOCKET=" + ts_socket + " TS_SAVELIST=" + ts_savelist + " "
ts_env=" TS_SAVELIST=" + ts_savelist + " "



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

    # Folders

    if 'DOCKER_STARK_MAIN_FOLDER' in os.environ:
        docker_stark_main_folder=os.environ['DOCKER_STARK_MAIN_FOLDER']
    else:
        docker_stark_main_folder=""

    if 'DOCKER_STARK_SERVICE_LISTENER_FOLDER_LOG' in os.environ:
        docker_stark_service_listener_folder_log=os.environ['DOCKER_STARK_SERVICE_LISTENER_FOLDER_LOG']
    else:
        docker_stark_service_listener_folder_log="analyses/stark-services/listener"

    if 'DOCKER_STARK_INNER_FOLDER_ANALYSES' in os.environ:
        docker_stark_inner_folder_analyses=os.environ['DOCKER_STARK_INNER_FOLDER_ANALYSES']
    else:
        docker_stark_inner_folder_analyses="/STARK/analyses"

    if 'DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_LAUNCHER' in os.environ:
        docker_stark_service_data_subfolder_services_launcher=os.environ['DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_LAUNCHER']
    else:
        docker_stark_service_data_subfolder_services_launcher="stark-services/launcher"



    # Find analysisID from RUN info in JSON
    runID=""
    if 'run' in json_input.keys():
        #print json_input['run']
        if len(json_input['run'].split(",")) == 1:
            #for i, name in enumerate(json_input['run'].split(":")):
            #runID=os.path.basename(name)
            name=json_input['run'].split(":")[0]
            runID=os.path.basename(name)
            myCmd ="find " + name + " -maxdepth 1 -type f -print0 | xargs -0 sha1sum | cut -b-40 | sha1sum | awk '{print $1}'"
            runMD5 = subprocess.check_output(myCmd, shell=True).strip();

    if runID != "":
        #analysesNAME=runID
        analysesNAME="ID-" + runMD5 + "-NAME-" + runID

    # Find analysisID from ANALYSIS_NAME info in JSON
    analysis_name=""
    if 'analysis_name' in json_input.keys():
        #print json_input['run']
        analysis_name=json_input['analysis_name']

    if analysis_name != "":
        analysesNAME=analysis_name

    #print "analysesID:"+analysesID
    # MD5
    #MD5=$(find $IFA -maxdepth 1 -xtype f -print0 | xargs -0 sha1sum | cut -b-40 | sha1sum | awk '{print $1}')
    #MD5=$(find $IFA -maxdepth 1 -type f -print0 | xargs -0 sha1sum | cut -b-40 | sha1sum | awk '{print $1}')
    #RUN_NAME=$(basename $IFA)
    #ID="ID-$MD5-NAME-$RUN_NAME"


    # DOCKER MOUNT PARAMETERS
    docker_mount=""
    docker_stark_folder_pattern="DOCKER_STARK_FOLDER_"
    for docker_stark_folder in os.environ:
        if docker_stark_folder.startswith(docker_stark_folder_pattern):
            #print "Found",docker_stark_folder
            docker_stark_folder_basename=re.sub('_', '/', re.sub(docker_stark_folder_pattern, '', docker_stark_folder)).lower()
            #print docker_stark_folder_basename
            docker_mount+=" -v "+docker_stark_main_folder+"/"+os.environ[docker_stark_folder]+":"+"/STARK/"+docker_stark_folder_basename
        #print docker_stark_folder
    if 'DOCKER_STARK_FOLDER_MOUNT' in os.environ:
        docker_mount+=" " + DOCKER_STARK_FOLDER_MOUNT + " "
    #print docker_mount



    docker_env=""
    if 'DOCKER_STARK_ENV' in os.environ:
        docker_env=os.environ['DOCKER_STARK_ENV']


    analysisIDNAME='STARK.' + analysesID + '.' + analysesNAME
    PostanalysisIDNAME='STARK-POSTANALYSIS.' + analysesID + '.' + analysesNAME

    # Docker name
    docker_name=" --name " + analysisIDNAME + " "

    # docker parameters
    docker_parameters=' --rm ' + docker_mount + ' ' + docker_env + ' ' + docker_name + ' '


    # ${DOCKER_STARK_INNER_FOLDER_ANALYSES}/${DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_LAUNCHER}
    analysisFOLDER=docker_stark_inner_folder_analyses+"/"+docker_stark_service_data_subfolder_services_launcher
    analysisFILE=analysisFOLDER+"/"+analysisIDNAME+".json"
    analysisLOG=analysisFOLDER+"/"+analysisIDNAME+".log"
    analysisERR=analysisFOLDER+"/"+analysisIDNAME+".err"
    analysisINFO=analysisFOLDER+"/"+analysisIDNAME+".info"
    analysisLOGERR_PARAM=' 1>' + analysisLOG + ' 2>' + analysisERR
    f = open(analysisFILE, "w")
    f.write(json_dump)
    f.close()


    # TS
    if ts != "":
        ts_cmd=ts_env + ts + ' -L ' + analysisIDNAME
        #ts_cmd= ts + ' -L ' + analysisIDNAME
    else:
        ts_cmd=""

    # Create a Notification by a command ?

    # Prepare and Launch TS Docker run command
    #myCmd = ts_cmd + ' docker run ' + docker_parameters + ' ' + docker_stark + ' --analysis_name=' + analysesNAME + ' --analysis=' + analysisFILE #+ analysisLOGERR_PARAM
    myCmd = ts_cmd + ' docker run ' + docker_parameters + ' ' + docker_stark + ' --analysis_name=' + analysesNAME + ' --analysis=' + analysisFILE + ' '#+ analysisLOGERR_PARAM
    #print(myCmd);
    getCmd = subprocess.check_output(myCmd, shell=True);
    STARKCmdID=getCmd.strip();

    # Retrive INFO from command
    # Example :
    # FIRST_TASKID=`ts ash -c "sleep 10; echo hi"`
    # ts ash -c "ts -w $FIRST_TASKID && echo there"
    if ts != "":
        #myPostCmd=ts_env + " " + ts + ' ash -c "' + ts_env + ts + ' -w ' + STARKCmdID + ' && ' + ts_env + ts + ' -i ' + STARKCmdID + ' > ' + analysisINFO + '"'
        myPostCmd=ts_env + " " + ts + ' -L ' + PostanalysisIDNAME + ' ash -c "' + ts_env + ts + ' -i ' + STARKCmdID + ' > ' + analysisINFO + '"'
        #myPostCmd=" TS_SOCKET=/tmp/info TS_SAVELIST=/tmp " + " " + ts + ' ash -c "' + ts_env + ts + ' -i ' + STARKCmdID + ' > ' + analysisINFO + '"'
        getPostCmd = subprocess.check_output(myPostCmd, shell=True);


    # Notification or something else
    # command to launch : ( ts -w ; xmessage Finished! ) &

    #return "Hello world from Distelli & Docker!"+json+getCmd
    if getCmd:
        return analysisIDNAME, 200
    else:
        return "KO", 400



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
