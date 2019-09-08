
source stark-service.env
#echo $DOCKER_STARK_DATA

# DOCKER STARK FOLDERS
DOCKER_STARK_FOLDERS=""

[ ! -z $DOCKER_STARK_DATABASES ] && [ -d $DOCKER_STARK_DATABASES ] && DOCKER_STARK_FOLDERS=$DOCKER_STARK_FOLDERS" -v $DOCKER_STARK_DATABASES:/STARK/databases "
[ ! -z $DOCKER_STARK_INPUT ] && [ -d $DOCKER_STARK_INPUT ] && DOCKER_STARK_FOLDERS=$DOCKER_STARK_FOLDERS" -v $DOCKER_STARK_INPUT:/STARK/input "
[ ! -z $DOCKER_STARK_OUTPUT ] && [ -d $DOCKER_STARK_OUTPUT ] && DOCKER_STARK_FOLDERS=$DOCKER_STARK_FOLDERS" -v $DOCKER_STARK_OUTPUT:/STARK/output "
[ ! -z $DOCKER_STARK_OUTPUT_TMP ] && [ -d $DOCKER_STARK_OUTPUT_TMP ] && DOCKER_STARK_FOLDERS=$DOCKER_STARK_FOLDERS" -v $DOCKER_STARK_OUTPUT_TMP:/STARK/output/tmp "
[ ! -z $DOCKER_STARK_OUTPUT_RESULTS ] && [ -d $DOCKER_STARK_OUTPUT_RESULTS ] && DOCKER_STARK_FOLDERS=$DOCKER_STARK_FOLDERS" -v $DOCKER_STARK_OUTPUT_RESULTS:/STARK/output/results "
[ ! -z $DOCKER_STARK_OUTPUT_LOG ] && [ -d $DOCKER_STARK_OUTPUT_LOG ] && DOCKER_STARK_FOLDERS=$DOCKER_STARK_FOLDERS" -v $DOCKER_STARK_OUTPUT_LOG:/STARK/output/log "
[ ! -z $DOCKER_STARK_REPOSITORY ] && [ -d $DOCKER_STARK_REPOSITORY ] && DOCKER_STARK_FOLDERS=$DOCKER_STARK_FOLDERS" -v $DOCKER_STARK_REPOSITORY:/STARK/repository "
[ ! -z $DOCKER_STARK_ARCHIVE ] && [ -d $DOCKER_STARK_ARCHIVE ] && DOCKER_STARK_FOLDERS=$DOCKER_STARK_FOLDERS" -v $DOCKER_STARK_ARCHIVE:/STARK/archive "
[ ! -z $DOCKER_STARK_DATA ] && [ -d $DOCKER_STARK_DATA ] && DOCKER_STARK_FOLDERS=$DOCKER_STARK_FOLDERS" -v $DOCKER_STARK_DATA:/STARK/data "

echo $DOCKER_STARK_FOLDERS


docker run --rm -d \
$DOCKER_STARK_FOLDERS \
-v /var/run/docker.sock:/var/run/docker.sock \
--env-file stark-service.env \
-p $DOCKER_STARK_PORT:80 \
--name stark-service-8 \
stark-service:0.9.18b-10


#--entrypoint=docker-entrypoint.sh \
# -v /Users/lebechea/Documents/NGS/DATA/RAW/MANIFESTS/:/STARK/input/manifests \
#-v /Users/lebechea/Documents/NGS/DATA/REPOSITORY:/STARK/repository \
#-v /Users/lebechea/Documents/NGS/DATA/ARCHIVE:/STARK/archive \
