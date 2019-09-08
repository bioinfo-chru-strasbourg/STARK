# Docker run Example


# Docker for development
docker run \
	-v /Users/lebechea/Documents/NGS/DATABASES_ARCHIVES/1.4:/STARK/databases \
	-v /Users/lebechea/Documents/NGS/GIT/Strasbourg/STARK:/STARK/tools/stark/current \
	-v /Users/lebechea/Documents/NGS:/STARK/data \
	-v /Users/lebechea/Documents/NGS/DOCKER_RESULTS/RESULTS:/STARK/results \
	-v /Users/lebechea/Documents/NGS/DOCKER_RESULTS/REPOSITORY:/STARK/repository \
	-v /Users/lebechea/Documents/NGS/DOCKER_RESULTS/ARCHIVE:/STARK/archive \
	--entrypoint=bash \
	--name STARK \
	--rm \
	-ti stark:0.9.18b-6

# Enter to the docker
docker exec -ti STARK bash
