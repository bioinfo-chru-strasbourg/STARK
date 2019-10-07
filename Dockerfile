
##############################################################
# Dockerfile Version:   1.0
# Software:             STARK
# Software Version:     0.9.18b
# Software Website:     none
# Licence:              GNU Affero General Public License (AGPL)
# Description:          STARK
# Usage:                docker run [-v [DATA FOLDER]:/data -v [DATABASE_FOLDER]:/databases] stark-base:version
##############################################################

##########
# README #
##########

# Config parameters
#    identify yum packages for installation
#    identify yum packages to remove
#
# Dependecies installation
#    identify tools dependences
#    config each tool
#    write installation procedure for each tools
#
# Tool
#    configure tool
#    write isntallation procedure for the tool
#    add link to current and root tool folder
#
# Workdir / Entrypoint / Cmd
#    configure workdir, endpoint and command
#    /!\ no variables in endpoint



########
# ARGS #
#######

ARG DOCKER_STARK_IMAGE_BASE=stark-base:latest



########
# FROM #
########

FROM $DOCKER_STARK_IMAGE_BASE
MAINTAINER Antony Le Bechec <antony.lebechec@gmail.com>
LABEL Software="STARK" \
	Version="0.9.18b" \
	Website="none" \
	Description="STARK" \
	License="GNU Affero General Public License (AGPL)" \
	Usage="docker run [-v [DATA FOLDER]:/STARK/data -v [DATABASE FOLDER]:/STARK/databases -v [RESULTS FOLDER]:/STARK/results -v [RUNS FOLDER]:/STARK/runs -v [MANIFESTS FOLDER]:/STARK/manifests] stark:version"



##############
# PARAMETERS #
##############

ENV STARK_FOLDER=/STARK
ENV TOOLS=$STARK_FOLDER/tools
ENV DATA=$STARK_FOLDER/data
ENV TOOL=/tool
ENV MYAPPS=$STARK_FOLDER/myapps
ENV DATABASES=$STARK_FOLDER/databases
ENV YUM_INSTALL="bc htop"
ENV YUM_REMOVE="null"




###############
# YUM INSTALL #
###############

RUN yum install -y $YUM_INSTALL ;



################
# DEPENDENCIES #
################



#########
# STARK #
#########

ENV TOOL_NAME=stark
ENV TOOL_VERSION=0.9.18d
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

COPY . $TOOLS/$TOOL_NAME/$TOOL_VERSION/

RUN mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
		ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ; \
		ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ $TOOL ; \
		ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/config/apps/myapps/ $MYAPPS ;


		#mkdir -p $DATABASES; \
		#ln -s /home /home1 ; \
		#ln -s $DATABASES $TOOL_DATABASE_FOLDER ;



######################
# YUM REMOVE & CLEAR #
######################

RUN yum erase -y $YUM_REMOVE ; yum clean all ;



##############################
# WORKDIR / ENTRYPOINT / CMD #
##############################


WORKDIR "/STARK/data"

ENTRYPOINT [ "/tool/bin/STARK" ]
# find /Users/lebechea/Documents/NGS/Dockerfile -maxdepth 2 -mindepth 2 -type d
#CMD [ "ls", "-l", "/home/TOOLS/tools" ]
#CMD [ "find", "/home/TOOLS/tools", "-maxdepth", "2", "-mindepth", "2", "-type", "d" ]
