
##############################################################
# Dockerfile Version:   1.1
# Software:             STARK
# Software Version:     0.9.18d
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
# FROM #
########

FROM stark
MAINTAINER Antony Le Bechec <antony.lebechec@gmail.com>
LABEL Software="STARK-SERVICE" \
	Version="0.9.0" \
	Website="none" \
	Description="STARK" \
	License="GNU Affero General Public License (AGPL)" \
	Usage="docker run -d [-v [DATA FOLDER]:/STARK/data -v [DATABASE FOLDER]:/STARK/databases -v [RESULTS FOLDER]:/STARK/results -v [RUNS FOLDER]:/STARK/runs -v [MANIFESTS FOLDER]:/STARK/manifests] stark-service:version"



##############
# PARAMETERS #
##############

ENV STARK_FOLDER=/STARK
ENV TOOLS=$STARK_FOLDER/tools
ENV DATA=$STARK_FOLDER/data
ENV TOOL=$STARK_FOLDER/tool
ENV DATABASES=$STARK_FOLDER/databases
ENV YUM_INSTALL="httpd httpd-tools php55w php55w-gd php55w-ldap php55w-mysql php55w-pdo php55w-pecl-apcu php55w-pgsql"
ENV YUM_REMOVE="gcc"






################
# DEPENDENCIES #
################



#########
# STARK #
#########

#ENV TOOL_NAME=stark
#ENV TOOL_VERSION=0.9.18d
#ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
#ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

#COPY . $TOOLS/$TOOL_NAME/$TOOL_VERSION/

#RUN cp /tool/config/httpd/webtatic.repo /etc/yum.repos.d/webtatic.repo ;




##########
# DOCKER #
##########


#RUN yum-config-manager --add-repo https://download.docker.com/linux/centos/docker-ce.repo && \
#	yum install docker-ce docker-ce-cli containerd.io -y && \
#	yum clean all ;

RUN yum install docker -y && yum clean all


###############
# YUM INSTALL #
###############

#RUN yum install -y $YUM_INSTALL ;

# add webtatic for php55
#COPY conf/webtatic.repo /etc/yum.repos.d/webtatic.repo

# install necessary packages
RUN cp /tool/config/httpd/webtatic.repo /etc/yum.repos.d/webtatic.repo && \
	yum upgrade -y && \
	yum updateinfo -y && \
	yum install -y $YUM_INSTALL && \
    yum clean all




################
# TASK SPOOLER #
################

ENV TOOL_NAME=task-spooler
ENV TOOL_VERSION=current
ENV GIT=https://github.com/thomaspreece/task-spooler.git
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

RUN git clone $GIT ; \
    cd $TOOL_NAME ; \
    make -j ; \
    mkdir -p $DEST/bin ; \
    cp ts $DEST/bin ; \
    rm -rf .git ; \
    cd ../ ;  \
    rm -rf $TOOL_NAME ; \
	mkdir /tmp/ts-tmp; \
	chmod 0777 /tmp/ts-tmp; \
	ln -s /tmp/ts-tmp/ /ts-tmp; \
	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;




#################
# STARK Sercice #
#################


# INSTALL
RUN cp /tool/config/httpd/* /var/www/html/ && \
	cp /tool/config/httpd/docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh && \
	cp /tool/config/httpd/httpd.conf /etc/httpd/conf/httpd.conf && \
	cp /tool/config/httpd/php.ini /etc/php.ini && \
	mkdir -p /var/www/html/tmp && \
	mkdir /STARK/output && chmod -R 0777 /STARK/output && \
	chmod -R 0777 /STARK/output && \
	chmod -R 0777 /var/www/html;





##############################
# WORKDIR / ENTRYPOINT / CMD #
##############################

EXPOSE 80

WORKDIR /var/www/html

ENTRYPOINT [ "docker-entrypoint.sh" ]
#CMD [ "docker-entrypoint.sh" ]
