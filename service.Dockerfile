
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
ENV YUM_INSTALL="gcc bc make wget perl-Switch perl-Digest-MD5 perl-Data-Dumper which zlib-devel zlib 	zlib2-devel zlib2 	bzip2-devel bzip2 	lzma-devel lzma 	xz-devel xz 	ncurses-devel 	unzip crontabs at"
ENV YUM_REMOVE="zlib-devel zlib2-devel bzip2-devel lzma-devel xz-devel ncurses-devel unzip gcc"






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
		cp $TOOLS/$TOOL_NAME/$TOOL_VERSION/config/httpd/webtatic.repo /etc/yum.repos.d/webtatic.repo && \
		ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ /tool ;




###############
# YUM INSTALL #
###############

#RUN yum install -y $YUM_INSTALL ;

# add webtatic for php55
#COPY conf/webtatic.repo /etc/yum.repos.d/webtatic.repo

# install necessary packages
RUN yum upgrade -y && \
	yum updateinfo -y && \
	yum install -y \
	autoconf \
	curl \
	cyrus-sasl-lib \
	gcc \
	g++ \
	httpd \
	libc-dev \
	libedit2 \
	libxml2 \
	make \
	php55w \
	php55w-gd \
	php55w-ldap \
	php55w-mysql \
	php55w-pdo \
	php55w-pecl-apcu \
	php55w-pgsql \
	pkg-config \
	which \
	xz-utils \
	openssl-devel \
	psmisc \
	postfix \
	ImageMagick \
	ImageMagick-devel \
	bc \
	make \
	wget \
	ghostscript \
	libgomp \
	java \
	enscript \
	perl-Switch \
	perl-Digest-MD5 \
	perl-Data-Dumper \
	which \
	rsync \
	wine \
	zlib-devel zlib \
    zlib2-devel zlib2 \
    bzip2-devel bzip2 \
    lzma-devel lzma \
    xz-devel xz \
    ncurses-devel \
	$YUM_INSTALL \
    && \
    yum clean all


#################
# JARVIS/VISION #
#################


ENV TOOL_NAME=stark
ENV TOOL_VERSION=0.9.18d
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

# INSTALL
RUN cp $TOOLS/$TOOL_NAME/$TOOL_VERSION/config/httpd/* /var/www/html/ && \
	cp $TOOLS/$TOOL_NAME/$TOOL_VERSION/config/httpd/docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh && \
	cp $TOOLS/$TOOL_NAME/$TOOL_VERSION/config/httpd/httpd.conf /etc/httpd/conf/httpd.conf && \
	cp $TOOLS/$TOOL_NAME/$TOOL_VERSION/config/httpd/php.ini /etc/php.ini && \
	mkdir -p /var/www/html/tmp && \
	chmod -R 0777 /var/www/html && \
	chmod -R 0777 $TOOLS





##############################
# WORKDIR / ENTRYPOINT / CMD #
##############################

EXPOSE 80

WORKDIR /var/www/html

ENTRYPOINT [ "docker-entrypoint.sh" ]
#CMD [ "docker-entrypoint.sh" ]
