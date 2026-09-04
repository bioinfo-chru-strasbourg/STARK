
##############################################################
# Dockerfile Version:   1.3
# Software:             STARK
# Software Version:     3
# Software Website:     https://github.com/bioinfo-chru-strasbourg/STARK
# Licence:              GNU Affero General Public License (AGPL)
# Description:          STARK
# Usage:                docker run -ti [-v [DATA FOLDER]:/data -v [DATABASE_FOLDER]:/databases] stark:version
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

ARG FROM_IMAGE="almalinux:9"

FROM $FROM_IMAGE
LABEL Software="STARK" \
	Version="19.0.0-devel" \
	Website="https://gitlab.bioinfo-diag.fr/Strasbourg/STARK" \
	maintainer="Antony Le Bechec <antony.lebechec@gmail.com>" \
	Description="STARK" \
	License="GNU Affero General Public License (AGPL)" \
	Usage="docker run [-v [DATA FOLDER]:/STARK/data -v [DATABASE FOLDER]:/STARK/databases -v [RESULTS FOLDER]:/STARK/output/results -v [RUNS FOLDER]:/STARK/input/runs -v [MANIFESTS FOLDER]:/STARK/input/manifests] stark:version"



########
# ARGS #
########

# Repeat ARG because of loading image
ARG FROM_IMAGE="almalinux:9"
# Threads
ARG THREADS="1"
# REPO from GIT
ARG REPO_SOURCES="https://gitlab.bioinfo-diag.fr/Strasbourg/STARK-repo/raw/master/"
# REPO from internal HTTP server
#ARG REPO="http://192.168.1.14:8080/"
# REPO null
#ARG REPO=""
ARG REMOVE_SOURCES="1"



##############
# PARAMETERS #
##############

ENV STARK_FOLDER="/STARK"
ENV TOOLS="$STARK_FOLDER/tools"
ENV DATA="$STARK_FOLDER/data"
ENV TOOL="$STARK_FOLDER/tool"
ENV REPO="$REPO_SOURCES"
ENV SOURCES_FOLDER="sources"
ENV SOURCES="$STARK_FOLDER/$SOURCES_FOLDER"
ENV DATABASES="$STARK_FOLDER/databases"
ENV CONFIGS="$STARK_FOLDER/config"
ENV GENOMES="$DATABASES/genomes/current"
ENV WORKDIR="/tmp"
ENV YUM_PARAM=" "
ENV WGET_PARAM=" "
ENV TAR_PARAM=" "
ENV ZIP_PARAM=" "
ENV MAKE_PARAM=" "



###########
# SOURCES #
###########

# Copy sources packages, scripts and tools
ADD ./$SOURCES_FOLDER $SOURCES


###########
# WORKDIR #
###########

WORKDIR $WORKDIR



##########
# HEADER #
##########

# RUN echo "#[INFO] SYSTEM Configuration" && \
# 	echo "#[INFO] STARK_FOLDER=$STARK_FOLDER" && \
# 	echo "#[INFO] TOOLS=$TOOLS" && \
# 	echo "#[INFO] DATA=$DATA" && \
# 	echo "#[INFO] TOOL=$TOOL" && \
# 	echo "#[INFO] SOURCES_FOLDER=$SOURCES_FOLDER" && \
# 	echo "#[INFO] SOURCES=$SOURCES" && \
# 	echo "#[INFO] DATABASES=$DATABASES" && \
# 	echo "#[INFO] WORKDIR=$WORKDIR" && \
# 	echo "#[INFO] REMOVE_SOURCES=$REMOVE_SOURCES" && \
# 	echo "#[INFO] THREADS=$THREADS" && \
# 	echo "#[INFO] REPO=$REPO" && \
# 	echo "#";



##################
# SYSTEM INSTALL #
##################
# This will install system packages, python packages and scripts to install tools

#ENV YUM_INSTALL="autoconf automake htop bc bzip2 bzip2-devel curl gcc gcc-c++ git make mlocate ncurses-devel tbb-devel unzip rsync wget which xz xz-devel zlib zlib-devel docker java-17 java-1.8.0 curl-devel openssl-devel htslib diffutils parallel aria2 jq"
ENV PACKAGES_INSTALL="autoconf automake htop tree bc bzip2 bzip2-devel pigz curl gcc gcc-c++ git make mlocate ncurses-devel tbb-devel unzip rsync wget which xz xz-devel zlib zlib-devel java-21 java-1.8.0 curl-devel openssl-devel diffutils parallel aria2 jq"

ENV PYTHON_MODULE=" pathos numpy scipy argparse pandas bx-python requests"
ENV PERL_INSTALL=" perl perl-Switch perl-Time-HiRes perl-Data-Dumper perl-Digest-MD5 perl-Tk perl-devel"

ENV REPO_SYSTEM_GIT="$REPO/sources.system.tar.gz?path=sources/system"
ENV REPO_SYSTEM_HTTP="$REPO/sources/system/"

ENV GET_TOOL_SOURCE=$SOURCES/get_tool_source.sh
ENV TOOL_INIT=$SOURCES/tool_init.sh
ENV TOOL_CHECK=$SOURCES/tool_check.sh


# System installation
RUN echo "#[INFO] SYSTEM Packages installation" && \
	${SOURCES}/install_system.sh --yum_install="${PACKAGES_INSTALL}" --yum_param="${YUM_PARAM}"


#############
# MINIFORGE #
#############

ENV TOOL_NAME=miniforge
ENV TOOL_VERSION=25.1.1-0
ENV TARBALL_LOCATION=https://github.com/conda-forge/miniforge/releases/download/$TOOL_VERSION
ENV TARBALL=Miniforge-pypy3.sh
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV MAMBA=$DEST/bin/mamba
ENV CONDA=$DEST/bin/conda
ENV PIP=$DEST/bin/pip
ENV PATH=$TOOLS/$TOOL_NAME/current/bin:$PATH

# INSTALL
RUN echo "#[INFO] SYSTEM Mamba installation '$TOOL_NAME:$TOOL_VERSION'" && \
    wget $TARBALL_LOCATION/Miniforge3-$TOOL_VERSION-$(uname)-$(uname -m).sh -O $TARBALL && \
    bash $TARBALL -b -p $DEST && \
    rm -f $TARBALL && \
	find ${DEST} -follow -type f -name '*.a' -or -name '*.pyc' -delete && \
	$MAMBA clean --force-pkgs-dirs --all --yes && \
	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;


##########
# PYTHON #
##########

ENV TOOL_NAME=python
ENV PATH=$TOOLS/$TOOL_NAME/current/bin:$PATH

# PYTHON 3.10 - current
ENV TOOL_NAME=python
ENV TOOL_VERSION=3.10
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PYTHON_ENV=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PYTHON=$PYTHON_ENV/bin/python


# INSTALL
RUN echo "#[INFO] SYSTEM Python installation '$TOOL_NAME:$TOOL_VERSION'" && \
    $MAMBA create -y python=$TOOL_VERSION -p ${DEST} && \
    $MAMBA clean -y --all && \
	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
	$PYTHON -m pip install $PYTHON_MODULE && \
	find ${DEST} -follow -ignore_readdir_race \( -name '*.a' -o -name '*.pyc' -o -name '*.txt' -o -name '*.md' -o -name '*.pdf' -o  -name '__pycache__' \) -exec rm -rf {} + || true




########
# PERL #
########



# PERL installation
RUN	echo "#[INFO] SYSTEM Perl installation - download from yum" && \
	mkdir -p $SOURCES/$SOURCES_FOLDER/perl/build/install && \
	echo "FROM_IMAGE=$FROM_IMAGE" && \
	yum $YUM_PARAM install -y --downloadonly --downloaddir=$SOURCES/$SOURCES_FOLDER/perl/build/install $PERL_INSTALL && \
	yum $YUM_PARAM localinstall -y --nogpgcheck $SOURCES/$SOURCES_FOLDER/perl/build/install/*.rpm && \
	rsync -auczqAXhi --no-links --no-perms --no-owner --no-group --ignore-missing-args $SOURCES/$SOURCES_FOLDER/perl/build/install/*rpm $SOURCES/$SOURCES_FOLDER/system/ && \
	rm -rf $SOURCES/$SOURCES_FOLDER/perl/build && \
	yum clean -y all && \
	rm -rf /var/cache/yum && \
	echo "#[INFO] System Clean" && \
	echo "#";


##########
# DOCKER #
##########

ENV TOOL_NAME="docker"
ENV TOOL_VERSION="29.2.1"
ENV TOOL_TARBALL="docker-$TOOL_VERSION.tgz"
ENV TOOL_SOURCE_EXTERNAL="https://download.docker.com/linux/static/stable/x86_64/$TOOL_TARBALL"
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# TOOL PARAMETERS
ENV DOCKER_HOST=unix:///var/run/docker.sock

# TOOL INSTALLATION
RUN echo "#[INFO] SYSTEM Docker installation '$TOOL_NAME:$TOOL_VERSION'" && \
	source $TOOL_INIT && \
	tar -xvf $TOOL_SOURCE -C $TOOL_DEST/ && \
	pwd $TOOL_DEST && \
	ls -lah $TOOL_DEST/* && \
	mv $TOOL_DEST/docker/* $TOOL_DEST/bin/ && \
	ln -s $TOOL_DEST/bin/docker /usr/local/bin/docker && \
	$TOOL_CHECK ;


##########
# JAVA7 #
##########

ENV TOOL_NAME="java"
ENV TOOL_VERSION="1.7.0"
ENV TOOL_TARBALL="openjdk-7u75-b13-linux-x64-18_dec_2014.tar.gz"
ENV TOOL_SOURCE_EXTERNAL="https://download.java.net/openjdk/jdk7u75/ri/$TOOL_TARBALL"
RUN echo "#[INFO] SYSTEM Java installation '$TOOL_NAME:$TOOL_VERSION'" && \
	source $TOOL_INIT && \
	tar -xvzf $TOOL_SOURCE -C $TOOL_DEST/ --strip-components=1 && \
	$TOOL_CHECK ;


##########
# JAVA8 #
##########

ENV TOOL_NAME="java"
ENV TOOL_VERSION="1.8.0"
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
RUN echo "#[INFO] SYSTEM Java installation '$TOOL_NAME:$TOOL_VERSION'" && \
	mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
	ln -s /usr/lib/jvm/jre-1.8.0/bin/java $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/java ;


########
# JAVA #
########

ENV TOOL_NAME="java"
ENV TOOL_VERSION="21"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
RUN echo "#[INFO] SYSTEM Java installation '$TOOL_NAME:$TOOL_VERSION'" && \
	mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
	ln -snf /usr/lib/jvm/jre-21/bin/java $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/java && \
	ln -snf $TOOL_VERSION/ $TOOLS/$TOOL_NAME/current ;



################
# DEPENDENCIES #
################



### TOOL INSTALLATION CODE FORMAT
# All variables <VARIABLE> must be changed

# ##########
# # <TOOL> #
# ##########
#
# # TOOL INFO
# ENV TOOL_NAME="<TOOL_NAME>"								# tool name
# ENV TOOL_VERSION="<TOOL_RELEASE>"							# tool release
# ENV TOOL_TARBALL="<TOOL_TARBALL>"							# filename of the tarball
# ENV TOOL_SOURCE_EXTERNAL="<TOOL_SOURCE_EXTERNAL>"			# Extenal tarball source of the tool
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH		# Add tool bin to the PATH
# # TOOL PARAMETERS
# ENV TOOL_PARAM_<VARIABLE1>="<VALUE1>"						# Parameter 1
# ENV TOOL_PARAM_<VARIABLE2>="<VALUE2>"						# Parameter 2
#
# # TOOL INSTALLATION
# RUN source $TOOL_INIT && \								# Init tool variables && tarball source && ...
# 	echo "#[INFO] TOOL installation" && \					# head
# 	<Installation CMD such as tar... && cp... && ...> && \	# Installation commands using $TOOL_SOURCE $TOOL_SOURCE_BUILD $TOOL_DEST...
# 	# Example: tar && make install  \						# Example of command (tar && make install)
# 	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
# 	make install -j $THREADS -C $(ls -d $TOOL_SOURCE_BUILD/*) prefix=$TOOL_DEST && \
# 	# Example: unzip && rpm  \								# Example of command (unzip && rpm)
# 	unzip -q $TOOL_SOURCE -d $TOOL_SOURCE_BUILD && \
# 	rpm -ih $TOOL_SOURCE_BUILD/*.rpm --excludedocs --prefix=$TOOL_DEST && \
# 	$TOOL_CHECK ;											# Check tool, rm source build folder...



################
# TOOLS SYSTEM #
################


################
# TOOLS EXTERN #
################


################
# SHARED TOOLS #
################

# SHARED TOOLS - current
ENV TOOL_NAME=shared
ENV TOOL_VERSION=current
ENV DEST=${TOOLS}/${TOOL_NAME}/${TOOL_VERSION}


# # Add tools and install script
# ADD tools.json ${TOOLS}/tools.json
# ADD install_tools.sh ${TOOLS}/install_tools.sh

# Install tools
RUN echo "#[INFO] TOOLS installation shared with Mamba" && \
	${SOURCES}/install_tools.sh --list_tools=${SOURCES}/tools.json --folder_tools=${TOOLS} --shared_tools=${DEST} --mamba=${MAMBA}




###########
# ANNOVAR #
###########

# TOOL INFO
ENV TOOL_NAME="annovar"
ENV TOOL_VERSION="2020Jun08"
ENV TOOL_VERSION="2025May02"
ENV TOOL_TARBALL="$TOOL_NAME.latest.tar.gz"
ENV TOOL_SOURCE_EXTERNAL="http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/$TOOL_TARBALL"
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# TOOL PARAMETERS
ENV TOOL_PARAM_TARBALL_FOLDER=$TOOL_NAME
ENV TOOL_PARAM_DATABASE_FOLDER_LINK=$DATABASES/annovar/current
ENV TOOL_PARAM_DATABASE_FOLDER=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/databases/

# TOOL INSTALLATION
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	source $TOOL_INIT && \
	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
	cp $TOOL_SOURCE_BUILD/*/*.pl $TOOL_DEST/bin/ -R && \
	echo "#[INFO] TOOL databases configuration" && \
	mkdir -p $TOOL_PARAM_DATABASE_FOLDER_LINK && \
	mkdir -p $TOOL_PARAM_DATABASE_FOLDER && \
	ln -s $TOOL_PARAM_DATABASE_FOLDER_LINK $TOOL_PARAM_DATABASE_FOLDER && \
	$TOOL_CHECK ;




#############
# BCL2FASTQ #
#############

# TOOL INFO
ENV TOOL_NAME="bcl2fastq"
ENV TOOL_VERSION="2.20.0"
ENV TOOL_TARBALL=$TOOL_NAME"2-v2-20-0-linux-x86-64.zip"
ENV TOOL_SOURCE_EXTERNAL="https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/$TOOL_TARBALL"
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	source $TOOL_INIT && \
	unzip -q $TOOL_SOURCE -d $TOOL_SOURCE_BUILD && \
	rpm -ih $TOOL_SOURCE_BUILD/*.rpm --excludedocs --prefix=$TOOL_DEST && \
	$TOOL_CHECK ;




#######
# CAP #
#######

# TOOL INFO
ENV TOOL_NAME="cap"
ENV TOOL_VERSION="0.9.13"
ENV TOOL_TARBALL="$TOOL_VERSION.tar.gz"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/bioinfo-chru-strasbourg/CAP/archive/refs/heads/$TOOL_TARBALL"
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	source $TOOL_INIT && \
	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
	cp -R $TOOL_SOURCE_BUILD/*/* $TOOL_DEST/ && \
	chmod a+x $TOOL_DEST/bin/* && \
	$TOOL_CHECK ;



#########
# GATK4 #
#########

# TOOL INFO
ENV TOOL_NAME="gatk4"
ENV TOOL_VERSION="4.6.1.0-0"
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	ln -s $(ls $TOOLS/$TOOL_NAME/$TOOL_VERSION/share/gatk4-$TOOL_VERSION/gatk-package-*-local.jar) $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/GenomeAnalysisTK4.jar



########
# GATK #
########

# TOOL INFO
ENV TOOL_NAME="gatk"
ENV TOOL_VERSION="3.8"
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	ln -s $(ls $TOOLS/$TOOL_NAME/$TOOL_VERSION/opt/*/GenomeAnalysisTK*.jar) $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/GenomeAnalysisTK.jar




###########
# ITDSEEK #
###########

# TOOL INFO
ENV TOOL_NAME="itdseek"
ENV TOOL_VERSION="1.2-2"
ENV TOOL_TARBALL="$TOOL_NAME-$TOOL_VERSION.zip"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/tommyau/itdseek/zipball/master/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	source $TOOL_INIT && \
	unzip -q $TOOL_SOURCE -d $TOOL_SOURCE_BUILD && \
	cp -R $TOOL_SOURCE_BUILD/*/* $TOOL_DEST/bin/ && \
	$TOOL_CHECK ;

##########
# MUTECT #
##########


# TOOL INFO
ENV TOOL_NAME="mutect"
ENV TOOL_VERSION="1.1.6"
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	ln -s $(ls $TOOLS/$TOOL_NAME/$TOOL_VERSION*/share/$TOOL_NAME*$TOOL_VERSION*/muTect*jar) $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/muTect.jar




############
# OUTLYZER #
############

# TOOL INFO
ENV TOOL_NAME="outlyzer"
ENV TOOL_VERSION="3.2"
ENV TOOL_TARBALL="outLyzer_V$TOOL_VERSION.py"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/EtieM/outLyzer/releases/download/$TOOL_VERSION/$TOOL_TARBALL"
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	source $TOOL_INIT && \
	cp $TOOL_SOURCE -d $TOOL_DEST/bin/ && \
	chmod a+x $TOOL_DEST/bin/*.py && \
	$TOOL_CHECK ;



##########
# PICARD #
##########


# TOOL INFO
ENV TOOL_NAME="picard"
ENV TOOL_VERSION="3.4.0-0"
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	ln -s $(ls $TOOLS/$TOOL_NAME/$TOOL_VERSION/share/$TOOL_NAME*$TOOL_VERSION*/picard*jar) $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/picard.jar




##########
# SNPEFF #
##########

# TOOL INFO
ENV TOOL_NAME="snpeff"
ENV TOOL_VERSION="5.4.0a-0"
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	ln -s $(ls $TOOLS/$TOOL_NAME/$TOOL_VERSION/share/$TOOL_NAME*$TOOL_VERSION*/snpEff.jar) $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/snpEff.jar




###########
# VARSCAN #
###########


# TOOL INFO
ENV TOOL_NAME="varscan"
ENV TOOL_VERSION="2.4.6-0"
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	ln -s $(ls $TOOLS/$TOOL_NAME/$TOOL_VERSION/share/*/VarScan*jar) $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/VarScan.jar




##################
# variantconvert #
##################

# TOOL INFO
ENV TOOL_NAME="variantconvert"
ENV TOOL_VERSION="2.0.1"
ENV TOOL_TARBALL=$TOOL_VERSION".tar.gz"  
ENV TOOL_SOURCE_EXTERNAL="https://github.com/SamuelNicaise/$TOOL_NAME/archive/refs/tags/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS
ENV CONFIGS_VARIANTCONVERT_FOLDER=$TOOLS/$TOOL_NAME/$TOOL_VERSION/configs

# TOOL INSTALLATION
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION PYTHON=3.10 && \
	source $TOOL_INIT && \
	tar -xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
	ls -l $TOOL_SOURCE_BUILD && \
	cp -R $TOOL_SOURCE_BUILD/*/* $TOOL_DEST/ && \
	$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/python -m pip install -e $TOOL_DEST && \
	$TOOL_DEST/bin/variantconvert init -d $CONFIGS_VARIANTCONVERT_FOLDER && \
	for variantconvert_assembly in $CONFIGS_VARIANTCONVERT_FOLDER/*; do \
		$TOOL_DEST/bin/variantconvert config -c $variantconvert_assembly/* --set GENOME.path=$GENOMES/$(basename $variantconvert_assembly)/$(basename $variantconvert_assembly).fa; \
	done && \
	$MAMBA clean -y --all && \
	$TOOL_CHECK ;




##########
# HOWARD #
##########

# https://github.com/bioinfo-chru-strasbourg/howard/archive/refs/heads/devel.zip

# TOOL INFO
ENV TOOL_NAME="howard"
ENV TOOL_VERSION="devel"
ENV TOOL_TARBALL="$TOOL_VERSION.zip"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/bioinfo-chru-strasbourg/howard/archive/refs/heads/$TOOL_TARBALL"
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# TOOL PARAMETERS


RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION PYTHON=3.10 && \
	source $TOOL_INIT && \
	unzip -q $TOOL_SOURCE -d $TOOL_SOURCE_BUILD && \
	cp -R $TOOL_SOURCE_BUILD/*/* $TOOL_DEST/ && \
	$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/python -m pip install -e $TOOL_DEST && \
	$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/python -m pip install polars-lts-cpu && \
	$MAMBA clean -y --all && \
	howard query --input=$TOOL_DEST/tests/data/example.vcf --query="SELECT 1" && \
	$TOOL_CHECK ;
	

##########
# MKDOCS #
##########


# TOOL INFO
ENV TOOL_NAME="mkdocs"
ENV TOOL_VERSION="1.6.1"
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# TOOL PARAMETERS

RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION mkdocs=$TOOL_VERSION python=3.10 && \
	$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/python -m pip install mkdocs mkdocs-material pymdown-extensions plotly mkdocs-macros-plugin mkdocs-include-dir-to-nav mkdocs-include-markdown-plugin "markdown-exec[ansi]" && \
	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
	$MAMBA clean -y --all



###########
# RNASeQC #
###########


# TOOL INFO
ENV TOOL_NAME="rnaseqc"
ENV TOOL_VERSION="2.4.2"
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# TOOL PARAMETERS

RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION  -c bioconda rna-seqc=$TOOL_VERSION bx-python=0.14.0 pandas=2.3.3 numpy=2.2.6 ucsc-genepredtogtf=482-0 python=3.10 && \
	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
	$MAMBA clean -y --all



#########
# STARK #
#########


# TOOL INFO
ENV TOOL_NAME="stark"
ENV TOOL_VERSION="19.0.1-devel"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS
ENV TOOL="/tool"
ENV CONFIG_MYAPPS_FOLDER="$CONFIGS/myapps"
ENV CONFIG_HOWARD_FOLDER="$CONFIGS/howard"
ENV CONFIG_VARIANTCONVERT_FOLDER="$CONFIGS/variantconvert"


ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

COPY bin $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
COPY config $TOOLS/$TOOL_NAME/$TOOL_VERSION/config
COPY docs $TOOLS/$TOOL_NAME/$TOOL_VERSION/docs
COPY toolbox $TOOLS/$TOOL_NAME/$TOOL_VERSION/toolbox
COPY .env $TOOLS/$TOOL_NAME/$TOOL_VERSION/
COPY docker-compose.yml $TOOLS/$TOOL_NAME/$TOOL_VERSION/
COPY Dockerfile $TOOLS/$TOOL_NAME/$TOOL_VERSION/
COPY mkdocs.yml $TOOLS/$TOOL_NAME/$TOOL_VERSION/

# TOOL INSTALLATION
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
	ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ $TOOL && \
	# MYAPPS CONFIG FOLDER \
	mkdir -p $CONFIG_MYAPPS_FOLDER && \
	cp -R $TOOLS/$TOOL_NAME/$TOOL_VERSION/config/apps/myapps/* $CONFIG_MYAPPS_FOLDER/ && \
	rm -rf $TOOLS/$TOOL_NAME/$TOOL_VERSION/config/apps/myapps && \
	ln -sf $CONFIG_MYAPPS_FOLDER/ $TOOLS/$TOOL_NAME/$TOOL_VERSION/config/apps/myapps && \
	# HOWARD CONFIG FOLDER \
	mkdir -p $CONFIG_HOWARD_FOLDER && \
	cp -R $TOOLS/$TOOL_NAME/$TOOL_VERSION/config/howard/* $CONFIG_HOWARD_FOLDER && \
	rm -rf $TOOLS/$TOOL_NAME/$TOOL_VERSION/config/howard && \
	ln -sf $CONFIG_HOWARD_FOLDER $TOOLS/$TOOL_NAME/$TOOL_VERSION/config/howard && \
	# VARIANTCONVERT CONFIG FOLDER \
	mkdir -p $CONFIG_VARIANTCONVERT_FOLDER && \
	cp -R $CONFIGS_VARIANTCONVERT_FOLDER/* $CONFIG_VARIANTCONVERT_FOLDER/ ;



######################
# YUM REMOVE & CLEAR #
######################

# RUN echo "#[INFO] Cleaning" && \
# 	yum erase -y $YUM_REMOVE && \
# 	yum clean all && \
# 	rm -rf /var/cache/yum && \
# 	rm -rf $WORKDIR/* && \
# 	rm -rf /tmp/* && \
# 	if (($REMOVE_SOURCES)); then rm -rf $SOURCES; fi;



##############################
# WORKDIR / ENTRYPOINT / CMD #
##############################


WORKDIR "/STARK/data"

ENTRYPOINT [ "/tool/bin/STARK" ]
