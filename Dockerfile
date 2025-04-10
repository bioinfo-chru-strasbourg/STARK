
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
ENV WORKDIR="/tmp"
#ENV YUM_PARAM=" -q -e 0 "
ENV YUM_PARAM=" "
ENV WGET_PARAM=" "
ENV TAR_PARAM=" "
ENV ZIP_PARAM=" "
ENV MAKE_PARAM=" "



###########
# SOURCES #
###########

# Copy sources packages, scripts and tools
#ADD ./$SOURCES_FOLDER $SOURCES/sources
ADD ./$SOURCES_FOLDER $SOURCES
# COPY tools.json ${TOOLS}/tools.json
# COPY install_tools.sh ${TOOLS}/install_tools.sh


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
ENV PACKAGES_INSTALL="autoconf automake htop tree bc bzip2 bzip2-devel curl gcc gcc-c++ git make mlocate ncurses-devel tbb-devel unzip rsync wget which xz xz-devel zlib zlib-devel docker java-17 java-1.8.0 curl-devel openssl-devel diffutils parallel aria2 jq"
#ENV YUM_REMOVE="autoconf automake bzip2-devel lzma-devel ncurses-devel tbb-devel xz-devel zlib-devel zlib2-devel python3-devel curl-devel openssl-devel"

ENV PYTHON_MODULE=" pathos numpy scipy argparse"
ENV PERL_INSTALL=" perl perl-Switch perl-Time-HiRes perl-Data-Dumper perl-Digest-MD5 perl-Tk perl-devel"

ENV REPO_SYSTEM_GIT="$REPO/sources.system.tar.gz?path=sources/system"
ENV REPO_SYSTEM_HTTP="$REPO/sources/system/"

# ENV GET_TOOL_SOURCE=$SOURCES/$SOURCES_FOLDER/get_tool_source.sh
# ENV TOOL_INIT=$SOURCES/$SOURCES_FOLDER/tool_init.sh
# ENV TOOL_CHECK=$SOURCES/$SOURCES_FOLDER/tool_check.sh

ENV GET_TOOL_SOURCE=$SOURCES/get_tool_source.sh
ENV TOOL_INIT=$SOURCES/tool_init.sh
ENV TOOL_CHECK=$SOURCES/tool_check.sh


# RUN echo "#[INFO] SYSTEM Sources scripts" && \
# 	if [ -e $GET_TOOL_SOURCE ]; then \
# 	echo "#[INFO] GET TOOL SOURCE script exists" ; \
# 	elif $(wget --no-cache --progress=bar:force -nv --quiet "$REPO/$SOURCES_FOLDER/$(basename $GET_TOOL_SOURCE)" -O $GET_TOOL_SOURCE); then \
# 	echo "#[INFO] GET TOOL SOURCE script downloaded from REPO '$REPO/$SOURCES_FOLDER/$(basename $GET_TOOL_SOURCE)'" ; \
# 	else \
# 	mkdir -p $(dirname $GET_TOOL_SOURCE) ; \
# 	echo 'echo "#[INFO] TOOL source ($TOOL_SOURCE)" && \
# 	mkdir -p $(dirname $TOOL_SOURCE) && \
# 	if [ -e $TOOL_SOURCE ]; then \
# 	echo "#[INFO] TOOL TARBALL already in $TOOL_SOURCE"; \
# 	elif $(wget --no-cache --progress=bar:force "$TOOL_SOURCE_REPO" -O $TOOL_SOURCE); then \
# 	echo "#[INFO] TOOL TARBALL downloaded from STARK REPO $TOOL_SOURCE_REPO"; \
# 	if $(wget --no-cache --progress=bar:force -nv --quiet "$(dirname $TOOL_SOURCE_REPO)/source.info" -O $(dirname $TOOL_SOURCE)/source.info); then \
# 	echo "#[INFO] TOOL TARBALL external source information downloaded from STARK REPO $TOOL_SOURCE_REPO " ; \
# 	fi ; \
# 	elif $(wget --no-cache --progress=bar:force "$TOOL_SOURCE_EXTERNAL" -O $TOOL_SOURCE); then \
# 	echo "#[INFO] TOOL TARBALL downloaded from EXTERNAL SOURCE $TOOL_SOURCE_EXTERNAL"; \
# 	echo "$TOOL_SOURCE_EXTERNAL" > $(dirname $TOOL_SOURCE)/source.info; \
# 	else \
# 	echo "#[ERROR] TOOL TARBALL NOT FOUND"; \
# 	exit 1; \
# 	fi && \
# 	if [ -e $(dirname $TOOL_SOURCE)/source.info ]; then \
# 	echo "#[INFO] TOOL TARBALL external source: "$(cat $(dirname $TOOL_SOURCE)/source.info) ; \
# 	fi && \
# 	exit 0;' > $GET_TOOL_SOURCE ; \
# 	echo "#[INFO] GET TOOL SOURCE script written" ; \
# 	fi && \
# 	chmod u+x $GET_TOOL_SOURCE && \
# 	if [ -e $TOOL_INIT ]; then \
# 	echo "#[INFO] TOOL INIT script exists" ; \
# 	elif $(wget --no-cache --progress=bar:force -nv --quiet "$REPO/$SOURCES_FOLDER/$(basename $TOOL_INIT)" -O $TOOL_INIT); then \
# 	echo "#[INFO] TOOLS INIT script downloaded from REPO '$REPO/$SOURCES_FOLDER/$(basename $TOOL_INIT)'" ; \
# 	else \
# 	mkdir -p $(dirname $TOOL_INIT) ; \
# 	echo 'echo "#[INFO] TOOL $TOOL_NAME/$TOOL_VERSION" && \
# 	export TOOL_SOURCE=$SOURCES/$SOURCES_FOLDER/tools/$TOOL_NAME/$TOOL_VERSION/$TOOL_TARBALL && \
# 	export TOOL_SOURCE_REPO=$REPO/$SOURCES_FOLDER/tools/$TOOL_NAME/$TOOL_VERSION/$TOOL_TARBALL && \
# 	export TOOL_SOURCE_BUILD=$SOURCES/$SOURCES_FOLDER/tools/$TOOL_NAME/$TOOL_VERSION/build && \
# 	export TOOL_DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION && \
# 	export PATH=$TOOL_DEST/bin:$PATH && \
# 	# Get TOOL SOURCE && \
# 	$GET_TOOL_SOURCE $TOOL_SOURCE $TOOL_SOURCE_REPO $TOOL_SOURCE_EXTERNAL && \
# 	# TOOL folder preparation \
# 	echo "#[INFO] TOOL preparation" && \
# 	mkdir -p $TOOL_SOURCE_BUILD && \
# 	mkdir -p $TOOL_DEST/bin && \
# 	echo "#[INFO] TOOL release as current (forced)" && \
# 	ln -snf $TOOL_VERSION/ $TOOLS/$TOOL_NAME/previous && \
# 	if [ -e $TOOLS/$TOOL_NAME/current ]; then ln -snf $(basename $(realpath $TOOLS/$TOOL_NAME/current))/ $TOOLS/$TOOL_NAME/previous; fi && \
# 	ln -snf $TOOL_VERSION/ $TOOLS/$TOOL_NAME/current && \
# 	ln -snf $TOOL_VERSION/ $TOOLS/$TOOL_NAME/latest' > $TOOL_INIT ; \
# 	echo "#[INFO] TOOLS INIT script written" ; \
# 	fi && \
# 	chmod u+x $TOOL_INIT && \
# 	if [ -e $TOOL_CHECK ]; then \
# 	echo "#[INFO] TOOLS CHECK script exists" ; \
# 	elif $(wget --no-cache --progress=bar:force -nv --quiet "$REPO/$SOURCES_FOLDER/$(basename $TOOL_CHECK)" -O $TOOL_CHECK); then \
# 	echo "#[INFO] TOOLS CHECK script downloaded from REPO '$REPO/$SOURCES_FOLDER/$(basename $TOOL_CHECK)'" ; \
# 	else \
# 	mkdir -p $(dirname $TOOL_CHECK) ; \
# 	echo 'echo "#[INFO] TOOL cleaning" && \
# 	rm -rf $TOOL_SOURCE_BUILD && \
# 	if (($REMOVE_SOURCES)); then rm -rf $SOURCES/$SOURCES_FOLDER/tools/$TOOL_NAME/$TOOL_VERSION; fi && \
# 	echo "#[INFO] TOOL $TOOL_NAME/$TOOL_VERSION installed" ;' > $TOOL_CHECK ; \
# 	echo "#[INFO] TOOLS CHECK script written" ; \
# 	fi && \
# 	chmod u+x $TOOL_CHECK && \
# 	echo "#";

# Add tools and install script
#ADD install_system.sh ${TOOLS}/install_system.sh

# System installation
RUN echo "#[INFO] SYSTEM Packages installation" && \
	${SOURCES}/install_system.sh --yum_install="${PACKAGES_INSTALL}" --yum_param="${YUM_PARAM}"

# RUN echo "#[INFO] SYSTEM YUM installation - and download" && \
# 	# Create system repository \
# 	mkdir -p $SOURCES/$SOURCES_FOLDER/system && \
# 	mkdir -p $SOURCES/$SOURCES_FOLDER/system/$(uname -m) && \
# 	# INSTALL WGET \
# 	echo "#[INFO] System install wget package" && \
# 	#ls $SOURCES/$SOURCES_FOLDER/system/*.rpm && \
# 	if ! ls $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/wget-*.rpm 1> /dev/null 2>&1; then \
# 	echo "#[INFO] System wget package not locally available"; \
# 	yum $YUM_PARAM install -y --nogpgcheck --downloadonly --downloaddir=$SOURCES/$SOURCES_FOLDER/system/$(uname -m)/ wget; \
# 	echo "#[INFO] System wget package downloaded from YUM Repository"; \
# 	fi && \
# 	echo "#[INFO] System install rsync package" && \
# 	if ! ls $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/rsync-*.rpm 1> /dev/null 2>&1; then \
# 	echo "#[INFO] System rsync package not locally available"; \
# 	yum $YUM_PARAM install -y --nogpgcheck --downloadonly --downloaddir=$SOURCES/$SOURCES_FOLDER/system/$(uname -m)/ rsync; \
# 	echo "#[INFO] System rsync package downloaded from YUM Repository"; \
# 	fi && \
# 	# Install packages locally \
# 	echo "#[INFO] System packages installation locally" && \
# 	yum $YUM_PARAM localinstall -y --allowerasing --nogpgcheck $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/wget-*.rpm $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/rsync-*.rpm && \
# 	# Test WGET installation \
# 	if ! command -v wget 1>/dev/null 2>/dev/null; then \
# 	echo "#[ERROR] System wget package not installed (Please open Internet connexion or provide WGET rpm in sources/system folder)"; \
# 	exit 1; \
# 	fi && \
# 	if ! command -v rsync 1>/dev/null 2>/dev/null; then \
# 	echo "#[ERROR] System rsync package not installed (Please open Internet connexion or provide RSYNC rpm in sources/system folder)"; \
# 	exit 1; \
# 	fi && \
# 	# DOWNLOAD packages from repository \
# 	echo "#[INFO] System packages download from REPO '$REPO'"; \
# 	mkdir -p $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build && \
# 	# in GIT mode
# 	if wget -q --progress=bar:force --tries=3 $REPO_SYSTEM_GIT -O $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/STARK-repo.sources.system.tar.gz; then \
# 	if tar xf $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/STARK-repo.sources.system.tar.gz -C $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/; then \
# 	rsync -auczqAXhi --no-links --no-perms --no-owner --no-group --ignore-missing-args $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/STARK-repo.sources.system*/sources/system/*rpm $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/; \
# 	echo "#[INFO] System packages downloaded from REPO '$REPO' (GIT)"; \
# 	else \
# 	echo "#[WARNING] System fail to uncompress packages from REPO '$REPO'"; \
# 	fi; \
# 	# in HTTP mode
# 	elif wget -q --progress=bar:force --tries=3 -r --no-parent $REPO_SYSTEM_HTTP -x --directory-prefix=$SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/STARK-repo.sources.system/; then \
# 	rsync -auczqAXhi --no-links --no-perms --no-owner --no-group --ignore-missing-args $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/STARK-repo.sources.system/*/sources/system/*rpm $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/; \
# 	echo "#[INFO] System packages downloaded from REPO '$REPO' (FTP/HTTP)"; \
# 	else \
# 	echo "#[WARNING] System fail packages download from REPO '$REPO'"; \
# 	fi && \
# 	rm -rf $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build && \
# 	# Install packages locally \
# 	echo "#[INFO] System packages installation locally" && \
# 	if ! ls $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/*.rpm 1> /dev/null 2>&1; then \
# 	yum $YUM_PARAM localinstall -y --nogpgcheck $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/*.rpm; \
# 	echo "#[INFO] System packages installation locally done."; \
# 	fi && \
# 	#yum $YUM_PARAM localinstall -y --nogpgcheck $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/*.rpm && \
# 	# Install EPEL Repository \
# 	echo "#[INFO] System EPEL Repository package" && \
# 	if ! ls $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/epel-release-*.rpm 1> /dev/null 2>&1; then \
# 	yum $YUM_PARAM install -y --nogpgcheck --downloadonly --downloaddir=$SOURCES/$SOURCES_FOLDER/system/$(uname -m)/ epel-release; \
# 	echo "#[INFO] System EPEL Repository package downloaded from YUM repository"; \
# 	fi && \
# 	if ls $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/epel-release-*.rpm 1> /dev/null 2>&1; then \
# 	yum $YUM_PARAM localinstall -y --nogpgcheck $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/epel-release-*.rpm; \
# 	echo "#[INFO] System EPEL Repository package enabled"; \
# 	else \
# 	echo "#[WARNING] System fail enable EPEL Repository"; \
# 	fi && \
# 	# Update YUM \
# 	echo "#[INFO] System packages update from YUM Repository" && \
# 	mkdir -p $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/update && \
# 	yum $YUM_PARAM update -y --downloadonly --downloaddir=$SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/update && \
# 	yum $YUM_PARAM localinstall -y --nogpgcheck $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/update/*.rpm && \
# 	rsync -auczqAXhi --no-links --no-perms --no-owner --no-group --ignore-missing-args $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/update/*rpm $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/ && \
# 	echo "#[INFO] System packages downloaded & updated from YUM Repository" && \
# 	echo "#[INFO] System packages install from YUM Repository" && \
# 	mkdir -p $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/install && \
# 	yum $YUM_PARAM install -y --downloadonly --allowerasing --downloaddir=$SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/install/ $YUM_INSTALL && \
# 	ls -lah $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/install/ && \
# 	yum $YUM_PARAM localinstall -y --nogpgcheck --allowerasing $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/install/*.rpm && \
# 	rsync -auczqAXhi --no-links --no-perms --no-owner --no-group --ignore-missing-args $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/install/*rpm $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/ && \
# 	echo "#[INFO] System packages downloaded & installed from YUM Repository" && \
# 	rm -rf $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build && \
# 	yum clean -y all && \
# 	rm -rf /var/cache/yum && \
# 	echo "#[INFO] System Clean" && \
# 	echo "#[INFO] SYSTEM Bashrc" && \
# 	echo "alias ll='ll -lah'" >> ~/.bashrc && \
# 	echo "#"
	


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
    #wget $TARBALL_LOCATION/Miniforge-pypy3-$TOOL_VERSION-$(uname)-$(uname -m).sh -O $TARBALL && \
	echo $TARBALL_LOCATION/Miniforge3-$TOOL_VERSION-$(uname)-$(uname -m).sh && \
    wget $TARBALL_LOCATION/Miniforge3-$TOOL_VERSION-$(uname)-$(uname -m).sh -O $TARBALL && \
    bash $TARBALL -b -p $DEST && \
    rm -f $TARBALL && \
	find ${DEST} -follow -type f -name '*.a' -or -name '*.pyc' -delete && \
	$MAMBA clean --force-pkgs-dirs --all --yes && \
	#find ${DEST} -follow -ignore_readdir_race \( -name '*.a' -o -name '*.pyc' -o -name '*.txt' -o -name '*.md' -o -name '*.pdf' -o  -name '__pycache__' \) -exec rm -rf {} + && \
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
	find ${DEST} -follow -ignore_readdir_race \( -name '*.a' -o -name '*.pyc' -o -name '*.txt' -o -name '*.md' -o -name '*.pdf' -o  -name '__pycache__' \) -exec rm -rf {} +
	



########
# PERL #
########



# PERL installation
RUN	echo "#[INFO] SYSTEM Perl installation - download from yum" && \
	mkdir -p $SOURCES/$SOURCES_FOLDER/perl/build/install && \
	echo "FROM_IMAGE=$FROM_IMAGE" && \
	# if [ "$FROM_IMAGE" == "almalinux:9" ]; then \
	# 	yum config-manager --set-enabled crb; \
	# else \
	# 	dnf config-manager --set-enabled powertools; \
	# fi && \
	#dnf config-manager --set-enabled powertools && \
	#yum config-manager --set-enabled crb && \
	#yum $YUM_PARAM install -y --downloadonly --downloaddir=$SOURCES/$SOURCES_FOLDER/perl/build/install --enablerepo=powertools $PERL_INSTALL && \
	yum $YUM_PARAM install -y --downloadonly --downloaddir=$SOURCES/$SOURCES_FOLDER/perl/build/install $PERL_INSTALL && \
	yum $YUM_PARAM localinstall -y --nogpgcheck $SOURCES/$SOURCES_FOLDER/perl/build/install/*.rpm && \
	rsync -auczqAXhi --no-links --no-perms --no-owner --no-group --ignore-missing-args $SOURCES/$SOURCES_FOLDER/perl/build/install/*rpm $SOURCES/$SOURCES_FOLDER/system/ && \
	rm -rf $SOURCES/$SOURCES_FOLDER/perl/build && \
	yum clean -y all && \
	rm -rf /var/cache/yum && \
	echo "#[INFO] System Clean" && \
	echo "#";


# RUN	echo "#[INFO] SYSTEM Perl packages installation - download from yum" && \
# 	mkdir -p $SOURCES/$SOURCES_FOLDER/perl/build/install && \
# 	dnf install -y dnf-plugins-core

# # PERL installation
# RUN	echo "#[INFO] SYSTEM Perl packages installation - download from yum" && \
# 	mkdir -p $SOURCES/$SOURCES_FOLDER/perl/build/install && \
# 	yum $YUM_PARAM install -y --downloadonly --downloaddir=$SOURCES/$SOURCES_FOLDER/perl/build/install $PERL_INSTALL && \
# 	yum $YUM_PARAM localinstall -y --nogpgcheck $SOURCES/$SOURCES_FOLDER/perl/build/install/*.rpm && \
# 	rsync -auczqAXhi --no-links --no-perms --no-owner --no-group --ignore-missing-args $SOURCES/$SOURCES_FOLDER/perl/build/install/*rpm $SOURCES/$SOURCES_FOLDER/system/ && \
# 	rm -rf $SOURCES/$SOURCES_FOLDER/perl/build && \
# 	yum clean -y all && \
# 	rm -rf /var/cache/yum && \
# 	echo "#[INFO] System Clean" && \
# 	echo "#";


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
ENV TOOL_VERSION="17"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
RUN echo "#[INFO] SYSTEM Java installation '$TOOL_NAME:$TOOL_VERSION'" && \
	mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
	ln -s /usr/lib/jvm/jre-17/bin/java $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/java && \
	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;



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



##########
# JAVA11 #
##########

# ENV TOOL_NAME="java"
# ENV TOOL_VERSION="11"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
# 	ln -s /usr/lib/jvm/jre-11/bin/java $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/java ;




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



# RUN echo "#[INFO] TOOLS installation shared with Mamba" && \
# 	$MAMBA create -y -p $DEST \
# 	-c bioconda -c compbiocore -c conda-forge \
#     htslib=1.21 \
#     bcftools=1.21 \
#     bedtools=2.31.1 \
#     bowtie2=2.5.1 \
#     bwa=0.7.18 \
#     bwa-mem2=2.2.1 \
#     fastp=0.23.2 \
#     gatk4=4.6.1.0 \
#     gatk=3.8 \
#     igvtools=2.17.3 \
#     mutect=1.1.6 \
#     picard=3.3.0 \
#     samtools=1.21 \
#     snpeff=5.1d \
#     star=2.7.11b \
#     #star-fusion=1.14.0 \
#     arriba=2.4.0 \
#     varscan=2.4.6 \
#     fgbio=2.4.0 && \
# 	$MAMBA clean -y --all && \
# 	find ${DEST} -follow -ignore_readdir_race \( -name '*.a' -o -name '*.pyc' -o -name '*.txt' -o -name '*.md' -o -name '*.pdf' -o  -name '__pycache__' \) -exec rm -rf {} +


#RUN du -h -d1 $DEST && du -h -d1 $DEST/share && truc


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




##########
# HTSLIB #
##########

# # TOOL INFO
# ENV TOOL_NAME="htslib"
# ENV TOOL_VERSION="1.18"
# ENV TOOL_TARBALL="$TOOL_NAME-$TOOL_VERSION.tar.bz2"
# ENV TOOL_SOURCE_EXTERNAL="https://github.com/samtools/$TOOL_NAME/releases/download/$TOOL_VERSION/$TOOL_TARBALL"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS

# # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	source $TOOL_INIT && \
# 	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
# 	make install --quiet -j $THREADS -C $(ls -d $TOOL_SOURCE_BUILD/*) prefix=$TOOL_DEST && \
# 	$TOOL_CHECK ;


# ENV TOOL_NAME="htslib"
# ENV TOOL_VERSION="1.21"
# ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
# ENV PATH=$PATH:$DEST/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p ${DEST} -c bioconda $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	find ${DEST} -follow -type f -name '*.a' -or -name '*.pyc' -delete && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;




############
# BCFTOOLS #
############

# TOOL INFO
#ENV TOOL_NAME="bcftools"
#ENV TOOL_VERSION="1.15.1"
#ENV TOOL_TARBALL="$TOOL_NAME-$TOOL_VERSION.tar.bz2"
#ENV TOOL_SOURCE_EXTERNAL="https://github.com/samtools/$TOOL_NAME/releases/download/$TOOL_VERSION/$TOOL_TARBALL"
#ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS

# TOOL INSTALLATION
#RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
#	source $TOOL_INIT && \
#	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
#	make install --quiet -j $THREADS -C $(ls -d $TOOL_SOURCE_BUILD/*) prefix=$TOOL_DEST && \
#	$TOOL_CHECK ;


# ENV TOOL_NAME="bcftools"
# ENV TOOL_VERSION="1.21"
# ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
# ENV PATH=$PATH:$DEST/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p ${DEST} -c bioconda $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	find ${DEST} -follow -type f -name '*.a' -or -name '*.pyc' -delete && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;


# find .  -follow \( -name '*.a' -o  -name '*.pyc' -o -name '__pycache__' -o -name '*.txt' -o -name '*.md' -o -name '*.pdf' \) -exec rm -rf {} +

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



############
# BEDTOOLS #
############

# # TOOL INFO
# ENV TOOL_NAME="bedtools"
# ENV TOOL_VERSION="2.31.0"
# ENV TOOL_TARBALL="$TOOL_NAME-$TOOL_VERSION.tar.gz"
# ENV TOOL_SOURCE_EXTERNAL="https://github.com/arq5x/bedtools2/releases/download/v$TOOL_VERSION/$TOOL_TARBALL"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS

# # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	source $TOOL_INIT && \
# 	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
# 	make install --quiet -j $THREADS -C $(ls -d $TOOL_SOURCE_BUILD/*) prefix=$TOOL_DEST && \
# 	$TOOL_CHECK ;


# ENV TOOL_NAME="bedtools"
# ENV TOOL_VERSION="2.31.1"
# ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION -c bioconda $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;


############
# BOWTIE2 #
############

# # TOOL INFO
# ENV TOOL_NAME="bowtie2"
# ENV TOOL_VERSION="2.5.1"
# ENV TOOL_TARBALL="$TOOL_NAME-$TOOL_VERSION-linux-x86_64.zip"
# ENV TOOL_SOURCE_EXTERNAL="https://github.com/BenLangmead/bowtie2/releases/download/v$TOOL_VERSION/$TOOL_TARBALL"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS

# # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	source $TOOL_INIT && \
# 	unzip -q $TOOL_SOURCE -d $TOOL_SOURCE_BUILD && \
# 	cp $TOOL_SOURCE_BUILD/*/bowtie2* $TOOL_DEST/bin/ && \
# 	rm -f $TOOL_DEST/bin/*debug && \
# 	$TOOL_CHECK ;


# ENV TOOL_NAME="bowtie2"
# ENV TOOL_VERSION="2.5.1"
# ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION -c bioconda $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;


#######
# BWA #
#######

# # TOOL INFO
# ENV TOOL_NAME="bwa"
# ENV TOOL_VERSION="0.7.17"
# ENV TOOL_TARBALL="$TOOL_NAME-$TOOL_VERSION.tar.bz2"
# ENV TOOL_SOURCE_EXTERNAL="https://github.com/lh3/bwa/releases/download/v$TOOL_VERSION/bwa-$TOOL_VERSION.tar.bz2"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS

# # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	source $TOOL_INIT && \
# 	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
# 	make --quiet -j $THREADS -C $(ls -d $TOOL_SOURCE_BUILD/*) && \
# 	cp $TOOL_SOURCE_BUILD/*/bwa $TOOL_DEST/bin/ && \
# 	$TOOL_CHECK ;


# # TOOL INFO
# ENV TOOL_NAME="bwa"
# ENV TOOL_VERSION="0.7.18"
# ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION -c bioconda $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;


########
# BWA2 #
########

# # TOOL INFO
# ENV TOOL_NAME="bwa"
# ENV TOOL_VERSION="2.2.1"
# ENV TOOL_TARBALL=$TOOL_NAME"-mem2-"$TOOL_VERSION"_x64-linux.tar.bz2"
# ENV TOOL_SOURCE_EXTERNAL="https://github.com/bwa-mem2/bwa-mem2/releases/download/v$TOOL_VERSION/$TOOL_TARBALL"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS

# # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	source $TOOL_INIT && \
# 	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
# 	cp $TOOL_SOURCE_BUILD/*/bwa-mem2* $TOOL_DEST/bin/ && \
# 	$TOOL_CHECK ;


# # TOOL INFO
# ENV TOOL_NAME="bwa"
# ENV TOOL_VERSION="2.2.1"
# ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION -c bioconda bwa-mem2~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all;
# 	#ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;



#########
# FASTP #
#########

# TOOL INFO
#ENV TOOL_NAME="fastp"
#ENV TOOL_VERSION="0.23.2"
#ENV TOOL_TARBALL="$TOOL_NAME"
#ENV TOOL_SOURCE_EXTERNAL="http://opengene.org/$TOOL_NAME/$TOOL_NAME.$TOOL_VERSION"
#ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS

# TOOL INSTALLATION
#RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
#	source $TOOL_INIT && \
#	cp $TOOL_SOURCE $TOOL_DEST/bin/ && \
#	chmod a+x $TOOL_DEST/bin/* && \
#	$TOOL_CHECK ;

# # TOOL INFO
# ENV TOOL_NAME="fastp"
# ENV TOOL_VERSION="0.23.2"
# # ENV TOOL_TARBALL="$TOOL_NAME"
# # ENV TOOL_SOURCE_EXTERNAL="http://opengene.org/$TOOL_NAME/$TOOL_NAME.$TOOL_VERSION"
# ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION -c bioconda $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;



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

# # TOOL INFO
# ENV TOOL_NAME="gatk"
# ENV TOOL_VERSION="4.4.0.0"
# ENV TOOL_TARBALL="gatk-$TOOL_VERSION.zip"
# ENV TOOL_SOURCE_EXTERNAL="https://github.com/broadinstitute/gatk/releases/download/$TOOL_VERSION/$TOOL_TARBALL"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS
# ENV TOOL_JAR=gatk-package-$TOOL_VERSION-local.jar

# # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	source $TOOL_INIT && \
# 	unzip -q $TOOL_SOURCE -d $TOOL_SOURCE_BUILD && \
# 	cp -R $TOOL_SOURCE_BUILD/gatk-$TOOL_VERSION/* $TOOL_DEST/bin/ && \
# 	$TOOL_CHECK ;


# # TOOL INFO
# ENV TOOL_NAME="gatk4"
# ENV TOOL_VERSION="4.6.1.0"
# ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION -c bioconda $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;



########
# GATK #
########

# # TOOL INFO
# ENV TOOL_NAME="gatk"
# ENV TOOL_VERSION="3.8-1-0"
# ENV TOOL_TARBALL="GenomeAnalysisTK-$TOOL_VERSION.tar.bz2"
# ENV TOOL_SOURCE_EXTERNAL="https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=$TOOL_VERSION-gf15c1c3ef"
# https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS
# ENV TOOL_JAR=GenomeAnalysisTK.jar

# # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	source $TOOL_INIT && \
# 	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
# 	cp -R $TOOL_SOURCE_BUILD/*/$TOOL_JAR $TOOL_DEST/bin/ && \
# 	$TOOL_CHECK ;


# # TOOL INFO
# ENV TOOL_NAME="gatk"
# ENV TOOL_VERSION="3.8"
# ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION -c bioconda $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;



##########
# HOWARD #
##########

# TOOL INFO
ENV TOOL_NAME="howard"
ENV TOOL_VERSION="0.9.15.6"
ENV TOOL_TARBALL="$TOOL_VERSION.tar.gz"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/bioinfo-chru-strasbourg/howard/archive/refs/heads/$TOOL_TARBALL"
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# TOOL PARAMETERS
ENV TOOL_PARAM_DATABASE_FOLDER_LINK=$DATABASES
ENV TOOL_PARAM_DATABASE_FOLDER=/databases


# TOOL INSTALLATION
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	source $TOOL_INIT && \
	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
	cp -R $TOOL_SOURCE_BUILD/*/* $TOOL_DEST/ && \
	chmod a+x $TOOL_DEST/* -R && \
	mkdir -p $TOOL_PARAM_DATABASE_FOLDER_LINK && \
	mkdir -p $TOOL_PARAM_DATABASE_FOLDER && \
	ln -s $DATABASES $TOOL_DATABASE_FOLDER && \
	$TOOL_CHECK ;



############
# HOWARD 2 #
############

# https://github.com/bioinfo-chru-strasbourg/howard/archive/refs/heads/devel.zip

# TOOL INFO
ENV TOOL_NAME="howard"
ENV TOOL_VERSION="devel"
ENV TOOL_TARBALL="$TOOL_VERSION.zip"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/bioinfo-chru-strasbourg/howard/archive/refs/heads/$TOOL_TARBALL"
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# TOOL PARAMETERS

# # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	source $TOOL_INIT && \
# 	unzip -q $TOOL_SOURCE -d $TOOL_SOURCE_BUILD && \
# 	cp -R $TOOL_SOURCE_BUILD/*/* $TOOL_DEST/ && \
# 	cd $TOOL_DEST/ && \
# 	$PYTHON -m pip install -e . && \
#     $TOOL_CHECK ;


RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION PYTHON=3.10 && \
	source $TOOL_INIT && \
	unzip -q $TOOL_SOURCE -d $TOOL_SOURCE_BUILD && \
	cp -R $TOOL_SOURCE_BUILD/*/* $TOOL_DEST/ && \
	$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/python -m pip install -e $TOOL_DEST && \
	$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/python -m pip install polars-lts-cpu && \
	#ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
	$MAMBA clean -y --all
	# echo ls -lah $TOOLS/$TOOL_NAME/ && \
	# ls -lah $TOOLS/$TOOL_NAME/ && \
	# echo ls -lah $TOOLS/$TOOL_NAME/current/ && \
	# ls -lah $TOOLS/$TOOL_NAME/current/ && \
	# echo ls -lah $TOOLS/$TOOL_NAME/current/bin/ && \
	# ls -lah $TOOLS/$TOOL_NAME/current/bin/ && \
	# whereis howard && \
	# howard --help
	
############
# IGVTOOLS #
############

# # TOOL INFO
# ENV TOOL_NAME="igvtools"
# ENV TOOL_VERSION="2.16.1"
# ENV TOOL_TARBALL="IGV_$TOOL_VERSION.zip"
# ENV TOOL_VERSION_MAIN="2.16"
# ENV TOOL_SOURCE_EXTERNAL="https://data.broadinstitute.org/igv/projects/downloads/$TOOL_VERSION_MAIN/$TOOL_TARBALL"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS

# # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	source $TOOL_INIT && \
# 	unzip -q $TOOL_SOURCE -d $TOOL_SOURCE_BUILD && \
# 	cp -R $TOOL_SOURCE_BUILD/*/* $TOOL_DEST/bin/ && \
# 	$TOOL_CHECK ;


# # TOOL INFO
# ENV TOOL_NAME="igvtools"
# ENV TOOL_VERSION="2.17.3"
# ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION -c bioconda $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;



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

# # TOOL INFO
# ENV TOOL_NAME="mutect"
# ENV TOOL_VERSION="1.1.7"
# ENV TOOL_TARBALL="$TOOL_NAME-$TOOL_VERSION.jar.zip"
# ENV TOOL_SOURCE_EXTERNAL="https://software.broadinstitute.org/gatk/download/auth?package=M1"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS
# ENV TARBALL_JAR=mutect-$TOOL_VERSION.jar
# ENV TOOL_JAR=mutect.jar

# # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	source $TOOL_INIT && \
# 	unzip -q $TOOL_SOURCE -d $TOOL_DEST/bin/ && \
# 	mv $TOOL_DEST/bin/$TARBALL_JAR $TOOL_DEST/bin/$TOOL_JAR && \
# 	$TOOL_CHECK ;


# # TOOL INFO - deprecated!!!
# ENV TOOL_NAME="mutect"
# ENV TOOL_VERSION="1.1.6"
# ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION -c compbiocore $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;



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

# # TOOL INFO
# ENV TOOL_NAME="picard"
# ENV TOOL_VERSION="3.0.0"
# ENV TOOL_TARBALL="picard.jar"
# ENV TOOL_SOURCE_EXTERNAL="https://github.com/broadinstitute/picard/releases/download/$TOOL_VERSION/picard.jar"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS

# # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	source $TOOL_INIT && \
# 	cp $TOOL_SOURCE $TOOL_DEST/bin/ && \
# 	$TOOL_CHECK ;

# # TOOL INFO
# ENV TOOL_NAME="picard"
# ENV TOOL_VERSION="3.3.0"
# ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION -c bioconda $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;




############
# SAMTOOLS #
############

# # TOOL INFO
# ENV TOOL_NAME="samtools"
# ENV TOOL_VERSION="1.18"
# ENV TOOL_TARBALL="$TOOL_NAME-$TOOL_VERSION.tar.bz2"
# ENV TOOL_SOURCE_EXTERNAL="https://github.com/samtools/$TOOL_NAME/releases/download/$TOOL_VERSION/$TOOL_TARBALL"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS

# # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	source $TOOL_INIT && \
# 	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
# 	make install --quiet -j $THREADS -C $(ls -d $TOOL_SOURCE_BUILD/*) prefix=$TOOL_DEST && \
# 	$TOOL_CHECK ;


# ENV TOOL_NAME="samtools"
# ENV TOOL_VERSION="1.21"
# ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION -c bioconda $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;


##########
# SNPEFF #
##########

# # TOOL INFO
# # Beware of TARBALL release
# ENV TOOL_NAME="snpeff"
# ENV TOOL_VERSION="5.1d"
# #ENV TOOL_TARBALL="snpEff_latest_core.zip"
# ENV TOOL_TARBALL="snpEff_v5_1d_core.zip"
# ENV TOOL_SOURCE_EXTERNAL="https://snpeff.blob.core.windows.net/versions/$TOOL_TARBALL"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS
# ENV TOOL_PARAM_DATABASE_FOLDER_LINK=$DATABASES/snpeff/$TOOL_VERSION
# ENV TOOL_PARAM_DATABASE_FOLDER=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/data

# # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	source $TOOL_INIT && \
# 	unzip -q $TOOL_SOURCE -d $TOOL_SOURCE_BUILD && \
# 	cp $TOOL_SOURCE_BUILD/*/*jar $TOOL_DEST/bin/ && \
# 	cp $TOOL_SOURCE_BUILD/*/*config $TOOL_DEST/bin/ && \
# 	mkdir -p $TOOL_PARAM_DATABASE_FOLDER_LINK && \
# 	ln -snf $TOOL_PARAM_DATABASE_FOLDER_LINK/ $TOOL_PARAM_DATABASE_FOLDER && \
# 	$TOOL_CHECK ;


# ENV TOOL_NAME="snpeff"
# ENV TOOL_VERSION="5.1d"
# ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION -c bioconda $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
# 	echo ln -s  $(find $TOOLS/$TOOL_NAME/$TOOL_VERSION -name "snpEff.jar" | head -n 1) $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/$(basename $(find $TOOLS/$TOOL_NAME/$TOOL_VERSION -name "snpEff.jar" | head -n 1))



#############
# GENCORE   #
#############
# https://github.com/OpenGene/gencore

# # TOOL INFO
# ENV TOOL_NAME="gencore"
# ENV TOOL_VERSION="0.17.2"
# ENV TOOL_TARBALL="$TOOL_NAME"
# ENV TOOL_SOURCE_EXTERNAL="http://opengene.org/$TOOL_NAME/$TOOL_NAME"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS

# # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	source $TOOL_INIT && \
# 	cp $TOOL_SOURCE $TOOL_DEST/bin/ && \
# 	chmod a+x $TOOL_DEST/bin/* && \
# 	$TOOL_CHECK ;

# ENV TOOL_NAME="gencore"
# ENV TOOL_VERSION="0.17.2"
# ENV PATH=$PYTHON_ENV/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA install -y -p $PYTHON_ENV -c bioconda $TOOL_NAME=$TOOL_VERSION




###########
# VARSCAN #
###########

# # TOOL INFO
# ENV TOOL_NAME="varscan"
# ENV TOOL_VERSION="2.4.6"
# ENV TOOL_TARBALL="VarScan.v$TOOL_VERSION.jar"
# ENV TOOL_SOURCE_EXTERNAL="https://github.com/dkoboldt/varscan/raw/master/$TOOL_TARBALL"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS
# ENV TOOL_PARAM_JAR_NAME=VarScan.jar

# # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	source $TOOL_INIT && \
# 	cp $TOOL_SOURCE $TOOL_DEST/bin/ && \
# 	ln -s $(basename $TOOL_SOURCE) $TOOL_DEST/bin/$TOOL_PARAM_JAR_NAME && \
# 	$TOOL_CHECK ;


# ENV TOOL_NAME="varscan"
# ENV TOOL_VERSION="2.4.6"
# ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION -c bioconda $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;




##########
# FGBIO #
##########

# # # TOOL INFO
#  ENV TOOL_NAME="fgbio"
#  ENV TOOL_VERSION="2.1.0"
#  ENV TOOL_TARBALL="fgbio.jar"
#  ENV TOOL_SOURCE_EXTERNAL="https://github.com/fulcrumgenomics/fgbio/releases/download/$TOOL_VERSION/fgbio-$TOOL_VERSION.jar"
#  ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # # TOOL PARAMETERS

# # # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
#  	source $TOOL_INIT && \
#  	cp $TOOL_SOURCE $TOOL_DEST/bin/ && \
#  	$TOOL_CHECK ;


# ENV TOOL_NAME="fgbio"
# ENV TOOL_VERSION="2.4.0"
# ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION -c bioconda $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;




########
# STAR #
########

# TOOL INFO
#ENV TOOL_NAME="STAR"
#ENV TOOL_VERSION="2.7.8a"
#ENV TOOL_TARBALL="$TOOL_VERSION.zip"
#ENV TOOL_SOURCE_EXTERNAL="https://github.com/alexdobin/$TOOL_NAME/archive/refs/tags/$TOOL_TARBALL"
#ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS
# TOOL INSTALLATION
#RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
#	source $TOOL_INIT && \
#	unzip -q $TOOL_SOURCE -d $TOOL_SOURCE_BUILD && \
#	ls && \
#	cd $TOOL_SOURCE_BUILD/$TOOL_NAME-$TOOL_VERSION/source && \
#	make STAR && \
#	cp STAR $TOOL_DEST/bin/ && \
#	$TOOL_CHECK ;


# ENV TOOL_NAME="star"
# ENV TOOL_VERSION="2.7.11b"
# ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION -c bioconda $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;




###############
# STAR FUSION #
###############
# Depends on STAR v2.7.8a

# # TOOL INFO
# ENV TOOL_NAME="STAR-Fusion"
# ENV TOOL_VERSION="1.12.0"
# ENV TOOL_TARBALL="$TOOL_NAME-v$TOOL_VERSION.FULL.tar.gz"
# ENV TOOL_SOURCE_EXTERNAL="https://github.com/$TOOL_NAME/$TOOL_NAME/releases/download/$TOOL_NAME-v$TOOL_VERSION/$TOOL_TARBALL"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS
# # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	source $TOOL_INIT && \
# 	tar xvf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
# 	cd  $TOOL_SOURCE_BUILD/$TOOL_NAME-v$TOOL_VERSION/ && \
# 	make && \
# 	cp -r * $TOOL_DEST/bin/ && \
# 	$TOOL_CHECK ;

# FusionInspector & Interval Tree included with STAR-Fusion

# ENV TOOL_NAME="star-fusion"
# ENV TOOL_VERSION="1.14.0"
# ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION -c bioconda $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;



##########
# ARRIBA #
##########

# # TOOL INFO
# ENV TOOL_NAME="arriba"
# ENV TOOL_VERSION="2.4.0"
# ENV TOOL_TARBALL=$TOOL_NAME"_v"$TOOL_VERSION".tar.gz"  
# ENV TOOL_SOURCE_EXTERNAL="https://github.com/suhrig/$TOOL_NAME/releases/download/v$TOOL_VERSION/$TOOL_TARBALL"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS

# # TOOL INSTALLATION
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	source $TOOL_INIT && \
# 	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
# 	cp $TOOL_SOURCE_BUILD/$TOOL_NAME"_v"$TOOL_VERSION/$TOOL_NAME $TOOL_DEST/bin/ && \
# 	chmod a+x $TOOL_DEST/bin/* && \
# 	mkdir -p $TOOL_DEST/scripts/ && \
# 	cp $TOOL_SOURCE_BUILD/$TOOL_NAME"_v"$TOOL_VERSION/scripts/*.sh $TOOL_DEST/scripts/ && \
# 	chmod a+x $TOOL_DEST/scripts/*.sh && \
# 	$TOOL_CHECK ;



# ENV TOOL_NAME="arriba"
# ENV TOOL_VERSION="2.4.0"
# ENV PATH=$PATH:$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION -c bioconda $TOOL_NAME~=$TOOL_VERSION && \
# 	$MAMBA clean -y --all && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;





##################
# variantconvert #
##################

# # INFO
# ENV TOOL_NAME="variantconvert"
# ENV TOOL_VERSION="pypi"
# # INSTALLATION
# ENV INSTALL_DIR=/usr/lib
# WORKDIR $INSTALL_DIR
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	git clone https://github.com/SamuelNicaise/variantconvert.git && cd $TOOL_NAME && git fetch && git checkout $TOOL_VERSION && $PIP install -e . && $PIP cache purge && \
# 	$PYTHON $INSTALL_DIR/$TOOL_NAME/src/$TOOL_NAME/__main__.py init && \
# 	$PYTHON $INSTALL_DIR/$TOOL_NAME/src/$TOOL_NAME/__main__.py config --set GENOME.path=/STARK/databases/genomes/current/hg19/hg19.fa --configFiles hg19/*.json
# WORKDIR $WORKDIR


##################
# variantconvert #
##################

# # INFO
# ENV TOOL_NAME="variantconvert"
# ENV TOOL_VERSION="pypi"
# #https://github.com/SamuelNicaise/variantconvert/archive/refs/tags/2.0.1.tar.gz
# # INSTALLATION
# ENV INSTALL_DIR=/usr/lib
# WORKDIR $INSTALL_DIR
# RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
# 	git clone https://github.com/SamuelNicaise/variantconvert.git && cd $TOOL_NAME && git fetch && git checkout $TOOL_VERSION && $PIP install -e . && $PIP cache purge && \
# 	$PYTHON $INSTALL_DIR/$TOOL_NAME/src/$TOOL_NAME/__main__.py init && \
# 	$PYTHON $INSTALL_DIR/$TOOL_NAME/src/$TOOL_NAME/__main__.py config --set GENOME.path=/STARK/databases/genomes/current/hg19/hg19.fa --configFiles hg19/*.json
# WORKDIR $WORKDIR

# TOOL INFO
ENV TOOL_NAME="variantconvert"
ENV TOOL_VERSION="2.0.1"
ENV TOOL_TARBALL=$TOOL_VERSION".tar.gz"  
ENV TOOL_SOURCE_EXTERNAL="https://github.com/SamuelNicaise/$TOOL_NAME/archive/refs/tags/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	$MAMBA create -y -p $TOOLS/$TOOL_NAME/$TOOL_VERSION PYTHON=3.10 && \
	source $TOOL_INIT && \
	tar -xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
	ls -l $TOOL_SOURCE_BUILD && \
	cp -R $TOOL_SOURCE_BUILD/*/* $TOOL_DEST/ && \
	$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/python -m pip install -e $TOOL_DEST && \
	#ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
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
ENV CONFIG_MYAPPS_FOLDER="$STARK_FOLDER/config/myapps"
ENV CONFIG_HOWARD_FOLDER="$STARK_FOLDER/config/howard"


ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

COPY bin $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
COPY config $TOOLS/$TOOL_NAME/$TOOL_VERSION/config
COPY docs $TOOLS/$TOOL_NAME/$TOOL_VERSION/docs
COPY toolbox $TOOLS/$TOOL_NAME/$TOOL_VERSION/toolbox
COPY .env $TOOLS/$TOOL_NAME/$TOOL_VERSION/
COPY docker-compose.yml $TOOLS/$TOOL_NAME/$TOOL_VERSION/
COPY Dockerfile $TOOLS/$TOOL_NAME/$TOOL_VERSION/

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
	ln -sf $CONFIG_HOWARD_FOLDER $TOOLS/$TOOL_NAME/$TOOL_VERSION/config/howard ;



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
