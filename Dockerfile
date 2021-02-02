
##############################################################
# Dockerfile Version:   1.2
# Software:             STARK
# Software Version:     0.9.18.2
# Software Website:     https://gitlab.bioinfo-diag.fr/Strasbourg/STARK
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

FROM centos:7
MAINTAINER Antony Le Bechec <antony.lebechec@gmail.com>
LABEL Software="STARK" \
	Version="0.9.18.2" \
	Website="https://gitlab.bioinfo-diag.fr/Strasbourg/STARK" \
	Description="STARK" \
	License="GNU Affero General Public License (AGPL)" \
	Usage="docker run [-v [DATA FOLDER]:/STARK/data -v [DATABASE FOLDER]:/STARK/databases -v [RESULTS FOLDER]:/STARK/output/results -v [RUNS FOLDER]:/STARK/input/runs -v [MANIFESTS FOLDER]:/STARK/input/manifests] stark:version"



########
# ARGS #
#######

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


##########
# HEADER #
##########

RUN echo "###########################" && \
	echo "# STARK BASE INSTALLATION #" && \
	echo "###########################" && \
	echo "#" && \
	echo "#### PARAMETERS" && \
	echo "#[INFO] STARK_FOLDER=$STARK_FOLDER" && \
	echo "#[INFO] TOOLS=$TOOLS" && \
	echo "#[INFO] DATA=$DATA" && \
	echo "#[INFO] TOOL=$TOOL" && \
	echo "#[INFO] SOURCES_FOLDER=$SOURCES_FOLDER" && \
	echo "#[INFO] SOURCES=$SOURCES" && \
	echo "#[INFO] DATABASES=$DATABASES" && \
	echo "#[INFO] WORKDIR=$WORKDIR" && \
	echo "#[INFO] REMOVE_SOURCES=$REMOVE_SOURCES" && \
	echo "#[INFO] THREADS=$THREADS" && \
	echo "#[INFO] REPO=$REPO" && \
	echo "#";




###########
# SOURCES #
###########

ADD . $SOURCES



###########
# WORKDIR #
###########

WORKDIR $WORKDIR



###############
# YUM INSTALL #
###############

ENV YUM_INSTALL="autoconf automake htop bc bzip2 bzip2-devel curl gcc gcc-c++ git lzma lzma-devel make ncurses-devel perl perl-Data-Dumper perl-Digest-MD5 perl-Switch perl-devel perl-Tk tbb-devel unzip rsync wget which xz xz-devel zlib zlib-devel zlib2 zlib2-devel docker java-1.7.0 java-1.8.0 python2 python2-pip python3 python3-pip curl-devel"
#ENV YUM_INSTALL=" wget rsync python2 python2-pip python3 python3-pip"
ENV YUM_REMOVE="autoconf automake bzip2-devel lzma-devel ncurses-devel perl-devel tbb-devel xz-devel zlib-devel zlib2-devel python3-devel curl-devel"
ENV PYTHON_MODULE=""
ENV PYTHON2_MODULE=$PYTHON_MODULE" pathos numpy scipy argparse multiprocess" 
ENV PYTHON3_MODULE=$PYTHON_MODULE" pathos numpy scipy argparse"

#pathos

ENV REPO_SYSTEM_GIT="$REPO/sources.system.tar.gz?path=sources/system"
ENV REPO_SYSTEM_HTTP="$REPO/sources/system/"
ENV REPO_PYTHON_GIT="$REPO/sources.python.tar.gz?path=sources/python"
ENV REPO_PYTHON_HTTP="$REPO/sources/python/"


RUN echo "#### SYSTEM INSTALLATION" && \
	# Create system repository \
	mkdir -p $SOURCES/$SOURCES_FOLDER/system && \
	# INSTALL WGET \
	echo "#[INFO] System install wget package" && \
	#ls $SOURCES/$SOURCES_FOLDER/system/*.rpm && \
	if ! ls $SOURCES/$SOURCES_FOLDER/system/wget-*.rpm 1> /dev/null 2>&1; then \
		echo "#[INFO] System wget package not locally available"; \
		yum $YUM_PARAM install -y --nogpgcheck --downloadonly --downloaddir=$SOURCES/$SOURCES_FOLDER/system/ wget; \
		echo "#[INFO] System wget package downloaded from YUM Repository"; \
	fi && \
	echo "#[INFO] System install rsync package" && \
	if ! ls $SOURCES/$SOURCES_FOLDER/system/rsync-*.rpm 1> /dev/null 2>&1; then \
		echo "#[INFO] System rsync package not locally available"; \
		yum $YUM_PARAM install -y --nogpgcheck --downloadonly --downloaddir=$SOURCES/$SOURCES_FOLDER/system/ rsync; \
		echo "#[INFO] System rsync package downloaded from YUM Repository"; \
	fi && \
	# Install packages locally \
	echo "#[INFO] System packages installation locally" && \
	yum $YUM_PARAM localinstall -y --nogpgcheck $SOURCES/$SOURCES_FOLDER/system/wget-*.rpm $SOURCES/$SOURCES_FOLDER/system/rsync-*.rpm && \
	# Test WGET installation \
	if ! command -v wget 1>/dev/null 2>/dev/null; then \
		echo "#[ERROR] System wget package not installed (Please open Internet connexion or provide WGET rpm in sources/system folder)"; \
		exit 1; \
	fi && \
	if ! command -v rsync 1>/dev/null 2>/dev/null; then \
		echo "#[ERROR] System rsync package not installed (Please open Internet connexion or provide RSYNC rpm in sources/system folder)"; \
		exit 1; \
	fi && \
	# DOWNLOAD packages from repository \
	echo "#[INFO] System packages download from REPO '$REPO'"; \
	mkdir -p $SOURCES/$SOURCES_FOLDER/system/build && \
	# in GIT mode
	if wget -q --progress=bar:force --tries=3 $REPO_SYSTEM_GIT -O $SOURCES/$SOURCES_FOLDER/system/build/STARK-repo.sources.system.tar.gz; then \
		if tar xf $SOURCES/$SOURCES_FOLDER/system/build/STARK-repo.sources.system.tar.gz -C $SOURCES/$SOURCES_FOLDER/system/build/; then \
			rsync -auczqAXhi --no-links --no-perms --no-owner --no-group $SOURCES/$SOURCES_FOLDER/system/build/STARK-repo.sources.system*/sources/system/*rpm $SOURCES/$SOURCES_FOLDER/system/; \
			echo "#[INFO] System packages downloaded from REPO '$REPO' (GIT)"; \
		else \
			echo "#[WARNING] System fail to uncompress packages from REPO '$REPO'"; \
		fi; \
	# in HTTP mode
	elif wget -q --progress=bar:force --tries=3 -r --no-parent $REPO_SYSTEM_HTTP -x --directory-prefix=$SOURCES/$SOURCES_FOLDER/system/build/STARK-repo.sources.system/; then \
		rsync -auczqAXhi --no-links --no-perms --no-owner --no-group $SOURCES/$SOURCES_FOLDER/system/build/STARK-repo.sources.system/*/sources/system/*rpm $SOURCES/$SOURCES_FOLDER/system/; \
		echo "#[INFO] System packages downloaded from REPO '$REPO' (FTP/HTTP)"; \
	else \
		echo "#[WARNING] System fail packages download from REPO '$REPO'"; \
	fi && \
	rm -rf $SOURCES/$SOURCES_FOLDER/system/build && \
	# Install packages locally \
	echo "#[INFO] System packages installation locally" && \
	yum $YUM_PARAM localinstall -y --nogpgcheck $SOURCES/$SOURCES_FOLDER/system/*.rpm && \
	# Install EPEL Repository \
	echo "#[INFO] System EPEL Repository package" && \
	if ! ls $SOURCES/$SOURCES_FOLDER/system/epel-release-*.rpm 1> /dev/null 2>&1; then \
		yum $YUM_PARAM install -y --nogpgcheck --downloadonly --downloaddir=$SOURCES/$SOURCES_FOLDER/system/ epel-release; \
		echo "#[INFO] System EPEL Repository package downloaded from YUM repository"; \
	fi && \
	if ls $SOURCES/$SOURCES_FOLDER/system/epel-release-*.rpm 1> /dev/null 2>&1; then \
		yum $YUM_PARAM localinstall -y --nogpgcheck $SOURCES/$SOURCES_FOLDER/system/epel-release-*.rpm; \
		echo "#[INFO] System EPEL Repository package enabled"; \
	else \
		echo "#[WARNING] System fail enable EPEL Repository"; \
	fi && \
	# Update YUM \
	echo "#[INFO] System packages update from YUM Repository" && \
	mkdir -p $SOURCES/$SOURCES_FOLDER/system/build/update && \
	yum $YUM_PARAM update -y --downloadonly --downloaddir=$SOURCES/$SOURCES_FOLDER/system/build/update && \
	yum $YUM_PARAM update -y && \
	#yum $YUM_PARAM localinstall -y --nogpgcheck $SOURCES/$SOURCES_FOLDER/system/build/update/*.rpm && \
	rsync -auczqAXhi --no-links --no-perms --no-owner --no-group $SOURCES/$SOURCES_FOLDER/system/build/update/*rpm $SOURCES/$SOURCES_FOLDER/system/ && \
	echo "#[INFO] System packages downloaded & updated from YUM Repository"; \
	echo "#[INFO] System packages install from YUM Repository" && \
	mkdir -p $SOURCES/$SOURCES_FOLDER/system/build/install && \
	yum $YUM_PARAM install -y --downloadonly --downloaddir=$SOURCES/$SOURCES_FOLDER/system/build/install/ $YUM_INSTALL && \
	yum $YUM_PARAM install -y $YUM_INSTALL && \
	#yum $YUM_PARAM localinstall -y --nogpgcheck $SOURCES/$SOURCES_FOLDER/system/build/install/*.rpm && \
	rsync -auczqAXhi --no-links --no-perms --no-owner --no-group $SOURCES/$SOURCES_FOLDER/system/build/install/*rpm $SOURCES/$SOURCES_FOLDER/system/ && \
	echo "#[INFO] System packages downloaded & installed from YUM Repository"; \
	rm -rf $SOURCES/$SOURCES_FOLDER/system/build && \
	# PYTHON
	# Install whl
	echo "#[INFO] System Python packages download from REPO '$REPO'"; \
	mkdir -p $SOURCES/$SOURCES_FOLDER/python/build && \
	# in GIT mode
	if wget --progress=bar:force --tries=3 $REPO_PYTHON_GIT -O $SOURCES/$SOURCES_FOLDER/python/build/STARK-repo.sources.python.tar.gz; then \
		if tar xf $SOURCES/$SOURCES_FOLDER/python/build/STARK-repo.sources.python.tar.gz -C $SOURCES/$SOURCES_FOLDER/python/build/; then \
			rsync -auczqAXhi --no-links --no-perms --no-owner --no-group $SOURCES/$SOURCES_FOLDER/python/build/STARK-repo.sources.python*/sources/system/* $SOURCES/$SOURCES_FOLDER/python/; \
			echo "#[INFO] System Python packages downloaded from REPO '$REPO' (GIT)"; \
		else \
			echo "#[WARNING] System fail to uncompress Python packages from REPO '$REPO'"; \
		fi; \
	# in HTTP mode
	elif wget --progress=bar:force --tries=3 -r --no-parent $REPO_PYTHON_HTTP -x --directory-prefix=$SOURCES/$SOURCES_FOLDER/python/build/STARK-repo.sources.system/; then \
		rsync -auczqAXhi --no-links --no-perms --no-owner --no-group $SOURCES/$SOURCES_FOLDER/system/build/STARK-repo.sources.system*/sources/python/* $SOURCES/$SOURCES_FOLDER/python/; \
		echo "#[INFO] System Python packages downloaded from REPO '$REPO' (FTP/HTTP)"; \
	else \
		echo "#[WARNING] System fail Python packages download from REPO '$REPO'"; \
	fi && \
	rm -rf $SOURCES/$SOURCES_FOLDER/python/build && \
	# Install Python packages locally \
	echo "#[INFO] System Python packages installation locally" && \
	if ls $SOURCES/$SOURCES_FOLDER/python/2/*whl 1> /dev/null 2>&1; then \
		pip2 --no-cache-dir install $SOURCES/$SOURCES_FOLDER/python/2/*whl ; \
	fi && \
	if ls $SOURCES/$SOURCES_FOLDER/python/3/*whl 1> /dev/null 2>&1; then \
		pip3 --no-cache-dir install $SOURCES/$SOURCES_FOLDER/python/3/*whl ; \
	fi && \
	# Update PIP \
	echo "#[INFO] System Python packages update from PIP Repository" && \
	echo "#[INFO] System Python packages update from PIP Repository - Python2" && \
	mkdir -p $TOOLS/python/2/bin && \
	curl https://bootstrap.pypa.io/2.7/get-pip.py --output $TOOLS/python/2/bin/get-pip.py && \
	chmod u+x $TOOLS/python/2/bin/get-pip.py && \
	$TOOLS/python/2/bin/get-pip.py && \
	python2 -m pip --no-cache-dir install --upgrade pip && \
	echo "#[INFO] System Python packages update from PIP Repository - Python3" && \
	python3 -m pip  --no-cache-dir install --upgrade pip && \
	echo "#[INFO] System Python packages downloaded & updated from PIP Repository" && \
	echo "#[INFO] System Python packages install from PIP Repository" && \
	if (( $(echo $PYTHON2_MODULE | wc -w | tr -d " ") )); then \
		python2 -m pip --no-cache-dir download $PYTHON2_MODULE --dest $SOURCES/$SOURCES_FOLDER/python/2/ ; \
		python2 -m pip --no-cache-dir install $SOURCES/$SOURCES_FOLDER/python/2/*whl ; \
	fi && \
	if (( $(echo $PYTHON3_MODULE | wc -w | tr -d " ") )); then \
		python3 -m pip --no-cache-dir download $PYTHON3_MODULE --dest $SOURCES/$SOURCES_FOLDER/python/3/ ; \
		python3 -m pip --no-cache-dir install $SOURCES/$SOURCES_FOLDER/python/3/*whl ; \
	fi && \
	echo "#[INFO] System Python packages downloaded & installed from PIP Repository" && \
	# CLEAN
	echo "#[INFO] System Cleaning" && \
	if (($REMOVE_SOURCES)); then \
		rm -rf $SOURCES/$SOURCES_FOLDER/* ; \
		echo "#[INFO] System Remove Sources" ; \
	fi && \
	yum clean all && \
	rm -rf /var/cache/yum && \
	echo "#[INFO] System Clean" && \
	echo "#";



##################
# SOURCE SCRIPTS #
##################

ENV GET_TOOL_SOURCE=$SOURCES/$SOURCES_FOLDER/get_tool_source.sh
ENV TOOL_INIT=$SOURCES/$SOURCES_FOLDER/tool_init.sh
ENV TOOL_CHECK=$SOURCES/$SOURCES_FOLDER/tool_check.sh

RUN echo "#### SOURCES SCRIPTS" && \
	if [ -e $GET_TOOL_SOURCE ]; then \
		echo "#[INFO] GET TOOL SOURCE script exists" ; \
	elif $(wget --no-cache --progress=bar:force -nv --quiet "$REPO/$SOURCES_FOLDER/$(basename $GET_TOOL_SOURCE)" -O $GET_TOOL_SOURCE); then \
		echo "#[INFO] GET TOOL SOURCE script downloaded from REPO '$REPO/$SOURCES_FOLDER/$(basename $GET_TOOL_SOURCE)'" ; \
	else \
		mkdir -p $(dirname $GET_TOOL_SOURCE) ; \
		echo 'echo "#[INFO] TOOL source" && \
		mkdir -p $(dirname $TOOL_SOURCE) && \
		if [ -e $TOOL_SOURCE ]; then \
			echo "#[INFO] TOOL TARBALL already in $TOOL_SOURCE"; \
		elif $(wget --no-cache --progress=bar:force "$TOOL_SOURCE_REPO" -O $TOOL_SOURCE); then \
			echo "#[INFO] TOOL TARBALL downloaded from STARK REPO $TOOL_SOURCE_REPO"; \
			if $(wget --no-cache --progress=bar:force -nv --quiet "$(dirname $TOOL_SOURCE_REPO)/source.info" -O $(dirname $TOOL_SOURCE)/source.info); then \
				echo "#[INFO] TOOL TARBALL external source information downloaded from STARK REPO $TOOL_SOURCE_REPO " ; \
			fi ; \
		elif $(wget --no-cache --progress=bar:force "$TOOL_SOURCE_EXTERNAL" -O $TOOL_SOURCE); then \
			echo "#[INFO] TOOL TARBALL downloaded from EXTERNAL SOURCE $TOOL_SOURCE_EXTERNAL"; \
			echo "$TOOL_SOURCE_EXTERNAL" > $(dirname $TOOL_SOURCE)/source.info; \
		else \
			echo "#[ERROR] TOOL TARBALL NOT FOUND"; \
			exit 1; \
		fi && \
		if [ -e $(dirname $TOOL_SOURCE)/source.info ]; then \
			echo "#[INFO] TOOL TARBALL external source: "$(cat $(dirname $TOOL_SOURCE)/source.info) ; \
		fi && \
		exit 0;' > $GET_TOOL_SOURCE ; \
		echo "#[INFO] GET TOOL SOURCE script written" ; \
	fi && \
	chmod u+x $GET_TOOL_SOURCE && \
	if [ -e $TOOL_INIT ]; then \
		echo "#[INFO] INIT TOOL script exists" ; \
	elif $(wget --no-cache --progress=bar:force -nv --quiet "$REPO/$SOURCES_FOLDER/$(basename $TOOL_INIT)" -O $TOOL_INIT); then \
		echo "#[INFO] INIT TOOL script downloaded from REPO '$REPO/$SOURCES_FOLDER/$(basename $TOOL_INIT)'" ; \
	else \
		mkdir -p $(dirname $TOOL_INIT) ; \
		echo 'echo "#[INFO] TOOL $TOOL_NAME/$TOOL_VERSION" && \
		export TOOL_SOURCE=$SOURCES/$SOURCES_FOLDER/tools/$TOOL_NAME/$TOOL_VERSION/$TOOL_TARBALL && \
		export TOOL_SOURCE_REPO=$REPO/$SOURCES_FOLDER/tools/$TOOL_NAME/$TOOL_VERSION/$TOOL_TARBALL && \
		export TOOL_SOURCE_BUILD=$SOURCES/$SOURCES_FOLDER/tools/$TOOL_NAME/$TOOL_VERSION/build && \
		export TOOL_DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION && \
		export PATH=$TOOL_DEST/bin:$PATH && \
		# Get TOOL SOURCE && \
		$GET_TOOL_SOURCE $TOOL_SOURCE $TOOL_SOURCE_REPO $TOOL_SOURCE_EXTERNAL && \
		# TOOL folder preparation \
		echo "#[INFO] TOOL preparation" && \
		mkdir -p $TOOL_SOURCE_BUILD && \
		mkdir -p $TOOL_DEST/bin && \
		ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current' > $TOOL_INIT ; \
		echo "#[INFO] INIT TOOL script written" ; \
	fi && \
	chmod u+x $TOOL_INIT && \
	if [ -e $TOOL_CHECK ]; then \
		echo "#[INFO] INIT TOOL script exists" ; \
	elif $(wget --no-cache --progress=bar:force -nv --quiet "$REPO/$SOURCES_FOLDER/$(basename $TOOL_CHECK)" -O $TOOL_CHECK); then \
		echo "#[INFO] INIT TOOL script downloaded from REPO '$REPO/$SOURCES_FOLDER/$(basename $TOOL_CHECK)'" ; \
	else \
		mkdir -p $(dirname $TOOL_CHECK) ; \
		echo 'echo "#[INFO] TOOL cleaning" && \
		rm -rf $TOOL_SOURCE_BUILD && \
		if (($REMOVE_SOURCES)); then rm -rf $SOURCES/$SOURCES_FOLDER/tools/$TOOL_NAME; fi && \
		echo "#[INFO] TOOL $TOOL_NAME/$TOOL_VERSION installed" ;' > $TOOL_CHECK ; \
		echo "#[INFO] INIT TOOL script written" ; \
	fi && \
	chmod u+x $TOOL_CHECK && \
	echo "#";



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



#########
# JAVA7 #
#########

ENV TOOL_NAME="java"
ENV TOOL_VERSION="1.7.0"
RUN mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
	ln -s /usr/lib/jvm/jre-1.7.0/bin/java $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/java ;



#########
# JAVA8 #
#########

ENV TOOL_NAME="java"
ENV TOOL_VERSION="1.8.0"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
RUN mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
	ln -s /usr/lib/jvm/jre-1.8.0/bin/java $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/java && \
	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;



###########
# PYTHON2 #
###########

ENV TOOL_NAME="python"
ENV TOOL_VERSION="2"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
RUN mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
	ln -s /usr/bin/python2 $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/python2 && \
	ln -s python2 $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/python && \
	ln -s /usr/bin/pip2 $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/pip2 && \
	ln -s pip2 $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/pip ;



###########
# PYTHON3 #
###########

ENV TOOL_NAME="python"
ENV TOOL_VERSION="3"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
RUN mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
	ln -s /usr/bin/python3 $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/python3 && \
	ln -s python3 $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/python && \
	ln -s /usr/bin/pip3 $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/pip3 && \
	ln -s pip3 $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/pip && \
	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;



################
# TOOLS EXTERN #
################



##########
# PATHOS #
##########

# # TOOL INFO
# ENV TOOL_NAME="pathos"
# ENV TOOL_VERSION="master"
# ENV TOOL_TARBALL="$TOOL_VERSION.tar.gz"
# ENV TOOL_SOURCE_EXTERNAL="https://github.com/uqfoundation/pathos/archive/$TOOL_TARBALL"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS
# ENV TOOL_TARBALL_FOLDER="$TOOL_NAME-$TOOL_VERSION"

# # TOOL INSTALLATION
# RUN source $TOOL_INIT && \
# 	echo "#[INFO] TOOL installation" && \
# 	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
# 	echo "#[INFO] TOOL installation - python2 setup.py build" && \
# 	(cd $TOOL_SOURCE_BUILD/$TOOL_TARBALL_FOLDER/; python2 setup.py build) && \
# 	echo "#[INFO] TOOL installation - python2 setup.py install" && \
# 	(cd $TOOL_SOURCE_BUILD/$TOOL_TARBALL_FOLDER/; python2 setup.py install) && \
#     $TOOL_CHECK ;


#echo "#[INFO] TOOL installation - pip2 install" && \
#	pip2 install pathos==0.2.6 && \
	


###########
# ANNOVAR #
###########

# TOOL INFO
ENV TOOL_NAME="annovar"
ENV TOOL_VERSION="2019Oct24"
ENV TOOL_TARBALL="$TOOL_NAME.latest.tar.gz"
ENV TOOL_SOURCE_EXTERNAL="http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS
ENV TOOL_PARAM_TARBALL_FOLDER=$TOOL_NAME
ENV TOOL_PARAM_DATABASE_FOLDER_LINK=$DATABASES/annovar/current
ENV TOOL_PARAM_DATABASE_FOLDER=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/databases/


# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
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

# TOOL INFO
ENV TOOL_NAME="htslib"
ENV TOOL_VERSION="1.11"
ENV TOOL_TARBALL="$TOOL_NAME-$TOOL_VERSION.tar.bz2"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/samtools/$TOOL_NAME/releases/download/$TOOL_VERSION/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
	make install --quiet -j $THREADS -C $(ls -d $TOOL_SOURCE_BUILD/*) prefix=$TOOL_DEST && \
    $TOOL_CHECK ;



############
# BCFTOOLS #
############

# TOOL INFO
ENV TOOL_NAME="bcftools"
ENV TOOL_VERSION="1.11"
ENV TOOL_TARBALL="$TOOL_NAME-$TOOL_VERSION.tar.bz2"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/samtools/$TOOL_NAME/releases/download/$TOOL_VERSION/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
	make install --quiet -j $THREADS -C $(ls -d $TOOL_SOURCE_BUILD/*) prefix=$TOOL_DEST && \
	$TOOL_CHECK ;



#############
# BCL2FASTQ #
#############

# TOOL INFO
ENV TOOL_NAME="bcl2fastq"
ENV TOOL_VERSION="2.20.0"
ENV TOOL_TARBALL=$TOOL_NAME"2-v2-20-0-linux-x86-64.zip"
ENV TOOL_SOURCE_EXTERNAL="https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	unzip -q $TOOL_SOURCE -d $TOOL_SOURCE_BUILD && \
    rpm -ih $TOOL_SOURCE_BUILD/*.rpm --excludedocs --prefix=$TOOL_DEST && \
    $TOOL_CHECK ;



############
# BEDTOOLS #
############

# TOOL INFO
ENV TOOL_NAME="bedtools"
ENV TOOL_VERSION="2.29.2"
ENV TOOL_TARBALL="$TOOL_NAME-$TOOL_VERSION.tar.gz"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/arq5x/bedtools2/releases/download/v$TOOL_VERSION/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
	make install --quiet -j $THREADS -C $(ls -d $TOOL_SOURCE_BUILD/*) prefix=$TOOL_DEST && \
    $TOOL_CHECK ;



############
# BOWTIE2 #
############

# TOOL INFO
ENV TOOL_NAME="bowtie2"
ENV TOOL_VERSION="2.4.2"
ENV TOOL_TARBALL="$TOOL_NAME-$TOOL_VERSION-linux-x86_64.zip"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/BenLangmead/bowtie2/releases/download/v$TOOL_VERSION/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	unzip -q $TOOL_SOURCE -d $TOOL_SOURCE_BUILD && \
	cp $TOOL_SOURCE_BUILD/*/bowtie2* $TOOL_DEST/bin/ && \
	rm -f $TOOL_DEST/bin/*debug && \
    $TOOL_CHECK ;



#######
# BWA #
#######

# TOOL INFO
ENV TOOL_NAME="bwa"
ENV TOOL_VERSION="0.7.17"
ENV TOOL_TARBALL="$TOOL_NAME-$TOOL_VERSION.tar.bz2"
ENV TOOL_SOURCE_EXTERNAL="https://sourceforge.net/projects/bio-bwa/files/$TOOL_NAME-$TOOL_VERSION.tar.bz2/download"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
	make --quiet -j $THREADS -C $(ls -d $TOOL_SOURCE_BUILD/*) && \
	cp $TOOL_SOURCE_BUILD/*/bwa $TOOL_DEST/bin/ && \
    $TOOL_CHECK ;



#########
# FASTP #
#########

# TOOL INFO
ENV TOOL_NAME="fastp"
ENV TOOL_VERSION="0.20.0"
ENV TOOL_TARBALL="$TOOL_NAME"
ENV TOOL_SOURCE_EXTERNAL="http://opengene.org/$TOOL_NAME/$TOOL_NAME"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	cp $TOOL_SOURCE $TOOL_DEST/bin/ && \
	chmod a+x $TOOL_DEST/bin/* && \
    $TOOL_CHECK ;



#######
# CAP #
#######

# TOOL INFO
ENV TOOL_NAME="cap"
ENV TOOL_VERSION="0.9.13"
ENV TOOL_TARBALL="archive.tar.gz"
ENV TOOL_SOURCE_EXTERNAL="https://gitlab.bioinfo-diag.fr/Strasbourg/CAP/repository/$TOOL_VERSION/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
	cp -R $TOOL_SOURCE_BUILD/*/* $TOOL_DEST/ && \
	chmod a+x $TOOL_DEST/bin/* && \
    $TOOL_CHECK ;



########
# GATK #
########

# TOOL INFO
ENV TOOL_NAME="gatk"
ENV TOOL_VERSION="3.8-1-0"
ENV TOOL_TARBALL="GenomeAnalysisTK-$TOOL_VERSION.tar.bz2"
ENV TOOL_SOURCE_EXTERNAL="https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=$TOOL_VERSION-gf15c1c3ef"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS
ENV TOOL_JAR=GenomeAnalysisTK.jar

# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
	cp -R $TOOL_SOURCE_BUILD/*/$TOOL_JAR $TOOL_DEST/bin/ && \
    $TOOL_CHECK ;



#########
# GATK4 #
#########

# TOOL INFO
ENV TOOL_NAME="gatk"
ENV TOOL_VERSION="4.1.9.0"
ENV TOOL_TARBALL="gatk-$TOOL_VERSION.zip"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/broadinstitute/gatk/releases/download/$TOOL_VERSION/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS
ENV TOOL_JAR=gatk-package-$TOOL_VERSION-local.jar

# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	unzip -q $TOOL_SOURCE -d $TOOL_SOURCE_BUILD && \
	cp -R $TOOL_SOURCE_BUILD/*/$TOOL_JAR $TOOL_DEST/bin/ && \
    $TOOL_CHECK ;



##########
# HOWARD #
##########

# TOOL INFO
ENV TOOL_NAME="howard"
ENV TOOL_VERSION="0.9.15.3"
ENV TOOL_TARBALL="archive.tar.gz"
ENV TOOL_SOURCE_EXTERNAL="https://gitlab.bioinfo-diag.fr/Strasbourg/HOWARD/repository/$TOOL_VERSION/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS
#ENV TOOL_DATABASE_FOLDER=/databases
ENV TOOL_PARAM_DATABASE_FOLDER_LINK=$DATABASES
ENV TOOL_PARAM_DATABASE_FOLDER=/databases


# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
	cp -R $TOOL_SOURCE_BUILD/*/* $TOOL_DEST/ && \
	chmod a+x $TOOL_DEST/* -R && \
	mkdir -p $TOOL_PARAM_DATABASE_FOLDER_LINK && \
	mkdir -p $TOOL_PARAM_DATABASE_FOLDER && \
	ln -s $DATABASES $TOOL_DATABASE_FOLDER && \
    $TOOL_CHECK ;



############
# IGVTOOLS #
############

# TOOL INFO
ENV TOOL_NAME="igvtools"
ENV TOOL_VERSION="2.4.19"
ENV TOOL_TARBALL="igvtools_$TOOL_VERSION.zip"
ENV TOOL_VERSION_MAIN="2.4"
ENV TOOL_SOURCE_EXTERNAL="http://data.broadinstitute.org/igv/projects/downloads/$TOOL_VERSION_MAIN/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	unzip -q $TOOL_SOURCE -d $TOOL_SOURCE_BUILD && \
	cp -R $TOOL_SOURCE_BUILD/*/* $TOOL_DEST/bin/ && \
    $TOOL_CHECK ;



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
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	unzip -q $TOOL_SOURCE -d $TOOL_SOURCE_BUILD && \
	cp -R $TOOL_SOURCE_BUILD/*/* $TOOL_DEST/bin/ && \
    $TOOL_CHECK ;



##########
# MUTECT #
##########

# TOOL INFO
ENV TOOL_NAME="mutect"
ENV TOOL_VERSION="1.1.7"
ENV TOOL_TARBALL="$TOOL_NAME-$TOOL_VERSION.jar.zip"
ENV TOOL_SOURCE_EXTERNAL="https://software.broadinstitute.org/gatk/download/auth?package=M1"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS
ENV TARBALL_JAR=mutect-$TOOL_VERSION.jar
ENV TOOL_JAR=mutect.jar

# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	unzip -q $TOOL_SOURCE -d $TOOL_DEST/bin/ && \
	mv $TOOL_DEST/bin/$TARBALL_JAR $TOOL_DEST/bin/$TOOL_JAR && \
    $TOOL_CHECK ;



############
# OUTLYZER #
############

# # TOOL INFO
# ENV TOOL_NAME="outlyzer"
# ENV TOOL_VERSION="2"
# ENV TOOL_TARBALL="outLyzer.py"
# ENV TOOL_SOURCE_EXTERNAL="https://raw.githubusercontent.com/EtieM/outLyzer/master/$TOOL_TARBALL"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS

# # TOOL INSTALLATION
# RUN source $TOOL_INIT && \
# 	echo "#[INFO] TOOL installation" && \
# 	cp $TOOL_SOURCE -d $TOOL_DEST/bin/ && \
# 	chmod a+x $TOOL_DEST/bin/*.py && \
#     $TOOL_CHECK ;


# TOOL INFO
ENV TOOL_NAME="outlyzer"
ENV TOOL_VERSION="3"
ENV TOOL_TARBALL="outLyzer.py"
ENV TOOL_SOURCE_EXTERNAL="https://raw.githubusercontent.com/EtieM/outLyzer/master/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	cp $TOOL_SOURCE -d $TOOL_DEST/bin/ && \
	chmod a+x $TOOL_DEST/bin/*.py && \
    $TOOL_CHECK ;




##########
# PATHOS #
##########

# # TOOL INFO
# ENV TOOL_NAME="pathos"
# ENV TOOL_VERSION="master"
# ENV TOOL_TARBALL="$TOOL_VERSION.tar.gz"
# ENV TOOL_SOURCE_EXTERNAL="https://github.com/uqfoundation/pathos/archive/$TOOL_TARBALL"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS
# ENV TOOL_TARBALL_FOLDER="$TOOL_NAME-$TOOL_VERSION"

# # TOOL INSTALLATION
# RUN source $TOOL_INIT && \
# 	echo "#[INFO] TOOL installation" && \
# 	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
# 	echo "#[INFO] TOOL installation - pip2 install" && \
# 	pip2 install pathos==0.2.6 && \
# 	echo "#[INFO] TOOL installation - python2 setup.py build" && \
# 	(cd $TOOL_SOURCE_BUILD/$TOOL_TARBALL_FOLDER/; python2 setup.py build) && \
# 	echo "#[INFO] TOOL installation - python2 setup.py install" && \
# 	(cd $TOOL_SOURCE_BUILD/$TOOL_TARBALL_FOLDER/; python2 setup.py install) && \
#     $TOOL_CHECK ;



##########
# PICARD #
##########

# TOOL INFO
ENV TOOL_NAME="picard"
ENV TOOL_VERSION="2.23.7"
ENV TOOL_TARBALL="picard.jar"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/broadinstitute/picard/releases/download/$TOOL_VERSION/picard.jar"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	cp $TOOL_SOURCE $TOOL_DEST/bin/ && \
    $TOOL_CHECK ;



############
# SAMTOOLS #
############

# TOOL INFO
ENV TOOL_NAME="samtools"
ENV TOOL_VERSION="1.11"
ENV TOOL_TARBALL="$TOOL_NAME-$TOOL_VERSION.tar.bz2"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/samtools/$TOOL_NAME/releases/download/$TOOL_VERSION/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	tar xf $TOOL_SOURCE -C $TOOL_SOURCE_BUILD && \
	make install --quiet -j $THREADS -C $(ls -d $TOOL_SOURCE_BUILD/*) prefix=$TOOL_DEST && \
	$TOOL_CHECK ;



##########
# SNPEFF #
##########

# TOOL INFO
ENV TOOL_NAME="snpeff"
ENV TOOL_VERSION="4.3t"
ENV TOOL_TARBALL="snpEff_v4_3t_core.zip"
ENV TOOL_SOURCE_EXTERNAL="https://sourceforge.net/projects/snpeff/files/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS
#ENV TOOL_PARAM_DATABASE_FOLDER_LINK=$DATABASES/snpeff/4.3t
ENV TOOL_PARAM_DATABASE_FOLDER_LINK=$DATABASES/snpeff/current
ENV TOOL_PARAM_DATABASE_FOLDER=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/data

# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	unzip -q $TOOL_SOURCE -d $TOOL_SOURCE_BUILD && \
	cp $TOOL_SOURCE_BUILD/*/*jar $TOOL_DEST/bin/ && \
	cp $TOOL_SOURCE_BUILD/*/*config $TOOL_DEST/bin/ && \
	mkdir -p $TOOL_PARAM_DATABASE_FOLDER_LINK && \
	mkdir -p $TOOL_PARAM_DATABASE_FOLDER && \
	ln -s $TOOL_PARAM_DATABASE_FOLDER_LINK $TOOL_PARAM_DATABASE_FOLDER && \
    $TOOL_CHECK ;



###########
# VARSCAN #
###########

# TOOL INFO
ENV TOOL_NAME="varscan"
ENV TOOL_VERSION="2.4.4"
ENV TOOL_TARBALL="VarScan.v$TOOL_VERSION.jar"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/dkoboldt/varscan/raw/master/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS
ENV TOOL_PARAM_JAR_NAME=VarScan.jar

# TOOL INSTALLATION
RUN source $TOOL_INIT && \
	echo "#[INFO] TOOL installation" && \
	cp $TOOL_SOURCE $TOOL_DEST/bin/ && \
	ln -s $(basename $TOOL_SOURCE) $TOOL_DEST/bin/$TOOL_PARAM_JAR_NAME && \
    $TOOL_CHECK ;



#########
# STARK #
#########


# TOOL INFO
ENV TOOL_NAME="stark"
ENV TOOL_VERSION="0.9.18.2"
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
RUN mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
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

RUN yum erase -y $YUM_REMOVE && \
	yum clean all && \
    rm -rf /var/cache/yum && \
	rm -rf $WORKDIR/* && \
	rm -rf /tmp/* && \
	if (($REMOVE_SOURCES)); then rm -rf $SOURCES; fi;



##############################
# WORKDIR / ENTRYPOINT / CMD #
##############################


WORKDIR "/STARK/data"

ENTRYPOINT [ "/tool/bin/STARK" ]
