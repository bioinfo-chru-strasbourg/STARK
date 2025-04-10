#!/bin/bash
#################################
##
## STARK Install System 
##
#################################

SCRIPT_NAME="STARKInstallSystem"
SCRIPT_DESCRIPTION="STARK Install System"
SCRIPT_RELEASE="0.9.0"
SCRIPT_DATE="05/03/2025"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="HUS/CPS"
SCRIPT_LICENCE="GNU GPLA V3"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.0-05/03/2025: Script creation\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


# Header
function header () {
	echo "#######################################";
	echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
	echo "# $SCRIPT_DESCRIPTION ";
	echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © $SCRIPT_LICENCE";
	echo "#######################################";
}

# Release
function release () {
	echo "# RELEASE NOTES:";
	echo -e $RELEASE_NOTES
}


# Usage
function usage {
	echo "# USAGE: $(basename $0) [--help] [options...]";
	echo "#";
	echo "### This script manages STARK modules and services.";
	echo "#";
	echo "# --yum_install=<LIST>         List of yum packages to install";
	echo "#                              Default: ''";
  	echo "# --yum_param=<LIST>           List of yum parameters";
	echo "#                              Default: ''";
  	echo "# -v|--verbose                 Verbose mode";
	echo "# -d|--debug                   Debug mode";
	echo "# -n|--release                 Script Release";
	echo "# -h|--help                    Help message";
	echo "#";

}




####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "vdnh" --long "yum_install:,yum_param:,verbose,debug,release,help" -- "$@" 2> /dev/null)
if [ $? -ne 0 ]; then
	:
	echo "#[ERROR] Error in the argument list:";
	echo "#[ERROR] $@"
	echo ""
	usage;
	exit;
fi;


PARAM=$@
DEBUG=0
VERBOSE=0

eval set -- "$ARGS"
while true
do
	case "$1" in
        --yum_install)
			YUM_INSTALL=$(echo "$2" | tr "," " ")
			shift 2
			;;
    	--yum_param)
			YUM_PARAM=$(echo "$2" | tr "," " ")
			shift 2
			;;
    	-h|--help)
			usage
			exit 0
			;;
		-v|--verbose)
			VERBOSE=1
			shift 1
			;;
		-d|--debug)
			DEBUG=1
			shift 1
			;;
		--) shift
			break
			;;
		*) 	echo "Option $1 is not recognized. " "Use -h or --help to display the help."

			exit 1
			;;
	esac
done

# header
(($NO_HEADER)) || header;



####################################################################################################################################
# Checking the input parameter
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# if [ -z "$SOURCE_RUNS" ] && [ -z "$DEST_RUNS" ] && ((!$DEBUG)); then
# 	echo "#[ERROR] Required parameter: --sources and --dest. Use --help to display the help." && echo "" && usage && exit 1;
# fi
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


(($DEBUG)) && VERBOSE=1

echo "#[INFO] SYSTEM YUM installation"

# # Create system repository
# mkdir -p $SOURCES/$SOURCES_FOLDER/system
# mkdir -p $SOURCES/$SOURCES_FOLDER/system/$(uname -m)


# echo "#[INFO] System install wget package"
# if ! ls $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/wget-*.rpm 1> /dev/null 2>&1; then \
# 	echo "#[INFO] System wget package not locally available"; \
# 	yum $YUM_PARAM install -y --nogpgcheck --downloadonly --downloaddir=$SOURCES/$SOURCES_FOLDER/system/$(uname -m)/ wget; \
# 	echo "#[INFO] System wget package downloaded from YUM Repository"; \
# fi


# echo "#[INFO] System install rsync package"
# if ! ls $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/rsync-*.rpm 1> /dev/null 2>&1; then \
# 	echo "#[INFO] System rsync package not locally available"; \
# 	yum $YUM_PARAM install -y --nogpgcheck --downloadonly --downloaddir=$SOURCES/$SOURCES_FOLDER/system/$(uname -m)/ rsync; \
# 	echo "#[INFO] System rsync package downloaded from YUM Repository"; \
# fi


# echo "#[INFO] System packages installation locally"
# yum $YUM_PARAM localinstall -y --allowerasing --nogpgcheck $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/wget-*.rpm $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/rsync-*.rpm
# if ! command -v wget 1>/dev/null 2>/dev/null; then \
# 	echo "#[ERROR] System wget package not installed (Please open Internet connexion or provide WGET rpm in sources/system folder)"; \
# 	exit 1; \
# fi
# if ! command -v rsync 1>/dev/null 2>/dev/null; then \
# 	echo "#[ERROR] System rsync package not installed (Please open Internet connexion or provide RSYNC rpm in sources/system folder)"; \
# 	exit 1; \
# fi


# echo "#[INFO] System packages download from REPO '$REPO'"; \
# mkdir -p $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build
# if wget -q --progress=bar:force --tries=3 $REPO_SYSTEM_GIT -O $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/STARK-repo.sources.system.tar.gz; then \
# 	if tar xf $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/STARK-repo.sources.system.tar.gz -C $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/; then \
# 		rsync -auczqAXhi --no-links --no-perms --no-owner --no-group --ignore-missing-args $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/STARK-repo.sources.system*/sources/system/*rpm $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/; \
# 		echo "#[INFO] System packages downloaded from REPO '$REPO' (GIT)"; \
# 	else \
# 		echo "#[WARNING] System fail to uncompress packages from REPO '$REPO'"; \
# 	fi; \
# elif wget -q --progress=bar:force --tries=3 -r --no-parent $REPO_SYSTEM_HTTP -x --directory-prefix=$SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/STARK-repo.sources.system/; then \
# 	rsync -auczqAXhi --no-links --no-perms --no-owner --no-group --ignore-missing-args $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/STARK-repo.sources.system/*/sources/system/*rpm $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/; \
# 	echo "#[INFO] System packages downloaded from REPO '$REPO' (FTP/HTTP)"; \
# else \
# 	echo "#[WARNING] System fail packages download from REPO '$REPO'"; \
# fi
# rm -rf $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build


# echo "#[INFO] System packages installation locally"
# if ! ls $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/*.rpm 1> /dev/null 2>&1; then \
# 	yum $YUM_PARAM localinstall -y --nogpgcheck $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/*.rpm; \
# 	echo "#[INFO] System packages installation locally done."; \
# fi


# echo "#[INFO] System EPEL Repository package"
# if ! ls $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/epel-release-*.rpm 1> /dev/null 2>&1; then \
# 	yum $YUM_PARAM install -y --nogpgcheck --downloadonly --downloaddir=$SOURCES/$SOURCES_FOLDER/system/$(uname -m)/ epel-release; \
# 	echo "#[INFO] System EPEL Repository package downloaded from YUM repository"; \
# fi
# if ls $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/epel-release-*.rpm 1> /dev/null 2>&1; then \
# 	yum $YUM_PARAM localinstall -y --nogpgcheck $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/epel-release-*.rpm; \
# 	echo "#[INFO] System EPEL Repository package enabled"; \
# else \
# 	echo "#[WARNING] System fail enable EPEL Repository"; \
# fi

# echo "#[INFO] System packages update from YUM Repository"
# mkdir -p $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/update
# yum $YUM_PARAM update -y --downloadonly --downloaddir=$SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/update
# yum $YUM_PARAM localinstall -y --nogpgcheck $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/update/*.rpm
# rsync -auczqAXhi --no-links --no-perms --no-owner --no-group --ignore-missing-args $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/update/*rpm $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/

# echo "#[INFO] System packages downloaded & updated from YUM Repository"

# echo "#[INFO] System packages install from YUM Repository"
# mkdir -p $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/install
# yum $YUM_PARAM install -y --downloadonly --allowerasing --downloaddir=$SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/install/ $YUM_INSTALL
# ls -lah $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/install/
# yum $YUM_PARAM localinstall -y --nogpgcheck --allowerasing $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/install/*.rpm
# rsync -auczqAXhi --no-links --no-perms --no-owner --no-group --ignore-missing-args $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build/install/*rpm $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/

# echo "#[INFO] System packages downloaded & installed from YUM Repository"
# rm -rf $SOURCES/$SOURCES_FOLDER/system/$(uname -m)/build
# yum clean -y all
# rm -rf /var/cache/yum

# echo "#[INFO] System Clean"

# YUM Install
echo "#[INFO] System packages Repository"
dnf install epel-release -y
dnf $YUM_PARAM config-manager --set-enabled crb 
# yum $YUM_PARAM install -y --allowerasing epel-release 
# if [ "$FROM_IMAGE" == "almalinux:9" ]; then
# 		yum config-manager --set-enabled crb;
# 	else \
# 		dnf config-manager --set-enabled powertools;
# 	fi

echo "#[INFO] System packages install"
dnf $YUM_PARAM install -y --allowerasing $YUM_INSTALL

echo "#[INFO] System packages install clean"
dnf clean -y all
rm -rf /var/cache/yum

echo "#[INFO] SYSTEM Bashrc"
echo "alias ll='ls -lah'" >> ~/.bashrc

echo "#"