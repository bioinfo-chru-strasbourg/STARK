#!/bin/bash
#################################
##
## STARK Install Tools
##
#################################

SCRIPT_NAME="STARKInstallTools"
SCRIPT_DESCRIPTION="STARK Install Tools"
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
  	echo "# --list_tools=<FILE>          STARK Tools list file";
	echo "#                              List of tools in JSON format";
	echo "#                              Default: 'tools.json'";
  	echo "# --folder_tools=<FOLDER>      STARK tools folder ";
	echo "#                              Default: '/STARK/tools'";
	echo "# --shared_tools=<FOLDER>      STARK shared tools folder ";
	echo "#                              Default: '<folder_tools>/share/current'";
	echo "# --mamba=<PATH>               Mamba binary ";
	echo "#                              Default: 'mamba'";
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
ARGS=$(getopt -o "vdnh" --long "list_tools:,folder_tools:,shared_tools:,mamba:,verbose,debug,release,help" -- "$@" 2> /dev/null)
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
		--list_tools)
			TOOLS_JSON=$(echo "$2" | tr "," " ")
			shift 2
			;;
    	--folder_tools)
			TOOLS_DIR="$2"
			shift 2
			;;
  		--shared_tools)
            SHARED_TOOLS="$2"
            shift 2
            ;;
  		--mamba)
			MAMBA="$2"
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
		*) 	echo "Option $1 is not recognized. " "Use -h or --help to display the help." && \
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

(($DEBUG)) && echo "DEBUG mode"

if [ "$TOOLS_JSON" == "" ]; then
    TOOLS_JSON="tools.json"
fi

if [ "$TOOLS_DIR" == "" ]; then
    TOOLS_DIR="/STARK/tools"
fi

if [ "$SHARED_TOOLS" == "" ]; then
    SHARED_TOOLS="$TOOLS_DIR/shared/current"
fi

if [ "$MAMBA" == "" ]; then
    MAMBA="mamba"
fi

echo "TOOLS_JSON=$TOOLS_JSON"
echo "TOOLS_DIR=$TOOLS_DIR"
echo "SHARED_TOOLS=$SHARED_TOOLS"
echo "MAMBA=$MAMBA"

#TOOLS_DIR="/STARK/tools"

# Lire la liste des outils depuis le fichier JSON
TOOLS=$(jq -c '.tools[]' $TOOLS_JSON)

# echo "TOOLS=$TOOLS"

# TOOLS_LIST_FOR_MAMBA=$(for TOOL in $TOOLS; do NAME=$(echo $TOOL | jq -r '.name'); VERSION=$(echo $TOOL | jq -r '.version'); echo "$NAME=$VERSION"; done)
# echo "TOOLS_LIST_FOR_MAMBA=$TOOLS_LIST_FOR_MAMBA"


if ((1)); then

    # Shared tools
    TOOLS_LIST_FOR_MAMBA=""
    PATH_TOOLS=""

    # List of tools
    for TOOL in $TOOLS; do 
	#jq -c '.tools[]' $TOOLS_JSON | while read -r TOOL; do
	#echo $TOOLS | while read -r TOOL; do

        # Info
        NAME=$(echo $TOOL | jq -r '.name');
        VERSION=$(echo $TOOL | jq -r '.version');
        TYPE=$(echo $TOOL | jq -r '.type');
        DEST="$TOOLS_DIR/$NAME/$VERSION"

		# Clean version for Mamba
		MAMBA_VERSION=$VERSION

        echo "Install $NAME=$VERSION ($TYPE)..."

        # Shared tools
        if [ "$TYPE" == "shared" ]; then
			echo "Shared install"
            TOOLS_LIST_FOR_MAMBA=$TOOLS_LIST_FOR_MAMBA" $NAME=$MAMBA_VERSION"
            #DEST=$SHARED_TOOLS
        elif [ "$TYPE" == "mamba" ]; then
			echo "Mamba install"
			if ! $MAMBA create -y -p $TOOLS_DIR/$NAME/$VERSION -c bioconda -c compbiocore -c conda-forge $NAME=$MAMBA_VERSION; then
				echo "#[ERROR] Mamba installation failed."
				exit 1
			fi
            $MAMBA clean -y --all
            find $TOOLS_DIR/$NAME/$VERSION -follow -ignore_readdir_race \( -name '*.a' -o -name '*.pyc' -o -name '*.txt' -o -name '*.md' -o -name '*.pdf' -o  -name '__pycache__' \) -exec rm -rf {} +
        elif [ "$TYPE" == "custom" ]; then
            echo "Custom install"

        else
            echo "WARNING: Unknown type $TYPE"
        fi

        # Linking
        
        mkdir -p $TOOLS_DIR/$NAME
        echo "Linking $NAME/$VERSION and $NAME/current ($TYPE)"
        echo ln -s $SHARED_TOOLS/ $DEST
        echo ln -s $DEST/ $TOOLS_DIR/$NAME/current

        if [ "$TYPE" == "shared" ]; then
            ln -s $SHARED_TOOLS/ $DEST
        fi;
        ln -s $DEST/ $TOOLS_DIR/$NAME/current

        PATH_TOOLS=$PATH_TOOLS:$TOOLS_DIR/$NAME/$VERSION/bin

    done;

	echo "TOOLS_LIST_FOR_MAMBA=$TOOLS_LIST_FOR_MAMBA"

    if [ "$TOOLS_LIST_FOR_MAMBA" != "" ] && ((1)); then
        echo $MAMBA create -y -p $SHARED_TOOLS -c bioconda -c compbiocore -c conda-forge $TOOLS_LIST_FOR_MAMBA
		if ! $MAMBA create -y -p $SHARED_TOOLS -c bioconda -c compbiocore -c conda-forge $TOOLS_LIST_FOR_MAMBA; then
			echo "#[ERROR] Mamba installation failed."
			exit 1
		fi
        $MAMBA clean -y --all
        find $SHARED_TOOLS -follow -ignore_readdir_race \( -name '*.a' -o -name '*.pyc' -o -name '*.txt' -o -name '*.md' -o -name '*.pdf' -o  -name '__pycache__' \) -exec rm -rf {} +
		
    fi;
    
fi;

# jq -c '.tools[]' $TOOLS_JSON | while read -r TOOL; do
#   NAME=$(echo $TOOL | jq -r '.name')
#   VERSION=$(echo $TOOL | jq -r '.version')
#   TYPE=$(echo $TOOL | jq -r '.type')
#   COMMENT=$(echo $TOOL | jq -r '.comment // empty')
  
#   echo "Name: $NAME, Version: $VERSION, Type: $TYPE, Comment: $COMMENT"
  
#   # Vous pouvez ajouter ici des actions spécifiques pour chaque outil
# done

# # Installer chaque outil et créer des liens symboliques
# for TOOL in $TOOLS; do
# 	echo $TOOL
#   NAME=$(echo $TOOL | jq -r '.name')
#   VERSION=$(echo "$TOOL" | jq -r '.version')
#   TYPE=$(echo "$TOOL" | jq -r '.type')
#   DEST="$TOOLS_DIR/$NAME/$VERSION"
#   echo "DEST=$DEST"
# #   mkdir -p $TOOLS_DIR/$NAME
# #   echo "Linking $NAME/$VERSION and $NAME/current ($TYPE)"
# #   ln -s $TOOLS_DIR/share/current $DEST
# #   ln -s $DEST $TOOLS_DIR/$NAME/current
#   #mkdir -p $DEST

  
# #   echo "[INFO] Installing $NAME version $VERSION"
# #   mamba create -y -p $DEST -c bioconda -c conda-forge $NAME=$VERSION
# #   mamba clean -y --all
  
# #   echo "[INFO] Creating symbolic link for $NAME"
# #   ln -s $VERSION $TOOLS_DIR/$NAME/current
# done

if ((1)); then
    echo "PATH_TOOLS=$PATH_TOOLS"
    echo "" >> ~/.bashrc
    echo "# PATH for tools" >> ~/.bashrc
    echo "export PATH=\$PATH$PATH_TOOLS" >> ~/.bashrc
fi;