#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARK_DATABASES_UPDATE_ANNOVAR"
SCRIPT_DESCRIPTION="STARK DATABASES UPDATE ANNOVAR databases"
SCRIPT_RELEASE="0.9.0.0"
SCRIPT_DATE="07/04/2021"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="HUS"
SCRIPT_LICENCE="GNU-AGPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.0.0-07/06/2021: Script creation\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Configuration
ENV_CONFIG=$(find -L $SCRIPT_DIR/.. -name config.app)

source $ENV_CONFIG 1>/dev/null 2>/dev/null

ENV_TOOLS=$(find -L $SCRIPT_DIR/.. -name tools.app)

source $ENV_TOOLS 1>/dev/null 2>/dev/null


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
	echo "# USAGE: $(basename $0) --config_annotation=<FILE> [options...]";
	# echo "# --application=<STRING|FILE>                   APP name or APP file configuration of the APPLICATION.";
	# echo "#                                               Must be in the STARK APPS folder if relative path";
	# echo "#                                               Default: 'default.app'";
	echo "# --config_annotation=<FILE>                    HOWARD configuration annotation ini file";
	echo "# --annotation=<STRING>                         ANNOVAR annotations";
	echo "#                                               Format: 'ann1,ann2,...'";
	echo "#                                               Default: 'ALL'";
	echo "# --annovar_new_release=<STRING>                ANNOVAR new release";
	echo "#                                               Default: 'date +%Y%m%d-%H%M%S'";
	echo "# --stark_main_folder=<FOLDER>                  STARK main folder";
	echo "#                                               Default: '$HOME/STARK'";
	echo "# --check_update                                 Check if update needed";	
	echo "# --verbose                                     VERBOSE";
	echo "# --debug                                       DEBUG";
	echo "# --release                                     RELEASE";
	echo "# --help                                        HELP";
	echo "#";
}

# header
header;


####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "vdnh" --long "env:,app:,application:,config_annotation:,annotation:,annovar_new_release:,stark_main_folder:,check_update,verbose,debug,release,help" -- "$@" 2> /dev/null)

eval set -- "$ARGS"
while true
do
	case "$1" in
		--config_annotation)
			CONFIG_ANNOTATION="$2"
			shift 2
			;;
		--annovar_new_release)
			ANNOVAR_NEW_RELEASE="$2";
			shift 2
			;;
		--annotation)
			ANNOTATION="$2";
			shift 2
			;;
		--stark_main_folder)
			STARK_MAIN_FOLDER=$2;
			shift 2
			;;
		--check_update)
			CHECK_UPDATE=1
			shift 1
			;;
		-v|--verbose)
			VERBOSE=1
			shift 1
			;;
		-d|--debug)
			VERBOSE=1
			DEBUG=1
			shift 1
			;;
		-n|--release)
			release;
			exit 0
			;;
		-h|--help)
			usage
			exit 0
			;;
		--) shift
			break
			;;
		*) 	echo "# Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done




####################################################################################################################################
# Checking the input parameter
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if [ -z "$CONFIG_ANNOTATION" ] && ! (($CHECK_UPDATE)); then
	echo "#[ERROR] Required parameter: --config_annotation or --check_update. Use --help to display the help." && echo "" && usage && exit 1;
fi
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


## PARAMETERS
##############


# CONFIG_ANNOTATION
if [ "$CONFIG_ANNOTATION" == "" ] || [ ! -e "$CONFIG_ANNOTATION" ] && ! (($CHECK_UPDATE)); then
    echo "#[ERROR] No Config Annotation file '$CONFIG_ANNOTATION'"
    exit 0
fi;
echo "#[INFO] Config Annotation file '$CONFIG_ANNOTATION'"


# ANNOTATION
if [ "$ANNOTATION" == "" ]; then
    ANNOTATION="ALL"
fi;
echo "#[INFO] Annotation '$ANNOTATION'"


# STARK_MAIN_ANNOTATION
if [ "$STARK_MAIN_FOLDER" == "" ]; then
	STARK_MAIN_FOLDER=$HOME/STARK
fi;
if [ "$STARK_MAIN_FOLDER" == "" ] || [ ! -d "$STARK_MAIN_FOLDER" ]; then
	echo "#[ERROR] No STARK main folder '$STARK_MAIN_FOLDER'"
    exit 0
fi;
echo "#[INFO] STARK main folder '$STARK_MAIN_FOLDER'"


# RELEASE
if [ "$ANNOVAR_NEW_RELEASE" == "" ]; then
	ANNOVAR_NEW_RELEASE=$(date '+%Y%m%d-%H%M%S')
fi;
[ "$ANNOVAR_NEW_RELEASE" == ""  ] && ANNOVAR_NEW_RELEASE=$(date '+%Y%m%d-%H%M%S') 
echo "#[INFO] Release '$ANNOVAR_NEW_RELEASE'"


# RELEASE LATEST
if [ -e $STARK_MAIN_FOLDER/databases/annovar/latest/ ]; then
	ANNOVAR_LATEST_RELEASE=$(basename $(realpath $STARK_MAIN_FOLDER/databases/annovar/latest/)) 
else
	echo "#[INFO] No latest release "
fi;


# FOLDERS
TMP_ID=$RANDOM$ANDOM

# HOST
TMP_FOLDER=$STARK_MAIN_FOLDER/output/tmp/$TMP_ID
HOWARD_UPDATE_RELEASE_FOLDER=$TMP_FOLDER/$ANNOVAR_NEW_RELEASE
HOWARD_UPDATE_RELEASE_FOLDER_DATABASES=$HOWARD_UPDATE_RELEASE_FOLDER/databases
HOWARD_UPDATE_RELEASE_FOLDER_TMP=$HOWARD_UPDATE_RELEASE_FOLDER/tmp

# INNER
STARK_MAIN_FOLDER_INNER=/STARK
TMP_FOLDER_INNER=$STARK_MAIN_FOLDER_INNER/output/tmp/$TMP_ID
HOWARD_UPDATE_RELEASE_FOLDER_INNER=$TMP_FOLDER_INNER/$ANNOVAR_NEW_RELEASE
HOWARD_UPDATE_RELEASE_FOLDER_DATABASES_INNER=$HOWARD_UPDATE_RELEASE_FOLDER_INNER/databases
HOWARD_UPDATE_RELEASE_FOLDER_TMP_INNER=$HOWARD_UPDATE_RELEASE_FOLDER_INNER/tmp


# DEBUG
(($DEBUG)) && echo "#[INFO] ANNOVAR_LATEST_RELEASE=$ANNOVAR_LATEST_RELEASE"
(($DEBUG)) && echo "#[INFO] HOWARD_UPDATE_RELEASE_FOLDER=$HOWARD_UPDATE_RELEASE_FOLDER"
(($DEBUG)) && echo "#[INFO] HOWARD_UPDATE_RELEASE_FOLDER_DATABASES=$HOWARD_UPDATE_RELEASE_FOLDER_DATABASES"
(($DEBUG)) && echo "#[INFO] HOWARD_UPDATE_RELEASE_FOLDER_TMP=$HOWARD_UPDATE_RELEASE_FOLDER_TMP"


# CHECK UPDATE
if (($CHECK_UPDATE)); then

	if [ -e $STARK_MAIN_FOLDER/databases/annovar/latest/hg19_avdblist.txt ]; then

		mkdir -p $TMP_FOLDER

		echo "#[INFO] Check for update..."
		docker exec -ti stark-module-stark-submodule-stark-service-cli bash -c 'source /tool/config/config.app && perl $ANNOVAR/annotate_variation.pl -webfrom annovar -downdb avdblist -buildver hg19 /tmp && grep "^NOTICE" -v /tmp/hg19_avdblist.txt ' > $TMP_FOLDER/hg19_avdblist.txt
		checksum_new=$(grep "^NOTICE" -v  $TMP_FOLDER/hg19_avdblist.txt | sort -u | sha256sum | cut -d" " -f1 )
		checksum_latest=$(grep "^NOTICE" -v  $STARK_MAIN_FOLDER/databases/annovar/latest/hg19_avdblist.txt | sort -u | sha256sum | cut -d" " -f1 )
		if [ "$checksum_new" != "$checksum_latest" ]; then
			echo "#[INFO] Update needed"
			diff -E <(grep "^NOTICE" -v $TMP_FOLDER/hg19_avdblist.txt | sort -u) <(grep "^NOTICE" -v $STARK_MAIN_FOLDER/databases/annovar/latest/hg19_avdblist.txt | sort -u) | grep "^<" | column -t
		else
			echo "#[INFO] NO update needed"
		fi
	
		rm -rf $TMP_FOLDER

	else 

		echo "#[INFO] No latest release. Can not check update "

	fi

	exit 0

fi


if ((1)); then

	mkdir -p $HOWARD_UPDATE_RELEASE_FOLDER $HOWARD_UPDATE_RELEASE_FOLDER_DATABASES $HOWARD_UPDATE_RELEASE_FOLDER_TMP

	# Fichier de configuration d'HOWARD
	# Créer et copier le fichier de configuration config.annotation.ini
	cp $CONFIG_ANNOTATION $HOWARD_UPDATE_RELEASE_FOLDER_DATABASES/config.annotation.ini


	# Téléchargement (par HOWARD dans STARK CLI)
	# Dans STARK CLI (pour l'acces a HOWARD et aux paramétrages de STARK-tools)
	docker exec -ti stark-module-stark-submodule-stark-service-cli bash -c 'cd /STARK && source /tool/config/config.app && $HOWARD --input=$HOWARD_FOLDER/docs/example.vcf --output='$HOWARD_UPDATE_RELEASE_FOLDER_DATABASES_INNER'/HOWARD.download.vcf --env=/tool/config/tools.app --annotation='$ANNOTATION' --annovar_folder=$ANNOVAR --annovar_databases='$HOWARD_UPDATE_RELEASE_FOLDER_DATABASES_INNER' --config_annotation='$HOWARD_UPDATE_RELEASE_FOLDER_DATABASES_INNER'/config.annotation.ini --snpeff_jar=$SNPEFF --snpeff_threads=1 --tmp='$HOWARD_UPDATE_RELEASE_FOLDER_TMP_INNER' --assembly=hg19 --java=$JAVA --threads=1 --verbose' 1>$HOWARD_UPDATE_RELEASE_FOLDER_DATABASES/HOWARD.download.log 2>$HOWARD_UPDATE_RELEASE_FOLDER_DATABASES/HOWARD.download.err

	# Telechargement de la liste des bases de données 'hg19_avdblist.txt'
	#docker exec -ti stark-module-stark-submodule-stark-service-cli bash -c 'cd /STARK && source /tool/config/config.app && perl $ANNOVAR/annotate_variation.pl -webfrom annovar -downdb avdblist -buildver hg19 '$HOWARD_UPDATE_RELEASE_FOLDER_DATABASES_INNER' '
	docker exec -ti stark-module-stark-submodule-stark-service-cli bash -c 'source /tool/config/config.app && perl $ANNOVAR/annotate_variation.pl -webfrom annovar -downdb avdblist -buildver hg19 /tmp && grep "^NOTICE" -v /tmp/hg19_avdblist.txt ' > $HOWARD_UPDATE_RELEASE_FOLDER_DATABASES/hg19_avdblist.txt
	

	# Update (sur le serveur)
	cd $STARK_MAIN_FOLDER
	# créer la version d'annovar
	mkdir -p $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_NEW_RELEASE
	# Copier la nouvelle version d'ANNOVAR dans database
	rsync -rauz $HOWARD_UPDATE_RELEASE_FOLDER_DATABASES/* $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_NEW_RELEASE/
	# Ajouter les fichiers manquants de la précédente version d'ANNOVAR dans database
	if [ "$ANNOVAR_LATEST_RELEASE" != "" ]; then
		rsync -rauz $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_LATEST_RELEASE/* $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_NEW_RELEASE/
	fi;
	# modifier le fichier STARK.database.release au besoin ("download")
	#mv $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_NEW_RELEASE/STARK.database.release $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_NEW_RELEASE/STARK.database.release.previous
	#sed "s/$ANNOVAR_LATEST_RELEASE/$ANNOVAR_NEW_RELEASE/gi" $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_NEW_RELEASE/STARK.database.release.previous > $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_NEW_RELEASE/STARK.database.release
	echo '{
			"release": "'$ANNOVAR_NEW_RELEASE'",
			"date": "'$ANNOVAR_NEW_RELEASE'",
			"files": [ "config.annotation.ini" ],
			"assembly": [ "hg19" ],
			"download": {
					"methode": "'$SCRIPT_DESCRIPTION' ['$SCRIPT_RELEASE'-'$SCRIPT_DATE']",
					"info": "From ANNOVAR binary",
					"databases": "'$ANNOTATION'",
					"configuration": "config.annotation.ini"
			}
	}' > $STARK_MAIN_FOLDER/databases/annovar/$ANNOVAR_NEW_RELEASE/STARK.database.release


	# Mise a disposition
	# Copier le fichier de configuration dans
	cp $HOWARD_UPDATE_RELEASE_FOLDER_DATABASES/config.annotation.ini $STARK_MAIN_FOLDER/config/howard/config.annotation.$ANNOVAR_NEW_RELEASE.ini
	# Changer latest
	if [ "$ANNOVAR_LATEST_RELEASE" != "" ]; then
		rm -f $STARK_MAIN_FOLDER/databases/annovar/latest
		ln -snf $ANNOVAR_NEW_RELEASE/ $STARK_MAIN_FOLDER/databases/annovar/latest
	fi


fi


### Production
echo "
### Manually execute these command to push ANNOVAR database into production
### !!! possibly destructive !!!
# change current release
rm -f $STARK_MAIN_FOLDER/databases/annovar/current
ln -snf $ANNOVAR_NEW_RELEASE $STARK_MAIN_FOLDER/databases/annovar/current
# change HOWARD configuration
cp $STARK_MAIN_FOLDER/config/howard/config.annotation.ini $STARK_MAIN_FOLDER/config/howard/config.annotation.ini.from_update_$ANNOVAR_NEW_RELEASE
cp $STARK_MAIN_FOLDER/config/howard/config.annotation.$ANNOVAR_NEW_RELEASE.ini $STARK_MAIN_FOLDER/config/howard/config.annotation.ini
"
if [ "$ANNOVAR_LATEST_RELEASE" != "" ]; then
echo "
# Compress previous release
(cd databases/annovar && tar -vczf $ANNOVAR_LATEST_RELEASE.tar.gz $ANNOVAR_LATEST_RELEASE/)
rm -rf databases/annovar/$ANNOVAR_LATEST_RELEASE
"
fi;


rm -rf $TMP_FOLDER

exit 0

