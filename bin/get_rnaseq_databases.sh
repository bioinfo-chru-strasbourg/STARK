#!/bin/bash
### Parameters ###
DATABASES_FOLDER=/home1/BAS/nicaises/Tests/deleteme/databases
TMP_DIR=/tmp
CTAT_TAR="GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz"
CTAT_URL="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/"$CTAT_TAR
ARRIBA_TAR="arriba_v2.4.0.tar.gz"
ARRIBA_URL="https://github.com/suhrig/arriba/releases/download/v2.4.0/"$ARRIBA_TAR
MAKE_CURRENT=1
#################


DATE_VERSION=$(date '+%Y%m%d')

# STAR-Fusion
CTAT_ROOT=$DATABASES_FOLDER"/ctat"
CTAT_FOLDER=$CTAT_ROOT/$DATE_VERSION

if [ ! -d $CTAT_FOLDER ]; then
    mkdir -p $CTAT_FOLDER
    wget $CTAT_URL -P $TMP_DIR
    tar -xzf  $TMP_DIR/$CTAT_TAR -C $CTAT_FOLDER --strip-components=1
    rm -rf $TMP_DIR/$CTAT_TAR
else
    echo "[WARNING] "$CTAT_FOLDER" already exists. Please delete it manually if you need a fresh download"
fi;

if [ $MAKE_CURRENT ]; then
    [ ! -e $CTAT_ROOT/current ] || unlink $CTAT_ROOT/current
    ln -s $CTAT_FOLDER $CTAT_ROOT/current
fi;

# Arriba
ARRIBA_ROOT=$DATABASES_FOLDER"/arriba"
ARRIBA_FOLDER=$ARRIBA_ROOT/$DATE_VERSION
UNPACK_DIR=$(mktemp -d -p $TMP_DIR)

if [ ! -d $ARRIBA_FOLDER ]; then
    mkdir -p $ARRIBA_FOLDER
    wget $ARRIBA_URL -P $TMP_DIR
    tar -xzf  $TMP_DIR/$ARRIBA_TAR -C $UNPACK_DIR --strip-components=1
    mv $UNPACK_DIR/database/* $ARRIBA_FOLDER
    rm -rf $TMP_DIR/$ARRIBA_TAR $UNPACK_DIR
else
    echo "[WARNING] "$ARRIBA_FOLDER" already exists. Please delete it manually if you need a fresh download"
fi;

if [ $MAKE_CURRENT ]; then
    [ ! -e $ARRIBA_ROOT/current ] || unlink $ARRIBA_ROOT/current
    ln -s $ARRIBA_FOLDER $ARRIBA_ROOT/current 
fi;