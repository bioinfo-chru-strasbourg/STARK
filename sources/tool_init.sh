#!/bin/bash
echo "#[INFO] TOOL $TOOL_NAME/$TOOL_VERSION" && \
export TOOL_SOURCE=$SOURCES/$SOURCES_FOLDER/tools/$TOOL_NAME/$TOOL_VERSION/$TOOL_TARBALL && \
export TOOL_SOURCE_REPO=$REPO/$SOURCES_FOLDER/tools/$TOOL_NAME/$TOOL_VERSION/$TOOL_TARBALL && \
export TOOL_SOURCE_BUILD=$SOURCES/$SOURCES_FOLDER/tools/$TOOL_NAME/$TOOL_VERSION/build && \
export TOOL_DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION && \
export PATH=$TOOL_DEST/bin:$PATH && \
# Get TOOL SOURCE && \
$GET_TOOL_SOURCE && \
# TOOL folder preparation \
echo "#[INFO] TOOL preparation" && \
mkdir -p $TOOL_SOURCE_BUILD && \
mkdir -p $TOOL_DEST/bin && \
ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current
