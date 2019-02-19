#!/bin/bash
## STARK env for HUSHEMATO TSOMYELOID RUNs

# DEFAULT ENV
######################
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env.HUSHEMATO.sh

# APPLICATION INFOS
#####################

APP_NAME="HUSHEMATO.TSOMYELOID"
APP_RELEASE="1.0"


GROUP="HUSHEMATO"
PROJECT="TSOMYELOID"


# ANALYSIS PARAMETERS
#######################



# FOOTER
#############
#source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env_footer.sh
source $STARK_FOLDER/env_footer.sh




