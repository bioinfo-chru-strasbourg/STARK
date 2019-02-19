#!/bin/bash
## STARK env for HUSTUMSOL.MTPTHS RUNs

# DEFAULT ENV
######################
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env.HUSTUMSOL.sh

# APPLICATION INFOS
#####################

APP_NAME="HUSTUMSOL.MTPTHS"
APP_RELEASE="1.0"

GROUP="HUSTUMSOL"
PROJECT="MTPTHS"


# ANALYSIS PARAMETERS
#######################



# FOOTER
###########
#source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env_footer.sh
source $STARK_FOLDER/env_footer.sh



