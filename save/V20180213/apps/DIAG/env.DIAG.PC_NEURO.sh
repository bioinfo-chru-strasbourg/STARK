	#!/bin/bash

# DEFAULT ENV
######################
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env.DIAG.sh

# APPLICATION INFOS
#####################

# Name of the APP. usefull to include specific rules if exists (in $STARK/$APP_NAME.rules.mk/*rules.mk)
# AUTO detect: $(basename ${BASH_SOURCE[0]} | sed 's/^env//' | sed 's/\.sh$//gi' | cut -d. -f2-)
APP_NAME="DIAG.PC_NEURO"

GROUP="DIAG"
PROJECT="PC_NEURO"


# ANALYSIS PARAMETERS
#######################


# FOOTER
###########
#source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env_footer.sh
source $STARK_FOLDER/env_footer.sh





