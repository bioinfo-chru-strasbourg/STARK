#!/bin/bash
## STARK env for SOLIDTUMOR

# DEFAULT ENV
######################
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env.sh

# APPLICATION INFOS
#####################
APP_NAME=SOLIDTUMOR
GROUP=SOLIDTUMOR
PROJECT=SOLIDTUMOR

# ANALYSIS PARAMETERS
#######################

# ALIGNERS CALLERS ANNOTATORS
ALIGNERS="bwamem"
CALLERS="gatkUG_SOLIDTUMOR VarScan_SOLIDTUMOR"
ANNOTATORS="howard"
# PIPELINES
PIPELINES="bwamem.gatkUG_SOLIDTUMOR.howard bwamem.VarScan_SOLIDTUMOR.howard"

# NB_BASES_AROUND
NB_BASES_AROUND=0

# HOWARD FILTER
HOWARD_FILTER=SOLIDTUMOR

# FOOTER
###################
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/env_footer.sh


