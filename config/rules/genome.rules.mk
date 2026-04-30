############################
# Main Rules 
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="1.0.0"
MK_DATE="30/04/2026"

# Release note
# 10/03/2015: New file. Generates a file containing reference genome location
# 23/09/2016: Cleaning
# 30/04/2026: Keep only .genome to provide genome path for retro compatibility on old rules


## Genome REF
###############


# write the genome file in a file
%.genome: #%.assembly
	echo $(GENOME) > $@;


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# GENOME '$(MK_RELEASE)': Write .genome file with reference genome path (deprecated)."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


	