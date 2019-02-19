############################
# Makefile Functions
# Release: 0.9.0
# Date: 20/03/2014
# Author: Antony Le Bechec
############################


# Retrieves the RUN ID
# Param:
#   1. String to parse in form 'RUN:SAMPLE:...'.
run = $(firstword $(subst $(SEP), ,$1))

# Retrieves the SAMPLE ID
# Param:
#   1. String to parse in form 'RUN:SAMPLE:...'.
sample = $(word 2,$(subst $(SEP), ,$1))

# Retrieves the ALIGNER
# Param:
#   1. String to parse in form 'ALIGNER.CALLER.ANNOTATOR'.
aligner = $(firstword $(subst ., ,$1))

# Retrieves the CALLER
# Param:
#   1. String to parse in form 'ALIGNER.CALLER.ANNOTATOR'.
caller = $(word 2,$(subst ., ,$1))

# Retrieves the ANNOTATOR
# Param:
#   1. String to parse in form 'ALIGNER.CALLER.ANNOTATOR'.
annotator = $(word 3,$(subst ., ,$1))

# Retrieves the SAMPLE ID
# Param:
#   1. String to parse in form 'RUN:SAMPLE:...'.
#project = $(firstword $(subst -, ,$(word 3,$(subst $(SEP), ,$1))))

# Retrieves the SAMPLE ID
# Param:
#   1. String to parse in form 'RUN.SAMPLE.ALIGNER.CALLER.ANNOTATOR'.
#fastq = $(shell ls ./RUN*)
#fastq = $1

