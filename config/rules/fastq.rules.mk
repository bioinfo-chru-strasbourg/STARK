############################
# Unaligned BAM Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.8.0"
MK_DATE="13/04/2021"

# Release note
# 10/03/2015: change genome reference location, in the file %.genome
# 22/05/2015: Deal with empty FASTQ files
# 04/12/2015: Change rules for bcl2fastq demultiplexing parameters
# 25/08/2016: Change fastq file generation to include R1 and R2. Simplification of code. Add $PICARD_FLAGS with compression level 0. Test number of reads using samtools. Add compression 9 with samtools
# 23/09/2016: Change PICARD LIB to BIN due to new release and PICARD use (unique JAR file)
# 0.9.8.0-13/04/2021: clarify code and removing unused rules


## FASTQ  ##


%.R1.fastq.gz: $(NEEDED)
	### Rule needed
	### code below depreciated: old code looking for FASTQ Reads1 file from INPUTDIR (demultiplexing folder)
	# Create directory
	#-mkdir -p $(@D)
	# Contatenate all fastq.gz files
	#-cat $(INPUTDIR)/$$(echo $$(basename $$(dirname $(@D))))/$(*F)_S*_R1_*.fastq.gz $(INPUTDIR)/$$(echo $$(basename $$(dirname $(@D))))/*/$(*F)_S*_R1_*.fastq.gz > $@;
	


%.R2.fastq.gz: $(NEEDED)
	### Rule needed
	### code below depreciated: old code looking for FASTQ Reads2 file from INPUTDIR (demultiplexing folder)
	# Create directory
	#-mkdir -p $(@D)
	# Contatenate all fastq.gz files
	#-cat $(INPUTDIR)/$$(echo $$(basename $$(dirname $(@D))))/$(*F)_S*_R2_*.fastq.gz $(INPUTDIR)/$$(echo $$(basename $$(dirname $(@D))))/*/$(*F)_S*_R2_*.fastq.gz > $@;



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# FASTQ '$(MK_RELEASE)': Generate uniq FASTQ from two FASTQ Read1 and Read2 files."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )
