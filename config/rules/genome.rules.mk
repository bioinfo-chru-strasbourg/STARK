############################
# Main Rules 
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.1b"
MK_DATE="23/09/2016"

# Release note
# 10/03/2015: New file. Generates a file containing reference genome location
# 23/09/2016: Cleaning



## Genome REF
###############



%.assembly: %.manifest %.bed
	# Find assembly from bed
	if [ -s $*.bed ]; then grep -Po "\thg[0-9]*\t|\\\\hg[0-9]*\\\\|\\hg[0-9]*\." $*.bed | tr -d "\t\\\." | uniq > $@; fi;
	# Find assembly from bed
	if [ ! -s $@ ] && [ -s $*.manifest ]; then grep -Po "\thg[0-9]*\t|\\\\hg[0-9]*\\\\|\\hg[0-9]*\." $*.manifest | tr -d "\t\\\." | uniq > $@; fi;
	# if no assembly on manifest, default (from config or parameter)
	if [ ! -s $@ ]; then echo $(ASSEMBLY) > $@; fi;
	echo "# ASSEMBLY: "`cat $@`;

	
# write the genome file in a file
%.genome: %.assembly
	if [ -s $(REF) ]; then \
		echo $(REF) > $@; \
		echo "# GENOME FILE: "$(REF); \
	elif [ -s $(GENOMES)/`cat $<`/`cat $<`.fa ]; then \
		echo $(GENOMES)/`cat $<`/`cat $<`.fa > $@; \
		echo "# GENOME FILE: "$(GENOMES)/`cat $<`/`cat $<`.fa; \
	elif [ -s $(GENOMES)/current/`cat $<`.fa ]; then \
		echo $(GENOMES)/current/`cat $<`.fa > $@; \
		echo "# GENOME FILE: "$(GENOMES)/current/`cat $<`.fa; \
	else \
		echo "# NO GENOME FILE '"$(GENOMES)/`cat $<`/`cat $<`.fa"'"; \
	fi;
	
	# index the reference for bwa
	if [ ! -s `cat $@`.bwt ] \
	|| [ ! -s `cat $@`.ann ] \
	|| [ ! -s `cat $@`.amb ] \
	|| [ ! -s `cat $@`.pac ] \
	|| [ ! -s `cat $@`.sa ]; then \
		$(BWA) index -a bwtsw `cat $@`; \
	fi;
	
	# index the reference for samtools
	if [ ! -s `cat $@`.fai ]; then \
		$(SAMTOOLS) faidx `cat $@`; \
	fi;
	

%.dict: %.genome
	cat $< | sed 's/\.fa$$/.dict/g' > $@
	if [ ! -s $$(cat $@) ]; then \
		$(JAVA) -jar $(PICARD) CreateSequenceDictionary REFERENCE=$$(cat $<) OUTPUT=$$(cat $@); \
	fi;



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# GENOME '$(MK_RELEASE)': Check for genome assembly from parameters, either bed file, manifest file, or option."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


	
