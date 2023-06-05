############################
# MUTECT Calling Rules
# Release: 0.9.2
# Date: 03/02/2023
# Author: Antony Le Bechec
############################


##########
# MUTECT #
##########

# Parameters
MUTECT_INTERVAL_PADDING?=0
DPMIN_MUTECT?=30


%.MuTect$(POST_CALLING).vcf: %.bam %.bam.bai %.empty.vcf %.genome %.design.bed.interval_listt
	# Calling
	$(JAVA) -jar $(MUTECT) \
		--analysis_type MuTect \
		--reference_sequence $$(cat $*.genome) \
		$$(if [ "`grep ^ -c $*.design.bed.interval_list`" == "0" ]; then echo ""; else echo "--intervals $*.design.bed.interval_list"; fi;) \
		--input_file:tumor $< \
		--vcf $@.tmp1 \
		-ip $(MUTECT_INTERVAL_PADDING);
	# Normalize
	grep "^##" $@.tmp1 > $@.tmp
	grep "^##" -v $@.tmp1 | grep "REJECT" -v | cut -f1-10 >> $@.tmp
	$(BCFTOOLS) view  -i 'FORMAT/DP>=$(DPMIN_MUTECT)' $@.tmp > $@.tmp2;
	# Empty
	if [ ! -e $@.tmp ]; then cp $*.empty.vcf $@.tmp2; fi;
	# Copy
	if [ ! -e $@ ]; then cp $@.tmp2 $@; fi;
	# Clean
	rm -f $@.tmp* $@.idx



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# CALLING: MUTECT to identify somatic variants and generate *.MuTect.vcf files. Parameteres MUTECT_INTERVAL_PADDING='$(MUTECT_INTERVAL_PADDING)' "
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:MuTect:MuTect - by default, no filtering. MUTECT_INTERVAL_PADDING as parameter"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
