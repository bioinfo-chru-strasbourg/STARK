############################
# HOWARD Annotation Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.4.2"
MK_DATE="13/04/2021"

## Release note
# 10/07/2015-V0.9b: Create HOWARD annotation and VCF translation
# 24/11/2015-V0.9.1b: Bug correction
# 22/04/2016-V0.9.2b: HOWARD and snpEff
# 10/05/2016-V0.9.3b: Add CORE annotation only and Minimal Annotation, and empty.vcf
# 02/10/2018-V0.9.4b: Modification of the HOWARD annotation
# 27/09/2019-V0.9.4.1b: Add HOWARD NOMEN field option
# 13/04/2021-V0.9.4.2: Add HOWARD_CONFIG_OPTIONS


# HOWARD Variables
####################


HOWARD_ANNOTATION?="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs"
HOWARD_CALCULATION?=VAF,NOMEN,VAF_STATS,DP_STATS,VARTYPE
HOWARD_NOMEN_FIELDS?="hgvs"


# RULES
########


# HOWARD ANNOTATION
%.howard$(POST_ANNOTATION).vcf: %.vcf %.empty.vcf %.transcripts %.genome
	# Prevent comma in description in vcf header;
	$(STARK_FOLDER_BIN)/fix_vcf_header.sh $< $@.tmp0 "$(THREADS_BY_CALLER)" "$(BCFTOOLS)" "$(FIX_VCF_HEADER_REFORMAT)"
	# Annotation step DEJAVU (deprecated)
	# +if [ "$(HOWARD_DEJAVU_ANNOTATION)" != "" ]; then \
	# 	$(HOWARD) $(HOWARD_DEJAVU_CONFIG_OPTIONS) --input=$@.tmp0 --output=$@.tmp1 --annotation=$(HOWARD_DEJAVU_ANNOTATION) --norm=$$(cat $*.genome); \
	#	$(STARK_FOLDER_BIN)/fix_vcf_header.sh $@.tmp1 "$(THREADS_BY_CALLER)" "$(BCFTOOLS)"; \
	#	mv $@.tmp1 $@.tmp0; \
	# fi;
	# Annotation calculation step HOWARD
	+$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.tmp0 --output=$@ --annotation=$(HOWARD_ANNOTATION) --calculation=$(HOWARD_CALCULATION) --transcripts=$*.transcripts --nomen_fields=$(HOWARD_NOMEN_FIELDS) --norm=$$(cat $*.genome);
	# Prevent comma in description in vcf header
	$(STARK_FOLDER_BIN)/fix_vcf_header.sh $@ $@ "$(THREADS_BY_CALLER)" "$(BCFTOOLS)" "$(FIX_VCF_HEADER_REFORMAT)"
	# Clear
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	# Downgrading VCF format 4.2 to 4.1
	-cat $@ | sed "s/##fileformat=VCFv4.2/##fileformat=VCFv4.1/" > $@.dowgrade.4.2.to.4.1.tmp;
	-rm -f $@;
	-mv $@.dowgrade.4.2.to.4.1.tmp $@;
	# clean
	rm -rf $@.tmp*



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# HOWARD ANNOTATION '$(MK_RELEASE)': HOWARD annotates and prioritizes variants on a VCF, and generates *.howard.vcf file. Releases: '$(HOWARD_VERSION)'. Options: , HOWARD_ANNOTATION='$(HOWARD_ANNOTATION)', HOWARD_CALCULATION='$(HOWARD_CALCULATION)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "ANNOTATOR:howard:HOWARD annotates and prioritizes variants:HOWARD_ANNOTATION='$(HOWARD_ANNOTATION)', HOWARD_CALCULATION='$(HOWARD_CALCULATION)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)', HOWARD_NOMEN_FIELDS='$(HOWARD_NOMEN_FIELDS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

