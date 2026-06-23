############################
# HOWARD Annotation Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.5.0"
MK_DATE="25/06/2025"

## Release note
# 10/07/2015-V0.9b: Create HOWARD annotation and VCF translation
# 24/11/2015-V0.9.1b: Bug correction
# 22/04/2016-V0.9.2b: HOWARD and snpEff
# 10/05/2016-V0.9.3b: Add CORE annotation only and Minimal Annotation, and empty.vcf
# 02/10/2018-V0.9.4b: Modification of the HOWARD annotation
# 27/09/2019-V0.9.4.1b: Add HOWARD NOMEN field option
# 13/04/2021-V0.9.4.2: Add HOWARD_CONFIG_OPTIONS
# 25/06/2025-V0.9.5.0: STARK release 19 compatibility


# HOWARD Variables
####################

HOWARD_ANNOTATION?="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs"
HOWARD_CALCULATION?=VAF,NOMEN,VAF_STATS,DP_STATS,VARTYPE
HOWARD_NOMEN_FIELDS?="hgvs"

HOWARD_ANNOTATION?="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs"
HOWARD_CALCULATION?=VAF,NOMEN,VAF_STATS,DP_STATS,VARTYPE
HOWARD_NOMEN_FIELDS?="hgvs"


# RULES
########

# HOWARD ANNOTATION
%.howard$(POST_ANNOTATION).vcf: %.vcf %.empty.vcf %.transcripts 
	# Annotation calculation step HOWARD
	if (( $$(grep -v "^#" $< | head -n1 | wc -l) )); then \
		$(HOWARD) process $(HOWARD_CONFIG_OPTIONS) $(HOWARD_PRIORITIZATION_CONFIG_OPTIONS) $(HOWARD_CALCULATION_CONFIG_OPTIONS) --input=$< --output=$@ --param=$(HOWARD_PARAM); \
	else \
		cp $< $@; \
	fi;
	# Clean INFO spaces
	$(STARK_FOLDER_BIN)/clean_vcf_info_spaces.sh --input=$@
	# Empty
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# HOWARD ANNOTATION '$(MK_RELEASE)': HOWARD annotates and prioritizes variants on a VCF, and generates *.howard.vcf file. Releases: '$(HOWARD_VERSION)'. Options: , HOWARD_ANNOTATION='$(HOWARD_ANNOTATION)', HOWARD_CALCULATION='$(HOWARD_CALCULATION)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "ANNOTATOR:howard:HOWARD annotates and prioritizes variants:HOWARD_ANNOTATION='$(HOWARD_ANNOTATION)', HOWARD_CALCULATION='$(HOWARD_CALCULATION)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)', HOWARD_NOMEN_FIELDS='$(HOWARD_NOMEN_FIELDS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

