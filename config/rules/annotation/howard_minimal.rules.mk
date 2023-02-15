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


HOWARD_ANNOTATION_MINIMAL?="Symbol,location,outcome,hgvs"
HOWARD_CALCULATION?=VAF,NOMEN,VAF_STATS,DP_STATS,VARTYPE
HOWARD_NOMEN_FIELDS?="hgvs"


# RULES
########


# HOWARD MINIMAL ANNOTATION
%.howard_minimal$(POST_ANNOTATION).vcf: %.vcf %.empty.vcf %.transcripts %.genome
	# Prevent comma in description in vcf header
	$(STARK_FOLDER_BIN)/fix_vcf_header.sh --input=$< --output=$@.tmp --threads=$(THREADS_BY_CALLER) --bcftools=$(BCFTOOLS) $(FIX_VCF_HEADER_REFORMAT_option);
	# Annotation step
	+$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.tmp --output=$@ --annotation=$(HOWARD_ANNOTATION_MINIMAL) --calculation=$(HOWARD_CALCULATION_MINIMAL) --transcripts=$*.transcripts --nomen_fields=$(HOWARD_NOMEN_FIELDS)  --norm=$$(cat $*.genome);
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	# Downgrading VCF format 4.2 to 4.1
	-cat $@ | sed "s/##fileformat=VCFv4.2/##fileformat=VCFv4.1/" > $@.dowgrade.4.2.to.4.1.tmp;
	-rm -f $@ $@.tmp*;
	-mv $@.dowgrade.4.2.to.4.1.tmp $@;



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# HOWARD ANNOTATION '$(MK_RELEASE)': HOWARD MINIMAL annotates and prioritizes variants on a VCF, and generates *.howard_minimal.vcf file. Releases: '$(HOWARD_VERSION)'. Options: , HOWARD_ANNOTATION='$(HOWARD_ANNOTATION_MINIMAL)', HOWARD_CALCULATION='$(HOWARD_CALCULATION_MINIMAL)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "ANNOTATOR:howard_minimal:HOWARD MINIMAL annotates and prioritizes variants:HOWARD_ANNOTATION='$(HOWARD_ANNOTATION_MINIMAL)', HOWARD_CALCULATION='$(HOWARD_CALCULATION_MINIMAL)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)', HOWARD_NOMEN_FIELDS='$(HOWARD_NOMEN_FIELDS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
