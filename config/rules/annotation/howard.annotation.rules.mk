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


ANNOTATION_TYPE?="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs"
HOWARD_ANNOTATION?="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs"
HOWARD_ANNOTATION_MINIMAL?="Symbol,location,outcome,hgvs"
HOWARD_ANNOTATION_REPORT?="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs"
HOWARD_FILTER?="default"
HOWARD_PRIORITIZATION?="default"
HOWARD_NOMEN_FIELDS?="hgvs"
HOWARD_CALCULATION?=VAF,NOMEN,VAF_STATS,DP_STATS,VARTYPE


# RULES
########


# HOWARD ANNOTATION
%.howard$(POST_ANNOTATION).vcf: %.vcf %.empty.vcf %.transcripts %.genome
	# Annotation step
	#+$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$< --output=$@.tmp --annotation=$(HOWARD_ANNOTATION) --calculation=$(HOWARD_CALCULATION) --prioritization=$(HOWARD_PRIORITIZATION) --transcripts=$*.transcripts --nomen_fields=$(HOWARD_NOMEN_FIELDS) --norm=$$(cat $*.genome);
	+if [ "$(HOWARD_DEJAVU_ANNOTATION)" != "" ]; then \
		$(HOWARD) $(HOWARD_DEJAVU_CONFIG_OPTIONS) --input=$< --output=$@.tmp0 --annotation=$(HOWARD_DEJAVU_ANNOTATION) --norm=$$(cat $*.genome); \
	else \
		mv $< $@.tmp0; \
	fi;
	+$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.tmp0 --output=$@.tmp --annotation=$(HOWARD_ANNOTATION) --calculation=$(HOWARD_CALCULATION) --transcripts=$*.transcripts --nomen_fields=$(HOWARD_NOMEN_FIELDS) --norm=$$(cat $*.genome);
	+$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$@.tmp --output=$@ --prioritization=$(HOWARD_PRIORITIZATION_VARANK) --norm=$$(cat $*.genome);
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	# Downgrading VCF format 4.2 to 4.1
	-cat $@ | sed "s/##fileformat=VCFv4.2/##fileformat=VCFv4.1/" > $@.dowgrade.4.2.to.4.1.tmp;
	-rm -f $@;
	-mv $@.dowgrade.4.2.to.4.1.tmp $@;
	# clean
	rm -rf $@.tmp*

# HOWARD MINIMAL ANNOTATION
%.howard_minimal$(POST_ANNOTATION).vcf: %.vcf %.empty.vcf %.transcripts %.genome
	# Annotation step
	+$(HOWARD) $(HOWARD_CONFIG_OPTIONS) --input=$< --output=$@ --annotation=$(HOWARD_ANNOTATION_MINIMAL) --calculation=$(HOWARD_CALCULATION_MINIMAL) --prioritization=$(HOWARD_PRIORITIZATION_MINIMAL) --transcripts=$*.transcripts --nomen_fields=$(HOWARD_NOMEN_FIELDS)  --norm=$$(cat $*.genome);
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	# Downgrading VCF format 4.2 to 4.1
	-cat $@ | sed "s/##fileformat=VCFv4.2/##fileformat=VCFv4.1/" > $@.dowgrade.4.2.to.4.1.tmp;
	-rm -f $@;
	-mv $@.dowgrade.4.2.to.4.1.tmp $@;



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# HOWARD ANNOTATION '$(MK_RELEASE)': HOWARD annotates and prioritizes variants on a VCF, and generates *.howard.vcf file. Releases: '$(HOWARD_VERSION)'. Options: , HOWARD_ANNOTATION='$(HOWARD_ANNOTATION)', HOWARD_CALCULATION='$(HOWARD_CALCULATION)', HOWARD_PRIORITIZATION='$(HOWARD_PRIORITIZATION)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# HOWARD ANNOTATION '$(MK_RELEASE)': HOWARD MINIMAL annotates and prioritizes variants on a VCF, and generates *.howard_minimal.vcf file. Releases: '$(HOWARD_VERSION)'. Options: , HOWARD_ANNOTATION='$(HOWARD_ANNOTATION_MINIMAL)', HOWARD_CALCULATION='$(HOWARD_CALCULATION_MINIMAL)', HOWARD_PRIORITIZATION='$(HOWARD_PRIORITIZATION_MINIMAL)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "ANNOTATOR:howard:HOWARD annotates and prioritizes variants:HOWARD_ANNOTATION='$(HOWARD_ANNOTATION)', HOWARD_CALCULATION='$(HOWARD_CALCULATION)', HOWARD_PRIORITIZATION='$(HOWARD_PRIORITIZATION)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)', HOWARD_NOMEN_FIELDS='$(HOWARD_NOMEN_FIELDS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "ANNOTATOR:howard_minimal:HOWARD MINIMAL annotates and prioritizes variants:HOWARD_ANNOTATION='$(HOWARD_ANNOTATION_MINIMAL)', HOWARD_CALCULATION='$(HOWARD_CALCULATION_MINIMAL)', HOWARD_PRIORITIZATION='$(HOWARD_PRIORITIZATION_MINIMAL)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)', HOWARD_NOMEN_FIELDS='$(HOWARD_NOMEN_FIELDS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
