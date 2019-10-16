############################
# HOWARD Annotation Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.4.1b"
MK_DATE="27/09/2019"

## Release note
# 10/07/2015-V0.9b: Create HOWARD annotation and VCF translation
# 24/11/2015-V0.9.1b: Bug correction
# 22/04/2016-V0.9.2b: HOWARD and snpEff
# 10/05/2016-V0.9.3b: Add CORE annotation only and Minimal Annotation, and empty.vcf
# 02/10/2018-V0.9.4b: Modification of the HOWARD annotation
# 27/09/2019-V0.9.4.1b: Add HOWARD NOMEN field option


# HOWARD Variables
####################

HOWARD?=$(NGSscripts)
HOWARD_CONFIG?=$(HOWARD)/config.ini
HOWARD_CONFIG_PRIORITIZATION?="config.prioritization.ini"
HOWARD_CONFIG_ANNOTATION?="config.annotation.ini"
HOWARD_ANNOTATION?=$(HOWARD)/VCFannotation.pl
HOWARD_PRIORITIZATION?=$(HOWARD)/VCFprioritization.pl
HOWARD_TRANSLATION?=$(HOWARD)/VCFtranslation.pl
ANNOTATION_TYPE?="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs"
HOWARD_ANNOTATION?="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs"
HOWARD_ANNOTATION_MINIMAL?="Symbol,location,outcome,hgvs"
HOWARD_ANNOTATION_REPORT?="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs"
HOWARD_FILTER?="default"
HOWARD_PRIORITIZATION?="default"
HOWARD_NOMEN_FIELDS?="hgvs"
ANNOVAR?=$(NGSscripts)
BDFOLDER?=$(NGSscripts)
HOWARD_CALCULATION?=VAF,NOMEN,VAF_STATS,VARTYPE


# RULES
########


# HOWARD ANNOTATION
%.howard$(POST_ANNOTATION).vcf: %.vcf %.empty.vcf %.transcripts %.genome
	# Annotation step
	+$(HOWARD) --input=$< --output=$@ --transcripts=$*.transcripts --config=$(HOWARD_CONFIG) --config_annotation=$(HOWARD_CONFIG_ANNOTATION) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION) --annotation=$(HOWARD_ANNOTATION) --calculation=$(HOWARD_CALCULATION) --nomen_fields=$(HOWARD_NOMEN_FIELDS) --prioritization=$(HOWARD_PRIORITIZATION) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS) --norm=$$(cat $*.genome);
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	# Downgrading VCF format 4.2 to 4.1
	-cat $@ | sed "s/##fileformat=VCFv4.2/##fileformat=VCFv4.1/" > $@.dowgrade.4.2.to.4.1.tmp;
	-rm -f $@;
	-mv $@.dowgrade.4.2.to.4.1.tmp $@;

# HOWARD MINIMAL ANNOTATION
%.howard_minimal$(POST_ANNOTATION).vcf: %.vcf %.empty.vcf %.transcripts %.genome
	# Annotation step
	+$(HOWARD) --input=$< --output=$@ --transcripts=$*.transcripts --config=$(HOWARD_CONFIG) --config_annotation=$(HOWARD_CONFIG_PRIORITIZATION_MINIMAL) --config_prioritization=$(HOWARD_CONFIG_PRIORITIZATION_MINIMAL) --annotation=$(HOWARD_ANNOTATION_MINIMAL) --calculation=$(HOWARD_CALCULATION_MINIMAL) --nomen_fields=$(HOWARD_NOMEN_FIELDS) --prioritization=$(HOWARD_PRIORITIZATION_MINIMAL) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP) --env=$(CONFIG_TOOLS) --norm=$$(cat $*.genome);
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	# Downgrading VCF format 4.2 to 4.1
	-cat $@ | sed "s/##fileformat=VCFv4.2/##fileformat=VCFv4.1/" > $@.dowgrade.4.2.to.4.1.tmp;
	-rm -f $@;
	-mv $@.dowgrade.4.2.to.4.1.tmp $@;



# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# ANNOTATION '$(MK_RELEASE)': HOWARD annotates and prioritizes variants on a VCF, and generates *.howard.vcf file. Releases: '$(HOWARD_VERSION)'. Options: , HOWARD_ANNOTATION='$(HOWARD_ANNOTATION)', HOWARD_CALCULATION='$(HOWARD_CALCULATION)', HOWARD_PRIORITIZATION='$(HOWARD_PRIORITIZATION)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# ANNOTATION '$(MK_RELEASE)': HOWARD MINIMAL annotates and prioritizes variants on a VCF, and generates *.howard_minimal.vcf file. Releases: '$(HOWARD_VERSION)'. Options: , HOWARD_ANNOTATION='$(HOWARD_ANNOTATION_MINIMAL)', HOWARD_CALCULATION='$(HOWARD_CALCULATION_MINIMAL)', HOWARD_PRIORITIZATION='$(HOWARD_PRIORITIZATION_MINIMAL)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "ANNOTATOR:howard:HOWARD annotates and prioritizes variants:HOWARD_ANNOTATION='$(HOWARD_ANNOTATION)', HOWARD_CALCULATION='$(HOWARD_CALCULATION)', HOWARD_PRIORITIZATION='$(HOWARD_PRIORITIZATION)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)', HOWARD_NOMEN_FIELDS='$(HOWARD_NOMEN_FIELDS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "ANNOTATOR:howard_minimal:HOWARD MINIMAL annotates and prioritizes variants:HOWARD_ANNOTATION='$(HOWARD_ANNOTATION_MINIMAL)', HOWARD_CALCULATION='$(HOWARD_CALCULATION_MINIMAL)', HOWARD_PRIORITIZATION='$(HOWARD_PRIORITIZATION_MINIMAL)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)', HOWARD_NOMEN_FIELDS='$(HOWARD_NOMEN_FIELDS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
