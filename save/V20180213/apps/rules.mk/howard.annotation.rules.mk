############################
# HOWARD Annotation Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.4b"
MK_DATE="02/10/2018"

## Release note
# 10/07/2015-V0.9b: Create HOWARD annotation and VCF translation
# 24/11/2015-V0.9.1b: Bug correction
# 22/04/2016-V0.9.2b: HOWARD and snpEff
# 10/05/2016-V0.9.3b: Add CORE annotation only and Minimal Annotation, and empty.vcf 
# 02/10/2018-V0.9.4b: Modification of the HOWARD annotation


## DEBUG: Annotation only CORE TODO

#HOWARD
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
#ANNOTATIONS="PZScore,PZFlag,PZComment,Symbol,hgvs,location,outcome,AlleleFrequency,AD,dbSNP,dbSNPNonFlagged,DP,AF,VF,GQ,Ensembl,TI,FC,GWASCatalog,COSMIC,LocalDB,1000genomesALL,1000genomesEUR,6500NHLBIALL,6500NHLBIEUR,PolyPhen2HumanVarPred,PolyPhen2HumanDivPred,MutationTasterPred,MutationAssessorPred,LTRPred,IARCTP53,GGCPolymorphismsBRCA,SIFT,phastCons,PhyloP,SiPhy,FATHMM,LRT,GERP,PolyPhen2HumanVar,PolyPhen2HumanDiv,MutationTaster,MutationAssessor,TFBS,FilterComment,ALL"
ANNOVAR?=$(NGSscripts)
BDFOLDER?=$(NGSscripts)
HOWARD_CALCULATION?=VAF,NOMEN,VAF_STATS,VARTYPE

# dev:
# $VCFTOOLS/vcf-merge empty.ann1.vcf.gz empty.ann2.vcf.gz > empty.ann1ann1vcfmerge.vcf
# $VCFTOOLS/vcf-subset -c sample  empty.ann1ann1vcfmerge.vcv > empty.ann1ann2.vcf


# HOWARD ANNOTATION
%.howard.vcf: %.norm.vcf %.empty.vcf %.transcripts
	# Annotation step
	mkdir -p $@.metrics 
	+$(HOWARD) --input=$< --output=$@ --transcripts=$*.transcripts --config=$(HOWARD_CONFIG) --config_filter=$(HOWARD_CONFIG_PRIORITIZATION) --config_filter=$(HOWARD_CONFIG_ANNOTATION) --annotation=$(HOWARD_ANNOTATION) --calculation=$(HOWARD_CALCULATION) --filter=$(HOWARD_PRIORITIZATION) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP);
	#--snpeff_stats=$@.metrics/$(@F).snpeff.metrics.html
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	
	# STATS
	#-$(BCFTOOLS) stats $@ > $@.metrics/$(@F).bcftools.metrics

	# Downgrading VCF format 4.2 to 4.1
	-cat $@ | sed "s/##fileformat=VCFv4.2/##fileformat=VCFv4.1/" > $@.dowgrade.4.2.to.4.1.tmp;
	-rm -f $@;
	-mv $@.dowgrade.4.2.to.4.1.tmp $@;

# HOWARD MINIMAL ANNOTATION
%.howard_minimal.vcf: %.norm.vcf %.empty.vcf %.transcripts
	# Annotation step
	mkdir -p $@.metrics 
	+$(HOWARD) --input=$< --output=$@ --transcripts=$*.transcripts --config=$(HOWARD_CONFIG) --config_filter=$(HOWARD_CONFIG_PRIORITIZATION) --config_filter=$(HOWARD_CONFIG_ANNOTATION) --annotation=$(HOWARD_ANNOTATION_MINIMAL) --calculation=$(HOWARD_CALCULATION_MINIMAL) --filter=$(HOWARD_PRIORITIZATION_MINIMAL) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP);
	#--snpeff_stats=$@.metrics/$(@F).snpeff.metrics.html 
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	
	# STATS
	-$(BCFTOOLS) stats $@ > $@.metrics/$(@F).bcftools.metrics

	# Downgrading VCF format 4.2 to 4.1
	-cat $@ | sed "s/##fileformat=VCFv4.2/##fileformat=VCFv4.1/" > $@.dowgrade.4.2.to.4.1.tmp;
	-rm -f $@;
	-mv $@.dowgrade.4.2.to.4.1.tmp $@;

%.vcf: %.unannotatedCORE.vcf %.empty.vcf
	# Annotation CORE and snpEff HGVS
	#$(HOWARD_ANNOTATION) --input_file=$< --output_file=$@ --annotation=core,snpeff_hgvs --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --verbose;
	$(HOWARD) --input=$< --output=$@ --config=$(HOWARD_CONFIG) --config_filter=$(HOWARD_CONFIG_PRIORITIZATION) --config_filter=$(HOWARD_CONFIG_ANNOTATION) --annotation=core,snpeff_hgvs --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --tmp=$(TMP_FOLDER_TMP) --verbose;
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	#-rm $*.empty.vcf $@;


%.vcf: %.unMinimallyAnnotated.vcf %.empty.vcf
	# Annotation CORE and snpEff HGVS
	if [ "$(ANNOTATION_TYPE_MINIMAL)" != "" ]; then \
		#$(HOWARD_ANNOTATION) --input_file=$< --output_file=$@ --annotation=$(ANNOTATION_TYPE_MINIMAL) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --verbose; \
		$(HOWARD) --input=$< --output=$@ --config=$(HOWARD_CONFIG) --config_filter=$(HOWARD_CONFIG_PRIORITIZATION) --config_filter=$(HOWARD_CONFIG_ANNOTATION) --annotation=$(ANNOTATION_TYPE_MINIMAL) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --tmp=$(TMP_FOLDER_TMP) --verbose; \
	else \
		cp $< $@; \
	fi;
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	#-rm $*.empty.vcf $@;

# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# ANNOTATION '$(MK_RELEASE)': HOWARD annotates and prioritizes variants on a VCF, and generates *.howard.vcf file. Releases: '$(HOWARD_VERSION)'. Options: , HOWARD_ANNOTATION='$(HOWARD_ANNOTATION)', HOWARD_CALCULATION='$(HOWARD_CALCULATION)', HOWARD_PRIORITIZATION='$(HOWARD_PRIORITIZATION)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

RELEASE_COMMENT := "\#\# ANNOTATION '$(MK_RELEASE)': HOWARD MINIMAL annotates and prioritizes variants on a VCF, and generates *.howard_minimal.vcf file. Releases: '$(HOWARD_VERSION)'. Options: , HOWARD_ANNOTATION='$(HOWARD_ANNOTATION_MINIMAL)', HOWARD_CALCULATION='$(HOWARD_CALCULATION_MINIMAL)', HOWARD_PRIORITIZATION='$(HOWARD_PRIORITIZATION_MINIMAL)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


PIPELINES_COMMENT := "ANNOTATOR:howard:HOWARD annotates and prioritizes variants:HOWARD_ANNOTATION='$(HOWARD_ANNOTATION)', HOWARD_CALCULATION='$(HOWARD_CALCULATION)', HOWARD_PRIORITIZATION='$(HOWARD_PRIORITIZATION)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

PIPELINES_COMMENT := "ANNOTATOR:howard_minimal:HOWARD MINIMAL annotates and prioritizes variants:HOWARD_ANNOTATION='$(HOWARD_ANNOTATION_MINIMAL)', HOWARD_CALCULATION='$(HOWARD_CALCULATION)', HOWARD_PRIORITIZATION='$(HOWARD_PRIORITIZATION)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


