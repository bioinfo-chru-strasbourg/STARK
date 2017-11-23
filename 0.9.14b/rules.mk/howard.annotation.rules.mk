############################
# HOWARD Annotation Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE="0.9.3b"
MK_DATE="10/05/2016"

## Release note
# 10/07/2015-V0.9b: Create HOWARD annotation and VCF translation
# 24/11/2015-V0.9.1b: Bug correction
# 22/04/2016-V0.9.2b: HOWARD and snpEff
# 10/05/2016-V0.9.3b: Add CORE annotation only and Minimal Annotation, and empty.vcf 


## DEBUG: Annotation only CORE TODO

#HOWARD
HOWARD?=$(NGSscripts)
HOWARD_CONFIG?=$(HOWARD)/config.ini
HOWARD_ANNOTATION?=$(HOWARD)/VCFannotation.pl
HOWARD_PRIORITIZATION?=$(HOWARD)/VCFprioritization.pl
HOWARD_TRANSLATION?=$(HOWARD)/VCFtranslation.pl
ANNOTATION_TYPE?="core,frequency,score,annotation,prediction,snpeff,snpeff_hgvs"
ANNOTATION_TYPE_MINIMAL?="Symbol,location,outcome,hgvs"
#ANNOTATIONS="PZScore,PZFlag,PZComment,Symbol,hgvs,location,outcome,AlleleFrequency,AD,dbSNP,dbSNPNonFlagged,DP,AF,VF,GQ,Ensembl,TI,FC,GWASCatalog,COSMIC,LocalDB,1000genomesALL,1000genomesEUR,6500NHLBIALL,6500NHLBIEUR,PolyPhen2HumanVarPred,PolyPhen2HumanDivPred,MutationTasterPred,MutationAssessorPred,LTRPred,IARCTP53,GGCPolymorphismsBRCA,SIFT,phastCons,PhyloP,SiPhy,FATHMM,LRT,GERP,PolyPhen2HumanVar,PolyPhen2HumanDiv,MutationTaster,MutationAssessor,TFBS,FilterComment,ALL"
ANNOVAR?=$(NGSscripts)
BDFOLDER?=$(NGSscripts)
HOWARD_CALCULATION?=VAF,NOMEN,VAF_STATS

# dev:
# $VCFTOOLS/vcf-merge empty.ann1.vcf.gz empty.ann2.vcf.gz > empty.ann1ann1vcfmerge.vcf
# $VCFTOOLS/vcf-subset -c sample  empty.ann1ann1vcfmerge.vcv > empty.ann1ann2.vcf

#%.howard.unsorted.vcf: %.vcf %.empty.vcf
%.howard.vcf: %.vcf %.empty.vcf
	# Annotation step
	#if [ "`grep -c -v ^# $<`" == "0" ]; then \
	#	cp $< $@; \
	#else \
	#	$(HOWARD_ANNOTATION) --input_file=$< --output_file=$@ --annotation=$(ANNOTATION_TYPE); \
	#fi;
	#$(HOWARD_ANNOTATION) --input_file=$< --output_file=$@ --annotation=$(ANNOTATION_TYPE) --annovar_folder=$(ANNOVAR) --annovar_databases=$(BDFOLDER) --verbose;
	mkdir -p $@.metrics 
	#$(HOWARD_ANNOTATION) --input_file=$< --output_file=$@ --annotation=$(ANNOTATION_TYPE) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --snpeff_stats=$@.metrics/$(@F).snpeff.metrics.html --verbose;
	+$(HOWARD) --input=$< --output=$@ --annotation=$(ANNOTATION_TYPE) --calculation=$(HOWARD_CALCULATION) --filter=$(HOWARD_FILTER) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --snpeff_stats=$@.metrics/$(@F).snpeff.metrics.html --multithreading --threads=$(THREADS) --snpeff_threads=$(THREADS_BY_SAMPLE) --tmp=$(TMP_FOLDER_TMP);
	#+$(HOWARD) --input=$< --output=$@ --annotation=$(ANNOTATION_TYPE) --filter=$(HOWARD_FILTER) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --snpeff_stats=$@.metrics/$(@F).snpeff.metrics.html --multithreading --threads=$(THREADS) --tmp=$(TMP_FOLDER_TMP) --verbose;
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
	$(HOWARD) --input=$< --output=$@ --annotation=core,snpeff_hgvs --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --tmp=$(TMP_FOLDER_TMP) --verbose;
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	#-rm $*.empty.vcf $@;


%.vcf: %.unMinimallyAnnotated.vcf %.empty.vcf
	# Annotation CORE and snpEff HGVS
	if [ "$(ANNOTATION_TYPE_MINIMAL)" != "" ]; then \
		#$(HOWARD_ANNOTATION) --input_file=$< --output_file=$@ --annotation=$(ANNOTATION_TYPE_MINIMAL) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --verbose; \
		$(HOWARD) --input=$< --output=$@ --annotation=$(ANNOTATION_TYPE_MINIMAL) --annovar_folder=$(ANNOVAR) --annovar_databases=$(ANNOVAR_DATABASES) --snpeff_jar=$(SNPEFF) --snpeff_databases=$(SNPEFF_DATABASES) --tmp=$(TMP_FOLDER_TMP) --verbose; \
	else \
		cp $< $@; \
	fi;
	-if [ ! -e $@ ]; then cp $*.empty.vcf $@; fi;
	#-rm $*.empty.vcf $@;


#%.txt: %.vcf
#	# Prioritization step
#	$(HOWARD_PRIORITIZATION) --config=$(HOWARD_CONFIG) --input_file=$< --output_file=$*.prioritized.vcf --filter=$(HOWARD_FILTER);
#	# Translation step
#	#$(HOWARD_TRANSLATION) --input_file=$*.prioritized.vcf --output_file=$@ --format=tab --annotation="$(ANNOTATIONS)" --sort_by="$(SORT_BY)" --order_by="$(ORDER_BY)"  --annovar_folder=$(ANNOVAR) --annovar_databases=$(BDFOLDER); # --noremove_filtered --output_file=$@ --output_format=tab
#	$(HOWARD_TRANSLATION) --input_file=$*.prioritized.vcf --output_file=$@ --format=tab --annotation="$(ANNOTATIONS)" --sort_by="$(SORT_BY)" --order_by="$(ORDER_BY)"  ; # --noremove_filtered --output_file=$@ --output_format=tab
#	# Cleaning
#	-rm $*.prioritized.vcf

#%.hard.vcf: %.vcf
#	$(HOWARD_PRIORITIZATION) --input_file=$< --output_file=$@ --hard;


# HARD FILTERING TODO or not
#ANNOTATION_FIELDS="PZScore,PZFlag,genesymbol,hgvs,location,outcome,AlleleFrequency,AD,dbSNP,dbSNPNonFlagged,DP,AF,VF,GQ,Ensembl,TI,FC,GWASCatalog,COSMIC,LocalDB,1000genomesALL,1000genomesEUR,6500NHLBIALL,6500NHLBIEUR,PolyPhen2HumanVarPred,PolyPhen2HumanDivPred,MutationTasterPred,MutationAssessorPred,LTRPred,IARCTP53,GGCPolymorphismsBRCA,SIFT,phastCons,PhyloP,SiPhy,FATHMM,LRT,GERP,PolyPhen2HumanVar,PolyPhen2HumanDiv,MutationTaster,MutationAssessor,TFBS,FilterComment,ALL"
#FILTERS_DEFAULT=filter:!=PASS|\\.:f,filter:=REJECT:-20,qual:<=20:f,qual:<=90:2,qual:<=99:5,qual:<=100:10,qual:<100:-17,DP:<10:f,DP:<30:-10,DP:<50:-5,GQ:<=20:f,\
#,dbSNPNonFlagged:rs:f,1000genomesALL:>$(MAF):f,1000genomesEUR:>$(MAF):f,6500NHLBIALL:>$(MAF):f,6500NHLBIEUR:>$(MAF):f\
#,location:=intronic:f,location:=intergenic:f,location:=UTR.*:f,location:=.*stream:f,location:=ncRNA.*:f\
#,location:=exonic.*:10,location:=.*splic.*:20\
#,outcome:=nonsynonymous.*:10,outcome:=frame.*:15,outcome:=stop.*:15,outcome:=nonframe.*:10,outcome:=synonymous.*:0\
#,AF:<0.05:f,FA:<0.05:f,AlleleFrequency:<0.05:f,TotalDepth:<10:f,ALTDepth:<10:f\
#,COSMIC:ID:20\
#,SIFT:&lt;=0.05:1,PolyPhen2HumanVarPred:=.*D:1,MutationTasterPred:=.*D:1,FATHMM:&lt;=-1.5:1
#
#%.hard.txt: %.vcf
#	$(NGSscripts)/VCF.pl --input_file=$< --annotation_list="$(ANNOTATION_FIELDS)" --filters="$(FILTERS_DEFAULT)" --order_by=FilterScore,desc --noremove_filtered --output_file=$@ --output_format=tab --remove_filtered;
#	$(HOWARD_TRANSLATION) --input_file=$*.prioritized.vcf --output_file=$@ --format=tab --annotation="$(ANNOTATIONS)" --sort_by="$(SORT_BY)" --order_by="$(ORDER_BY)" # --noremove_filtered --output_file=$@ --output_format=tab
#	if [ ! -s $@ ]; then touch $@; fi;
#
#%.hard.vcf: %.vcf
#	$(NGSscripts)/VCF.pl --input_file=$< --annotation_list="$(ANNOTATION_FIELDS)" --filters="$(FILTERS_DEFAULT)" --order_by=FilterScore,desc --noremove_filtered --output_file=$@ --output_format=vcf --remove_filtered;
#	if [ ! -s $@ ]; then touch $@; fi;


# CONFIG/RELEASE
#HOWARD_ANNOTATION_RELEASE:=$(shell $(HOWARD_ANNOTATION) --release | grep release | awk '{print $$3}')
#HOWARD_PRIORITIZATION_RELEASE:=$(shell $(HOWARD_PRIORITIZATION) --release | grep release | awk '{print $$3}')
#HOWARD_TRANSLATION_RELEASE:=$(shell $(HOWARD_TRANSLATION) --release | grep release | awk '{print $$3}')
RELEASE_COMMENT := "\#\# ANNOTATION '$(MK_RELEASE)': HOWARD annotates and prioritizes variants on a VCF, and generates *.howard.vcf and *.howard.txt files. Releases: '$(HOWARD_VERSION)'. Options: HOWARD_FILTER='$(HOWARD_FILTER)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)', ANNOTATION_TYPE=$(ANNOTATION_TYPE), ANNOTATIONS='$(ANNOTATIONS)'"
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "ANNOTATOR:howard:HOWARD annotates and prioritizes variants:HOWARD_FILTER='$(HOWARD_FILTER)', SORT_BY='$(SORT_BY)', ORDER_BY='$(ORDER_BY)', ANNOTATION_TYPE=$(ANNOTATION_TYPE), ANNOTATIONS='$(ANNOTATIONS)'"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )


