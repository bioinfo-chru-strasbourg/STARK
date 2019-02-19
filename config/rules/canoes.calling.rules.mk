####################################################################
# CANOES Calling Rules
# Author: Vincent ZILLIOX, Sinthuja PACHCHEK, Antony LE BECHEC
####################################################################
# Release
MK_RELEASE="0.9"
MK_DATE="04/05/2018"

%.multicov: %.metrics.bed %.bam %.bam.bai
	mkdir -p $*.canoes
	$(BEDTOOLS)/bedtools multicov -bams $*.bam -bed $*.metrics.bed -q 20 > $@
	cp $*.metrics.bed $*.canoes
	cp $@ $*.canoes

%.multicovs.list: $(foreach RUN_SAMPLE,$(RUNS_SAMPLES),$(foreach ALIGNER,$(ALIGNERS),$(OUTDIR)/$(call run,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE))/$(call sample,$(RUN_SAMPLE)).$(ALIGNER).multicov ))
	mkdir -p $(@D)
	ls $^ > $@

%.canoes.vcf: %.empty.vcf %.beds.list %.bams.list %.SampleSheets.list %.multicovs.list
	mkdir -p $*.canoes
	$(STARK)/canoes_CNV.sh -o $*.canoes -s $*.SampleSheets.list -b $*.beds.list -r $(R) -t $(BEDTOOLS) -u $(GATK) -g /home1/TOOLS/genomes/hg19/hg19.fa -c $*.multicovs.list -a $(ANNOTCNV)
	cp $*.empty.vcf $@
	


RELEASE_CMD := $(shell echo "\#\# Canoes: CNV detection, generate *.canoe folder" >> $(RELEASE_INFOS) )

PIPELINES_COMMENT := "CALLER:canoes:CANOES - CNV caller"
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

