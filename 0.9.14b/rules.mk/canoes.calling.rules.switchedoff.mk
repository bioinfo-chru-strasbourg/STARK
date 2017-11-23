############################
# CANOES Calling Rules
# Author: Sinthuja PACHCHEK
############################
# Release
MK_RELEASE="0.9.12b"
MK_DATE="23/12/2016"


%.canoes.vcf: %.bam %.bam.bai %.from_manifest.intervals %.empty.vcf %.genome %.bams.list %.bed %.SampleSheet.csv
	mkdir -p $*.canoe ;
	#pour les tests : chmod 0777 $*.canoe ;
	sample=`basename $* | cut -d"." -f1`; \
	#Dir_path=$$(dirname $(@D))/$$sample ; \
	output_canoes='$*.canoe' ; \
	$(STARK)/canoes_CNV.sh -i $$sample -o $$output_canoes -s $*.SampleSheet.csv -b $*.bed -r $(R) -t $(BEDTOOLS) -u $(GATK) -g `cat $*.genome` -c $*.bams.list -a $(ANNOTCNV); \
	# if we use export blank in env
	#$(NGSscripts)/canoes_CNV.sh -i $$sample -o $$output_canoes -s $$SampleSheet_path -b $*.bed -r $(R) -t $(BEDTOOLS) -u $(GATK) -g `cat $*.genome` -c $*.bams.list -n $(BLANK) -a $(ANNOTCNV); \
	cp $*.empty.vcf $@