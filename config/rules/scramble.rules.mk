############################
# Scramble Rules
# Authoring : Thomas LAVAUX, Antony Le Béchec
############################
# Release

MK_RELEASE=0.1.1.0
MK_DATE=15/09/2021

## Release note
# 0.1.0.0-15/09/2021 : DEV version
# 0.1.1.0-19/09/2021 : DEV version, splitted rules


### Options
# Generate cluster txt file
# Options : -m min soft clipped bases to be put in a clsuter defaut 10 / -s min soft clipped reads to be a cluster default 5 / -r region default all
# Call Rscript cluster_analyser 
# MEI options
#-n or --nCluster default=5 : min cluster size to analyze
#-m or --mei-score default=50 : min MEI alignment score to call
#--poly-a-frac default=0.75 : fraction of clipped length for calling polyA tail in MEIs
#--poly-a-dist default=100 : how far from MEI to look for polyA tail
# Option  --eval-dels to eval deletions ; the ref need to have .nhr, .nsq and .ndb files for BLAST to work ; it can be generate with the BLAST following command : makeblastdb -in ref.fa -dbtype nucl)
# Deletions options
#--indel-score default=80 : min indel alignment score to call
#--pct-align default=90 : percent alignment of clipped read for calling deletions
#--min-del-len default=50 : minimum deletion length to call


# Generate cluster file
%.scramble_clusters : %.bam %.bam.bai %.genome
	$(SCRAMBLE_FOLDER)/cluster_identifier $< > $@ ;


# Call Rscript cluster_analyser - MEIS
%.scramble$(POST_CALLING).vcf: %.scramble_clusters %.empty.vcf
	$(RSCRIPT) --vanilla $(SCRAMBLE) --out-name $@.tmp.results.vcf --cluster-file $< --install-dir $(SCRAMBLE_FOLDER) --mei-refs $(SCRAMBLE_FOLDER)/MEI_consensus_seqs.fa --ref $$(cat $*.genome) --eval-meis ;
	+$(HOWARD) --input=$@.tmp.results.vcf --output=$@ --translation=VCF --verbose;
	#$@.tmp*


# Call Rscript cluster_analyser - MEIS and DELS
%.scramble_meis_dels$(POST_CALLING).vcf: %.scramble_clusters %.empty.vcf
	$(RSCRIPT) --vanilla $(SCRAMBLE) --out-name $@.tmp.results.vcf --cluster-file $< --install-dir $(SCRAMBLE_FOLDER) --mei-refs $(SCRAMBLE_FOLDER)/MEI_consensus_seqs.fa --ref $$(cat $*.genome) --eval-meis --eval-dels ;
	+$(HOWARD) --input=$@.tmp.results.vcf --output=$@ --translation=VCF --verbose;
	#$@.tmp*


# Call Rscript cluster_analyser - MEIS
%.scramble_meis$(POST_CALLING).vcf: %.scramble_clusters %.empty.vcf
	$(RSCRIPT) --vanilla $(SCRAMBLE) --out-name $@.tmp.results.vcf --cluster-file $< --install-dir $(SCRAMBLE_FOLDER) --mei-refs $(SCRAMBLE_FOLDER)/MEI_consensus_seqs.fa --ref $$(cat $*.genome) --eval-meis ;
	+$(HOWARD) --input=$@.tmp.results.vcf --output=$@ --translation=VCF --verbose;
	#$@.tmp*


# Call Rscript cluster_analyser - DELS
%.scramble_dels$(POST_CALLING).vcf: %.scramble_clusters %.empty.vcf
	$(RSCRIPT) --vanilla $(SCRAMBLE) --out-name $@.tmp.results.vcf --cluster-file $< --install-dir $(SCRAMBLE_FOLDER) --mei-refs $(SCRAMBLE_FOLDER)/MEI_consensus_seqs.fa --ref $$(cat $*.genome) --eval-dels ;
	+$(HOWARD) --input=$@.tmp.results.vcf --output=$@ --translation=VCF --verbose;
	#$@.tmp*



# RELEASE_COMMENT := "\#\# Scramble: MEIs and DELs identification"
# RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )

# PIPELINES_COMMENT := "CALLING:scramble: MEIs identification"
# PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

# PIPELINES_COMMENT := "CALLING:scramble_meis_dels: MEIs and DELs identification"
# PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

# PIPELINES_COMMENT := "CALLING:scramble_meis: MEIs identification"
# PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )

# PIPELINES_COMMENT := "CALLING:scramble_dels: DELs identification"
# PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )
