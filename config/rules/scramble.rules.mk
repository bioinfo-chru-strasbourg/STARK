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
SCRAMBLE_CLUSTER_OPTION?=
SCRAMBLE_OPTION?=


# Generate cluster file
%.scramble_clusters : %.bam %.bam.bai %.genome
	$(SCRAMBLE_FOLDER)/cluster_identifier $< $(SCRAMBLE_CLUSTER_OPTION) > $@ ;


# Call Rscript cluster_analyser - MEIS
%.scramble$(POST_CALLING).vcf: %.scramble_clusters %.empty.vcf
	-$(RSCRIPT) --vanilla $(SCRAMBLE) --out-name $@.tmp.results --cluster-file $< --install-dir $(SCRAMBLE_FOLDER) --mei-refs $(SCRAMBLE_FOLDER)/MEI_consensus_seqs.fa --ref $$(cat $*.genome) --eval-meis $(SCRAMBLE_OPTION);
	+if [ ! -e $@.tmp.results_MEIs.txt ]; then \
		exit_error; \
	else \
		if [ ! -e $@.tmp.results.vcf ] || ! (( $$(grep "^#" -v -c $@.tmp.results.vcf) )); then \
			cp $*.empty.vcf $@; \
		else \
			# add sample \
			(grep "^##" $@.tmp.results.vcf && echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' && grep "^#CHROM" $@.tmp.results.vcf | awk -v SAMPLE=$$(basename $(@D)) '{print $$0"\tFORMAT\t"SAMPLE}' && grep "^#" -v $@.tmp.results.vcf | awk '{print $$0"\tGT\t0/1"}' ) > $@.tmp.add_sample.vcf; \
			# normalize with HOWARD \
			$(HOWARD) --input=$@.tmp.add_sample.vcf --output=$@ --translation=VCF --verbose; \
		fi; \
	fi;
	rm -rf $@.tmp*


# Call Rscript cluster_analyser - MEIS
%.scramble_meis$(POST_CALLING).vcf: %.scramble_clusters %.empty.vcf
	-$(RSCRIPT) --vanilla $(SCRAMBLE) --out-name $@.tmp.results --cluster-file $< --install-dir $(SCRAMBLE_FOLDER) --mei-refs $(SCRAMBLE_FOLDER)/MEI_consensus_seqs.fa --ref $$(cat $*.genome) --eval-meis $(SCRAMBLE_OPTION);
	+if [ ! -e $@.tmp.results_MEIs.txt ]; then \
		exit_error; \
	else \
		if [ ! -e $@.tmp.results.vcf ] || ! (( $$(grep "^#" -v -c $@.tmp.results.vcf) )); then \
			cp $*.empty.vcf $@; \
		else \
			# add sample \
			(grep "^##" $@.tmp.results.vcf && echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' && grep "^#CHROM" $@.tmp.results.vcf | awk -v SAMPLE=$$(basename $(@D)) '{print $$0"\tFORMAT\t"SAMPLE}' && grep "^#" -v $@.tmp.results.vcf | awk '{print $$0"\tGT\t0/1"}' ) > $@.tmp.add_sample.vcf; \
			# normalize with HOWARD \
			$(HOWARD) --input=$@.tmp.add_sample.vcf --output=$@ --translation=VCF --verbose; \
		fi; \
	fi;
	rm -rf $@.tmp*


# # Call Rscript cluster_analyser - MEIS and DELS
# %.scramble_meis_dels$(POST_CALLING).vcf: %.scramble_clusters %.empty.vcf
#  	-$(RSCRIPT) --vanilla $(SCRAMBLE) --out-name $@.tmp.results --cluster-file $< --install-dir $(SCRAMBLE_FOLDER) --mei-refs $(SCRAMBLE_FOLDER)/MEI_consensus_seqs.fa --ref $$(cat $*.genome) --eval-meis --eval-dels $(SCRAMBLE_OPTION);
# 	+if [ ! -e $@.tmp.results_MEIs.txt ]; then \
# 		exit_error; \
# 	else \
# 		if [ ! -e $@.tmp.results.vcf ] || ! (( $$(grep "^#" -v -c $@.tmp.results.vcf) )); then \
# 			cp $*.empty.vcf $@; \
# 		else \
# 			# add sample \
# 			(grep "^##" $@.tmp.results.vcf && echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' && grep "^#CHROM" $@.tmp.results.vcf | awk '{print $$0"\tFORMAT\tHORIZON"}' && grep "^#" -v $@.tmp.results.vcf | awk '{print $0"\tGT\t0/1"}' ) > $@.tmp.add_sample.vcf; \
# 			# normalize with HOWARD \
# 			$(HOWARD) --input=$@.tmp.add_sample.vcf --output=$@ --translation=VCF --verbose; \
# 		fi; \
# 	fi;
# 	$@.tmp*


# # Call Rscript cluster_analyser - DELS
# %.scramble_dels$(POST_CALLING).vcf: %.scramble_clusters %.empty.vcf
#  	-$(RSCRIPT) --vanilla $(SCRAMBLE) --out-name $@.tmp.results --cluster-file $< --install-dir $(SCRAMBLE_FOLDER) --mei-refs $(SCRAMBLE_FOLDER)/MEI_consensus_seqs.fa --ref $$(cat $*.genome) --eval-dels $(SCRAMBLE_OPTION);
# 	+if [ ! -e $@.tmp.results_MEIs.txt ]; then \
# 		exit_error; \
# 	else \
# 		if [ ! -e $@.tmp.results.vcf ] || ! (( $$(grep "^#" -v -c $@.tmp.results.vcf) )); then \
# 			cp $*.empty.vcf $@; \
# 		else \
# 			# add sample \
# 			(grep "^##" $@.tmp.results.vcf && echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' && grep "^#CHROM" $@.tmp.results.vcf | awk '{print $$0"\tFORMAT\tHORIZON"}' && grep "^#" -v $@.tmp.results.vcf | awk '{print $0"\tGT\t0/1"}' ) > $@.tmp.add_sample.vcf; \
# 			# normalize with HOWARD \
# 			$(HOWARD) --input=$@.tmp.add_sample.vcf --output=$@ --translation=VCF --verbose; \
# 		fi; \
# 	fi;
# 	$@.tmp*



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
