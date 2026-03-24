############################
# STAR Aligner Rules
# Release: 1.0.0
# Date: 24/03/2026
# Author: Samuel Nicaise, Thomas Lavaux, Antony Le Béchec
############################
# Release
MK_RELEASE="1.0.0"
MK_DATE="24/03/2026"

# Release note
# 1.0.0-24/03/2026: Rule for STAR aligner with 1pass, with improved genome sharing, memory management and lockfile system, and raw bam file handling.


###################################
# STAR By Default (FROM FASTQ) #
###################################

# Option to keep raw BAM file from STAR (with splitNcigar for SNV calling). This option is useful to keep BAM without post-alignment steps (e.g. realignment, recalibration...).
STAR_KEEP_RAW_BAM?=0

# STAR Flags to set specific parameters (see STAR help)
STAR_FLAGS_1PASS?=--outSAMtype BAM SortedByCoordinate --chimOutJunctionFormat 1 --outSAMunmapped Within --outBAMcompression 0 --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 10 --chimOutType Junctions WithinBAM --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 --twopassMode None --quantMode TranscriptomeSAM GeneCounts --quantTranscriptomeSAMoutput BanSingleEnd --limitBAMsortRAM 1073741824


# Raw alignemnt with STAR
# File star_raw.bam is an intermediate file
%.star_1pass.star_raw.bam: %.R1$(POST_SEQUENCING).fastq.gz %.R2$(POST_SEQUENCING).fastq.gz $(GENOME_RNA).star.idx
	# Create metrics folder
	mkdir -p $*.star_1pass.bam.metrics/$(*F).star_raw.bam.metrics;
	# Load STAR genome into shared memory to speed up alignments. This rule is used as a prerequisite for STAR alignments to ensure that the genome is loaded into memory before launching STAR alignments.
	# Launch STAR genome load with lock to ensure that only one STAR genome load is done. Other STAR genome load will fail because of genome is already loaded, but this is not a problem because we just need to ensure that the genome is loaded into memory before launching STAR alignments
	$(PYTHON3) $(STARK_FOLDER_BIN)/concurrency.py launch \
		--cmd "$(STAR) --genomeDir $(GENOME_RNA).star.idx --genomeLoad LoadAndExit --outFileNamePrefix /dev/null" \
		--lockfile_prefix $$(echo $@ | xargs -0 dirname | xargs -0 dirname)/lockfile.load_star_1pass_genome. \
		--target $@ \
		--max_jobs 1;
	# Launch STAR alignment with lockfile to manage memory and avoid multiple STAR processes to run at the same time if not enough memory is available. The lockfile is created in the parent folder of the target (e.g. postalignment) to avoid multiple STAR processes to run at the same time if not enough memory is available. This allow also to not kill the rule if STAR alignment fail because of memory issue, create the endpoint file, to release memory. By default, the maximum number of concurrent STAR alignments is set to the number of threads available, which correspond to no limitation (same as without lockfile).
	$(PYTHON3) $(STARK_FOLDER_BIN)/concurrency.py launch \
		--cmd "$(STAR) --genomeDir $(GENOME_RNA).star.idx \
			--genomeLoad LoadAndKeep \
			--runThreadN $(THREADS_BY_ALIGNER) \
			--readFilesIn $*.R1$(POST_SEQUENCING).fastq.gz $*.R2$(POST_SEQUENCING).fastq.gz \
			--readFilesCommand zcat \
			--outFileNamePrefix $*.star_1pass.bam.metrics/$(*F).star_raw.bam.metrics/$(*F).star_raw. \
			--outSAMattrRGline ID:1 PL:ILLUMINA PU:PU LB:001 \"SM:$(*F)\" \
			$(STAR_FLAGS_1PASS)" \
		--lockfile_prefix $$(echo $@ | xargs -0 dirname | xargs -0 dirname)/lockfile.star_1pass.star_raw. \
		--target $@ \
		--max_jobs $(THREADS);
	# Create endppoint file to indicate that STAR alignment is done. This is used to manage the lockfile system and to ensure that the genome is unloaded from memory when all STAR alignments are done.
	touch $@.tmp.done; 
	# Catch output bam file
	if (( $(STAR_KEEP_RAW_BAM) )); then \
		cp $*.star_1pass.bam.metrics/$(*F).star_raw.bam.metrics/$(*F).star_raw.Aligned.sortedByCoord.out.bam $@; \
		cp $*.star_1pass.bam.metrics/$(*F).star_raw.bam.metrics/$(*F).star_raw.Chimeric.out.junction $*.star_1pass.junction; \
	else \
		mv $*.star_1pass.bam.metrics/$(*F).star_raw.bam.metrics/$(*F).star_raw.Aligned.sortedByCoord.out.bam $@; \
		cp $*.star_1pass.bam.metrics/$(*F).star_raw.bam.metrics/$(*F).star_raw.Chimeric.out.junction $*.star_1pass.junction; \
		rm -rf $*.star_1pass.bam.metrics/$(*F).star_raw.bam.metrics/*.bam; \
	fi;
	# Clean
	-rm -rf $*.star_1pass.bam.metrics/$(*F).star_raw.bam.metrics/$(*F).star_raw._STARgenome $*.star_1pass.bam.metrics/$(*F).star_raw.bam.metrics/$(*F).star_raw._STARpass1;
	# Check if all BAM files are generated, if yes, unload STAR genome from shared memory to free up memory for other processes. This is done to avoid keeping the STAR genome in memory if all alignments are done, and to free up memory for other processes that may need it.
	if (( $$(ls $$(echo $(BAM) | tr " " "\n" | grep "star_1pass.bam$$" | sort -u | sed 's/star_1pass.bam$$/star_1pass.star_raw.bam.tmp.done/gi') | wc -w) == $$(echo $$(echo $(BAM) | tr " " "\n" | grep "star_1pass.bam$$" | sort -u | sed 's/star_1pass.bam$$/star_1pass.star_raw.bam/gi') | wc -w) )); then \
		echo "All STAR BAM files are generated, unloading STAR genome from shared memory..."; \
		$(STAR) --genomeDir $(GENOME_RNA).star.idx --genomeLoad Remove --outFileNamePrefix /dev/null; \
		rm -rf $$(dirname $$(dirname $@))/*/*star_1pass.star_raw.bam.tmp.done; \
	else \
		echo "Not all STAR BAM files are generated yet, keeping STAR genome in shared memory..."; \
	fi;


# Alignement with STAR from raw alignement (with post alignment for SNV calling)
%.star_1pass$(POST_ALIGNMENT).bam: %.star_1pass.star_raw.bam
	ln -s $< $@


# DEVEL - TODO
# Ensure genome is really removed in shared memory in case of STAR commnand fail
# before=$(ipcs -m | awk '/0x/ {print $2}')
# STAR command...
# after=$(ipcs -m | awk '/0x/ {print $2}')
# ipcids=$(printf "%s\n" $after | grep -vxFf <(printf "%s\n" $before))
# for ipcid in $ipcids; do ipcrm -m $ipcid; done;


# CONFIG/RELEASE
RELEASE_COMMENT := "\#\# STAR ALIGNMENT '$(MK_RELEASE)': STAR generates an aligned BAM file from FASTQ files."
RELEASE_CMD := $(shell echo "$(RELEASE_COMMENT)" >> $(RELEASE_INFOS) )


# PIPELINES INFOS
PIPELINES_COMMENT := "ALIGNER:STAR:STAR - Excellent for RNA-Seq data. From FASTQ files."
PIPELINES_CMD := $(shell echo -e "$(PIPELINES_COMMENT)" >> $(PIPELINES_INFOS) )