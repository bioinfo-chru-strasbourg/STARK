

### Variables
##############


FASTQ_DEMULTIPLEXING_KEEP?=0
REMOVE_FASTQ_COMMENT?=0
FASTQ_COMPRESSION_LEVEL?=5
FASTP_PARAM?=
FASTP_THREADS_BY_SAMPLE?=1
FASTP_COMPRESSION_LEVEL?=1
ENABLE_ADAPTER_TRIMMING?=0
FASTQ_QUALITY_FILTERING?=0
POLY_G_MIN_LEN?=0
READ_LENGTH_REQUIRED?=0
UMI_LOC?=
UMI_RELOC?=
R1_RELOC?=
R2_RELOC?=
UMI_BARCODE_PATTERN?=
UMI_BARCODE_PATTERN_1?=
UMI_BARCODE_PATTERN_2?=
UMI_BARCODE_PATTERN_LENGTH?=0
UMI_BARCODE_PATTERN_1_LENGTH?=0
UMI_BARCODE_PATTERN_2_LENGTH?=0
FASTP?=$(FASTP)
GZ?=$(GZ)
UNGZ?=$(UNGZ)
FASTQ_CLEAN_HEADER?=$(FASTQ_CLEAN_HEADER)
RELOCATE_UMI?=$(RELOCATE_UMI)
FASTQ_PROCESSED_STEPS?=.fastq_reheader.sort.fastp.fastq_clean_header.compress

FASTQ_CLEAN_HEADER_PARAM=-v SAM_TAG=1 -v UMI_REFORMAT=1 -v UMI_TAG=1 -v UMI_LOC=$(UMI_LOC)

#.SECONDARY: %.input.R1.fastq.gz %.input.R2.fastq.gz %.compress.R1.fastq.gz %.compress.R2.fastq.gz %.fastp.R1.fastq.gz %.fastp.R2.fastq.gz %.fastq_clean_header.R1.fastq.gz %.fastq_clean_header.R2.fastq.gz %.relocate_umi.R1.fastq.gz %.relocate_umi.R2.fastq.gz

#.PRECIOUS: %.input.R1.fastq.gz %.input.R2.fastq.gz %.compress.R1.fastq.gz %.compress.R2.fastq.gz %.fastp.R1.fastq.gz %.fastp.R2.fastq.gz %.fastq_clean_header.R1.fastq.gz %.fastq_clean_header.R2.fastq.gz %.relocate_umi.R1.fastq.gz %.relocate_umi.R2.fastq.gz

### Mandatory rules
###################


### Demultiplexing keep

%.demultiplexing_keep.R1.fastq.gz: %.log
	if (($(FASTQ_DEMULTIPLEXING_KEEP))); then \
		for f in $(@D)/../*.R1.fastq.gz; do \
			cp -p $$f $$(echo $@ | sed s/demultiplexing_keep/demultiplexing/); \
		done; \
	fi;

%.demultiplexing_keep.R2.fastq.gz: %.log
	if (($(FASTQ_DEMULTIPLEXING_KEEP))); then \
		for f in $(@D)/../*.R2.fastq.gz; do \
			cp -p $$f $$(echo $@ | sed s/demultiplexing_keep/demultiplexing/); \
		done; \
	fi;

%.demultiplexing_keep.I1.fastq.gz: %.log
	if (($(FASTQ_DEMULTIPLEXING_KEEP))); then \
		for f in $(@D)/../*.I1.fastq.gz; do \
			cp -p $$f $$(echo $@ | sed s/demultiplexing_keep/demultiplexing/); \
		done; \
	fi;

%.demultiplexing_keep.I2.fastq.gz: %.log
	if (($(FASTQ_DEMULTIPLEXING_KEEP))); then \
		for f in $(@D)/../*.I2.fastq.gz; do \
			cp -p $$f $$(echo $@ | sed s/demultiplexing_keep/demultiplexing/); \
		done; \
	fi;

%.demultiplexing_keep.log: %.log %.demultiplexing_keep.R1.fastq.gz %.demultiplexing_keep.R2.fastq.gz %.demultiplexing_keep.I1.fastq.gz %.demultiplexing_keep.I2.fastq.gz
	echo 'demultiplexing_keep log' >> $@;



### Input

%.input.R1.fastq.gz: %.log
	if (($(REMOVE_FASTQ_COMMENT))); then \
		$(UNGZ) -c $(@D)/../$(*F).R1.fastq.gz | cut -d' ' -f1 | $(GZ) -1 -c > $@; \
	else \
		ln -s $(@D)/../$(*F).R1.fastq.gz $@; \
	fi;

%.input.R2.fastq.gz: %.log
	if (($(REMOVE_FASTQ_COMMENT))); then \
		$(UNGZ) -c $(@D)/../$(*F).R2.fastq.gz | cut -d' ' -f1 | $(GZ) -1 -c > $@; \
	else \
		ln -s $(@D)/../$(*F).R2.fastq.gz $@; \
	fi;

%.input.I1.fastq.gz: %.log
	if (($(REMOVE_FASTQ_COMMENT))); then \
		$(UNGZ) -c $(@D)/../$(*F).I1.fastq.gz | cut -d' ' -f1 | $(GZ) -1 -c > $@; \
	else \
		ln -s $(@D)/../$(*F).I1.fastq.gz $@; \
	fi;

%.input.I2.fastq.gz: %.log
	if (($(REMOVE_FASTQ_COMMENT))); then \
		$(UNGZ) -c $(@D)/../$(*F).I2.fastq.gz | cut -d' ' -f1 | $(GZ) -1 -c > $@; \
	else \
		ln -s $(@D)/../$(*F).I2.fastq.gz $@; \
	fi;

%.input.log: %.log %.input.R1.fastq.gz %.input.R2.fastq.gz %.input.I1.fastq.gz %.input.I2.fastq.gz
	echo 'input log !!!' > $@;



### sequencing

%.sequencing.R1.fastq.gz: %.input$(FASTQ_PROCESSED_STEPS).log %.demultiplexing_keep.log
	# Copy
	if ! [ $*.input$(FASTQ_PROCESSED_STEPS).R1.fastq.gz -ef $(@D)/../$(*F).R1.fastq.gz ]; then \
		rm -f $(@D)/../$(*F).R1.fastq.gz; \
		cp -p $*.input$(FASTQ_PROCESSED_STEPS).R1.fastq.gz $(@D)/../$(*F).R1.fastq.gz; \
	else \
		if [ -L $(@D)/../$(*F).R1.fastq.gz ]; then \
			cp -p $$(realpath $(@D)/../$(*F).R1.fastq.gz) $(@D)/../$(*F).R1.fastq.gz.tmp; \
			rm -f $(@D)/../$(*F).R1.fastq.gz; \
			mv $(@D)/../$(*F).R1.fastq.gz.tmp $(@D)/../$(*F).R1.fastq.gz; \
		fi; \
	fi;

%.sequencing.R2.fastq.gz: %.input$(FASTQ_PROCESSED_STEPS).log %.demultiplexing_keep.log
	# Copy
	if ! [ $*.input$(FASTQ_PROCESSED_STEPS).R2.fastq.gz -ef $(@D)/../$(*F).R2.fastq.gz ]; then \
		rm -f $(@D)/../$(*F).R2.fastq.gz; \
		cp -p $*.input$(FASTQ_PROCESSED_STEPS).R2.fastq.gz $(@D)/../$(*F).R2.fastq.gz; \
	else \
		if [ -L $(@D)/../$(*F).R2.fastq.gz ]; then \
			cp -p $$(realpath $(@D)/../$(*F).R2.fastq.gz) $(@D)/../$(*F).R2.fastq.gz.tmp; \
			rm -f $(@D)/../$(*F).R2.fastq.gz; \
			mv $(@D)/../$(*F).R2.fastq.gz.tmp $(@D)/../$(*F).R2.fastq.gz; \
		fi; \
	fi;

%.sequencing.I1.fastq.gz: %.input$(FASTQ_PROCESSED_STEPS).log %.demultiplexing_keep.log
	# Copy
	if ! [ $*.input$(FASTQ_PROCESSED_STEPS).I1.fastq.gz -ef $(@D)/../$(*F).I1.fastq.gz ]; then \
		rm -f $(@D)/../$(*F).I1.fastq.gz; \
		cp -p $*.input$(FASTQ_PROCESSED_STEPS).I1.fastq.gz $(@D)/../$(*F).I1.fastq.gz; \
	else \
		if [ -L $(@D)/../$(*F).I1.fastq.gz ]; then \
			cp -p $$(realpath $(@D)/../$(*F).I1.fastq.gz) $(@D)/../$(*F).I1.fastq.gz.tmp; \
			rm -f $(@D)/../$(*F).I1.fastq.gz; \
			mv $(@D)/../$(*F).I1.fastq.gz.tmp $(@D)/../$(*F).I1.fastq.gz; \
		fi; \
	fi;

%.sequencing.I2.fastq.gz: %.input$(FASTQ_PROCESSED_STEPS).log %.demultiplexing_keep.log
	# Copy
	if ! [ $*.input$(FASTQ_PROCESSED_STEPS).I2.fastq.gz -ef $(@D)/../$(*F).I2.fastq.gz ]; then \
		rm -f $(@D)/../$(*F).I2.fastq.gz; \
		cp -p $*.input$(FASTQ_PROCESSED_STEPS).I2.fastq.gz $(@D)/../$(*F).I2.fastq.gz; \
	else \
		if [ -L $(@D)/../$(*F).I2.fastq.gz ]; then \
			cp -p $$(realpath $(@D)/../$(*F).I2.fastq.gz) $(@D)/../$(*F).I2.fastq.gz.tmp; \
			rm -f $(@D)/../$(*F).I2.fastq.gz; \
			mv $(@D)/../$(*F).I2.fastq.gz.tmp $(@D)/../$(*F).I2.fastq.gz; \
		fi; \
	fi;

%.sequencing.log: %.sequencing.R1.fastq.gz %.sequencing.R2.fastq.gz %.sequencing.I1.fastq.gz %.sequencing.I2.fastq.gz %.demultiplexing_keep.log
	echo 'sequencing log with $<' > $@;
	-touch $$(ls $(@D)/../*bam 2>/dev/null) 2>/dev/null
	ls -l $(@D)/ >> $@;



### Optionnal rules
###################


### Compression

%.compress.R1.fastq.gz: %.R1.fastq.gz %.log 
	$(UNGZ) $< -c | $(GZ) -$(FASTQ_COMPRESSION_LEVEL) -c  > $@;

%.compress.R2.fastq.gz: %.R2.fastq.gz %.log 
	$(UNGZ) $< -c | $(GZ) -$(FASTQ_COMPRESSION_LEVEL) -c  > $@;

%.compress.I1.fastq.gz: %.I1.fastq.gz %.log
	$(UNGZ) $< -c | $(GZ) -$(FASTQ_COMPRESSION_LEVEL) -c  > $@;

%.compress.I2.fastq.gz: %.I2.fastq.gz %.log 
	$(UNGZ) $< -c | $(GZ) -$(FASTQ_COMPRESSION_LEVEL) -c  > $@;

%.compress.log: %.log %.compress.R1.fastq.gz %.compress.R2.fastq.gz %.compress.I1.fastq.gz %.compress.I2.fastq.gz
	#rm -rf $*.*.fastq.gz;
	echo 'compress log !!!' > $@;



### Compression

%.sort.R1.fastq.gz: %.R1.fastq.gz %.log
	mkdir -p $@.tmp.SORT;
	$(UNGZ) $< -c | tr "\t" " " | paste - - - - | sort -T $@.tmp.SORT -n -k1,1 -t " " | tr "\t" "\n" | $(GZ) -1 -c > $@;
	rm -rf $@.tmp*

%.sort.R2.fastq.gz: %.R2.fastq.gz %.log 
	mkdir -p $@.tmp.SORT;
	$(UNGZ) $< -c | tr "\t" " " | paste - - - - | sort -T $@.tmp.SORT -n -k1,1 -t " " | tr "\t" "\n" | $(GZ) -1 -c > $@;
	rm -rf $@.tmp*

%.sort.I1.fastq.gz: %.I1.fastq.gz %.log
	mkdir -p $@.tmp.SORT;
	$(UNGZ) $< -c | tr "\t" " " | paste - - - - | sort -T $@.tmp.SORT -n -k1,1 -t " " | tr "\t" "\n" | $(GZ) -1 -c > $@;
	rm -rf $@.tmp*

%.sort.I2.fastq.gz: %.I2.fastq.gz %.log 
	mkdir -p $@.tmp.SORT;
	$(UNGZ) $< -c | tr "\t" " " | paste - - - - | sort -T $@.tmp.SORT -n -k1,1 -t " " | tr "\t" "\n" | $(GZ) -1 -c > $@;
	rm -rf $@.tmp*

%.sort.log: %.log %.sort.R1.fastq.gz %.sort.R2.fastq.gz %.sort.I1.fastq.gz %.sort.I2.fastq.gz
	echo 'sort log !!!' > $@;




# zcat Sample3.R1.fastq.gz | tr "\t" " " | paste - - - - | sort -n -k1,1 -t " " | tr "\t" "\n" | gzip -1 -c > Sample3.R1.sorted.fastq.gz


### FASTP

%.fastp.R1.fastq.gz: %.log
	# FASTP parameters
	echo " $(FASTP_PARAM) " > $@.param;
	#echo " $(FASTP_PARAM) " > $@.param;
	#echo " --thread=$(FASTP_THREADS_BY_SAMPLE) " >> $@.param;
	echo " --thread=1 " >> $@.param;
	echo " --compression=$(FASTP_COMPRESSION_LEVEL) " >> $@.param;
	# UMI test
	#if [ "$(UMI_RELOC)" != "" ]; then
	if ! (( $$($(UNGZ) -c $*.R1.fastq.gz | head -n 1 | cut -d" " -f1 | awk -F: '{if ($$8!="") {print $$0}}' | wc -l) )); then \
		echo "UMI_RELOC "$(UMI_RELOC) >> $@.param.test; \
		$(UNGZ) -c $*.R1.fastq.gz | head -n 1 | cut -d" " -f2- >> $@.param.test; \
		$(UNGZ) -c $*.R1.fastq.gz | head -n 1 | cut -d" " -f2- | grep -P '[0-9]*:N:0:[^\t $$]*' >> $@.param.test; \
		if [[ "$(UMI_RELOC)" =~ .*index.* ]] && (( $$($(UNGZ) -c $*.R1.fastq.gz | head -n 1 | cut -d" " -f2- | grep -P '[0-9]*:N:0:[^\t $$]*' | wc -l) )); then \
			echo " --umi --umi_loc=$(UMI_RELOC) " >> $@.param; \
		fi; \
		if [[ "$(UMI_RELOC)" =~ .*read.* ]]; then \
			echo " --umi --umi_loc=$(UMI_RELOC) --umi_len=$(UMI_BARCODE_PATTERN_1_LENGTH) " >> $@.param; \
		fi; \
	fi;
	# Report
	echo " --html=$(@D)/$$(echo $(*F) | cut -d. -f1).fastp.html --json=$(@D)/$$(echo $(*F) | cut -d. -f1).fastp.json --report_title=$$(echo $(@D) | xargs dirname | xargs dirname | xargs basename)/$$(echo $(*F) | cut -d. -f1) " >> $@.param;
	# Paired-End or Single-End
	if (( $$($(UNGZ) -c $*.R2.fastq.gz | head -n 1 | wc -l) )); then \
		echo " --in1=$*.R1.fastq.gz --in2=$*.R2.fastq.gz " >> $@.param ; \
		echo " --out1=$*.fastp.R1.fastq.gz --out2=$*.fastp.R2.fastq.gz " >> $@.param ; \
		#if (($(DETECT_ADAPTER_FOR_PE))); then echo " --detect_adapter_for_pe " >> $@.param ; fi ; \
		#echo " --detect_adapter_for_pe " >> $@.param ; \
		if (($(ENABLE_ADAPTER_TRIMMING))); then echo " --detect_adapter_for_pe " >> $@.param ; fi ; \
	else \
		echo " --in1=$*.R1.fastq.gz " >> $@.param ; \
		echo " --out1=$*.fastp.R1.fastq.gz " >> $@.param ; \
	fi;
	# Read quality filtering
	if (($(FASTQ_QUALITY_FILTERING))); then  \
		echo " --cut_mean_quality=$(FASTQ_QUALITY_FILTERING) " >> $@.param; \
	else \
		echo " --disable_quality_filtering " >> $@.param; \
	fi;
	# POLY_G_MIN_LEN
	if (($(POLY_G_MIN_LEN))); then  \
		echo " --trim_poly_g --poly_g_min_len=$(POLY_G_MIN_LEN) " >> $@.param; \
	else \
		echo " --disable_trim_poly_g " >> $@.param; \
	fi;
	# READ_LENGTH_REQUIRED
	if (($(READ_LENGTH_REQUIRED))); then  \
		echo " --length_required=$(READ_LENGTH_REQUIRED) " >> $@.param; \
	else \
		echo " --disable_length_filtering " >> $@.param; \
	fi;
	# ENABLE_ADAPTER_TRIMMING
	if ! (($(ENABLE_ADAPTER_TRIMMING))); then  \
		echo " --disable_adapter_trimming " >> $@.param; \
	fi;
	# FASTP Process
	$(FASTP) $$(cat $@.param) 1>\$@.log 2>\$@.err
	if [ ! -e $*.fastp.R2.fastq.gz ]; then \
		ln -s $*.R2.fastq.gz $*.fastp.R2.fastq.gz; \
	fi;

%.fastp.R2.fastq.gz: %.log
	touch $@;

%.fastp.I1.fastq.gz: %.log
	ln -s $*.I1.fastq.gz $@;

%.fastp.I2.fastq.gz: %.log
	ln -s $*.I2.fastq.gz $@;

%.fastp.log: %.log %.fastp.R1.fastq.gz %.fastp.R2.fastq.gz %.fastp.I1.fastq.gz %.fastp.I2.fastq.gz
	echo 'fastp log' > $@;
	cat $*.fastp.R1.fastq.gz.log $*.fastp.R1.fastq.gz.err >> $@;



### UMITools

%.umi_tools.R1.fastq.gz: %.log
	if [ -z $(UMI_RELOC) ] || (( $$($(UNGZ) -c $*.R1.fastq.gz | head -n 1 | cut -d" " -f1 | awk -F: '{if ($$8!="") {print $$0}}' | wc -l) )); then \
		ln -s $*.R1.fastq.gz $@; \
	elif [ "$(UMI_LOC)" == "index1" ]; then \
		$(UMITOOLS) extract --bc-pattern=$(UMI_BARCODE_PATTERN_1) --stdin=$*.I1.fastq.gz --read2-in=$*.R1.fastq.gz --read2-stdout --compresslevel=1 --log=$@.log --error=$@.err | sed '/^@/ s/_\([A-Z]*\)/:\1/' | $(GZ) -1 -c > $@; \
	elif [ "$(UMI_LOC)" == "index2" ]; then \
		$(UMITOOLS) extract --bc-pattern=$(UMI_BARCODE_PATTERN_1) --stdin=$*.I2.fastq.gz --read2-in=$*.R1.fastq.gz --read2-stdout --compresslevel=1 --log=$@.log --error=$@.err | sed '/^@/ s/_\([A-Z]*\)/:\1/' | $(GZ) -1 -c > $@; \
	elif [ "$(UMI_LOC)" == "per_index" ]; then \
		$(UMITOOLS) extract --bc-pattern=$(UMI_BARCODE_PATTERN_1) --stdin=$*.I1.fastq.gz --read2-in=$*.R1.fastq.gz --read2-stdout --compresslevel=1 --log=$@.log --error=$@.err | sed '/^@/ s/_\([A-Z]*\)/:\1/' | $(GZ) -1 -c > $@; \
	elif [ "$(UMI_LOC)" == "read1" ]; then \
		$(UMITOOLS) extract --bc-pattern=$(UMI_BARCODE_PATTERN_1) --stdin=$*.R1.fastq.gz --read2-in=$*.R1.fastq.gz --read2-stdout --compresslevel=1 --log=$@.log --error=$@.err | sed '/^@/ s/_\([A-Z]*\)/:\1/' | $(GZ) -1 -c > $@; \
	elif [ "$(UMI_LOC)" == "read2" ]; then \
		$(UMITOOLS) extract --bc-pattern=$(UMI_BARCODE_PATTERN_1) --stdin=$*.R2.fastq.gz --read2-in=$*.R1.fastq.gz --read2-stdout --compresslevel=1 --log=$@.log --error=$@.err | sed '/^@/ s/_\([A-Z]*\)/:\1/' | $(GZ) -1 -c > $@; \
	elif [ "$(UMI_LOC)" == "per_read" ]; then \
		$(UMITOOLS) extract --extract-method=string --bc-pattern=$(UMI_BARCODE_PATTERN_1) --bc-pattern2=$(UMI_BARCODE_PATTERN_2) --stdin=$*.R1.fastq.gz --read2-in=$*.R2.fastq.gz --stdout=$@.tmp.R1.fastq.gz --read2-out=$@.tmp.R2.fastq.gz --compresslevel=1 --log=$@.log --error=$@.err; \
		$(UNGZ) -c $@.tmp.R1.fastq.gz | sed '/^@/ s/_\([A-Z]\{$(UMI_BARCODE_PATTERN_1_LENGTH)\}\)\([A-Z]\{$(UMI_BARCODE_PATTERN_2_LENGTH)\}\)/:\1-\2/' | $(GZ) -1 -c > $@; \
		$(UNGZ) -c $@.tmp.R2.fastq.gz | sed '/^@/ s/_\([A-Z]\{$(UMI_BARCODE_PATTERN_1_LENGTH)\}\)\([A-Z]\{$(UMI_BARCODE_PATTERN_2_LENGTH)\}\)/:\1-\2/' | $(GZ) -1 -c > $*.umi_tools.R2.fastq.gz; \
	else \
		echo '#[ERROR] UMI source unknown' > $@.err; \
		ln -s $*.R1.fastq.gz $@; \
	fi;
	rm -rf $@.tmp*

%.umi_tools.R2.fastq.gz: %.log
	if [ -z $(UMI_RELOC) ] || (( $$($(UNGZ) -c $*.R1.fastq.gz | head -n 1 | cut -d" " -f1 | awk -F: '{if ($$8!="") {print $$0}}' | wc -l) )); then \
		ln -s $*.R2.fastq.gz $@; \
	elif [ "$(UMI_LOC)" == "index1" ]; then \
		$(UMITOOLS) extract --bc-pattern=$(UMI_BARCODE_PATTERN_1) --stdin=$*.I1.fastq.gz --read2-in=$*.R2.fastq.gz --read2-stdout --compresslevel=1 --log=$@.log --error=$@.err | sed '/^@/ s/_\([A-Z]*\)/:\1/' | $(GZ) -1 -c > $@; \
	elif [ "$(UMI_LOC)" == "index2" ]; then \
		$(UMITOOLS) extract --bc-pattern=$(UMI_BARCODE_PATTERN_1) --stdin=$*.I2.fastq.gz --read2-in=$*.R2.fastq.gz --read2-stdout --compresslevel=1 --log=$@.log --error=$@.err | sed '/^@/ s/_\([A-Z]*\)/:\1/' | $(GZ) -1 -c > $@; \
	elif [ "$(UMI_LOC)" == "per_index" ]; then \
		$(UMITOOLS) extract --bc-pattern=$(UMI_BARCODE_PATTERN_1) --stdin=$*.I2.fastq.gz --read2-in=$*.R2.fastq.gz --read2-stdout --compresslevel=1 --log=$@.log --error=$@.err | sed '/^@/ s/_\([A-Z]*\)/:\1/' | $(GZ) -1 -c > $@; \
	elif [ "$(UMI_LOC)" == "read1" ]; then \
		$(UMITOOLS) extract --bc-pattern=$(UMI_BARCODE_PATTERN_1) --stdin=$*.R1.fastq.gz --read2-in=$*.R2.fastq.gz --read2-stdout --compresslevel=1 --log=$@.log --error=$@.err | sed '/^@/ s/_\([A-Z]*\)/:\1/' | $(GZ) -1 -c > $@; \
	elif [ "$(UMI_LOC)" == "read2" ]; then \
		$(UMITOOLS) extract --bc-pattern=$(UMI_BARCODE_PATTERN_1) --stdin=$*.R2.fastq.gz --read2-in=$*.R2.fastq.gz --read2-stdout --compresslevel=1 --log=$@.log --error=$@.err | sed '/^@/ s/_\([A-Z]*\)/:\1/' | $(GZ) -1 -c > $@; \
	elif [ "$(UMI_LOC)" == "per_read" ]; then \
		touch $@; \
	else \
		echo '#[ERROR] UMI source unknown' > $@.err; \
		ln -s $*.R2.fastq.gz $@; \
	fi;
	rm -rf $@.tmp*
	

%.umi_tools.I1.fastq.gz: %.log
	ln -s $*.I1.fastq.gz $@;

%.umi_tools.I2.fastq.gz: %.log
	ln -s $*.I2.fastq.gz $@;

%.umi_tools.log: %.log %.umi_tools.R1.fastq.gz %.umi_tools.R2.fastq.gz %.umi_tools.I1.fastq.gz %.umi_tools.I2.fastq.gz
	echo 'umi_tools log' > $@;
	


### Clean fastq header

%.fastq_clean_header.R1.fastq.gz: %.log
	$(UNGZ) -c $*.R1.fastq.gz | head -n1 > $@.tmp.source
	$(UNGZ) -c $*.R1.fastq.gz | head -n1 | awk $(FASTQ_CLEAN_HEADER_PARAM) -f $(FASTQ_CLEAN_HEADER) > $@.tmp.target
	if (( $$( diff $@.tmp.source $@.tmp.target | wc -l) )); then \
		$(UNGZ) -c $*.R1.fastq.gz | awk $(FASTQ_CLEAN_HEADER_PARAM) -f $(FASTQ_CLEAN_HEADER) | $(GZ) -1 -c > $@; \
	else \
		ln -s $*.R1.fastq.gz $@; \
	fi;
	rm -f $@.tmp.source $@.tmp.target 

%.fastq_clean_header.R2.fastq.gz: %.log
	$(UNGZ) -c $*.R2.fastq.gz | head -n1 > $@.tmp.source
	$(UNGZ) -c $*.R2.fastq.gz | head -n1 | awk $(FASTQ_CLEAN_HEADER_PARAM) -f $(FASTQ_CLEAN_HEADER) > $@.tmp.target
	if (( $$( diff $@.tmp.source $@.tmp.target | wc -l) )); then \
			$(UNGZ) -c $*.R2.fastq.gz | awk $(FASTQ_CLEAN_HEADER_PARAM) -f $(FASTQ_CLEAN_HEADER) | $(GZ) -1 -c > $@; \
	else \
		ln -s $*.R2.fastq.gz $@; \
	fi;
	rm -f $@.tmp.source $@.tmp.target 

%.fastq_clean_header.I1.fastq.gz: %.log
	$(UNGZ) -c $*.I1.fastq.gz | head -n1 > $@.tmp.source
	$(UNGZ) -c $*.I1.fastq.gz | head -n1 | awk $(FASTQ_CLEAN_HEADER_PARAM) -f $(FASTQ_CLEAN_HEADER) > $@.tmp.target
	if (( $$( diff $@.tmp.source $@.tmp.target | wc -l) )); then \
			$(UNGZ) -c $*.I1.fastq.gz | awk $(FASTQ_CLEAN_HEADER_PARAM) -f $(FASTQ_CLEAN_HEADER) | $(GZ) -1 -c > $@; \
	else \
		ln -s $*.I1.fastq.gz $@; \
	fi;
	rm -f $@.tmp.source $@.tmp.target 

%.fastq_clean_header.I2.fastq.gz: %.log
	$(UNGZ) -c $*.I2.fastq.gz | head -n1 > $@.tmp.source
	$(UNGZ) -c $*.I2.fastq.gz | head -n1 | awk $(FASTQ_CLEAN_HEADER_PARAM) -f $(FASTQ_CLEAN_HEADER) > $@.tmp.target
	if (( $$( diff $@.tmp.source $@.tmp.target | wc -l) )); then \
			$(UNGZ) -c $*.I2.fastq.gz | awk $(FASTQ_CLEAN_HEADER_PARAM) -f $(FASTQ_CLEAN_HEADER) | $(GZ) -1 -c > $@; \
	else \
		ln -s $*.I2.fastq.gz $@; \
	fi;
	rm -f $@.tmp.source $@.tmp.target 

%.fastq_clean_header.log: %.log %.fastq_clean_header.R1.fastq.gz %.fastq_clean_header.R2.fastq.gz %.fastq_clean_header.I1.fastq.gz %.fastq_clean_header.I2.fastq.gz
	echo 'fastq_clean_header log' > $@;
	



### Reheader
# if (( $(zcat /STARK/data/STARKData/RUN_TEST_UMI_index2/demultiplexing/RUN_TEST/Sample3UMI_S3_R1_001.fastq.gz  | head -n1 | grep -e "BC:" -e "RX" -c) )); then echo "exists BC or RX"; fi;

%.fastq_reheader.R1.fastq.gz: %.log
	#if (( $$($(UNGZ) -c $*.R1.fastq.gz | head -n1 | grep -e "BC:Z:" -e "RX:Z:" -c) )); then \
	if (( $$($(UNGZ) -c $*.R1.fastq.gz | head -n1 | wc -l) )) && ! (( $$($(UNGZ) -c $*.R1.fastq.gz | head -n1 | grep -e "BC:Z:" -e "RX:Z:" -c) )); then \
		echo "$(UNGZ) -c $*.R1.fastq.gz | head -n1 | tr '\t' ' ' > $@.tmp.read_head" > $@.tmp.read_head.cmd; \
		chmod u+x $@.tmp.read_head.cmd; \
		/bin/bash $@.tmp.read_head.cmd; \
		echo "paste <($(UNGZ) -c $*.R1.fastq.gz | head -n2 | tr '\t' ' ') <($(UNGZ) -c $*.I1.fastq.gz | head -n2) <($(UNGZ) -c $*.I2.fastq.gz | head -n2) | awk -F'\t' -v READ=1 -f $(FASTQ_REHEADER) | head -n1 > $@.tmp.read_rehead" > $@.tmp.read_rehead.cmd; \
		chmod u+x $@.tmp.read_rehead.cmd; \
		/bin/bash $@.tmp.read_rehead.cmd; \
		if [ "$$(cat $@.tmp.read_head)" != "$$(cat $@.tmp.read_rehead)" ]; then \
			echo 'paste <($(UNGZ) -c $*.R1.fastq.gz | tr "\t" " ") <($(UNGZ) -c $*.I1.fastq.gz) <($(UNGZ) -c $*.I2.fastq.gz) | awk -F"\t" -v READ=1 -f $(FASTQ_REHEADER) | $(GZ) -1 -c > $@' > $@.tmp.cmd; \
			chmod u+x $@.tmp.cmd; \
			/bin/bash $@.tmp.cmd; \
			#rm -f $@.tmp.cmd; \
		else \
			ln -s $*.R1.fastq.gz $@; \
		fi; \
	else \
		ln -s $*.R1.fastq.gz $@; \
	fi;
	#rm -rf $@.tmp*

%.fastq_reheader.R2.fastq.gz: %.log
	if (( $$($(UNGZ) -c $*.R2.fastq.gz | head -n1 | wc -l) )) && ! (( $$($(UNGZ) -c $*.R2.fastq.gz | head -n1 | grep -e "BC:Z:" -e "RX:Z:" -c) )); then \
		echo "$(UNGZ) -c $*.R2.fastq.gz | head -n1 | tr '\t' ' ' > $@.tmp.read_head" > $@.tmp.read_head.cmd; \
		chmod u+x $@.tmp.read_head.cmd; \
		/bin/bash $@.tmp.read_head.cmd; \
		echo "paste <($(UNGZ) -c $*.R2.fastq.gz | head -n2 | tr '\t' ' ') <($(UNGZ) -c $*.I1.fastq.gz | head -n2) <($(UNGZ) -c $*.I2.fastq.gz | head -n2) | awk -F'\t' -v READ=1 -f $(FASTQ_REHEADER) | head -n1 > $@.tmp.read_rehead" > $@.tmp.read_rehead.cmd; \
		chmod u+x $@.tmp.read_rehead.cmd; \
		/bin/bash $@.tmp.read_rehead.cmd; \
		if [ "$$(cat $@.tmp.read_head)" != "$$(cat $@.tmp.read_rehead)" ]; then \
			echo 'paste <($(UNGZ) -c $*.R2.fastq.gz | tr "\t" " ") <($(UNGZ) -c $*.I1.fastq.gz) <($(UNGZ) -c $*.I2.fastq.gz) | awk -F"\t" -v READ=1 -f $(FASTQ_REHEADER) | $(GZ) -1 -c > $@' > $@.tmp.cmd; \
			chmod u+x $@.tmp.cmd; \
			/bin/bash $@.tmp.cmd; \
			#rm -f $@.tmp.cmd; \
		else \
			ln -s $*.R2.fastq.gz $@; \
		fi; \
	else \
		ln -s $*.R2.fastq.gz $@; \
	fi;
	#rm -rf $@.tmp*

%.fastq_reheader.I1.fastq.gz: %.log
	if (( $$($(UNGZ) -c $*.I1.fastq.gz | head -n1 | wc -l) )) && ! (( $$($(UNGZ) -c $*.I1.fastq.gz | head -n1 | grep -e "BC:Z:" -e "RX:Z:" -c) )); then \
		echo "$(UNGZ) -c $*.I1.fastq.gz | head -n1 | tr '\t' ' ' > $@.tmp.read_head" > $@.tmp.read_head.cmd; \
		chmod u+x $@.tmp.read_head.cmd; \
		/bin/bash $@.tmp.read_head.cmd; \
		echo "paste <($(UNGZ) -c $*.I1.fastq.gz | head -n2 | tr '\t' ' ') <($(UNGZ) -c $*.I1.fastq.gz | head -n2) <($(UNGZ) -c $*.I2.fastq.gz | head -n2) | awk -F'\t' -v READ=1 -f $(FASTQ_REHEADER) | head -n1 > $@.tmp.read_rehead" > $@.tmp.read_rehead.cmd; \
		chmod u+x $@.tmp.read_rehead.cmd; \
		/bin/bash $@.tmp.read_rehead.cmd; \
		if [ "$$(cat $@.tmp.read_head)" != "$$(cat $@.tmp.read_rehead)" ]; then \
			echo 'paste <($(UNGZ) -c $*.I1.fastq.gz | tr "\t" " ") <($(UNGZ) -c $*.I1.fastq.gz) <($(UNGZ) -c $*.I2.fastq.gz) | awk -F"\t" -v READ=1 -f $(FASTQ_REHEADER) | $(GZ) -1 -c > $@' > $@.tmp.cmd; \
			chmod u+x $@.tmp.cmd; \
			/bin/bash $@.tmp.cmd; \
			#rm -f $@.tmp.cmd; \
		else \
			ln -s $*.I1.fastq.gz $@; \
		fi; \
	else \
		ln -s $*.I1.fastq.gz $@; \
	fi;
	#rm -rf $@.tmp*

%.fastq_reheader.I2.fastq.gz: %.log
	if (( $$($(UNGZ) -c $*.I2.fastq.gz | head -n1 | wc -l) )) && ! (( $$($(UNGZ) -c $*.I2.fastq.gz | head -n1 | grep -e "BC:Z:" -e "RX:Z:" -c) )); then \
		echo "$(UNGZ) -c $*.I2.fastq.gz | head -n1 | tr '\t' ' ' > $@.tmp.read_head" > $@.tmp.read_head.cmd; \
		chmod u+x $@.tmp.read_head.cmd; \
		/bin/bash $@.tmp.read_head.cmd; \
		echo "paste <($(UNGZ) -c $*.I2.fastq.gz | head -n2 | tr '\t' ' ') <($(UNGZ) -c $*.I1.fastq.gz | head -n2) <($(UNGZ) -c $*.I2.fastq.gz | head -n2) | awk -F'\t' -v READ=1 -f $(FASTQ_REHEADER) | head -n1 > $@.tmp.read_rehead" > $@.tmp.read_rehead.cmd; \
		chmod u+x $@.tmp.read_rehead.cmd; \
		/bin/bash $@.tmp.read_rehead.cmd; \
		if [ "$$(cat $@.tmp.read_head)" != "$$(cat $@.tmp.read_rehead)" ]; then \
			echo 'paste <($(UNGZ) -c $*.I2.fastq.gz | tr "\t" " ") <($(UNGZ) -c $*.I1.fastq.gz) <($(UNGZ) -c $*.I2.fastq.gz) | awk -F"\t" -v READ=1 -f $(FASTQ_REHEADER) | $(GZ) -1 -c > $@' > $@.tmp.cmd; \
			chmod u+x $@.tmp.cmd; \
			/bin/bash $@.tmp.cmd; \
			#rm -f $@.tmp.cmd; \
		else \
			ln -s $*.I2.fastq.gz $@; \
		fi; \
	else \
		ln -s $*.I2.fastq.gz $@; \
	fi;
	#rm -rf $@.tmp*

%.fastq_reheader.log: %.log %.fastq_reheader.R1.fastq.gz %.fastq_reheader.R2.fastq.gz %.fastq_reheader.I1.fastq.gz %.fastq_reheader.I2.fastq.gz
	echo 'fastq_reheader log' > $@;
	# for f in $**gz; do \
	# 	echo $$f >> $@; \
	# 	zcat $$f | grep -cvP '\S' >> $@; \
	# 	zcat $$f | grep -cP '\S' >> $@; \
	# 	zcat $$f | head -n4 >> $@; \
	# 	echo "empty:" >> $@; \
	# 	zcat $$f | grep -vP '\S' -B4 -A3 | head -n 4 >> $@; \
	# done;




