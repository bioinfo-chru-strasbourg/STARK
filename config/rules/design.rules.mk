############################
# Main Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE=1.0.0
MK_DATE=17/11/2026

# Release note
# 17/11/2026-1.0.0: Initial release of design rules (extracted from main.rules.mk)


# BED / INTERVALS
###################

# Manifest from SampleSheet to Manifest
%.manifest: %.manifest_from_samplesheet
	-mkdir -p $(@D);
	# Create MANIFEST file
	# 1. test if main manifest file (SAMPLE.manifest) exists (use for Sample analysis)
	# 2. test if a manifest is found on SampleSheet (use for Run analysis)
	# 3. test if a manifest is defined in the APP (as an environment variable)
	-if [ -e $(@D)/`echo $$(basename $(@D))`.manifest ]; then \
		echo "# MANIFEST for the sample '$(@D)/`echo $$(basename $(@D))`.manifest' exists" ; \
		if  [ "$(@D)/`echo $$(basename $(@D))`.manifest" != "$@" ]; then \
			echo "# MANIFEST for the sample '$(@D)/`echo $$(basename $(@D))`.manifest' exists and is different than the wanted MANIFEST file '$@'" ; \
			ln -s $(@D)/`echo $$(basename $(@D))`.manifest $@; \
			if [ ! -e $@ ]; then rm -f $@; cp $(@D)/`echo $$(basename $(@D))`.manifest $@; fi; \
		fi; \
	elif [ -s $*.manifest_from_samplesheet ]; then \
		echo "# MANIFEST for the SampleSheet found in the manifest" ; \
		cp $*.manifest_from_samplesheet $@; \
	elif [ -e "$(MANIFEST)" ] && [ "$(MANIFEST)" != "" ]; then \
		echo "# Input MANIFEST '$(MANIFEST)' exists" ; \
		cp $(MANIFEST) $@; \
	else \
		touch $@; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;


# Manifest from SampleSheet name to Manifest
%.manifest_name: %.manifest_from_samplesheet_name
	-mkdir -p $(@D)
	# Create MANIFEST file
	# 1. test if main manifest name file (SAMPLE.manifest_name) exists (use for Sample analysis)
	# 2. test if a manifest name is found on SampleSheet (use for Run analysis)
	-if [ -e $(@D)/`echo $$(basename $(@D))`.manifest_name ]; then \
		echo "# MANIFEST name for the sample '$(@D)/`echo $$(basename $(@D))`.manifest_name' exists" ; \
		if  [ "$(@D)/`echo $$(basename $(@D))`.manifest_name" != "$@" ]; then \
			echo "# MANIFEST name for the sample '$(@D)/`echo $$(basename $(@D))`.manifest_name' exists and is different than the wanted MANIFEST name file '$@'" ; \
			ln -s $(@D)/`echo $$(basename $(@D))`.manifest_name $@; \
			if [ ! -e $@ ]; then rm -f $@; cp $(@D)/`echo $$(basename $(@D))`.manifest_name $@; fi; \
		fi; \
	elif [ -e $(@D)/../`echo $$(basename $(@D))`.manifest_name ]; then \
		echo "# MANIFEST name for the sample '$(@D)/../`echo $$(basename $(@D))`.manifest_name' exists" ; \
		if  [ "$(@D)/../`echo $$(basename $(@D))`.manifest_name" != "$@" ]; then \
			echo "# MANIFEST name for the sample '$(@D)/../`echo $$(basename $(@D))`.manifest_name' exists and is different than the wanted MANIFEST name file '$@'" ; \
			ln -s $(@D)/../`echo $$(basename $(@D))`.manifest_name $@; \
			if [ ! -e $@ ]; then rm -f $@; cp $(@D)/../`echo $$(basename $(@D))`.manifest_name $@; fi; \
		fi; \
	elif [ -s $*.manifest_from_samplesheet_name ]; then \
		echo "# MANIFEST name for the SampleSheet found in the manifest" ; \
		cp $*.manifest_from_samplesheet_name $@; \
	elif [ -e "$(MANIFEST)" ] && [ "$(MANIFEST)" != "" ]; then \
		echo "# Input MANIFEST '$(MANIFEST)' exists" ; \
		echo -e "$(MANIFEST)\tfrom APP" > $@; \
	else \
		touch $@; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;


# SampleSheet to Manifest from SampleSheet
%.manifest_from_samplesheet: %.SampleSheet.csv
	# 1. Create the sample Manifest
	# Find the index of Manifest field
	# Find the Manifest index (A,B...)
	# Find the Manifest (blabla.manifest)
	# copy the manifest into the sample folder
	# clean
	# Create directory if needed
	-mkdir -p $(@D);
	# SampleSheet empty
	if [ ! -e $< ]; then touch $<; fi;
	# Manifest List
	-if [ -s $(INPUTDIR)/`echo $$(basename $$(dirname $(@D)))`/manifests_list.txt ]; then \
		cp -p $(INPUTDIR)/`echo $$(basename $$(dirname $(@D)))`/manifests_list.txt $*.manifests_list.txt; \
	else \
		touch $*.manifests_list.txt; \
	fi;
	if [ ! -e $@ ]; then touch $*.manifests_list.txt; fi;
	# Column number of sample manifest
	-grep -i ^Sample_ID $< | tr -d '\r\n' | sed -e 's/,/\n/g' | grep -i -n Manifest | cut -d \: -f 1 > $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I;
	-if [ ! -s $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I ]; then echo "1" > $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I; fi;
	# Manifest letter for the index
	-grep -e ^$$(basename $(@D)), $< | tail -n 1 | tr -d '\r\n' | cut -d \, -f `cat $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I` > $(@D)/$(*F).SAMPLE_MANIFEST_I;
	# Manifest name for the sample
	-grep ^`cat $(@D)/$(*F).SAMPLE_MANIFEST_I`, $*.manifests_list.txt | cut -d \, -f 2 > $(@D)/$(*F).SAMPLE_MANIFEST;
	# Found manifest, or Default manifest is the first in the list
	if [ "`cat $(@D)/$(*F).SAMPLE_MANIFEST`" != "" ] && [ -e "$(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST`" ]; then \
		echo "# TEST Found Manifest '"`cat $(@D)/$(*F).SAMPLE_MANIFEST`"' in '$(MANIFEST_FOLDER)' for sample `echo $$(basename $$(dirname $(@D)))`/$(*F)"; \
		cat "$(MANIFEST_FOLDER)/`cat $(@D)/$(*F).SAMPLE_MANIFEST`" > $@; \
	elif [ "`cut -d, -f2 $*.manifests_list.txt`" != "" ] && [ -e "$(MANIFEST_FOLDER)/`cut -d, -f2 $*.manifests_list.txt`" ]; then \
		echo "# TEST Manifest NOT found in '$(MANIFEST_FOLDER)'. Default Manifest used '"`cut -d, -f2 $*.manifests_list.txt`"' for sample `echo $$(basename $$(dirname $(@D)))`/$(*F)"; \
		cat "$(MANIFEST_FOLDER)/`cut -d, -f2 $*.manifests_list.txt`" > $@; \
	fi;
	#
	# IF root manifest exists
	if [ -s $(@D)/`echo $$(basename $(@D))`.manifest ] && [ "$(@D)/`echo $$(basename $(@D))`.manifest" != "$@" ]; then \
		cp -p $(@D)/`echo $$(basename $(@D))`.manifest $@; \
	fi;
	# Empty manifest if failed!
	if [ ! -e $@ ]; then touch $@; fi;
	# Touch manifest to use time of SampleSheet
	if [ -e $< ]; then \
		touch $@ -r $<; \
	fi;
	for manifest in `cut -d"," -f1 $*.manifests_list.txt`; do grep ",$$manifest," $< | cut -d"," -f1 | sed "s/$$/,$$manifest/g" >> $*.manifest_by_samples.txt; echo "-----" >> $*.manifest_by_samples.txt; done; \
	# Clean
	-rm -f $*.manifests_list.txt;
	# Remove intermediate files
	-rm -f $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I $(@D)/$(*F).SAMPLE_MANIFEST_I $(@D)/$(*F).SAMPLE_MANIFEST;


# SampleSheet to Manifest from SampleSheet Name
%.manifest_from_samplesheet_name: %.SampleSheet.csv
	# 1. Create the sample Manifest
	# Find the index of Manifest field
	# Find the Manifest index (A,B...)
	# Find the Manifest (blabla.manifest)
	# copy the manifest into the sample folder
	# clean
	# Create directory if needed
	-mkdir -p $(@D);
	# SampleSheet empty
	if [ ! -e $< ]; then touch $<; fi;
	# Manifest List
	-if [ -s $(INPUTDIR)/`echo $$(basename $$(dirname $(@D)))`/manifests_list.txt ]; then \
		cp -p $(INPUTDIR)/`echo $$(basename $$(dirname $(@D)))`/manifests_list.txt $*.manifests_list_name.txt; \
	else \
		touch $*.manifests_list_name.txt; \
	fi;
	if [ ! -e $@ ]; then touch $*.manifests_list_name.txt; fi;
	# Column number of sample manifest
	-grep -i ^Sample_ID $< | tr -d '\r\n' | sed -e 's/,/\n/g' | grep -i -n Manifest | cut -d \: -f 1 > $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I_name;
	-if [ ! -s $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I_name ]; then echo "1" > $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I_name; fi;
	# Manifest letter for the index
	-grep -e ^$$(basename $(@D)), $< | tail -n 1 | tr -d '\r\n' | cut -d \, -f `cat $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I_name` > $(@D)/$(*F).SAMPLE_MANIFEST_I_name;
	# Manifest name for the sample
	-grep ^`cat $(@D)/$(*F).SAMPLE_MANIFEST_I_name`, $*.manifests_list_name.txt | cut -d \, -f 2 > $(@D)/$(*F).SAMPLE_MANIFEST_name;
	-if [ ! -s $(@D)/$(*F).SAMPLE_MANIFEST_name ]; then \
		> $@; \
	else \
		echo -e $$(cat $(@D)/$(*F).SAMPLE_MANIFEST_name)"\tfrom SampleSheet" > $@; \
	fi;
	# Empty manifest if failed!
	if [ ! -e $@ ]; then touch $@; fi;
	# Clean
	-rm -f $*.manifests_list_name.txt;
	# Remove intermediate files
	-rm -f $(@D)/$(*F).INDEX_SAMPLEMANIFEST_I_name $(@D)/$(*F).SAMPLE_MANIFEST_I_name $(@D)/$(*F).SAMPLE_MANIFEST_name;


# Interval from bed from manifest?
%.from_manifest.interval_list: %.bed %.bam %.bam.bai %.bam.bed
	# try to extract from the bam if exists, in order to not call in the whome genome
	-+if [ -s $< ]; then \
		cp $< $@.bed; \
	elif [ -s $*.bam ]; then \
		echo "# Generation of BED from BAM because bed/manifest is empty"; \
		cp $*.bam.bed $@.bed ; \
	else \
		echo "[ERROR] No intervals generated '$@'" ; \
	fi;
	# INTERVAL WITH PICARD
	if [ -s $@.bed ]; then \
		echo "[INFO] Generate $@ from $@.bed with PICARD BedToIntervalList" ; \
		awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t"$$4}' $@.bed > $@.bed.4fields ; \
		$(JAVA) -jar $(PICARD) BedToIntervalList -I $@.bed.4fields -O $@ -SD $(DICT) ; \
		rm $@.bed.4fields ; \
	fi;
	# If error, try intervals with GREP/SED/AWK
	if [ ! -s $@ ] && [ -s $@.bed ]; then \
		echo "[INFO] Generate $@ from $@.bed with GREP/SED" ; \
		grep -v ^@ $@.bed  | tr -d '\r' | sed -e "s/^M//" | awk -F"\t" '{print $$1":"$$2"-"$$3}' | sed s/:0-/:1-/gi > $@; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;

# Interval from BED
%.bed.interval_list: %.bed
	# INTERVAL WITH PICARD
	if [ -s $< ]; then \
		awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t"$$4}' $< > $@.4fields ; \
		$(JAVA) -jar $(PICARD) BedToIntervalList -I $@.4fields -O $@ -SD $(DICT) ; \
		rm $@.4fields ; \
	fi;
	# If error, try intervals with GREP/SED/AWK
	if [ ! -s $@ ] && [ -s $< ]; then \
		grep -v ^@ $<  | tr -d '\r' | sed -e "s/^M//" | awk -F"\t" '{print $$1":"$$2"-"$$3}' | sed s/:0-/:1-/gi > $@; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;


# interval_list to intervals
%.intervals: %.interval_list
	touch $@;
	-grep -v ^@ $< > $@;


# BED primer file from a Manifest
%.primers.bed: %.manifest
	# Create BED file
	-if [ -e "$(PRIMER_BED)" ]; then \
		cp $(PRIMER_BED) $@; \
	elif [ -s $< ]; then \
		$(CAP_ManifestToBED) --input=$< --output=$@ --output_type=primer; \
	else \
		touch $@; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;


# BED file from a Manifest
%.bed: %.manifest
	# Create BED file
	# 1. test if main bed file (SAMPLE.bed) exists (use for Sample analysis)
	# 2. test if a bed can be generated from the manifest
	# 3. test if a manifest is defined in the APP (as an environment variable)
	-rm -f $@
	-if [ -e $(@D)/`echo $$(basename $(@D))`.bed ]; then \
		echo "# BED for the sample '$(@D)/`echo $$(basename $(@D))`.bed' exists" ; \
		if  [ "$(@D)/`echo $$(basename $(@D))`.bed" != "$@" ]; then \
			echo "# BED for the sample '$(@D)/`echo $$(basename $(@D))`.bed' exists and is different than the wanted BED file '$@'" ; \
			ln -s $(@D)/`echo $$(basename $(@D))`.bed $@; \
			if [ ! -e $@ ]; then rm -f $@; cp $(@D)/`echo $$(basename $(@D))`.bed $@; fi; \
		fi; \
	elif [ -s $< ]; then \
		echo "# BED for the sample generated from the manifest '$<'" ; \
		rm -f $@.tmp $@.sorted.tmp; \
		$(CAP_ManifestToBED) --input=$< --output=$@.tmp --output_type=region_clipped; \
		if (( $$(grep ^ $@.tmp -c) )); then \
			echo "region_clipped in Manifest is NOT empty "; \
			$(BEDTOOLS) sort -i $@.tmp | $(BEDTOOLS) merge -c 4,5 -o first,distinct  | awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t"$$5"\t0\t+"}' > $@; \
		else \
			echo "region_clipped in Manifest is empty "; \
			touch $@; \
		fi; \
		rm -f $@.tmp; \
	elif [ -e "$(BED)" ] && [ "$(BED)" != "" ]; then \
		echo "# Input BED '$(BED)' exists" ; \
		cp $(BED) $@; \
	else \
		touch $@; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;

# BED file from a Manifest
%.bed_name: %.manifest_name
	# Create BED file
	# 1. test if main bed file (SAMPLE.bed) exists (use for Sample analysis)
	# 2. test if a bed can be generated from the manifest
	# 3. test if a manifest is defined in the APP (as an environment variable)s
	#
	-if [ -e $(@D)/`echo $$(basename $(@D))`.bed_name ]; then \
		echo "# BED name for the sample '$(@D)/`echo $$(basename $(@D))`.bed_name' exists" ; \
		if  [ "$(@D)/`echo $$(basename $(@D))`.bed_name" != "$@" ]; then \
			echo "# BED name for the sample '$(@D)/`echo $$(basename $(@D))`.bed_name' exists and is different than the wanted BED name file '$@'" ; \
			ln -s $(@D)/`echo $$(basename $(@D))`.bed_name $@ || cp $(@D)/`echo $$(basename $(@D))`.bed_name $@; \
		fi; \
	elif [ -e $(@D)/../`echo $$(basename $(@D))`.bed_name ]; then \
		echo "# BED name for the sample '$(@D)/../`echo $$(basename $(@D))`.bed_name' exists" ; \
		if  [ "$(@D)/../`echo $$(basename $(@D))`.bed_name" != "$@" ]; then \
			echo "# BED name for the sample '$(@D)/../`echo $$(basename $(@D))`.bed_name' exists and is different than the wanted BED name file '$@'" ; \
			ln -s $(@D)/../`echo $$(basename $(@D))`.bed_name $@ || cp $(@D)/../`echo $$(basename $(@D))`.bed_name $@; \
		fi; \
	elif [ -s $< ]; then \
		echo "# BED name for the sample from the manifest '$<'" ; \
		echo -e "-\tfrom Manifest" > $@; \
	elif [ -e "$(BED)" ] && [ "$(BED)" != "" ]; then \
		echo "# Input BED '$(BED)' exists" ; \
		echo -e "$(BED)\tfrom APP" > $@; \
	else \
		touch $@; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;


# BED file from a Manifest
%.transcripts: %.manifest
	# Create BED file
	# Merge bacause of PCR manifest defined only Amplicons
	-if [ -e "$(TRANSCRIPTS)" ] && [ "$(TRANSCRIPTS)" != "" ]; then \
		echo "# Input TRANSCRIPTS '$(TRANSCRIPTS)' exists" ; \
		cp $(TRANSCRIPTS) $@; \
	elif [ -e $(@D)/`echo $$(basename $(@D))`.transcripts ]; then \
		echo "# TRANSCRIPTS for the sample '$(@D)/`echo $$(basename $(@D))`.transcripts' exists" ; \
		if  [ "$(@D)/`echo $$(basename $(@D))`.transcripts" != "$@" ]; then \
			echo "# TRANSCRIPTS for the sample '$(@D)/`echo $$(basename $(@D))`.transcripts' exists and is different than the wanted TRANSCRIPTS file '$@'" ; \
			ln -s $(@D)/`echo $$(basename $(@D))`.transcripts $@ || cp $(@D)/`echo $$(basename $(@D))`.transcripts $@; \
		fi; \
	elif [ -e $(@D)/../`echo $$(basename $$(dirname $(@D)))`.transcripts ]; then \
		echo "# TRANSCRIPTS for the sample '$(@D)/../`echo $$(basename $$(dirname $(@D)))`.transcripts' exists" ; \
		if  [ "$(@D)/../`echo $$(basename $$(dirname $(@D)))`.transcripts" != "$@" ]; then \
			echo "# TRANSCRIPTS for the sample '$(@D)/../`echo $$(basename $$(dirname $(@D)))`.transcripts' exists and is different than the wanted TRANSCRIPTS file '$@'" ; \
			ln -s $(@D)/../`echo $$(basename $$(dirname $(@D)))`.transcripts $@ || cp $(@D)/../`echo $$(basename $$(dirname $(@D)))`.transcripts $@; \
		fi; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;

# BED clipped region file from a Manifest
%.region_clipped.bed: %.manifest %.bed
	# Create BED file
	# Merge bacause of PCR manifest defined only Amplicons
	-if [ -e "$(REGIONCLIPPEDBED)" ]; then \
		echo "# Input Region Clipped BED '$(REGIONCLIPPEDBED)' exists" ; \
		cp $(REGIONCLIPPEDBED) $@; \
	elif [ -e $(@D)/`echo $$(basename $(@D))`.region_clipped.bed ]; then \
		echo "# Region Clipped BED for the sample '$(@D)/`echo $$(basename $(@D))`.region_clipped.bed' exists" ; \
		if  [ "$(@D)/`echo $$(basename $(@D))`.region_clipped.bed" != "$@" ]; then \
			echo "# Region Clipped BED for the sample '$(@D)/`echo $$(basename $(@D))`.region_clipped.bed' exists and is different than the wanted BED file '$@'" ; \
			cp $(@D)/`echo $$(basename $(@D))`.region_clipped.bed $@; \
		fi; \
	elif [ -s $< ]; then \
		echo "# Region Clipped BED for the sample generated from the manifest '$<'" ; \
		rm -f $@.tmp $@.sorted.tmp; \
		$(CAP_ManifestToBED) --input=$< --output=$@.tmp --output_type=region_clipped; \
		if (( $$(grep ^ $@.tmp -c) )); then \
			echo "region_clipped in Manifest is NOT empty "; \
			$(BEDTOOLS) sort -i $@.tmp | $(BEDTOOLS) merge -c 4,5 -o first,distinct | awk -F"\t" '{print $$1"\t"$$2"\t"$$3"\t"$$5"\t0\t"$$4}' > $@; \
		else \
			echo "region_clipped in Manifest is empty "; \
			touch $@; \
		fi; \
		rm -f $@.tmp; \
	else \
		echo "# Error creating region BED"; \
	fi;
	if [ "`grep ^ -c $@`" == "0" ]; then \
		echo "# Region Clipped BED generated is empty. BED will be used as Region Clipped BED" ; \
		cp $*.bed $@; \
	fi;
	if [ ! -e $@ ]; then touch $@; fi;

