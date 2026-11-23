############################
# SnpEff Annotation Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE=0.9.1.0
MK_DATE=25/06/2025

## Release note
# 21/04/2016-V0.9b: Create SNPEFF
# 25/06/2025-V0.9.1.0: STARK release 19 compatibility


##########
# SNPEFF #
##########

# SNPEFF Flags
SNPEFF_FLAGS?=

# SNPEFF Stats 
SNPEFF_STATS?=0
# Generate command option if SNPEFF_STATS is 1
SNPEFF_STATS_OPTION:=$(if $(filter 1,$(SNPEFF_STATS)),-stats $@.stats/$(@F).stats.html,)


# Rule to generate SNPEFF annotated VCF files
%.snpeff$(POST_ANNOTATION).vcf: %.vcf  %.empty.vcf
	# Create stats directory
	if [ $(SNPEFF_STATS) -eq 1 ]; then \
		mkdir -p $@.stats; \
	fi;
	# Annotation step
	$(JAVA) -jar $(SNPEFF) $(SNPEFF_FLAGS) $(ASSEMBLY) $< -v $(SNPEFF_STATS_OPTION) > $@
	# Empty
	if [ ! -e $@ ]; then \
		echo "Creating empty VCF file $@ as it does not exist"; \
		cp $*.empty.vcf $@; \
	fi;


# Documentation
$(file >>$(RELEASE_INFOS),## ANNOTATOR snpeff '$(MK_RELEASE)': snpEff annotates variants, and generates *.snpeff.vcf file and *.snpeff.vcf.stats folder, if any. Parameters SNPEFF_FLAGS='$(SNPEFF_FLAGS)', SNPEFF_STATS='$(SNPEFF_STATS)')
$(file >>$(PIPELINES_INFOS),ANNOTATOR:snpeff:snpEff annotation:SNPEFF_FLAGS='$(SNPEFF_FLAGS)', SNPEFF_STATS='$(SNPEFF_STATS)')
