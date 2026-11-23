############################
# NO ANNOTATION Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE=1.0.0
MK_DATE=12/11/2026

# Release note
# 12/11/2026-1.0.0: Initial release with NO_ANNOTATION rule


#################
# NO_ANNOTATION #
#################

%.no_annotation.vcf: %.vcf
	# Annotation - Same VCF
	cp $< $@
	

# Documentation
$(file >>$(RELEASE_INFOS),## ANNOTATOR no_annotation '$(MK_RELEASE)': NO_ANNOTATION generates a VCF file without performing any annotation.)
$(file >>$(PIPELINES_INFOS),ANNOTATOR:no_annotation:NO_ANNOTATION - Generates VCF without annotation)