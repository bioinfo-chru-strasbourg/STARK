############################
# NO CALLING Rules
# Author: Antony Le Bechec
############################
# Release
MK_RELEASE=1.0.0
MK_DATE=12/11/2026

# Release note
# 12/11/2026-1.0.0: Initial release with NO_CALLING rule


##############
# NO_CALLING #
##############

%.no_calling.vcf: %.bam %.bam.bai %.empty.vcf
	# Calling - Empty VCF
	cp $*.empty.vcf $@


# Documentation
$(file >>$(RELEASE_INFOS),## CALLER no_calling '$(MK_RELEASE)': NO_CALLING generates a VCF file without performing any variant calling. Parameters BAM_COMPRESSION='$(BAM_COMPRESSION)')
$(file >>$(PIPELINES_INFOS),CALLER:no_calling:NO_CALLING - Generates VCF without calling variants:BAM_COMPRESSION='$(BAM_COMPRESSION)')
