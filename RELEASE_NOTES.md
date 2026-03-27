# STARK Release notes

## Release 19.0.0

### New

- RNASeq
  - Alignement with STAR
  - Caller STAR Fusion and Arriba
- Docker in Docker (DinD)
  - Rules can use docker image
- Alignment
  - STAR for RNASeq
  - BWA2 for DNASeq
- Annotation
  - New HOWARD release
- New services:
  - STARK XTerm
  - STARK Docs

### Improvements

- New releases: Java:17(LTS), Picard:3.0.0, Java:17(LTS), Samtools:1.17, Bedtools:2.31.0, IGVTools:2.16.1
- Python:3 as default release (available: $PYTHON)
- Java:17 as default release (available: $JAVA8 and $JAVA17)
- GATK:4 as default release (available: $GATK3 and GATK4)
- Manage high consuming rules by sequencializing (e.g. aligners such as bwamem2 and star)

## Fixes

- Removed: Java:7, Java:11, Python:2
- Many fixes
