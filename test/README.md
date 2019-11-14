picard_bcl_to_fastq
=======================

### What is this?
This is a wrapper script for [Picard's IlluminaBasecallsToFastq](https://broadinstitute.github.io/picard/command-line-overview.html#IlluminaBasecallsToFastq). 

### What does it do?
It sniffs the run parameters (barcodes, read lengths, etc) from various files in an Illumina GAIIx, MiSeq, Nextseq or HiSeq run folder and automates the Picard tools to produce fastq files. 

### Why not use the default fastq generation on the instrument?
- In case fastq generation on the instrument fails
- Most instruments don't do fastq generation
- You really can't tolerate reads from one barcode showing up in another (this is a surprisingly frequent occurence)

### Why not just run the picard tools by hand?
- You don't want to manually create the files needed for those tools (and possibly introduce mistakes)

### To run it:
- [Configure](#configuration) the path to Picard and your mismatch tolerances at the top of the file
- Feed this script the path to a valid sample sheet file (e.g. ./bcl2fastq.sh SampleSheet.csv) 
- Wait a bit (~10 minutes for a typical miseq run)
- Enjoy your freshly demultiplexed fastq files, conveniently placed in a new folder called "fastq"

### Configuration:
You will need to modify the path to picard variable to fit your environment.
The defaults are for very strict barcode splitting (you may want to increase MAX_MISMATCHES and MAX_NO_CALLS).
This script is configured to use 200G of RAM. If your machine is smaller, set the Xmx value appropriately.

This script works by attempting to sniff run format, barcodes and platform from the files in the run directory.
It works for all the GAII-x, HiSeq, Nextseq, and MiSeq runs I've encountered, but YMMV.

### Dependencies:
- [Picard tools](https://broadinstitute.github.io/picard/) and Java runtime
- perl
- xmllint
- bash shell

Patches and bug reports are welcome.
