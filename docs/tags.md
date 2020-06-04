STARK
============
STARK is a Next-Generation Sequencing data analysis pipeline for clinical diagnosis
* Stellar Tools for variants Analysis and RanKing
* Author: Antony Le Béchec
* Copyright: HUS/CPS
* License: GNU GPLA V3
* Release : 0.9.18.1
* Date : 20200602




Tags
-----------


---
**Description**

Tags are used to add meta information to Runs/Analyses or Samples.

Usually to better follow samples important information in order to help biological interpretation (e.g "#CANCER", "#MALE!#BLOOD", "AGE#42!PATHOLOGY#lymphoma").

Also used to identify application to trigger. The tag "APPLICATION#<your_application>" is used to associate run/analysis and samples to a specific STARK application. A tag like "CNV#CANOE" may be used to identify a specific service triggered after STARK analysis (currently not available within STARK).

Tags for analysis will be applied to all samples within the analysis.

---
**Format**

Tags format is basically "TYPE#TAG".

```
[TYPE]#TAG[#TAG]{!}[[TYPE]#TAG[#TAG]]
```

TYPE can be null (e.g. "#TAG"). TAG null will not be considered (e.g. "TYPE#"). TAG are separated by "!" (e.g. "TYPE1#TAG1!TYPE2#TAG"). Words without "#" will not be considered as TAG (e.g. "information not a tag"). TAG can be cumulative for a TYPE (e.g. "TYPE#TAG1#TAG2", same as "TYPE#TAG1!TYPE#TAG2").


Example:

```
TYPE1#TAG1!TAG_TO_FIND#TAG_FOUND!TAG_TO_FIND#TAG_FOUND2!TAG_TO_FIND2#TAG_FOUND3#TAG_FOUND4!#TAG2!TAG_NOT_CONSIDERED# Other word not considered as TAG
```

Correspond to (TYPE: TAGS...):

- TYPE1: TAG1

- TAG_TO_FIND: TAG_FOUND TAG_FOUND2

- TAG_TO_FIND2: TAG_FOUND3 TAG_FOUND4

- null: TAG2


---
**How to use**

Tags can be used within STARK command line (see STARK --help):

```
$ STARK --help
...
# --sample_tag=<STRING1,STRING2...>        List of corresponding SAMPLE Tags.
...
# --analysis_tag=<STRING>                  List of ANALYSIS Tags.
...
```

Tags can be setup in Run SampleSheet, for all the run or for each sample:

```
[Header]
...
Description,APPLICATION#EXOME!#ONCOLOGY MTP-Nano500V2
...

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Manifest,GenomeFolder,Sample_Project,Description
Sample1,,,B01,p7-02,TGATAGTG,p5-01,GTCAATAC,A,hg19,,SEX#F!PATHOLOGY#lymphoma
Sample2,,,B02,p7-08,TTGGCCAC,p5-02,ACTTGAGT,A,hg19,,SEX#M!#LUNG#CANCER#RELAPSE
Control,,,H01,p7-05,GTGGCTAC,p5-08,TGTTGGGT,A,hg19,,SAMPLE_TYPE#internal_control
...
```

---
**Where to find**

Tags are available in sample result folder, within HTML reports, and in "<sample>.tag" (for sample) and "<sample.analysis.tag" (for analysis).


---
**External service**

By searching tags within ".tag" files, triggers can launch external services to complement main STARK analysis.
