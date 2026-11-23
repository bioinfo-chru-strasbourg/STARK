#!/bin/bash
## STARK application RNASEQ

# DEFAULT ENV
######################
source_app $CONFIG_DEFAULT_APP


# APPLICATION INFOS
#####################
APP_NAME="RNASEQ"
APP_RELEASE="1.0"
APP_DESCRIPTION="Application to detect mutations in RNA-Seq data"
APP_GROUP="UNKNOWN"
APP_PROJECT="UNKNOWN"


# PIPELINES
#############

# PIPELINES
# Use STAR to align RNA-Seq data and then different callers to detect mutations (Fusions and SNVs/Indels)
PIPELINES="star.Arriba.howard star.STARFusion.howard star.deepvariant.howard star.gatkHC.howard"

# POST ALIGNMENT STEPS 
# No need to realign because of STAR alignemnt
# Need splitncigar for SNV calling
POST_ALIGNMENT_STEPS="sorting splitncigar markduplicates recalibration compress"

# Option to keep raw BAM file from STAR (with splitNcigar for SNV calling). This option is useful to keep BAM without post-alignment steps (e.g. realignment, recalibration...).
# Keep default value (0)
STAR_KEEP_RAW_BAM=0


# DATABASES
#############

# CTAT_DATABASES_GENE_SOURCE can be "refgene" or "gencode". If "refgene" the gene annotation used for CTAT lib will be based on refGene (NCBI RefSeq) annotation. If "gencode" the gene annotation used for CTAT lib will be based on gencode annotation. The choice of gene annotation can impact the fusion detection results, especially for fusions involving genes that are not well annotated in refGene but are annotated in gencode. By default, CTAT_DATABASES_GENE_SOURCE is set to "refgene" because refGene annotation is more conservative and less likely to include dubious gene models that could lead to false positive fusion calls. However, if you want to use gencode annotation for CTAT lib, you can set CTAT_DATABASES_GENE_SOURCE to "gencode". Note that if you change the gene source for CTAT lib, you should also change the gene source for STAR fusion annotation (STAR_FUSION_GENE_SOURCE) to the same value to ensure consistency between fusion detection and annotation.
# Either "refgene" or "gencode"
# Keep default value ("refgene")
CTAT_DATABASES_GENE_SOURCE="refgene"
