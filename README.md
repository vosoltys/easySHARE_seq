# easySHARE_seq
This repository contains scripts used in https://doi.org/10.1101/2024.02.26.581705

All data including processed single-cell matrices and accompanying files can be accessed under GEO accession GSE256434.

Communications should be sent to volker.soltys@gmx.de

Analysis steps in order:
- Full_demultiplex_pipeline in combination with SHARE_demult_fastq.c
- scRNA-seq_Seurat_Analysis.r
- scATACseq_Signac_Analysis.r
- scRNA_ATAC_integration.r
- Cell_Type_Classification.r
- Peak_Gene_associations.r
- Zonation_LSECs.r



Software Versions:

bwa v.0.7.17 
samtools v.1.16.1 
featureCounts v.2.0.1 
UMItools v.1.1.2 
STAR v.2.7.9a 
cutadapt v.4.0   
PicardTools v.2.18.25 
macs2 v.2.2.7.1 

R (v4.3.1) package Versions:

Seurat v. 5.0.1 
Signac v.1.12.9 
chromVAR v.1.12 
GenomicRanges v.1.52.1 
BSgenome.Mmusculus.UCSC.mm10 v.1.4.0 
ChIPpeakAnno v.3.34.1 
monocle v.2.22 
ggplot2 v3.4.4 
reshape2 v1.4.4 
dplyr v1.1.4 
clusterProfiler v3.18.1 
