library(ChIPpeakAnno)
library(reshape2)
library(Matrix)
library(Seurat)
library(gdata)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(Signac)
library(GenomicRanges)
library(rtracklayer)
library(SeuratWrappers)
library(patchwork)
library(gdata)
library(tidyverse)

load("RNA_Liver_All_Integrated_Object.RData")
liver_RNA <- liver.combined.sct
load("ATAC_BL6_PWD_All_Integrated_Object_ATAC.RData")
liver_ATAC <- Liver_integrated

#how many cells are in both assays? 
summary(colnames(liver_RNA) %in% colnames(liver_ATAC))

#filter for common cells 
rna_cells <- colnames(liver_RNA)
atac_cells <- colnames(liver_ATAC)
common <- intersect(rna_cells,atac_cells)
liver_RNA[["CellName"]] <- colnames(liver_RNA)
liver_ATAC[["CellName"]] <- colnames(liver_ATAC)
liver_RNA_WNN <- subset(liver_RNA, subset = CellName %in% common)
liver_ATAC_WNN <- subset(liver_ATAC, subset = CellName %in% common)

#Integrate
Liver_WNN <- liver_RNA_WNN
Liver_WNN[["integrated_lsi"]] <- liver_ATAC_WNN[["integrated_lsi"]]
Liver_WNN[["ATAC"]] <- liver_ATAC_WNN[["ATAC"]]

Liver_WNN <- FindMultiModalNeighbors(Liver_WNN, reduction.list = list("pca", "integrated_lsi"), dims.list = list(1:30, 2:50))
Liver_WNN <- RunUMAP(Liver_WNN, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Liver_WNN <- FindClusters(Liver_WNN, graph.name = "wsnn", algorithm = 3, verbose = FALSE,resolution = 1)

Idents(Liver_WNN) <- "seurat_clusters"

DimPlot(Liver_WNN,label=FALSE, shuffle=TRUE) 
