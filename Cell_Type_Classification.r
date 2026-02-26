library(Seurat)
library(Signac)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(SeuratWrappers)
library(Matrix)
library(ggplot2)
library(patchwork)
library(gdata)
library(tidyverse)
library(future)
library(gprofiler2)
library(grid.extra)

load("Liver_WNN.RData") 

DefaultAssay(Liver_WNN) <- "RNA"

#CellCycle scoring
#using a pre-loaded set of genes for humans but converting it to mouse gene names
mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
#Scoring the cell cycle 
Liver_WNN <- CellCycleScoring(Liver_WNN, s.features = mmus_s, g2m.features = mmus_g2m, set.ident = TRUE)
head(Liver_WNN[[]])

#Cluster 13 and 15 consist of cells exclusively in G2M
Liver_WNN <- subset(Liver_WNN, invert=TRUE, subset= seurat_clusters == 13)
Liver_WNN <- subset(Liver_WNN, invert=TRUE, subset= seurat_clusters == 15)

#Hepatocyte marker
FeaturePlot(Liver_WNN, features = c("Alb","Cyp3a25","Cyp2e1","Ahsg","Serpina3m"),ncol=2)
VlnPlot(Liver_WNN, features = c("Alb","Cyp3a25","Cyp2e1","Ahsg","Serpina3m"),ncol=2)

#Kupffer Cell Markers
FeaturePlot(Liver_WNN, features = c("Clec4f","Vsig4","Adgre1"),ncol=2,order=T)
VlnPlot(Liver_WNN, features = c("Clec4f","Vsig4","Adgre1"),ncol=2)

#LSEC markers
FeaturePlot(Liver_WNN, features = c("Kdr","Pecam1","Stab2","Cleg4g","Meis2","Bmp6","Gata4"),ncol=2,order=T)
VlnPlot(Liver_WNN, features = c("Kdr","Pecam1","Stab2","Clec4g","Meis2","Bmp6"),ncol=2) 

#HSC Markers
FeaturePlot(Liver_WNN, features = c("Tgfbi","Dcn","Tmem56","Colec10","Reln"),ncol=2,order=T)
VlnPlot(Liver_WNN, features = c("Tgfbi","Dcn","Tmem56","Colec10","Reln"),ncol=2)

#BCell markers
FeaturePlot(Liver_WNN, features = c("Ptprc","Ebf1","Pax5"),ncol=2,order=T)
VlnPlot(Liver_WNN, features = c("Ptprc","Ebf1","Itga4","Grap2","Ikzf3"),ncol=2)

#Cholangiocyte markers
FeaturePlot(Liver_WNN, features = c("Spp1","Onecut2","Epcam","Tspan8"),ncol=2,order=T)
VlnPlot(Liver_WNN, features = c("Spp1","Onecut2","Epcam"),ncol=2)

#Monocyte / Macrophage Marker Cluster 9
FeaturePlot(Liver_WNN, features = c("Ly6e","Cd274","Oasl1"),ncol=2,order=T)
VlnPlot(Liver_WNN, features = c("Ly6e","Cd274","Oasl1"),ncol=2)

new.cluster.ids <- c("Hepatocytes","Hepatocytes","Hepatocytes","Hepatocytes", "Hepatocytes","LSECs","HSCs","Kupffer_Cells","LSECs","B_Cells","LSECs","Neurons","Monocytes","Cholangiocytes")
names(new.cluster.ids) <- levels(Liver_WNN)
Liver_WNN_broad <- RenameIdents(Liver_WNN, new.cluster.ids)

#Ambient RNA detection and decontamination using decontX
sce <- as.SingleCellExperiment(Liver_WNN_broad)
sce <- decontX(sce, z=Liver_WNN_broad$CellType)

#add back de-contaminated counts
Liver_WNN_broad[["decontXcounts"]] <- CreateAssayObject(counts = decontXcounts(sce))

#normalize them for plotting
Liver_WNN_broad <- NormalizeData(object =Liver_WNN_broad,assay = 'decontXcounts',normalization.method = 'CLR',scale.factor = median(Liver_WNN_broad$nCount_decontXcounts))
Liver_WNN_broad <- ScaleData(Liver_WNN_broad)
DefaultAssay(Liver_WNN_broad) <- 'decontXcounts'

#plot
my_levels <- c("Cholangiocytes","Neurons","Monocytes","B_Cells","Kupffer_Cells","HSCs","LSECs","Hepatocytes")
Idents(Liver_WNN_broad) <- factor(x = Idents(Liver_WNN_broad), levels = my_levels)
VlnPlot(Liver_WNN_broad, features=c("Alb","Cyp3a25","Clec4g","Stab2","Reln","Dcn","Clec4f","Vsig4","Ptprc","Itga4","Ly6e","Oasl1","Scd1","Kcnip1","Spp1"), pt.size=0, stack=TRUE, assay="RNA",fill.by="ident", cols=c("#2A91A2","#9B59B6","#3498DB","#F39C12","#34495E","#F1C40F","#2ECC71","#E74C3C"), same.y.lims=TRUE) + theme(legend.position="none")

DefaultAssay(Liver_WNN_broad) <- "ATAC"
p0 <- CoveragePlot(Liver_WNN_broad, assay="ATAC", region=c("chr12-103863000-103864200"),  peaks=FALSE, tile=FALSE,height=c(5,1)) + scale_color_manual(values=c("#E74C3C","#2ECC71","#F1C40F","#34495E","#F39C12","#3498DB","#9B59B6","#2A91A2"))#Serpina1a
p1 <- CoveragePlot(Liver_WNN_broad,assay="ATAC",region=c("chr8-3717500-3719000"), peaks=FALSE, tile=FALSE, height=c(5,1)) #Clec4g, LSEC
p2 <- CoveragePlot(Liver_WNN_broad, assay="ATAC", region=c("chr10-86989800-86994000"), peaks=FALSE,  tile=FALSE, height=c(5,1)) #stab2, LSEC
p3 <- CoveragePlot(Liver_WNN_broad, assay="ATAC", region=c("chr10-97479000-97484000"), peaks=FALSE, tile=FALSE,height=c(5,1)) #Dcn, HSCs
p4 <- CoveragePlot(Liver_WNN_broad, assay="ATAC", region=c("chr6-83655600-83656600"), peaks=FALSE, tile=FALSE,height=c(5,1)) #Clec4f, Kupffer
p5 <- CoveragePlot(Liver_WNN_broad, assay="ATAC", region=c("chr15-54409800-54411000"), peaks=FALSE, tile=FALSE,height=c(5,1)) #Colec10, HSCs,
p6 <- CoveragePlot(Liver_WNN_broad, assay="ATAC", region=c("chr1-138173500-138176000"), peaks=FALSE,  tile=FALSE,height=c(5,1)) #Ptprc, BCells & Kupffer
p7 <- CoveragePlot(Liver_WNN_broad, assay="ATAC", region=c("chr2-79254400-79256500"), peaks=FALSE,  tile=FALSE,height=c(5,1)) #Itga4, BCells
p8 <- CoveragePlot(Liver_WNN_broad, assay="ATAC", region=c("chr5-114923000-114923500"), peaks=FALSE,  tile=FALSE,height=c(5,1)) #Oasl1, Monocytes
p9 <- CoveragePlot(Liver_WNN_broad, assay="ATAC", region=c("chr5-104432000-104436000"), window=1000, peaks=FALSE,  tile=FALSE,height=c(5,1)) #Spp1, Chol
p10 <- CoveragePlot(Liver_WNN_broad, assay="ATAC", region=c("chr19-44393500-44395000"), window=200,peaks=FALSE,  tile=FALSE,height=c(5,1)) #Scd1, Neurons

cowplot::plot_grid(p0,p1,p2,p3,p5,p4,p6,p7,p8,p10,p9,ncol = 11)


DefaultAssay(Liver_WNN_broad) <- "RNA"
p1 <- FeaturePlot(Liver_WNN_broad, features=c("Alb"), order=TRUE, ncol=3, min.cutoff='q60', label=FALSE, repel=TRUE) + theme(axis.text.x=element_blank(),axis.text.y=element_blank()) 
p2 <- FeaturePlot(Liver_WNN_broad, features=c("Cyp3a25"), order=TRUE, ncol=3, min.cutoff='q60', label=FALSE, repel=TRUE) + theme(axis.text.x=element_blank(),axis.text.y=element_blank())
p3 <- FeaturePlot(Liver_WNN_broad, features=c("Clec4g"), order=TRUE, ncol=3, min.cutoff='q10', label=FALSE, repel=TRUE) + theme(axis.text.x=element_blank(),axis.text.y=element_blank())
p12 <- FeaturePlot(Liver_WNN_broad, features=c("Stab2"), order=TRUE, ncol=3, min.cutoff='q20', label=FALSE, repel=TRUE) + theme(axis.text.x=element_blank(),axis.text.y=element_blank())
p4 <- FeaturePlot(Liver_WNN_broad, features=c("Reln"), order=TRUE, ncol=3, min.cutoff='q10', label=FALSE, repel=TRUE) + theme(axis.text.x=element_blank(),axis.text.y=element_blank())
p5 <- FeaturePlot(Liver_WNN_broad, features=c("Dcn"), order=TRUE, ncol=3, min.cutoff='q10', label=FALSE, repel=TRUE) + theme(axis.text.x=element_blank(),axis.text.y=element_blank())
p6 <- FeaturePlot(Liver_WNN_broad, features=c("Clec4f"), order=TRUE, ncol=3, label=FALSE, repel=TRUE) + theme(axis.text.x=element_blank(),axis.text.y=element_blank())
p7 <- FeaturePlot(Liver_WNN_broad, features=c("Vsig4"), order=TRUE, ncol=3, label=FALSE, repel=TRUE) + theme(axis.text.x=element_blank(),axis.text.y=element_blank())
p8 <- FeaturePlot(Liver_WNN_broad, features=c("Ptprc"), order=TRUE, ncol=3, label=FALSE, repel=TRUE, min.cutoff='q15') + theme(axis.text.x=element_blank(),axis.text.y=element_blank())
p9 <- FeaturePlot(Liver_WNN_broad, features=c("Oasl1"), order=TRUE, ncol=3, label=FALSE, repel=TRUE, min.cutoff='q10') + 
theme(axis.text.x=element_blank(),axis.text.y=element_blank())
p10 <- FeaturePlot(Liver_WNN_broad, features=c("Kcnip1"), order=TRUE, ncol=3, label=FALSE, repel=TRUE, min.cutoff='q50') + theme(axis.text.x=element_blank(),axis.text.y=element_blank())
p11 <- FeaturePlot(Liver_WNN_broad, features=c("Spp1"), order=TRUE, ncol=3, label=FALSE, repel=TRUE) + theme(axis.text.x=element_blank(),axis.text.y=element_blank())

grid.arrange(p1,p2,p3,p12,p4,p5,p6,p7,p8,p9,p10,p11, ncol=4)

#Fraction of cell types
Liver_WNN_broad$Cell_Types <- Idents(Liver_WNN_broad)
library(data.table)
library(magrittr)
md <- Liver_WNN_broad@meta.data %>% as.data.table
md$Cell_Types <- Idents(Liver_WNN_bbroad)
Cell_dist <- md[, .N, by = c("Cell_Types")]
sum(Cell_dist$N)
Cell_dist$perc <- Cell_dist$N/19664

ggplot(Cell_dist, aes(fill=Cell_Types,y=perc,x="Liver")) + geom_bar(position="stack", stat="identity") + theme_classic() + ylab("Total Cells: 19.964") + xlab("") + scale_fill_manual(name = "Cell Types", values=c("#E74C3C","#2ECC71","#F1C40F","#34495E","#F39C12","#3498DB","#9B59B6","#2A91A2")) + coord_polar("y", start=0) + theme(axis.text.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.line=element_blank())


