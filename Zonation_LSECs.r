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
library(gridExtra)
library(tidyverse)
library(future)
library(gprofiler2)
library(gridExtra)
library(SummarizedExperiment)
library(chromVAR)
library(matrixStats)
library(ChIPpeakAnno)
library(BSDA)

load("Liver_WNN_broad_annotated.RObject") 
Liver_subset <- subset(Liver_WNN_broad, subset= CellType =="LSECs")


#Pulling out & binning RNA data
rna_coords <- as.data.frame(Liver_subset@reductions$wnn.umap@cell.embeddings)
data <- Liver_subset@assays$RNA@data

wnt2 <- as.numeric(unlist(data["Wnt2",]))
wnt9 <- as.numeric(unlist(data["Wnt9b",]))
efnb2 <- as.numeric(unlist(data["Efnb2",]))
ltbp4 <- as.numeric(unlist(data["Ltbp4",]))
rspo3 <- as.numeric(unlist(data["Rspo3",]))
lyve1 <- as.numeric(unlist(data["Lyve1",]))
kit <- as.numeric(unlist(data["Kit",]))
meis1 <- as.numeric(unlist(data["Meis1",]))
bmp2 <- as.numeric(unlist(data["Bmp2",]))
pear1 <- as.numeric(unlist(data["Pear1",]))
lama4 <- as.numeric(unlist(data["Lama4",]))
thbd <- as.numeric(unlist(data["Thbd",]))
dkk3 <- as.numeric(unlist(data["Dkk3",]))

rna_coords$Efnb2 <- efnb2
rna_coords$Ltbp4 <- ltbp4
rna_coords$Rspo3 <- rspo3
rna_coords$Lyve1 <- lyve1
rna_coords$wnt2 <- wnt2
rna_coords$Kit <- kit
rna_coords$Bmp2 <- bmp2
rna_coords$Meis1 <- meis1
rna_coords$Pear1 <- pear1
rna_coords$wnt9 <- wnt9
rna_coords$lama4 <- lama4
rna_coords$thbd <- thbd
rna_coords$dkk3 <- dkk3

rna_coords <- rna_coords[order(rna_coords$wnnUMAP_2),]

rna_efnb2 <- rna_coords[,1:3] %>% mutate(Bins = cut(wnnUMAP_2, breaks=10, label=FALSE)) %>% group_by(Bins) %>% dplyr::filter(Efnb2 > 0) %>% mutate(m_efnb2 = mean(Efnb2), sd_efnb2 = sd(Efnb2))
rna_ltbp4 <- rna_coords[,c(1,2,4)] %>% mutate(Bins = cut(wnnUMAP_2, breaks=10, label=FALSE)) %>% group_by(Bins) %>% dplyr::filter(Ltbp4 > 0) %>% mutate(m_Ltbp4 = mean(Ltbp4), sd_Ltbp4 = sd(Ltbp4))
rna_rspo3 <- rna_coords[,c(1,2,5)] %>% mutate(Bins = cut(wnnUMAP_2, breaks=10, label=FALSE)) %>% group_by(Bins) %>% dplyr::filter(Rspo3 > 0) %>% mutate(m_rspo3 = mean(Rspo3), sd_rspo3 = sd(Rspo3))
rna_lyve1 <- rna_coords[,c(1,2,6)] %>% mutate(Bins = cut(wnnUMAP_2, breaks=10, label=FALSE)) %>% group_by(Bins) %>% dplyr::filter(Lyve1 > 0) %>% mutate(m_Lyve1 = mean(Lyve1), sd_Lyve1 = sd(Lyve1))
rna_wnt2 <- rna_coords[,c(1,2,7)] %>% mutate(Bins = cut(wnnUMAP_2, breaks=10, label=FALSE)) %>% group_by(Bins) %>% dplyr::filter(wnt2 > 0) %>% mutate(m_wnt2 = mean(wnt2), sd_wnt2 = sd(wnt2))
rna_kit <- rna_coords[,c(1,2,8)] %>% mutate(Bins = cut(wnnUMAP_2, breaks=10, label=FALSE)) %>% group_by(Bins) %>% dplyr::filter(Kit > 0) %>% mutate(m_Kit = mean(Kit), sd_Kit = sd(Kit))
rna_bmp2 <- rna_coords[,c(1,2,9)] %>% mutate(Bins = cut(wnnUMAP_2, breaks=10, label=FALSE)) %>% group_by(Bins) %>% dplyr::filter(Bmp2 > 0) %>% mutate(m_Bmp2 = mean(Bmp2), sd_Bmp2 = sd(Bmp2))
rna_meis1 <- rna_coords[,c(1,2,10)] %>% mutate(Bins = cut(wnnUMAP_2, breaks=10, label=FALSE)) %>% group_by(Bins) %>% dplyr::filter(Meis1 > 0) %>% mutate(m_Meis1 = mean(Meis1), sd_Meis1 = sd(Meis1))
rna_pear1 <- rna_coords[,c(1,2,11)] %>% mutate(Bins = cut(wnnUMAP_2, breaks=10, label=FALSE)) %>% group_by(Bins) %>% dplyr::filter(Pear1 > 0) %>% mutate(m_Pear1 = mean(Pear1), sd_Pear1 = sd(Pear1))
rna_wnt9 <- rna_coords[,c(1,2,12)] %>% mutate(Bins = cut(wnnUMAP_2, breaks=10, label=FALSE)) %>% group_by(Bins) %>% dplyr::filter(wnt9 > 0) %>% mutate(m_wnt9 = mean(wnt9), sd_wnt9 = sd(wnt9))
rna_lama4 <- rna_coords[,c(1,2,13)] %>% mutate(Bins = cut(wnnUMAP_2, breaks=10, label=FALSE)) %>% group_by(Bins) %>% dplyr::filter(lama4 > 0) %>% mutate(m_lama4 = mean(lama4), sd_lama4 = sd(lama4))
rna_thbd <- rna_coords[,c(1,2,14)] %>% mutate(Bins = cut(wnnUMAP_2, breaks=10, label=FALSE)) %>% group_by(Bins) %>% dplyr::filter(thbd > 0) %>% mutate(m_thbd = mean(thbd), sd_thbd = sd(thbd))
rna_dkk3 <- rna_coords[,1:3] %>% mutate(Bins = cut(wnnUMAP_2, breaks=10, label=FALSE)) %>% group_by(Bins) %>% dplyr::filter(dkk3 > 0) %>% mutate(m_dkk3 = mean(dkk3), sd_dkk3 = sd(dkk3))

rna_binned_mean <- rna_coords %>% mutate(Bins = cut(wnnUMAP_2, breaks=10, label=FALSE)) %>% group_by(Bins) %>% mutate(m_efnb2 = mean(Efnb2), sd_efnb2 = sd(Efnb2), m_Ltbp4 = mean(Ltbp4), sd_Ltbp4 = sd(Ltbp4), m_rspo3 = mean(Rspo3), sd_rspo3 = sd(Rspo3), m_Lyve1 = mean(Lyve1), sd_Lyve1 = sd(Lyve1), m_wnt2 = mean(wnt2), sd_wnt2 = sd(wnt2), m_Kit = mean(Kit), sd_Kit = sd(Kit), m_Bmp2 = mean(Bmp2), sd_Bmp2 = sd(Bmp2), m_Meis1 = mean(Meis1), sd_Meis1 = sd(Meis1), m_Pear1 = mean(Pear1), sd_Pear1 = sd(Pear1), m_wnt9 = mean(wnt9), sd_wnt9 = sd(wnt9), m_lama4 = mean(lama4), sd_lama4 = sd(lama4),m_thbd = mean(thbd), sd_thbd = sd(thbd))


# Pulling out & binning ATAC-seq data
DefaultAssay(Liver_subset) <- "ATAC"

atac_coords <- as.data.frame(Liver_subset@reductions$wnn.umap@cell.embeddings)
data <- Liver_subset@assays$ATAC@data

wnt2_end <- as.numeric(unlist(data["chr6-17991462-17991911",]))
wnt2_prom_1 <- as.numeric(unlist(data["chr6-18030331-18030585",]))
wnt2_prom_2 <- as.numeric(unlist(data["chr6-18031261-18032326",]))
wnt9b_end <- as.numeric(unlist(data["chr11-103729308-103729554",]))
wnt9b_int <- as.numeric(unlist(data["chr11-103734570-103735917",]))
bmp2_up <- as.numeric(unlist(data["chr2-133550269-133550680",]))
bmp2_prom <- as.numeric(unlist(data["chr2-133552024-133552776",]))
efnb2_prom <- as.numeric(unlist(data["chr8-8660841-8662085",]))
lama4_prom <- as.numeric(unlist(data["chr10-38965164-38965990",]))
meis1_int <- as.numeric(unlist(data["chr11-19007769-19008171",]))
meis1_prom <- as.numeric(unlist(data["chr11-19018439-19020880",]))
lyve1_int <- as.numeric(unlist(data["chr7-110862706-110863248",]))
klf4_prom <- as.numeric(unlist(data["chr4-55532450-55533288",]))
jak1_prom <- as.numeric(unlist(data["chr4-101275456-101276779",]))
ltbp4_int <- as.numeric(unlist(data["chr7-27325703-27326369",]))
ltbp4_int2 <- as.numeric(unlist(data["chr7-27330856-27331432",]))
pear1_int1 <- as.numeric(unlist(data["chr3-87764893-87765440",]))
pear1_int2 <- as.numeric(unlist(data["chr3-87770496-87771366",]))
kit_prom <- as.numeric(unlist(data["chr5-75574893-75575195",]))
rspo3_prom <- as.numeric(unlist(data["chr10-29532964-29533569",]))
thbd_prom <- as.numeric(unlist(data["chr2-148409074-148409577",]))
dkk3_int <- as.numeric(unlist(data["chr7-112138555-112139266",]))

atac_coords$wnt2_end <- wnt2_end
atac_coords$wnt2_prom_1 <- wnt2_prom_1
atac_coords$wnt2_prom_2 <- wnt2_prom_2
atac_coords$wnt9b_end <- wnt9b_end
atac_coords$wnt9b_int <- wnt9b_int
atac_coords$bmp2_up <- bmp2_up
atac_coords$bmp2_prom <- bmp2_prom
atac_coords$efnb2_prom <- efnb2_prom
atac_coords$lama4_prom <- lama4_prom
atac_coords$meis1_int <- meis1_int
atac_coords$meis1_prom <- meis1_prom
atac_coords$lyve1_int <- lyve1_int
atac_coords$klf4_prom <- klf4_prom
atac_coords$jak1_prom <- jak1_prom
atac_coords$ltbp4_int <- ltbp4_int
atac_coords$ltbp4_int2 <- ltbp4_int2
atac_coords$pear1_int1 <- pear1_int1
atac_coords$pear1_int2 <- pear1_int2
atac_coords$kit_prom <-kit_prom
atac_coords$rspo3_prom <-rspo3_prom
atac_coords$thbd_prom <- thbd_prom
atac_coords$dkk3_int <- dkk3_int


atac_coords <- atac_coords %>% dplyr::filter(wnnUMAP_1 > 5 )
atac_coords <- atac_coords[order(atac_coords$wnnUMAP_2),]

#This line of code gets rid of all cells that contain 0 data across all peaks we look at, which is ~2/3rds of cells
#Not sure if this is fair - in ATACseq no data can mean closed chromatin
atac_coords <- atac_coords[rowSums(atac_coords[3:18])>0,]

atac_binned_mean <- atac_coords %>% mutate(Bins = cut(wnnUMAP_2, breaks=10, labels=FALSE)) %>% group_by(Bins) %>% mutate(m_wnt2_end = mean(wnt2_end), m_wnt2_prom_1 = mean(wnt2_prom_1),m_wnt2_prom_2 = mean(wnt2_prom_2),m_wnt9b_end = mean(wnt9b_end),m_wnt9b_int = mean(wnt9b_int),m_bmp2_up = mean(bmp2_up),m_bmp2_prom = mean(bmp2_prom),m_efnb2_prom = mean(efnb2_prom),m_lama4_prom = mean(lama4_prom), m_meis1_prom = mean(meis1_prom),m_meis1_int = mean(meis1_int),m_lyve1_int = mean(lyve1_int),m_klf4_prom = mean(klf4_prom),m_jak1_prom = mean(jak1_prom), m_ltbp4_int = mean(ltbp4_int),m_ltbp4_int2 = mean(ltbp4_int2),m_pear1_int1 = mean(pear1_int1),m_pear1_int2 = mean(pear1_int2),m_kit_prom = mean(kit_prom), m_rspo3_prom = mean(rspo3_prom),m_rspo3_prom = mean(rspo3_prom),m_thbd_prom = mean(thbd_prom),m_dkk3_int = mean(dkk3_int))   


#example plot
ggplot(rna_binned_mean,aes(x=as.integer(Bins), y=m_wnt2)) + geom_point() + theme_classic() + ggtitle("Wnt2") + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ylab("Mean Expression") + geom_line() + geom_smooth(aes(x=as.integer(Bins),y=m_wnt2,color=m_wnt2,fill=m_wnt2),method=loess) + theme(legend.position="none")
ggplot(atac_binned_mean,aes(x=Bins, y=m_wnt2_end)) + geom_point() + theme_classic() + ggtitle("Wnt2_int") + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ylab("Mean Accessibility") + geom_line() + geom_smooth(aes(x=as.integer(Bins),y=m_wnt2_end,color="red",fill=m_wnt2_end),method=loess)+ theme(legend.position="none")




## Pseudotime ordering & plotting
cds <- as.cell_data_set(Liver_subset, assay=NULL)
cds <- cluster_cells(cds, reduction_method=WNN.UMAP)
cds <- learn_graph(cds,use_partition=FALSE,learn_graph_control= list(nn.k=16))
#Set the root cells
max.wnt2 <- which.max(unlist(FetchData(Liver_subset, "Wnt2")))
max.wnt2 <- colnames(Liver_subset)[max.wnt2]
cds <- order_cells(cds, root_cells = max.wnt2)
mcols(cds@rowRanges)$gene_short_name <- rownames(cds)
traj.coord <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]

#add pseudotime back to Seurat
Liver_subset$pseudotime <- traj.coord
FeaturePlot(Liver_subset, "pseudotime", reduction='wnn.umap')

rna_coords <- as.data.frame(Liver_subset@reductions$wnn.umap@cell.embeddings)
rna_coords$pseudotime <- Liver_subset$pseudotime
data <- Liver_subset@assays$RNA@data
wnt2 <- as.numeric(unlist(data["Wnt2",]))
rna_coords$wnt2 <- wnt2

ggplot(rna_coords) + geom_point(aes(x=wnnUMAP_1, y=wnnUMAP_2, color=pseudotime), size=0.8) + scale_colour_viridis_c(option="plasma") + theme_classic()
ggplot(rna_coords,aes(x=pseudotime, y=wnt2)) +  geom_smooth(method=loess) + theme_classic() + xlab("Pseudotime") + ggtitle("Wnt2")+ ylab("")

