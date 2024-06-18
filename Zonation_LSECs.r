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
library(monocle3)

#object with RNA and ATAC-seq already integrated.
load("Liver_WNN_broad_annotated.RObject") 
Liver_subset <- subset(Liver_WNN_broad, subset= CellType =="LSECs")

#add coordinates as variables and remove a few outlier cells [61]
umapCoord <- as.data.frame(Embeddings(object = Liver_subset[["wnn.umap"]]))
Liver_subset$wnn_1 <- umapCoord$wnnUMAP_1
Liver_subset$wnn_2 <- umapCoord$wnnUMAP_2
Liver_subset$umap_2 <- as.data.frame(Embeddings(object = Liver_subset[["umap"]]))[,2]
Liver_subset <- subset(x = Liver_subset, subset = wnn_1 > 6)
UpdateSeuratObject(Liver_subset)


#Start Pseudotime analysis
cds <- as.cell_data_set(Liver_subset, assay=NULL)
cds <- cluster_cells(cds, reduction_method='UMAP')
cds <- learn_graph(cds,use_partition=FALSE,learn_graph_control= list(nn.k=16))
#set the root cells
max.wnt2 <- which.max(unlist(FetchData(Liver_subset, "Wnt2")))
max.wnt2 <- colnames(Liver_subset)[max.wnt2]
cds <- order_cells(cds, root_cells = max.wnt2)
mcols(cds@rowRanges)$gene_short_name <- rownames(cds)
traj.coord <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]


#add pseudotime back to Seurat
# 2 cells are outliers with pseudotime ~10, all others are below 5. Filter them out
Liver_subset$pseudotime <- traj.coord
Liver_subset <- subset(x = Liver_subset, subset =umap_2 < 3)
FeaturePlot(Liver_subset, "pseudotime", reduction='wnn.umap')
plot_cells(cds, show_trajectory_graph = TRUE)

#gene expression of individual markers along pseudotime
rna_coords <- as.data.frame(Liver_subset@reductions$wnn.umap@cell.embeddings)
rna_coords$pseudotime <- Liver_subset$pseudotime
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
fabp4 <- as.numeric(unlist(data["Fabp4",]))
jak1 <- as.numeric(unlist(data["Jak1",]))
nkd1 <- as.numeric(unlist(data["Nkd1",]))
fmn1 <- as.numeric(unlist(data["Fmn1",]))
flrt1 <- as.numeric(unlist(data["Flrt1",]))
olfm1 <- as.numeric(unlist(data["Olfm1",]))
tmed8 <- as.numeric(unlist(data["Tmed8",]))
itprip <- as.numeric(unlist(data["Ptgs1",]))

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
rna_coords$fabp4 <- fabp4
rna_coords$jak1 <- jak1
rna_coords$nkd1 <- nkd1
rna_coords$fmn1 <- fmn1
rna_coords$flrt1 <- flrt1
rna_coords$olfm1 <- olfm1
rna_coords$tmed8 <- tmed8
rna_coords$itprip <- itprip
rna_coords <- rna_coords[order(rna_coords$pseudotime),]

#accessibility of individual peaks along pseudotime
DefaultAssay(Liver_subset) <- "ATAC"
atac_coords <- as.data.frame(Liver_subset@reductions$wnn.umap@cell.embeddings)
atac_coords$pseudotime <- Liver_subset$pseudotime
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
lyve1_prom <- as.numeric(unlist(data["chr7-110862706-110863248",]))
lyve1_int <- as.numeric(unlist(data["chr7-110860211-110860407",]))
lyve1_int2 <- as.numeric(unlist(data["chr7-110853672-110855790",]))
klf4_prom <- as.numeric(unlist(data["chr4-55532450-55533288",]))
jak1_prom <- as.numeric(unlist(data["chr4-101275456-101276779",]))
ltbp4_int <- as.numeric(unlist(data["chr7-27325703-27326369",]))
ltbp4_int2 <- as.numeric(unlist(data["chr7-27330856-27331432",]))
pear1_int1 <- as.numeric(unlist(data["chr3-87764893-87765440",]))
pear1_int2 <- as.numeric(unlist(data["chr3-87770496-87771366",]))
kit_prom <- as.numeric(unlist(data["chr5-75574893-75575195",]))
kit_prom2 <- as.numeric(unlist(data["chr5-75575541-75575930",]))
rspo3_prom <- as.numeric(unlist(data["chr10-29532964-29533569",]))
thbd_prom <- as.numeric(unlist(data["chr2-148409074-148409577",]))
dkk3_int <- as.numeric(unlist(data["chr7-112138555-112139266",]))
dkk3_int2 <- as.numeric(unlist(data["chr7-112120964-112121300",]))
nkd1_int <- as.numeric(unlist(data["chr8-88561623-88561730",]))
fabp4_int <- as.numeric(unlist(data["chr3-10211267-10212179",]))
jak1_prom <- as.numeric(unlist(data["chr4-101264353-101265597",]))
jak1_int <- as.numeric(unlist(data["chr4-101175960-101176905",]))
olfm1_int <- as.numeric(unlist(data["chr2-28211865-28212033",]))
tmed8_int <- as.numeric(unlist(data["chr12-87199460-87200724",]))

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
atac_coords$lyve1_int2 <- lyve1_int2
atac_coords$lyve1_prom <- lyve1_prom
atac_coords$klf4_prom <- klf4_prom
atac_coords$jak1_prom <- jak1_prom
atac_coords$ltbp4_int <- ltbp4_int
atac_coords$ltbp4_int2 <- ltbp4_int2
atac_coords$pear1_int1 <- pear1_int1
atac_coords$pear1_int2 <- pear1_int2
atac_coords$kit_prom <-kit_prom
atac_coords$kit_prom2 <-kit_prom2
atac_coords$rspo3_prom <-rspo3_prom
atac_coords$thbd_prom <- thbd_prom
atac_coords$dkk3_int <- dkk3_int
atac_coords$dkk3_int2 <- dkk3_int2
atac_coords$nkd1_int <- nkd1_int
atac_coords$fabp4_int <- fabp4_int
atac_coords$jak1_prom <- jak1_prom
atac_coords$jak1_int <- jak1_int
atac_coords$olfm1_int <- olfm1_int
atac_coords$tmed8_int <- tmed8_int


## Identification of novel markers 

#RNA-seq
library(zoo)
rna_coords <- as.data.frame(Liver_subset@reductions$wnn.umap@cell.embeddings)
rna_coords$pseudotime <- Liver_subset$pseudotime
data <- Liver_subset@assays$RNA@data
data <- cbind(rna_coords, t(data))
data <- data[order(data$pseudotime),]
#bin data along pseudotime and calculate mean per bin across every gene
data <- data %>% mutate(Bins = cut(pseudotime, breaks=10, label=FALSE)) %>% group_by(Bins) %>% summarise(across(everything(), list(mean)))

## Pericentral
#filter out genes with low expression (maximum mean < 0.25)
index_keep <- c()
for (i in 1:length(data)) {
    if (max(data[,i]) < 0.25) {
        next
    } else {
        index_keep <- append(index_keep,i)
    }}
int_data <- data[,index_keep]

#require mean(bin9+10) < mean(bin1+2)
index_keep <- c()
for (i in 1:length(int_data)) {
    if (((int_data[1,i] + int_data[2,i])/2) > ((int_data[9,i] + int_data[10,i])/2)){
        index_keep <- c(index_keep,i)
    }}
int_data <- int_data[,index_keep]

#Calculate moving average & require decrease every time with 2 exceptions
mov_avg <- rollmeanr(int_data,3)
index_keep <- c()

for (i in 1:ncol(mov_avg)){
    for (x in 1:8){
        #if first bin, set avg value
        if (x == 1){ 
            avg <- mov_avg[x,i]
            exceptions <-  1
        if (exceptions < 0) {break}
        #if not, check if moving average decreases. If it doesn't, substract from exception    
        } else {
            if (avg > mov_avg[x,i]){
                avg <- mov_avg[x,i]
                next
            } else {
                avg <- mov_avg[x,i]
                exceptions <- exceptions -1
            }
        }
    }
    if (exceptions < 0) { 
        next
    } else {
        index_keep <- c(index_keep, i)
}}
int_data <- int_data[,index_keep]

#Divide each column by its maximum to have values from 1 to zero
for (i in 1:length(int_data)) {
    int_data[,i] <- int_data[,i] / max(int_data[,i])}
int_data$Bins <- c(1:10)
pericentral_data <- int_data 

## Periportal
#filter out genes with low expression (maximum mean < 0.25)
index_keep <- c()
for (i in 1:length(data)) {
    if (max(data[,i]) < 0.25) {
        next
    } else {
        index_keep <- append(index_keep,i)
    }}
int_data <- data[,index_keep]

#require mean(bin9+10) > mean(bin1+2)
index_keep <- c()
for (i in 1:length(int_data)) {
    if (((int_data[1,i] + int_data[2,i])/2) < ((int_data[9,i] + int_data[10,i])/2)){
        index_keep <- c(index_keep,i)
    }}
int_data <- int_data[,index_keep]
#Calculate moving average & require decrease every time with 2 exceptions
library(zoo)
mov_avg <- rollmeanr(int_data,3)
index_keep <- c()

for (i in 1:ncol(mov_avg)){
    for (x in 1:8){
        #if first bin, set avg value
        if (x == 1){ 
            avg <- mov_avg[x,i]
            exceptions <-  1
        if (exceptions < 0) {break}
        #if not, check if moving average increases If it doesn't, substract from exception    
        } else {
            if (avg < mov_avg[x,i]){
                avg <- mov_avg[x,i]
                next
            } else {
                avg <- mov_avg[x,i]
                exceptions <- exceptions -1
            }
        }
    }
    if (exceptions < 0) { 
        next
    } else {
        index_keep <- c(index_keep, i)
}}
int_data <- int_data[,index_keep]

#Divide each column by its maximum to have values from 1 to zero
for (i in 1:ncol(int_data)) {
    int_data[,i] <- int_data[,i] / max(int_data[,i])}
int_data <- int_data[-c(2)]
periportal_data <- int_data 

#ATAC-seq
DefaultAssay(Liver_subset) <- "ATAC"
atac_coords <- as.data.frame(Liver_subset@reductions$wnn.umap@cell.embeddings)
atac_coords$pseudotime <- Liver_subset$pseudotime
datax <- as.data.frame(Liver_subset@assays$ATAC@data)
datax <- cbind(atac_coords, t(datax))
datax <- datax[order(datax$pseudotime),]

#filter out peaks with no data
index_keep <- c()
for (i in 1:ncol(datax)){
    if (i ==2) {
        index_keep <- c(index_keep,i)
    } else if (sum(datax[,i]) > 25) {
        index_keep <- c(index_keep,i)
    }
}
datax <- datax[,index_keep]
datax <- datax[order(datax$pseudotime),]

#bin data and calculate mean per bin across every gene
datax$Bins <- as.numeric(cut_number(datax$pseudotime, 10))
#datax %>% mutate(Bins = cut(pseudotime, breaks=10, label=FALSE))
data <- datax %>% group_by(Bins) %>% summarise(across(everything(), list(mean)))

## Pericentral
#filter out peaks with low accessibility (maximum mean < 0.25)
index_keep <- c()
for (i in 1:length(data)) {
    if (max(data[,i]) < 0.25) {
        next
    } else {
        index_keep <- append(index_keep,i)
    }}
int_data <- data[,index_keep]

#require mean(bin9+10) < mean(bin1+2)
index_keep <- c()
for (i in 1:length(int_data)) {
    if (((int_data[1,i] + int_data[2,i])/2) > ((int_data[9,i] + int_data[10,i])/2)){
        index_keep <- c(index_keep,i)
    }}
int_data <- int_data[,index_keep]

#Calculate moving average & require decrease every time with 2 exceptions
mov_avg <- rollmeanr(int_data,3)
index_keep <- c()

for (i in 1:ncol(mov_avg)){
    for (x in 1:8){
        #if first bin, set avg value
        if (x == 1){ 
            avg <- mov_avg[x,i]
            exceptions <-  1
        if (exceptions < 0) {break}
        #if not, check if moving average decreases. If it doesn't, substract from exception    
        } else {
            if (avg > mov_avg[x,i]){
                avg <- mov_avg[x,i]
                next
            } else {
                avg <- mov_avg[x,i]
                exceptions <- exceptions -1
            }
        }
    }
    if (exceptions < 0) { 
        next
    } else {
        index_keep <- c(index_keep, i)
}}
int_data <- int_data[,index_keep]

#Divide each column by its maximum to have values from 1 to zero
for (i in 1:length(int_data)) {
    int_data[,i] <- int_data[,i] / max(int_data[,i])}
int_data$Bins <- c(1:10)
atac_pericentral_data <- int_data 

## Periportal
#filter out genes with low expression (maximum mean < 0.25)
index_keep <- c()
for (i in 1:length(data)) {
    if (max(data[,i]) < 0.25) {
        next
    } else {
        index_keep <- append(index_keep,i)
    }}
int_data <- data[,index_keep]

#require mean(bin9+10) > mean(bin1+2)
index_keep <- c()
for (i in 1:length(int_data)) {
    if (((int_data[1,i] + int_data[2,i])/2) < ((int_data[9,i] + int_data[10,i])/2)){
        index_keep <- c(index_keep,i)
    }}
int_data <- int_data[,index_keep]
#Calculate moving average & require decrease every time with 2 exceptions
library(zoo)
mov_avg <- rollmeanr(int_data,3)
index_keep <- c()

for (i in 1:ncol(mov_avg)){
    for (x in 1:8){
        #if first bin, set avg value
        if (x == 1){ 
            avg <- mov_avg[x,i]
            exceptions <-  1
        if (exceptions < 0) {break}
        #if not, check if moving average increases If it doesn't, substract from exception    
        } else {
            if (avg < mov_avg[x,i]){
                avg <- mov_avg[x,i]
                next
            } else {
                avg <- mov_avg[x,i]
                exceptions <- exceptions -1
            }
        }
    }
    if (exceptions < 0) { 
        next
    } else {
        index_keep <- c(index_keep, i)
}}
int_data <- int_data[,index_keep]

#Divide each column by its maximum to have values from 1 to zero
for (i in 1:ncol(int_data)) {
    int_data[,i] <- int_data[,i] / max(int_data[,i])}
int_data <- int_data[-c(2)]
atac_periportal_data <- int_data 

rna_periportal_data <- periportal_data
rna_pericentral_data <- pericentral_data


#Plotting
plot_atac_port <- melt(atac_periportal_data, id.vars="Bins")
plot_atac_cent <- melt(atac_pericentral_data, id.vars="Bins")
plot_rna_port <- melt(rna_periportal_data, id.vars="Bins")
plot_rna_cent <- melt(rna_pericentral_data, id.vars="Bins")

ggplot(plot_rna_cent, aes(x=Bins, y=reorder(variable, -value), fill=value))  + geom_tile() + scale_fill_viridis_c() + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
ggplot(plot_rna_port, aes(x=Bins, y=reorder(variable, value), fill=value))  + geom_tile() + scale_fill_viridis_c() + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
ggplot(plot_atac_cent, aes(x=Bins, y=reorder(variable, -value), fill=value))  + geom_tile() + scale_fill_viridis_c() + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
ggplot(plot_atac_port, aes(x=Bins, y=reorder(variable, value), fill=value))  + geom_tile() + scale_fill_viridis_c() + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
