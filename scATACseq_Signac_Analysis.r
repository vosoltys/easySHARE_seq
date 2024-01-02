library(ChIPpeakAnno)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(dplyr)
library(rtracklayer)
library(EnsDb.Mmusculus.v79)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(patchwork)
library(gdata)
library(tidyverse)
library(devtools)

#function for loading data
load_ATAC_matrices <- function(number, tissue){
    fragment_path <- paste(as.character(number),".txt.frags.sort.bed.gz",sep="") #fragment file from respective sample created with sinto
    cells_path <- paste(as.character(number),".stats.Cells",sep="") #to shorten loading time, only used top 5k barcodes from each sample
    cells <- read.table(cells_path,header=FALSE)
    fragments <- CreateFragmentObject(path = fragment_path, cells=cells$V1, verbose = TRUE, validate.fragments = FALSE)
    peak_matrix <- FeatureMatrix(fragments = fragments, features = tissue)
    total <- list(peak_matrix,fragments,cells)
    return(total)
}

#loading raw files
liv_peaks <- toGRanges("$PATH/B6_PWD_Liver.peaks_peaks.narrowPeak", format="narrowPeak", header=FALSE)
#This is a list of r object names, needs to be in the same order as the list created below (Liver_Samples)
liver.list <- read.table("liver.list",header=F)
liver.list <- list(liver.list$V1)

M046_2_PWD_M_Liv <- load_ATAC_matrices("D046_2", liv_peaks)
M047_2_PWD_M_Liv <- load_ATAC_matrices("D047_2", liv_peaks)
M030_2_PWD_F_Liv <- load_ATAC_matrices("D030_2", liv_peaks)
M029_2_PWD_F_Liv <- load_ATAC_matrices("D029_2", liv_peaks)
M075_1_BL6_F_Liv <- load_ATAC_matrices("D075_1", liv_peaks)
M076_3_BL6_F_Liv <- load_ATAC_matrices("D076_3", liv_peaks)
M090_1_BL6_M_Liv <- load_ATAC_matrices("D090_1", liv_peaks)
M095_1_BL6_M_Liv <- load_ATAC_matrices("D095_1", liv_peaks)

Liver_Samples <- list(M090_1_BL6_M_Liv,M075_1_BL6_F_Liv,M095_1_BL6_M_Liv,M076_3_BL6_F_Liv, M029_2_PWD_F_Liv, M030_2_PWD_F_Liv, M046_2_PWD_M_Liv, M047_2_PWD_M_Liv)
names(Liver_Samples) <- c("M090_1_BL6_M_Liv","M075_1_BL6_F_Liv","M095_1_BL6_M_Liv","M076_3_BL6_F_Liv","M029_2_PWD_F_Liv","M030_2_PWD_F_Liv","M046_2_PWD_M_Liv", "M047_2_PWD_M_Liv")

#make the objects
Liver_Samples <- lapply(X=seq_along(Liver_Samples), FUN = function(i){
    chrom_assay <- CreateChromatinAssay(counts = Liver_Samples[[i]][[1]],sep = c("-", "-"),genome = 'mm10',fragments = Liver_Samples[[i]][[2]], min.cells = 10,min.features = 100)
    y <- paste(liver.list[[1]][i])
    Liver_Samples[i] <- CreateSeuratObject(counts = chrom_assay,assay = "ATAC", project = y)        
})

#add sex and genotype information
Liver_Samples <- lapply(X = seq_along(Liver_Samples), FUN = function(i) {
    if (grepl("F", liver.list[[1]][i])){
        sex <- LETTERS[6][Idents(Liver_Samples[[i]])]
        names(sex) <- colnames(Liver_Samples[[i]])    
        Liver_Samples[[i]] <- AddMetaData(Liver_Samples[[i]],metadata = sex,col.name = 'Sex')
    } else {
      sex <- LETTERS[13][Idents(Liver_Samples[[i]])]
      names(sex) <- colnames(Liver_Samples[[i]])
      Liver_Samples[[i]] <- AddMetaData(Liver_Samples[[i]],metadata = sex,col.name = 'Sex')  
    }
}) 

Liver_Samples <- lapply(X = seq_along(Liver_Samples), FUN = function(i) {
    if (grepl("_BL6_", liver.list[[1]][i])){
        genotype <- paste("B",LETTERS[12][Idents(Liver_Samples[[i]])],"6",sep="")
        names(genotype) <- colnames(Liver_Samples[[i]])    
        Liver_Samples[[i]] <- AddMetaData(Liver_Samples[[i]],metadata = genotype,col.name = 'Genotype')
        } 
    else if (grepl("_PWD_", liver.list[[1]][i])){
        genotype <- paste("PWD",LETTERS[5][Idents(Liver_Samples[[i]])],"T",sep="")
        names(genotype) <- colnames(Liver_Samples[[i]])    
        Liver_Samples[[i]] <- AddMetaData(Liver_Samples[[i]],metadata = genotype,col.name = 'Genotype')
        } 
}) 

names(Liver_Samples) <- c("M090_1_BL6_M_Liv","M075_1_BL6_F_Liv","M095_1_BL6_M_Liv","M076_3_BL6_F_Liv","M029_2_PWD_F_Liv","M030_2_PWD_F_Liv","M046_2_PWD_M_Liv", "M047_2_PWD_M_Liv")

#FRiP Score
Liver_Samples <- lapply(X = seq_along(Liver_Samples), FUN = function(i) {
    one <- gsub("M","D", substr(liver.list[[1]][i],1,6))
    fragment_path <- paste(one,".txt.frags.sort.bed.gz",sep="")
    frags <- CountFragments(fragment_path)
    new_fragments <- dplyr::filter(frags,frags$CB %in% colnames(Liver_Samples[[i]])) 
    colOrd <- rownames(FetchData(Liver_Samples[[i]],"ident"))
    RcolData <- new_fragments %>% dplyr::slice(match(colOrd, CB))
    Liver_Samples[[i]] <- AddMetaData(object=Liver_Samples[[i]],  metadata=RcolData$frequency_count, col.name = 'fragments')
    Liver_Samples[[i]] <- FRiP(Liver_Samples[[i]], 'ATAC','fragments',col.name = "FRiP", verbose = TRUE)
})

names(Liver_Samples) <- c("M090_1_BL6_M_Liv","M075_1_BL6_F_Liv","M095_1_BL6_M_Liv","M076_3_BL6_F_Liv","M029_2_PWD_F_Liv","M030_2_PWD_F_Liv","M046_2_PWD_M_Liv", "M047_2_PWD_M_Liv")

#We need to change names so that they match between the RNA-seq and ATAC-seq when we integrate downstream. Currently, they have different indices (barcode combination is the same!)
#The uploaded single-cell data in GEO XXXX has its names already changed! 
Liver_Samples <- list(Liver_Samples$M090_1_BL6_M_Liv,Liver_Samples$M075_1_BL6_F_Liv,Liver_Samples$M095_1_BL6_M_Liv,Liver_Samples$M076_3_BL6_F_Liv,Liver_Samples$M029_2_PWD_F_Liv,Liver_Samples$M030_2_PWD_F_Liv,Liver_Samples$M046_2_PWD_M_Liv, Liver_Samples$M047_2_PWD_M_Liv)
names(Liver_Samples) <-c("M090_1_BL6_M_Liv","M075_1_BL6_F_Liv","M095_1_BL6_M_Liv","M076_3_BL6_F_Liv","M029_2_PWD_F_Liv","M030_2_PWD_F_Liv","M046_2_PWD_M_Liv", "M047_2_PWD_M_Liv")

#in R get a vector of cell names from each object 
lapply(X = seq_along(Liver_Samples), FUN = function(i){
    x <- paste(liver.list[[1]][i],".txt",sep="")
    write.table(colnames(Liver_Samples[[i]]), file=x)
})

#This line is in bash
#for sample in *Liv.txt; do tail -n+2 "$sample" | sed 's/D075/D083/g' | sed 's/D076/D084/g' | sed 's/D090/D002/g'| sed 's/D095/D007/g'| sed 's/D029/D021/g' | sed 's/D030/D022/g' | sed 's/D046/D038/g' | sed 's/D047/D039/g' > "$sample"_ATAC_new_cell_names; done
new_M029_2_PWD_F_Liv <- read.table("M029_2_PWD_F_Liv.txt_ATAC_new_cell_names")
names_M029_2_PWD_F_Liv <- new_M029_2_PWD_F_Liv$V2
Liver_Samples$M029_2_PWD_F_Liv <- RenameCells(Liver_Samples$M029_2_PWD_F_Liv,new.names = names_M029_2_PWD_F_Liv)

new_M030_2_PWD_F_Liv <- read.table("M030_2_PWD_F_Liv.txt_ATAC_new_cell_names")
names_M030_2_PWD_F_Liv <- new_M030_2_PWD_F_Liv$V2
Liver_Samples$M030_2_PWD_F_Liv <- RenameCells(Liver_Samples$M030_2_PWD_F_Liv,new.names = names_M030_2_PWD_F_Liv)

new_M046_2_PWD_M_Liv <- read.table("M046_2_PWD_M_Liv.txt_ATAC_new_cell_names")
names_M046_2_PWD_M_Liv <- new_M046_2_PWD_M_Liv$V2
Liver_Samples$M046_2_PWD_M_Liv <- RenameCells(Liver_Samples$M046_2_PWD_M_Liv,new.names = names_M046_2_PWD_M_Liv)

new_M047_2_PWD_M_Liv <- read.table("M047_2_PWD_M_Liv.txt_ATAC_new_cell_names")
names_M047_2_PWD_M_Liv <- new_M047_2_PWD_M_Liv$V2
Liver_Samples$M047_2_PWD_M_Liv <- RenameCells(Liver_Samples$M047_2_PWD_M_Liv,new.names = names_M047_2_PWD_M_Liv)

new_M075_1_BL6_F_Liv <- read.table("M075_1_BL6_F_Liv.txt_ATAC_new_cell_names")
names_M075_1_BL6_F_Liv <- new_M075_1_BL6_F_Liv$V2
Liver_Samples$M075_1_BL6_F_Liv <- RenameCells(Liver_Samples$M075_1_BL6_F_Liv,new.names = names_M075_1_BL6_F_Liv)

new_M076_3_BL6_F_Liv <- read.table("M076_3_BL6_F_Liv.txt_ATAC_new_cell_names")
names_M076_3_BL6_F_Liv <- new_M076_3_BL6_F_Liv$V2
Liver_Samples$M076_3_BL6_F_Liv <- RenameCells(Liver_Samples$M076_3_BL6_F_Liv,new.names = names_M076_3_BL6_F_Liv)

new_M090_1_BL6_M_Liv <- read.table("M090_1_BL6_M_Liv.txt_ATAC_new_cell_names")
names_M090_1_BL6_M_Liv <- new_M090_1_BL6_M_Liv$V2
Liver_Samples$M090_1_BL6_M_Liv <- RenameCells(Liver_Samples$M090_1_BL6_M_Liv,new.names = names_M090_1_BL6_M_Liv)

new_M095_1_BL6_M_Liv <- read.table("M095_1_BL6_M_Liv.txt_ATAC_new_cell_names")
names_M095_1_BL6_M_Liv <- new_M095_1_BL6_M_Liv$V2
Liver_Samples$M095_1_BL6_M_Liv <- RenameCells(Liver_Samples$M095_1_BL6_M_Liv,new.names = names_M095_1_BL6_M_Liv)


#Add Annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(Liver_Samples[[1]]) <- annotations
Annotation(Liver_Samples[[2]]) <- annotations
Annotation(Liver_Samples[[3]]) <- annotations
Annotation(Liver_Samples[[4]]) <- annotations
Annotation(Liver_Samples[[5]]) <- annotations
Annotation(Liver_Samples[[6]]) <- annotations
Annotation(Liver_Samples[[7]]) <- annotations
Annotation(Liver_Samples[[8]]) <- annotations

#QC
Liver_Samples <- lapply(X = Liver_Samples, FUN = function(x) {
    x <- UpdateSeuratObject(x)
})
Liver_Samples <- lapply(X = Liver_Samples, FUN = function(x) {
    x <- NucleosomeSignal(object=x)
    x$nucleosome_group <- ifelse(x$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
    x <- TSSEnrichment(object = x, fast = FALSE)
})

Liver_Samples[[1]]$high.tss <- ifelse(Liver_Samples[[1]]$TSS.enrichment > 3, 'High', 'Low')
Liver_Samples[[2]]$high.tss <- ifelse(Liver_Samples[[2]]$TSS.enrichment > 3, 'High', 'Low')
Liver_Samples[[3]]$high.tss <- ifelse(Liver_Samples[[3]]$TSS.enrichment > 3, 'High', 'Low')
Liver_Samples[[4]]$high.tss <- ifelse(Liver_Samples[[4]]$TSS.enrichment > 3, 'High', 'Low')
Liver_Samples[[5]]$high.tss <- ifelse(Liver_Samples[[5]]$TSS.enrichment > 3, 'High', 'Low')
Liver_Samples[[6]]$high.tss <- ifelse(Liver_Samples[[6]]$TSS.enrichment > 3, 'High', 'Low')
Liver_Samples[[7]]$high.tss <- ifelse(Liver_Samples[[7]]$TSS.enrichment > 3, 'High', 'Low')
Liver_Samples[[8]]$high.tss <- ifelse(Liver_Samples[[8]]$TSS.enrichment > 3, 'High', 'Low')
Liver_Samples[[1]]$nucleosome_group <- ifelse(Liver_Samples[[1]]$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
Liver_Samples[[2]]$nucleosome_group <- ifelse(Liver_Samples[[2]]$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
Liver_Samples[[3]]$nucleosome_group <- ifelse(Liver_Samples[[3]]$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
Liver_Samples[[4]]$nucleosome_group <- ifelse(Liver_Samples[[4]]$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
Liver_Samples[[5]]$nucleosome_group <- ifelse(Liver_Samples[[5]]$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
Liver_Samples[[6]]$nucleosome_group <- ifelse(Liver_Samples[[6]]$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
Liver_Samples[[7]]$nucleosome_group <- ifelse(Liver_Samples[[7]]$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
Liver_Samples[[8]]$nucleosome_group <- ifelse(Liver_Samples[[8]]$nucleosome_signal > 4, 'NS > 4', 'NS < 4')


#Filtering
Liver_Samples <- lapply(X =  Liver_Samples, FUN = function(x) {
    subset(x, subset = nCount_ATAC > 300 & nCount_ATAC < 15000 & FRiP > 0.3 & nucleosome_signal < 4 & TSS.enrichment > 2)
})

#Merge technical replicates
PWD.F_Liv <- merge(Liver_Samples$M029_2_PWD_F_Liv, y = Liver_Samples$M030_2_PWD_F_Liv, add.cell.ids = c("021", "022"))
PWD.M_Liv <- merge(Liver_Samples$M046_2_PWD_M_Liv, y = Liver_Samples$M047_2_PWD_M_Liv, add.cell.ids = c("038", "039"))
BL6.M_Liv <- merge(Liver_Samples$M090_1_BL6_M_Liv, y = Liver_Samples$M095_1_BL6_M_Liv, add.cell.ids = c("002", "007"))
BL6.F_Liv <- merge(Liver_Samples$M075_1_BL6_F_Liv, y = Liver_Samples$M076_3_BL6_F_Liv, add.cell.ids = c("083", "084"))

Liver_Reps <-list(BL6.M_Liv,BL6.F_Liv,PWD.F_Liv,PWD.M_Liv)
names(Liver_Reps) <-c("BL6.M_Liv","BL6.F_Liv","PWD.F_Liv","PWD.M_Liv")

#process to UMAP
Liver_Reps <- lapply(X = Liver_Reps, FUN = function(x) {
    x <- FindTopFeatures(x, min.cutoff = 'q0')
    x <- RunTFIDF(x)
    x <- RunSVD(x)
}) 

Liver_Reps <- lapply(X = Liver_Reps, FUN = function(x) {
    x <- RunUMAP(object = x, reduction = 'lsi', dims = 2:30)
    x <- FindNeighbors(object = x, reduction = 'lsi', dims = 2:30)
    x <- FindClusters(object = x, verbose = FALSE, algorithm = 3)
})

#Integrate all samples together
#make list
Liv_Comb <- list(Liver_Reps$PWD.M_Liv,Liver_Reps$PWD.F_Liv,Liver_Reps$BL6.M_Liv,Liver_Reps$BL6.F_Liv)
names(Liv_Comb) <- c("PWD.M","PWD.F","BL6.M_Liv","BL6.F_Liv")
liver.anchors <- FindIntegrationAnchors(Liv_Comb,anchor.features = rownames(Liv_Comb$BL6_Liv),reduction = "rlsi",dims = 2:30, l2.norm=TRUE, l.anchor=50)

temp.Liv <- merge(Liver_Reps$PWD.M_Liv,y=c(Liver_Reps$PWD.F_Liv,Liver_Reps$BL6.M_Liv,Liver_Reps$BL6.F_Liv))
temp.Liv <- FindTopFeatures(temp.Liv, min.cutoff = 'q0')
temp.Liv <- RunTFIDF(temp.Liv)
temp.Liv <- RunSVD(temp.Liv)
temp.Liv <- FindNeighbors(temp.Liv, reduction = 'lsi', dims = 1:20)
temp.Liv <- RunUMAP(object = temp.Liv, reduction = 'lsi', dims = 1:20)
temp.Liv <- FindClusters(temp.Liv, verbose = FALSE, algorithm = 3)

# integrate LSI embeddings
Liver_integrated <- IntegrateEmbeddings(anchorset = liver.anchors,reductions = temp.Liv[["lsi"]], new.reduction.name = "integrated_lsi", dims.to.integrate = 3:30, k.weight=50)

# create a new UMAP using the integrated embeddings
Liver_integrated <- FindNeighbors(Liver_integrated, reduction = 'integrated_lsi', dims =6:20)
Liver_integrated <- RunUMAP(Liver_integrated, reduction = "integrated_lsi", dims = 6:20)
Liver_integrated <- FindClusters(Liver_integrated, verbose = FALSE, algorithm = 3)

DimPlot(Liver_integrated)
DimPlot(Liver_integrated, group.by = "Genotype")
DimPlot(Liver_integrated, group.by = "orig.ident")

save(Liver_Integrated, file="ATAC_BL6_PWD_All_Integrated_Object_ATAC.RData")

