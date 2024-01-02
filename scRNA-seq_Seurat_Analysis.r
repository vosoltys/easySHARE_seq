library(Matrix)
library(data.table)  
library(Seurat)  
library(reshape2)
                             
load_Matrices <- function(path)
{
    print(path)
    mi <- read.table(path, header=TRUE)
    matrix_i <- acast(mi, gene ~ cell, value.var = "count", fill = 0)
    sparse_m <- as.sparse(matrix_i)
    return(sparse_m)
}                             

#BL6 samples
index.list <- c("002_1","083_1","007_1","084_3")
for (i in index.list){
    path <- paste("$PATH",i,".cells.combined.csv",sep="")
    sparse_m <- load_Matrices(path)
    filename <- paste("$PATH/M",i,".Robject",sep="")
    save(sparse_m, file=filename) 
}

#PWD samples, unfortunately they have wide format so slightly different code
index.list <- c("021_2","022_2","038_2","039_2")
for (i in index.list){
    path <- paste("$PATH",i,".cells.combined.csv",sep="")
    matrix <- read.table(path,header=T,sep="\t")
    matrix2 <- matrix[,-1]
    rownames(matrix2) <- matrix[,1]
    sparse_m <- as.sparse(matrix2)
    filename <- paste("$PATH/M",i,".Robject",sep="")
    save(sparse_m, file=filename)    
}

load("./M021_2.Robject")
M021_2_PWD_F_Liv <- sparse_m
load("./M022_2.Robject")
M022_2_PWD_F_Liv <- sparse_m
load("./RM038_2.Robject")
M038_2_PWD_M_Liv <- sparse_m
load("./M039_2.Robject")
M039_2_PWD_M_Liv <- sparse_m

load("./M002_1.Robject")
M002_1_BL6_M_Liv <- sparse_m
load("./M007_1.Robject")
M007_1_BL6_M_Liv <- sparse_m
load("./M083_1.Robject")
M083_1_BL6_F_Liv <- sparse_m
load("./M084_3.Robject")
M084_3_BL6_F_Liv <- sparse_m

#This is a vector of the R object names, e.g. M002_1_BL6_M_Liv, M007_1_BL6_M_Liv etc.
liver.list <- read.table("./RNA_files/liver.list",header=F)
liver.list <- list(liver.list$V1)

#making a Seurat Object for every single Dataset using the already made matrices
Liver_Samples <-list(M002_1_BL6_M_Liv,M007_1_BL6_M_Liv,M021_2_PWD_F_Liv,M022_2_PWD_F_Liv,M038_2_PWD_M_Liv,M039_2_PWD_M_Liv,M083_1_BL6_F_Liv,M084_3_BL6_F_Liv)
names(Liver_Samples) <-c("M002_1_BL6_M_Liv","M007_1_BL6_M_Liv","M021_2_PWD_F_Liv","M022_2_PWD_F_Liv","M038_2_PWD_M_Liv","M039_2_PWD_M_Liv","M083_1_BL6_F_Liv","M084_3_BL6_F_Liv")

#make all the Seurat Objects
Liver_Samples <- lapply(X = seq_along(Liver_Samples), FUN = function(i) {
    x <- paste(liver.list[[1]][i])
    Liver_Samples[i] <- CreateSeuratObject(counts = Liver_Samples[[i]], project = x, min.cells = 3, min.features = 100)
})  

#add sex 
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

#add genotype 
Liver_Samples <- lapply(X = seq_along(Liver_Samples), FUN = function(i) {
    if (grepl("_BL6_", liver.list[[1]][i])){
        genotype <- paste("B",LETTERS[12][Idents(Liver_Samples[[i]])],"6",sep="")
        names(genotype) <- colnames(Liver_Samples[[i]])    
        Liver_Samples[[i]] <- AddMetaData(Liver_Samples[[i]],metadata = genotype,col.name = 'Genotype')
        } 
    else if (grepl("_PWD_", liver.list[[1]][i])){
        genotype <- paste("PW",LETTERS[4][Idents(Liver_Samples[[i]])],sep="")
        names(genotype) <- colnames(Liver_Samples[[i]])    
        Liver_Samples[[i]] <- AddMetaData(Liver_Samples[[i]],metadata = genotype,col.name = 'Genotype')
        }    
})

names(Liver_Samples) <- c("M002_1_BL6_M_Liv","M007_1_BL6_M_Liv","M021_2_PWD_F_Liv","M022_2_PWD_F_Liv","M038_2_PWD_M_Liv","M039_2_PWD_M_Liv","M083_1_BL6_F_Liv","M084_3_BL6_F_Liv")

#For each object, plot counts and features!
lapply(X =  seq_along(Liver_Samples), FUN = function(i) {
    x <- paste("$PATH/Liver_UMIstats_",liver.list[[1]][i],".pdf",sep="")
    pdf(file= x,height=6,width=8) 
    print(VlnPlot(Liver_Samples[[i]], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2))
    print(FeatureScatter(Liver_Samples[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
    dev.off()
})

#Filtering
Liver_Samples <- lapply(X = Liver_Samples, FUN = function(i) {
    x <- subset(i, subset = nFeature_RNA > 200 & nCount_RNA > 500 & nCount_RNA < 20000)
})
save.image("RNA_All_Seurat_Objects_filtered_Liver.RData")

#Merge all the technical replicates
BL6.M_Liv <- merge(Liver_Samples$M002_1_BL6_M_Liv, y = Liver_Samples$M007_1_BL6_M_Liv, add.cell.ids = c("002", "007"), project = "BL6_M")
BL6.F_Liv <- merge(Liver_Samples$M083_1_BL6_F_Liv, y = Liver_Samples$M084_3_BL6_F_Liv, add.cell.ids = c("083", "084"), project = "BL6_F")
PWD.F_Liv <- merge(Liver_Samples$M021_2_PWD_F_Liv, y = Liver_Samples$M022_2_PWD_F_Liv, add.cell.ids = c("021", "022"), project = "PWD_F")
PWD.M_Liv <- merge(Liver_Samples$M038_2_PWD_M_Liv, y = Liver_Samples$M039_2_PWD_M_Liv, add.cell.ids = c("038", "039"), project = "PWD_M")

Liver_Reps <-list(BL6.M_Liv,BL6.F_Liv,PWD.F_Liv,PWD.M_Liv)
names(Liver_Reps) <-c("BL6.M_Liv","BL6.F_Liv","PWD.F_Liv","PWD.M_Liv")

Liver_Reps <- lapply(X = Liver_Reps, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(x)
    x <- ScaleData(x, features = all.genes)
    x <- RunPCA(x, features = VariableFeatures(object = x))  
}) 

Liver_Reps <- lapply(X = Liver_Reps, FUN = function(x) {
    x <- FindNeighbors(x, dims = 1:10)
    x <- FindClusters(x, resolution = 0.5)
    x <- RunUMAP(x, dims = 1:10) 
})

#Integrate the M and F from every species using SCT
BL6_Liver <-list(Liver_Reps$BL6.F_Liv, Liver_Reps$BL6.M_Liv)
names(BL6_Liver) <- c("BL6_F","BL6_M")
PWD_Liver <-list(Liver_Reps$PWD.F_Liv, Liver_Reps$PWD.M_Liv)
names(PWD_Liver) <- c("PWD_F","PWD_M")

#In order to sue this workflow, we need v3 assay objects
BL6_Liver$BL6_F[["RNA"]] <- as(object = BL6_Liver$BL6_F[["RNA"]], Class = "Assay")
BL6_Liver$BL6_M[["RNA"]] <- as(object = BL6_Liver$BL6_M[["RNA"]], Class = "Assay")
PWD_Liver$PWD_F[["RNA"]] <- as(object = PWD_Liver$PWD_F[["RNA"]], Class = "Assay")
PWD_Liver$PWD_M[["RNA"]] <- as(object = PWD_Liver$PWD_M[["RNA"]], Class = "Assay")

#do SCT
BL6_Liver <- lapply(X = BL6_Liver, FUN = SCTransform)
PWD_Liver <- lapply(X = PWD_Liver, FUN = SCTransform)

#prep integration
bl6.liver.features <- SelectIntegrationFeatures(object.list = BL6_Liver, nfeatures = 3000)
pwd.liver.features <- SelectIntegrationFeatures(object.list = PWD_Liver, nfeatures = 3000)
BL6_Liver <- PrepSCTIntegration(object.list = BL6_Liver, anchor.features = bl6.liver.features)
PWD_Liver <- PrepSCTIntegration(object.list = PWD_Liver, anchor.features = pwd.liver.features)

#find anchors
bl6.liver.anchors <- FindIntegrationAnchors(object.list = BL6_Liver, normalization.method = "SCT",anchor.features = bl6.liver.features)
pwd.liver.anchors <- FindIntegrationAnchors(object.list = PWD_Liver, normalization.method = "SCT",anchor.features = pwd.liver.features)
#integrate
bl6.liver.combined.sct <- IntegrateData(anchorset = bl6.liver.anchors, normalization.method = "SCT")
pwd.liver.combined.sct <- IntegrateData(anchorset = pwd.liver.anchors, normalization.method = "SCT")
#PCA etc.
bl6.liver.combined.sct <- RunPCA(bl6.liver.combined.sct, verbose = FALSE)
pwd.liver.combined.sct <- RunPCA(pwd.liver.combined.sct, verbose = FALSE)
bl6.liver.combined.sct <- RunUMAP(bl6.liver.combined.sct, reduction = "pca", dims = 1:30)
pwd.liver.combined.sct <- RunUMAP(pwd.liver.combined.sct, reduction = "pca", dims = 1:30)

#Lastly, integrate bl6 & pwd
Liver_Int <- list(pwd.liver.combined.sct,bl6.liver.combined.sct)
names(Liver_Int) <- c("pwd.liver.combined.sct","bl6.liver.combined.sct")
liver1.features <- SelectIntegrationFeatures(object.list = Liver_Int, nfeatures=1000) 
liver1.anchors <- FindIntegrationAnchors(object.list = Liver_Int, anchor.features = liver1.features) 
liver.combined.sct <- IntegrateData(anchorset = liver1.anchors)
liver.combined.sct <- ScaleData(liver.combined.sct, verbose=FALSE)
liver.combined.sct <- RunPCA(liver.combined.sct, npcs = 30, verbose = FALSE) 
liver.combined.sct <- FindNeighbors(liver.combined.sct, dims = 1:30)
liver.combined.sct <- FindClusters(liver.combined.sct, resolution = 0.2)
liver.combined.sct <- RunUMAP(liver.combined.sct, dims = 1:30) 

DefaultAssay(liver.combined.sct) <- "SCT"
DimPlot(liver.combined.sct, group.by="seurat_clusters", label=FALSE, shuffle=TRUE) + ggtitle("SCT Integrated Liver w/ SCT int. biol. reps")
DimPlot(liver.combined.sct, group.by="Genotype", label=FALSE, shuffle=TRUE) + ggtitle("SCT Integrated Liver w/ SCT int. biol. reps")
dev.off()

save(liver.combined.sct, file="RNA_Liver_All_Integrated_Object.RData")
