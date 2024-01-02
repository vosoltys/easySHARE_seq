library(SummarizedExperiment)
library(chromVAR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(matrixStats)
library(ChIPpeakAnno)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(dplyr)
library(rtracklayer)
library(SeuratWrappers)
library(Matrix)
library(patchwork)
library(gdata)
library(tidyverse)
library(devtools)
library(BSDA)

calcGCContent = function(query, ref) {
    .validateInputs(list(query=c("GRanges","GRangesList"),
                         ref="BSgenome"))
    if (is(query, "GRangesList")) {
        # Recurse over each GRanges object
        x = lapply(query, calcGCContent, ref)
        namelist = names(query)
        if (is.null(namelist)) {
            newnames = seq_along(query)
            namelist = newnames
            # Append names
            names(x) = namelist
        }
        return(x)
    }
    # Restrict the seqnames to known chromosomes
    query = GenomeInfoDb::keepStandardChromosomes(query, pruning.mode="coarse")
    v = IRanges::Views(ref, query)
    gcvec = apply(data.frame(Biostrings::alphabetFrequency(v)[,c("C","G")]),1,sum)/width(v)
    return(gcvec)
}

.validateInputs = function(checkList) {
    nms = names(checkList)
    for(i in seq_along(checkList)){
        fail = FALSE
        clss = checkList[[i]]
        x = get(nms[i], envir=parent.frame(1))
        for(cls in clss){
            if (is(x, cls)) fail = append(fail, TRUE)
        }
        if(!any(fail)) 
            stop(paste0(nms[i], " must be a ", paste(clss, collapse=" or "), 
                        ".  Got: ", class(x)))
    }
}


#Load SC dataset & process 
load("Liver_WNN_broad_annotated.RObject") #Liver_WNN_broad, RNA-seq
load("ATAC_BL6_PWD_All_Integrated_Object_ATAC.RData") #Liver_integrated, ATAC-seq
cell_types <- read.table("CellType_IDs", header=T, sep="\t")
rna_cells <- rownames(cell_types)
atac_cells <- colnames(Liver_integrated)
common <- intersect(rna_cells,atac_cells)
Liver_integrated[["CellName"]] <- colnames(Liver_integrated)
Liver_integrated <- subset(Liver_integrated, subset = CellName %in% common)
atac_cells <- colnames(Liver_integrated)
cell_types$ID <- rownames(cell_types)
cell_types <- cell_types %>% arrange(factor(ID, levels=atac_cells))
Liver_integrated[["CellType"]] <- cell_types$Idents.Liver_WNN_broad.
Liver_WNN_broad$CellType <- Idents(Liver_WNN_broad)
Liver_WNN_broad[["ATAC"]] <- Liver_integrated[["ATAC"]]
Liver_subset <- subset(Liver_WNN_broad, subset= CellType =="LSECs")

#load peakset (peaks.use is the subset used by LSECs)
peaks_new <- toGRanges("peaks.use", format="BED", header=FALSE)

#get background peaks per peak
coun <- AverageExpression(Liver_WNN_broad, assay="ATAC", slot="data", group.by="CellType")
counts <- data.frame(data.frame(coun)[,2]) #column 2 is LSECs
rownames(counts) <- rownames(data.frame(coun))
counts <- counts[order(row.names(counts)), ,drop=F]
names(counts) <- "count"
rownames(counts) <- names(peaks_new)
counts$count[counts$count == 0] <- 0.000001
counts <- as.matrix(counts)
pe_names <- read.table("peaks.use")
pe_names$name <- paste(pe_names$V1, pe_names$V2, pe_names$V3, sep="-")
rownames(counts) <- pe_names$name
names(peaks_new) <- pe_names$name
background_matrix <- getBackgroundPeaks(counts, bias=peaks_new$GC_Bias,niterations=100)
rownames(background_matrix) <- names(peaks_new)

#Make a list of expressed TSS
tss_pos <- read.table("TSS.positions", sep="\t")
expr_genes <- data.frame(data.frame(AverageExpression(Liver_WNN_broad, assay="RNA",group.by="CellType"))[,2])
rownames(expr_genes) <- rownames(data.frame(AverageExpression(Liver_WNN_broad, assay="RNA",group.by="CellType")))
colnames(expr_genes) <- "Expr"
expr_genes <- expr_genes%>% dplyr::filter(Expr > 0)
tss_pos <- tss_pos %>% dplyr::filter(V4 %in% rownames(expr_genes))

#Read in chromosome length so we don't overshoot when selecting peaks within range
lengths <- read.table("Chr_lengths.txt", row.names=1)

#Calculate Associations
#This takes a very long time - split into ~20 by adjusting the first for loop
peak_gene_lookup_list <- list()
for (z in 1:19615){  #number of expressed genes
    output_dataframe <- data.frame()
    gene <- tss_pos[z,4]
    TSS <- data.frame(tss_pos[z,1], tss_pos[z,2], tss_pos[z,3])
    names(TSS) <- c("chr","start","end")
    #This here determines interval we look at. Here: +- 500kb
    TSS$start <- as.numeric(TSS$start) - 500000
    if (TSS$start < 0){TSS$start=0}
    TSS$end <- as.numeric(TSS$end) + 500000
    if (TSS$end > lengths[TSS$chr,]) {TSS$end = lengths[TSS$chr,]}
    TSS <- makeGRangesFromDataFrame(TSS)
    peaks_tss <- subsetByOverlaps(peaks_new,TSS)
    if (length(peaks_tss) == 0) next
    #Get counts for the gene in question
    gene_count <- data.frame(Liver_subset@assays$RNA@counts[gene,])
    names(gene_count) <- "counts"
    #setup summary dataframe
    tss_cors <- data.frame()
    background_correlations <- data.frame()
    #Calculate focal and background correlation
    for (i in 1:length(peaks_tss)){
        #get name of focal peak
        chr <- peaks_tss@seqnames[i]@values
        start <- peaks_tss@ranges@start[i] -1 
        width <-  peaks_tss@ranges@width[i]
        end <- (start + width) 
        peak <- paste(chr,start,end,sep="-")
        #calculate soeamrna of peak-gene
        peak_count <- data.frame(Liver_subset@assays$ATAC@counts[peak,])
        names(peak_count) <- "counts"
        spm <- cor(peak_count, gene_count)
        colnames(spm) <- "spearman"
        rownames(spm) <- peak
        tss_cors <- rbind(tss_cors, spm)
        #Calculate background spearman
        background_indices <- background_matrix[names(peaks_tss[i]),]
        background_cors=data.frame()
        for (x in seq(1:100)){
            chr <- peaks_new@seqnames[background_indices[x]]@values
            start <- peaks_new@ranges@start[background_indices[x]] -1 
            width <-  peaks_new@ranges@width[background_indices[x]]
            end <- (start + width) 
            background_peak <- paste(chr,start,end,sep="-")
            #calculate spearman of peak-gene
            background_peak_count <- data.frame(Liver_subset@assays$ATAC@counts[background_peak,])
            names(background_peak_count) <- "counts"
            background_spm <- cor(background_peak_count, gene_count)
            colnames(background_spm) <- "spearman"
            rownames(background_spm) <- x
            background_cors <- rbind(background_cors, background_spm)
        } #end background loop
        
        colnames(background_cors) <- names(peaks_tss[i])
        background_correlations <- rbind(background_correlations, t(background_cors))
        ## tss_cors: dataframe of true Spearman correlations of every peak within X kb to the focal gene. rownames = peaks
        ## background_correlations: df of 100 background correlations for every peak in 500 kb to the focal gene.rownames = peaks, 100 cols
        ## Next step: statistical testing & making final output!
        pop_mean <- rowMeans(background_correlations[names(peaks_tss[i]),], na.rm=TRUE)
        pop_sd <- sd(t(background_correlations[names(peaks_tss[i]),])[,1], na.rm=TRUE)
        z_score <- (tss_cors[i,1] - pop_mean) / pop_sd
        pval <- pnorm(z_score, lower.tail=FALSE)
        focal_tss_cor <- tss_cors[i,1]
        temp <- data.frame(cbind(pval, focal_tss_cor,z_score, pop_mean, pop_sd))
        output_dataframe <- rbind(output_dataframe, temp)
    } #end gene loop. 
    peak_gene_lookup_list[[gene]] <- output_dataframe
    if (z %% 100 == 0 ) {
        save(peak_gene_lookup_list, file=("Peak_Gene_Cor_500kb")) #save intermediate results
    }
}

#Summarize results
sig_con_per_gene <- data.frame()
all_nonsig_peaks <- data.frame()
sig_peaks_names <- c()
big_summary <- data.frame()
tss_pos <- read.table("tss_positions")
colnames(tss_pos) <- c("chr","start","end","gene","strand")
#:length(LSEC_peak_gene_connections)
for (i in 13930:length(LSEC_peak_gene_connections)){
    if (i %% 100 == 0) {print(i)} #just for checking how fast it is
    sig_hits <- dplyr::filter(LSEC_peak_gene_connections[[i]], pval < 0.05)
    fin_sig_hits <- data.frame()
    gene <- names(LSEC_peak_gene_connections)[i]
    tss <- (tss_pos %>% dplyr::filter(gene == names(LSEC_peak_gene_connections)[i]))[1,2]  
    for (x in 1:nrow(sig_hits)){
        if (nrow(sig_hits) == 0){next}
        chr <- str_split(rownames(sig_hits)[x],"-")[[1]][1]
        start <- str_split(rownames(sig_hits)[x],"-")[[1]][2]
        end <- str_split(rownames(sig_hits)[x],"-")[[1]][3]
        peakname <- rownames(sig_hits)[x]
        midpoint <- as.numeric(start) + ((as.numeric(end)-as.numeric(start))/2)
        newline <- cbind(peakname,chr,start,end,midpoint,names(LSEC_peak_gene_connections)[i],tss,sig_hits[x,1],sig_hits[x,2],sig_hits[x,3],sig_hits[x,4],sig_hits[x,5])
        fin_sig_hits <- rbind(fin_sig_hits,newline)}
    sig_peaks_names <- c(sig_peaks_names,rownames(sig_hits))
    temp2 <- dplyr::filter(LSEC_peak_gene_connections[[i]], pval > 0.05)
    all_nonsig_peaks <- rbind(all_nonsig_peaks,temp2)
    Num_assoc <- nrow(sig_hits)
    gene <- names(LSEC_peak_gene_connections)[i]
    table <- data.frame(cbind(gene,Num_assoc))
    names(table) <- c("Gene","Num_Assoc")
    sig_con_per_gene <- rbind(sig_con_per_gene, table)
    big_summary <- rbind(big_summary, fin_sig_hits)
}
big_summary$midpoint <- round(as.numeric(big_summary$midpoint))
names(big_summary) <- c("peakname","chr","start","end","midpoint","gene","tss","pval","corr","zscore","pop_mean","pop_sd")
write.table(big_summary, file="LSEC_peak_gene_associations_all_500kb.txt", quote=F, sep="\t")

sig_sum <- sig_con_per_gene %>% group_by(Num_Assoc) %>% summarise(count=n())

#How many peaks were assigned multiple times?
genes_per_peak <- data.frame(table(big_summary$peakname))
line <- data.frame(Var1=c(1:79931), Freq = 0)
genes_per_peak <- rbind(genes_per_peak,line)

#rankorder
sig_con_per_gene$num <- as.numeric(sig_con_per_gene$Freq)
sig_con_per_gene <- sig_con_per_gene[order(sig_con_per_gene$num),]
sig_con_per_gene$rank <- seq(1:15498)


#Distance to tss
tss_positions <- read.table("tss_positions")
tss <- regioneR::toGRanges(tss_positions)

#all peaks 
all_peaks <- toGRanges("peaks.use", format="BED", header=FALSE)
nearest_tss <- GenomicRanges::nearest(all_peaks, tss, ignore.strand=FALSE)

#significant peaks
all_sig_peaks_dist <- data.frame()
for (i in 1:nrow(all_sig_peaks)){
    chr <- str_split(rownames(all_sig_peaks)[i],"-")[[1]][1]
    start <- str_split(rownames(all_sig_peaks)[i],"-")[[1]][2]
    end <- str_split(rownames(all_sig_peaks)[i],"-")[[1]][3]
    entry <- data.frame(chr,start,end)
    all_sig_peaks_dist <- rbind(all_sig_peaks_dist,entry)
}
all_peaks <- toGRanges(all_sig_peaks_dist)
nearest_tss <- GenomicRanges::nearest(all_peaks, tss, ignore.strand=FALSE)

all_peak_distances <- c()
for (i in 1:length(all_peaks)){
    peak_start <- all_peaks[i]@ranges@start
    peak_width <- all_peaks[i]@ranges@width
    tss_position <- tss[nearest_tss[i]]@ranges@start
    strand <- tss[nearest_tss[i]]$V5
    #First, if peak > tss and positive strand --> positive distance
    if ((peak_start > tss_position) & strand == "+")
        {
        distance = peak_start - tss_position
        all_peak_distances <- c(all_peak_distances, distance)
        next
    #Second, peak > tss and negative strand --> negative distance    
    } else if ((peak_start > tss_position) & strand == "-"){
        distance = -1 * (peak_start - tss_position)
        all_peak_distances <- c(all_peak_distances, distance)
        next
    # Third, peak overlaps tss 
    } else if ((peak_start < tss_position) & ((peak_start + peak_width) > tss_position)){
        distance = 0
        all_peak_distances <- c(all_peak_distances, distance)
        next
    # Fourth, peak_start < tss & positive strand --> negative distance    
    } else if ((peak_start < tss_position)  & strand =="+"){
        distance = -1 * (tss_position - peak_start)
        all_peak_distances <- c(all_peak_distances, distance)
        next
    # Lastly, peak start < tss & negative strand --> positive distance    
    } else if ((peak_start < tss_position)  & strand =="-"){
        distance = tss_position - peak_start
        all_peak_distances <- c(all_peak_distances, distance)
    } 
}

all_dist <- data.frame(all_peak_distances)
all_dist$set <- "All"
colnames(all_dist) <- c("Distance","Set")
sig_dist <- data.frame(all_sig_peak_distances)
sig_dist$set <- "Significant"
colnames(sig_dist) <- c("Distance","Set")
comb_dist <- rbind( all_dist,sig_dist)
comb_dist_plot <- melt(comb_dist, id.vars="Set")
comb_dist_plot <- comb_dist_plot[,c(1,3)]

#Plotting a showcase example
focal_links <- (dplyr::filter(big_summary, (chr == "chr14") & (start > 62830000) & (end < 63320000)))[,c(2,5,7,9,10,8,6,1)]
focal_links <-  transform(focal_links, corr = as.numeric(corr), zscore=as.numeric(zscore),pval = as.numeric(pval))
names(focal_links) <- c("chr","start","end","score","zscore","pvalue","gene","peak")
for (i in 1:nrow(focal_links)){
    if (as.numeric(focal_links[i,]$start) > as.numeric(focal_links[i,]$end)){
        te_va <- focal_links[i,]$start
        focal_links[i,]$start <- focal_links[i,]$end
        focal_links[i,]$end <- te_va}}
focal_links <- toGRanges(focal_links)
Links(Liver_subset[["ATAC"]]) <- focal_links

CoveragePlot(Liver_subset, region="chr14-62830000-63320000",extend.upstream = 0,extend.downstream = 0, show.bulk=TRUE,scale.linewidth=TRUE)
