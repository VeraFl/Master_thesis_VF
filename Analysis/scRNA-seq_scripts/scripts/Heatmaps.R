#Heatmaps

#load libraries
library(data.table)
library(VennDiagram)
library(RColorBrewer)
library(ComplexHeatmap)
library(dplyr)
library(SeuratData)
library(SeuratDisk)
library(Seurat)
library(patchwork)
library(readr)
library(gdata)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(tidyseurat)
library(grid)
library(viridis)
library(gridExtra)
library(ComplexHeatmap)
library(ggplot2)
library(ggcorrplot)
library(scales)
library(grDevices)
library(colorspace)
library(forcats)
library(glmGamPoi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(EnhancedVolcano)
library(ggrepel)
library(DESeq2, quietly=T)
library(ComplexHeatmap, quietly = T)
library(circlize, quietly = T)
library(DEGreport, quietly=T)
library(tidyverse, quietly = T)
library(data.table, quietly = T)
library(clusterProfiler, quietly = T)
library(enrichplot, quietly = T)
library(ReactomePA, quietly = T)
library(ggVennDiagram,quietly = T)
library(PCAtools, quietly = T)
library(gprofiler2)
library(biomaRt)
library(gprofiler2)
library(ggVennDiagram)
library(ggupset)
library(MAST)

memory.limit(24000)
##################################################################################
#Camel Heatmap
# Create color palettes 
conditionPalette <- structure(c("#F39C6B","#A5668B"), names = c("ACN4", "MERS"))
celltypePalette <- structure(c("#807EBF", "#7EA8BE", "#089392", "#40B48B", "#9DCD84", "#EAE29C"), 
                             names = c('Cluster1','Cluster2','Secretory', 'Ciliated', 'Club', 'Basal'))
statusPalette <- structure(c("#519872", "#CD4631"), names = c("Bystanders", "Infected"))

#Later for Llama
#(c("Llama Cluster 1", "Llama Cluster 2", "Secretory", "Ciliated", "Club", "Basal"))
#c("#8A6E7E", "#40679B", "#089392", "#40B48B", "#9DCD84", "#EAE29C")

# Import DEG files
#Loading filtered and merged Seurat Object into Script
#Infected and Mock Samples merged and cell clusters assigned

# set wd
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Ferus_both")
# import files
filelist = list.files(pattern="*.csv", full.names = TRUE)
datalist <- lapply(filelist, FUN= read.csv, na.strings="")
names(datalist) <- gsub(".txt", "", basename(filelist))

# Define viral genes
virus_genes <- c("mers-3UTR","mers-orf8b","mers-N","mers-M","mers-E","mers-orf5","mers-orf4b","mers-orf4a","mers-orf3","mers-S","mers-orf1ab","mers-5UTR",
                 "camel229E-3UTR","camel229E-ORF8","camel229E-N","camel229E-M","camel229E-E","camel229E-ORF4","camel229E-S","camel229E-ORF1a", "camel229E-ORF1ab","camel229E-5UTR")

# Modify datalist
new_colnames <- c("Gene", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")
names(datalist) <- gsub(".csv", "", basename(filelist))
datalist <- lapply(datalist, setNames, new_colnames) # Rename column names in list
datalist <- lapply(datalist, function(x) x[!x$Gene %in% virus_genes, ]) # Remove viral genes from files



##########################################################################################################

# Figure 5C
# Extract DEGs from datalist
gene_extract <- lapply(datalist, '[[', 1) # Extract Gene column from datalist
gene_extract <- unique(unlist(gene_extract, recursive = FALSE)) # Concatenate multiple lists to one

# Create count matrix
annGSEA <- data.frame(row.names=gene_extract)

# Score which genes are present, or not.
for (j in names(datalist)){
  cat("j =",j,"\n") 
  for (i in 1:length(gene_extract)){
    cat("i =",i,"\n") 
    ifelse(gene_extract[i] %in% datalist[[j]]$Gene,
           annGSEA[i,j] <- 1,
           annGSEA[i,j] <- 0)
  }
}

# Only keep Genes expressed in at least 2 conditions/cell types
annGSEA_fil = annGSEA[apply(annGSEA, 1, function(x) sum(x > 0)/length(x) > 0.1), ,drop = FALSE]

# Add avg_logFC values from datalist, based on presence of gene
for (j in colnames(annGSEA_fil)) {
  cat("j =",j,"\n") 
  for (i in row.names(annGSEA_fil)) {
    cat("i =",i,"\n") 
    if (annGSEA_fil[i,j] == 1){
      annGSEA_fil[i,j] <- lapply(datalist[j], function(x) x$avg_logFC[x$Gene==i])
    } else {
      annGSEA_fil[i,j] <- 0
    }
  }
}

# Add avg_logFC values from datalist, based on presence of gene to calculate top genes
for (j in colnames(annGSEA)) {
  cat("j =",j,"\n") 
  for (i in row.names(annGSEA)) {
    cat("i =",i,"\n") 
    if (annGSEA[i,j] == 1){
      annGSEA[i,j] <- lapply(datalist[j], function(x) x$avg_logFC[x$Gene==i])
    } else {
      annGSEA[i,j] <- 0
    }
  }
}

# Annotate top 5 up and down-regulated genes
no <- 5
topGeneAnno <- unique(unlist(lapply(colnames(annGSEA), function(x){as.vector(rbind(row.names(annGSEA)[order(annGSEA[,x], decreasing = TRUE)][1:no], 
                                                                                   row.names(annGSEA)[order(annGSEA[,x], decreasing = FALSE)][1:no]))})))

# Set text and figure dimensions
geneLab=20
termLab=8

# colnames(annGSEA) <- gsub("([:punct:]+)|Mock", " ", colnames(annGSEA)) # change colnames for heatmap

# Annotate heatmap
ha_column <- HeatmapAnnotation('Status' = unlist(lapply(strsplit(colnames(annGSEA_fil), "_"), "[[", 4)),
                               'Cell type' = unlist(lapply(strsplit(colnames(annGSEA_fil), "_"), "[[", 1)),
                               'Condition' = unlist(lapply(strsplit(colnames(annGSEA_fil), "_"), "[[", 3)),
                               col = list('Cell type' = celltypePalette,
                                          'Condition' = conditionPalette,
                                          'Status' = statusPalette),
                               annotation_legend_param = colnames(annGSEA_fil))
?enrichGO
# Generate heatmap
hmapGSEA <- Heatmap(annGSEA_fil, name = "Average logFC",
                    col=colorRamp2(c(-2, 0, 2), c("darkblue", "white", "red")),
                    #rect_gp=gpar(col="grey85"),
                    cluster_rows=T,
                    show_row_dend=T,
                    row_title="Differential expressed genes",
                    row_title_side="left",
                    row_title_gp=gpar(fontsize=12, fontface="bold"),
                    row_title_rot=90,
                    show_row_names=FALSE,
                    row_names_gp=gpar(fontsize=5, fontface="bold"),
                    row_names_side="left",
                    row_names_max_width=unit(5, "cm"),
                    row_dend_width=unit(10,"mm"),
                    
                    cluster_columns=T,
                    show_column_dend=T,
                    column_title="DEG among different cell types and conditions",
                    column_title_side="top",
                    column_title_gp=gpar(fontsize=12, fontface="bold"),
                    column_title_rot=0,
                    show_column_names=FALSE,
                    column_names_gp=gpar(fontsize=termLab, fontface="bold"),
                    column_names_max_height=unit(15, "cm"),
                    show_heatmap_legend=TRUE,
                    
                    width=unit(20, "cm"),
                    
                    clustering_distance_columns="euclidean",
                    clustering_method_columns="complete",
                    clustering_distance_rows="euclidean",
                    clustering_method_rows="complete",
                    
                    bottom_annotation = ha_column)

# Add row annotations
htmap = hmapGSEA + rowAnnotation(up = anno_mark(at = which(rownames(annGSEA_fil) %in% topGeneAnno), 
                                                labels = rownames(annGSEA_fil)[rownames(annGSEA_fil) %in% topGeneAnno], 
                                                labels_gp = gpar(fontsize=7)))

# Save heatmap                          

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/heatmap_camel.pdf",
    width=12, height=9)
draw(htmap)
dev.off()

##################################################################################
#Llama Heatmap
# Create color palettes 
conditionPalette <- structure(c("#F39C6B","#A5668B"), names = c("ACN4", "MERS"))
celltypePalette <- structure(c("#8A6E7E", "#40679B", "#089392", "#40B48B", "#9DCD84", "#EAE29C"), 
                             names = c("Cluster1", "Cluster2", "Secretory", "Ciliated", "Club", "Basal"))
statusPalette <- structure(c("#519872", "#CD4631"), names = c("Bystanders", "Infected"))


# Import DEG files
#Loading filtered and merged Seurat Object into Script
#Infected and Mock Samples merged and cell clusters assigned

# set wd
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_both")
# import files
filelist = list.files(pattern="*.csv", full.names = TRUE)
datalist <- lapply(filelist, FUN= read.csv, na.strings="")
names(datalist) <- gsub(".txt", "", basename(filelist))

# Define viral genes
virus_genes <- c("mers-3UTR","mers-orf8b","mers-N","mers-M","mers-E","mers-orf5","mers-orf4b","mers-orf4a","mers-orf3","mers-S","mers-orf1ab","mers-5UTR",
                 "camel229E-3UTR","camel229E-ORF8","camel229E-N","camel229E-M","camel229E-E","camel229E-ORF4","camel229E-S","camel229E-ORF1a", "camel229E-ORF1ab","camel229E-5UTR")

# Modify datalist
new_colnames <- c("Gene", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")
names(datalist) <- gsub(".csv", "", basename(filelist))
datalist <- lapply(datalist, setNames, new_colnames) # Rename column names in list
datalist <- lapply(datalist, function(x) x[!x$Gene %in% virus_genes, ]) # Remove viral genes from files



##########################################################################################################

# Figure 5C
# Extract DEGs from datalist
gene_extract <- lapply(datalist, '[[', 1) # Extract Gene column from datalist
gene_extract <- unique(unlist(gene_extract, recursive = FALSE)) # Concatenate multiple lists to one

# Create count matrix
annGSEA <- data.frame(row.names=gene_extract)

# Score which genes are present, or not.
for (j in names(datalist)){
  cat("j =",j,"\n") 
  for (i in 1:length(gene_extract)){
    cat("i =",i,"\n") 
    ifelse(gene_extract[i] %in% datalist[[j]]$Gene,
           annGSEA[i,j] <- 1,
           annGSEA[i,j] <- 0)
  }
}

# Only keep Genes expressed in at least 2 conditions/cell types
annGSEA_fil = annGSEA[apply(annGSEA, 1, function(x) sum(x > 0)/length(x) > 0.1), ,drop = FALSE]

# Add avg_logFC values from datalist, based on presence of gene
for (j in colnames(annGSEA_fil)) {
  cat("j =",j,"\n") 
  for (i in row.names(annGSEA_fil)) {
    cat("i =",i,"\n") 
    if (annGSEA_fil[i,j] == 1){
      annGSEA_fil[i,j] <- lapply(datalist[j], function(x) x$avg_logFC[x$Gene==i])
    } else {
      annGSEA_fil[i,j] <- 0
    }
  }
}

# Add avg_logFC values from datalist, based on presence of gene to calculate top genes
for (j in colnames(annGSEA)) {
  cat("j =",j,"\n") 
  for (i in row.names(annGSEA)) {
    cat("i =",i,"\n") 
    if (annGSEA[i,j] == 1){
      annGSEA[i,j] <- lapply(datalist[j], function(x) x$avg_logFC[x$Gene==i])
    } else {
      annGSEA[i,j] <- 0
    }
  }
}

# Annotate top 5 up and down-regulated genes
no <- 5
topGeneAnno <- unique(unlist(lapply(colnames(annGSEA), function(x){as.vector(rbind(row.names(annGSEA)[order(annGSEA[,x], decreasing = TRUE)][1:no], 
                                                                                   row.names(annGSEA)[order(annGSEA[,x], decreasing = FALSE)][1:no]))})))

# Set text and figure dimensions
geneLab=20
termLab=8

# colnames(annGSEA) <- gsub("([:punct:]+)|Mock", " ", colnames(annGSEA)) # change colnames for heatmap

# Annotate heatmap
ha_column <- HeatmapAnnotation('Status' = unlist(lapply(strsplit(colnames(annGSEA_fil), "_"), "[[", 4)),
                               'Cell type' = unlist(lapply(strsplit(colnames(annGSEA_fil), "_"), "[[", 1)),
                               'Condition' = unlist(lapply(strsplit(colnames(annGSEA_fil), "_"), "[[", 3)),
                               col = list('Cell type' = celltypePalette,
                                          'Condition' = conditionPalette,
                                          'Status' = statusPalette),
                               annotation_legend_param = colnames(annGSEA_fil))

# Generate heatmap
hmapGSEA <- Heatmap(annGSEA_fil, name = "Average logFC",
                    col=colorRamp2(c(-2, 0, 2), c("darkblue", "white", "red")),
                    #rect_gp=gpar(col="grey85"),
                    cluster_rows=T,
                    show_row_dend=T,
                    row_title="Differential expressed genes",
                    row_title_side="left",
                    row_title_gp=gpar(fontsize=12, fontface="bold"),
                    row_title_rot=90,
                    show_row_names=FALSE,
                    row_names_gp=gpar(fontsize=5, fontface="bold"),
                    row_names_side="left",
                    row_names_max_width=unit(5, "cm"),
                    row_dend_width=unit(10,"mm"),
                    
                    cluster_columns=T,
                    show_column_dend=T,
                    column_title="DEG among different cell types and conditions",
                    column_title_side="top",
                    column_title_gp=gpar(fontsize=12, fontface="bold"),
                    column_title_rot=0,
                    show_column_names=FALSE,
                    column_names_gp=gpar(fontsize=termLab, fontface="bold"),
                    column_names_max_height=unit(15, "cm"),
                    show_heatmap_legend=TRUE,
                    
                    width=unit(20, "cm"),
                    
                    clustering_distance_columns="euclidean",
                    clustering_method_columns="complete",
                    clustering_distance_rows="euclidean",
                    clustering_method_rows="complete",
                    
                    bottom_annotation = ha_column)

# Add row annotations
htmap = hmapGSEA + rowAnnotation(up = anno_mark(at = which(rownames(annGSEA_fil) %in% topGeneAnno), 
                                                labels = rownames(annGSEA_fil)[rownames(annGSEA_fil) %in% topGeneAnno], 
                                                labels_gp = gpar(fontsize=7)))


# Save heatmap
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/heatmap_llama.pdf",
    width=12, height=14)
draw(htmap)
dev.off()

