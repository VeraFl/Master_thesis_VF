# UMAP and Clustering of Single-cell Seq data
# 21.06.2022

#load libraries
library(dplyr)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(patchwork)
library(readr)
library(gdata)
library(Matrix)
library(writexl)
library(grDevices)

##########################################################################################################################
#Loading filtered and merged Seurat Object into Script
#Infected and Mock Samples merged

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/V.pacos")

camel229E.clustered <- LoadH5Seurat("camel229E.combined.h5seurat")
camel229E.clustered

##########################################################################################################################
#Analysis of the combined dataset
#cell cycle state genes from the Seurat package
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# 
# camel229E.clustered <- CellCycleScoring(camel229E.clustered, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


# Run the standard workflow for visualization and clustering
camel229E.clustered <- ScaleData(camel229E.clustered, vars.to.regress = c("Percent.Viral", "S.Score", "G2M.Score"), do.par = TRUE, num.cores = 16, verbose = FALSE)
camel229E.clustered <- RunPCA(camel229E.clustered, npcs = 30, verbose = FALSE)
camel229E.clustered <- FindNeighbors(camel229E.clustered, reduction = "pca", dims = 1:5)
camel229E.clustered <- FindClusters(camel229E.clustered, resolution = 0.2)
camel229E.clustered <- RunUMAP(camel229E.clustered, reduction = "pca", dims = 1:5)
# camel229E.clustered <- RunTSNE(camel229E.clustered, reduction = "pca", dims = 1:5)

DimPlot(camel229E.clustered, reduction = "umap", label.size = 6, repel = TRUE, label = TRUE, pt.size = 0.5)
# DimPlot(camel229E.clustered, reduction = "tsne")
DimPlot(camel229E.clustered, reduction = "umap", group.by = "Status")
DimPlot(camel229E.clustered, reduction = "umap", split.by = "Treat")

##################################################################################################################
# Markers to find cell type clusters
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
camel229E.clustered_markers <- FindAllMarkers(camel229E.clustered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
camel229E.clustered_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

Alpaca_camel229E_Markers_clean <- camel229E.clustered_markers %>% select(6:7,)

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_camel229E")
write.table(Alpaca_camel229E_Markers_clean,"Alpaca_camel229E_Markers_clean.txt", row.names = FALSE)
Alpaca_camel229E_Markers_clean_Cluster_0 <- Alpaca_camel229E_Markers_clean[Alpaca_camel229E_Markers_clean$cluster == 0, ]
Alpaca_camel229E_Markers_clean_Cluster_1 <- Alpaca_camel229E_Markers_clean[Alpaca_camel229E_Markers_clean$cluster == 1, ]
Alpaca_camel229E_Markers_clean_Cluster_2 <- Alpaca_camel229E_Markers_clean[Alpaca_camel229E_Markers_clean$cluster == 2, ]
Alpaca_camel229E_Markers_clean_Cluster_3 <- Alpaca_camel229E_Markers_clean[Alpaca_camel229E_Markers_clean$cluster == 3, ]
Alpaca_camel229E_Markers_clean_Cluster_4 <- Alpaca_camel229E_Markers_clean[Alpaca_camel229E_Markers_clean$cluster == 4, ]
Alpaca_camel229E_Markers_clean_Cluster_5 <- Alpaca_camel229E_Markers_clean[Alpaca_camel229E_Markers_clean$cluster == 5, ]
Alpaca_camel229E_Markers_clean_Cluster_6 <- Alpaca_camel229E_Markers_clean[Alpaca_camel229E_Markers_clean$cluster == 6, ]
# Alpaca_camel229E_Markers_clean_Cluster_7 <- Alpaca_camel229E_Markers_clean[Alpaca_camel229E_Markers_clean$cluster == 7, ]


Alpaca_camel229E_Markers_clean_Cluster_0 <- Alpaca_camel229E_Markers_clean_Cluster_0 %>% select(2,)
Alpaca_camel229E_Markers_clean_Cluster_1 <- Alpaca_camel229E_Markers_clean_Cluster_1 %>% select(2,)
Alpaca_camel229E_Markers_clean_Cluster_2 <- Alpaca_camel229E_Markers_clean_Cluster_2 %>% select(2,)
Alpaca_camel229E_Markers_clean_Cluster_3 <- Alpaca_camel229E_Markers_clean_Cluster_3 %>% select(2,)
Alpaca_camel229E_Markers_clean_Cluster_4 <- Alpaca_camel229E_Markers_clean_Cluster_4 %>% select(2,)
Alpaca_camel229E_Markers_clean_Cluster_5 <- Alpaca_camel229E_Markers_clean_Cluster_5 %>% select(2,)
Alpaca_camel229E_Markers_clean_Cluster_6 <- Alpaca_camel229E_Markers_clean_Cluster_6 %>% select(2,)
# Alpaca_camel229E_Markers_clean_Cluster_7 <- Alpaca_camel229E_Markers_clean_Cluster_7 %>% select(2,)





setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_camel229E")
write.table(camel229E.clustered_markers,"Alpaca_camel229E_Cluster_markers.xlsx")
write.table(Alpaca_camel229E_Markers_clean_Cluster_0,"Alpaca_camel229E_Cluster_0.txt", row.names = FALSE)
write.table(Alpaca_camel229E_Markers_clean_Cluster_1,"Alpaca_camel229E_Cluster_1.txt", row.names = FALSE)
write.table(Alpaca_camel229E_Markers_clean_Cluster_2,"Alpaca_camel229E_Cluster_2.txt", row.names = FALSE)
write.table(Alpaca_camel229E_Markers_clean_Cluster_3,"Alpaca_camel229E_Cluster_3.txt", row.names = FALSE)
write.table(Alpaca_camel229E_Markers_clean_Cluster_4,"Alpaca_camel229E_Cluster_4.txt", row.names = FALSE)
write.table(Alpaca_camel229E_Markers_clean_Cluster_5,"Alpaca_camel229E_Cluster_5.txt", row.names = FALSE)
write.table(Alpaca_camel229E_Markers_clean_Cluster_6,"Alpaca_camel229E_Cluster_6.txt", row.names = FALSE)
# write.table(Alpaca_camel229E_Markers_clean_Cluster_7,"Alpaca_camel229E_Cluster_7.txt", row.names = FALSE)


#Plotting expression of different markers (making the clusters)

####################################################################################################################
#Plotting expression of different markers (making the clusters)

#BASAL: KRT5, TP63, DLK2, KRT6A
#Suprabasal: KRT4, KRT13, KRT16, KRT23
#Club: SCGB1A1, KRT7, KRT19
#Deuterosomal: PLK4, CCNO, CEP78, DEUP1
#Secretory: MUC5B, MUC5AC
#Goblet: FOXQ1, SPDEF, MUC5AC
#Ciliated: FOXJ1, TPPP3, SNTN, TUBA1A, CETN2, SPEF2, PIFO, LRRC6
#Tuft: POUF2f3, TRPM5, LRMP, RGS13, HOXC5, HMX2, ANXA4
#PNEC: PSMD5, NGF, PCSK1N, SCGN, NEB, HOXB1, ASCL1/2, FOXA2
#Ionocytes: ASCL3, CFTR, FOXI1, DMRT2, V-ATPase


###########################################################################################################################
#Cluster0
FeaturePlot(camel229E.clustered, features = c("LOC116662762", "LOC102523613")) 
FeaturePlot(camel229E.clustered, features = "LOC116662762") 
FeaturePlot(camel229E.clustered, features = "LOC102523613") 

#Cluster1
FeaturePlot(camel229E.clustered, features = c("MUC21", "TIMP3"))
FeaturePlot(camel229E.clustered, features = "MUC5B")
FeaturePlot(camel229E.clustered, features = "TIMP3")

#Cluster2
FeaturePlot(camel229E.clustered, features = c("LOXL4", "FN1"))
FeaturePlot(camel229E.clustered, features = "LOXL4")
FeaturePlot(camel229E.clustered, features = "FN1")

#Cluster3
FeaturePlot(camel229E.clustered, features = c("TPPP3", "SEC14L3"))
FeaturePlot(camel229E.clustered, features = "TPPP3")
FeaturePlot(camel229E.clustered, features = "SEC14L3")

#Cluster4
FeaturePlot(camel229E.clustered, features = c("BPIFA1", "LOC102519002"))
FeaturePlot(camel229E.clustered, features = "BPIFA1")
FeaturePlot(camel229E.clustered, features = "LOC102519002")

#Cluster5
FeaturePlot(camel229E.clustered, features = c("CCDC17", "LOC102524239"))
FeaturePlot(camel229E.clustered, features = "CCDC17")
FeaturePlot(camel229E.clustered, features = "LOC102524239")



##################################################################################################################
#Renaming the found clusters to the respective cell types for downstream analysis differential expression analysis
# new.cluster.ids <- c("Secretory", "Ciliated", "Suprabasal", "Basal", "NA", "Deuterosomal")
# names(new.cluster.ids) <- levels(camel229E.clustered)
# camel229E.clustered <- RenameIdents(camel229E.clustered, new.cluster.ids)
# DimPlot(camel229E.clustered, reduction = "tsne", label.size = 6, repel = TRUE, label = TRUE, pt.size = 0.5) + NoLegend()


#Saving the seurat object to disk to access in other script for analysis
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/V.pacos")
SaveH5Seurat(camel229E.clustered, "camel229E.clustered", overwrite = TRUE)

