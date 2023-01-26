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

##########################################################################################################################
#Loading filtered and merged Seurat Object into Script
#Infected and Mock Samples merged

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.ferus")

MERS.clustered <- LoadH5Seurat("MERS.combined.h5seurat")
MERS.clustered

str(MERS.clustered)

##########################################################################################################################
#Analysis of the combined dataset
# Run the standard workflow for visualization and clustering
MERS.clustered <- ScaleData(MERS.clustered, vars.to.regress = c("Percent.Viral", "S.Score", "G2M.Score"), do.par = TRUE, num.cores = 16, verbose = FALSE)
MERS.clustered <- RunPCA(MERS.clustered, npcs = 30, verbose = FALSE)
MERS.clustered <- FindNeighbors(MERS.clustered, reduction = "pca", dims = 1:5)
MERS.clustered <- FindClusters(MERS.clustered, resolution = 0.1)
MERS.clustered <- RunUMAP(MERS.clustered, reduction = "pca", dims = 1:5)
#MERS.clustered <- RunTSNE(MERS.clustered, reduction = "pca", dims = 1:5)


DimPlot(MERS.clustered, reduction = "umap")
DimPlot(MERS.clustered, reduction = "umap", group.by = "Status")
DimPlot(MERS.clustered, reduction = "umap", split.by = "Treat")

##################################################################################################################
# Markers to find cell type clusters
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
MERS.clustered_markers <- FindAllMarkers(MERS.clustered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MERS.clustered_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


Camel_MERS_Markers_clean <- MERS.clustered_markers %>% select(7,)

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Ferus_camel_MERS")
write.table(Camel_MERS_Markers_clean,"Camel_Ferus_MERS_Markers_clean.txt", row.names = FALSE)


Camel_MERS_Markers_clean <- MERS.clustered_markers %>% select(6:7,)

Camel_MERS_Markers_clean_Cluster_0 <- Camel_MERS_Markers_clean[Camel_MERS_Markers_clean$cluster == 0, ]
Camel_MERS_Markers_clean_Cluster_1 <- Camel_MERS_Markers_clean[Camel_MERS_Markers_clean$cluster == 1, ]
Camel_MERS_Markers_clean_Cluster_2 <- Camel_MERS_Markers_clean[Camel_MERS_Markers_clean$cluster == 2, ]
Camel_MERS_Markers_clean_Cluster_3 <- Camel_MERS_Markers_clean[Camel_MERS_Markers_clean$cluster == 3, ]
Camel_MERS_Markers_clean_Cluster_4 <- Camel_MERS_Markers_clean[Camel_MERS_Markers_clean$cluster == 4, ]
Camel_MERS_Markers_clean_Cluster_5 <- Camel_MERS_Markers_clean[Camel_MERS_Markers_clean$cluster == 5, ]

Camel_MERS_Markers_clean_Cluster_0 <- Camel_MERS_Markers_clean_Cluster_0 %>% select(2,)
Camel_MERS_Markers_clean_Cluster_1 <- Camel_MERS_Markers_clean_Cluster_1 %>% select(2,)
Camel_MERS_Markers_clean_Cluster_2 <- Camel_MERS_Markers_clean_Cluster_2 %>% select(2,)
Camel_MERS_Markers_clean_Cluster_3 <- Camel_MERS_Markers_clean_Cluster_3 %>% select(2,)
Camel_MERS_Markers_clean_Cluster_4 <- Camel_MERS_Markers_clean_Cluster_4 %>% select(2,)
Camel_MERS_Markers_clean_Cluster_5 <- Camel_MERS_Markers_clean_Cluster_5 %>% select(2,)


setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Ferus_camel_MERS")
write.table(Camel_MERS_Markers_clean_Cluster_0,"Ferus_MERS_Cluster_0.txt", row.names = FALSE)
write.table(Camel_MERS_Markers_clean_Cluster_1,"Ferus_MERS_Cluster_1.txt", row.names = FALSE)
write.table(Camel_MERS_Markers_clean_Cluster_2,"Ferus_MERS_Cluster_2.txt", row.names = FALSE)
write.table(Camel_MERS_Markers_clean_Cluster_3,"Ferus_MERS_Cluster_3.txt", row.names = FALSE)
write.table(Camel_MERS_Markers_clean_Cluster_4,"Ferus_MERS_Cluster_4.txt", row.names = FALSE)
write.table(Camel_MERS_Markers_clean_Cluster_5,"Ferus_MERS_Cluster_5.txt", row.names = FALSE)




#############################################################################################################################################

#Cluster0
FeaturePlot(MERS.clustered, features = c("KRT5", "LOC102514448"))
FeaturePlot(MERS.clustered, features = "KRT5")
FeaturePlot(MERS.clustered, features = "LOC102514448")

#Cluster1
FeaturePlot(MERS.clustered, features = c("SCGB1A1", "LOC102509380"))
FeaturePlot(MERS.clustered, features = "SCGB1A1") #secretory
FeaturePlot(MERS.clustered, features = "LOC102509380") #secretory

#Cluster2
FeaturePlot(MERS.clustered, features = c("SPP1", "LGALS3"))
FeaturePlot(MERS.clustered, features = "SPP1")
FeaturePlot(MERS.clustered, features = "LGALS3")

#Cluster3
FeaturePlot(MERS.clustered, features = c("LOC102523613", "WFDC2"))
FeaturePlot(MERS.clustered, features = "LOC102523613")
FeaturePlot(MERS.clustered, features = "WFDC2")

#Cluster4
FeaturePlot(MERS.clustered, features = c("LOC102519239", "LOC102507471"))
FeaturePlot(MERS.clustered, features = "LOC102519239")
FeaturePlot(MERS.clustered, features = "LOC102507471")

#Cluster5
FeaturePlot(MERS.clustered, features = c("TPPP3", "TUBA1A"))
FeaturePlot(MERS.clustered, features = "TPPP3") #ciliated1
FeaturePlot(MERS.clustered, features = "TUBA1A")

#Cluster6
FeaturePlot(MERS.clustered, features = c("TOP2A", "CENPF"))
FeaturePlot(MERS.clustered, features = "TOP2A")
FeaturePlot(MERS.clustered, features = "CENPF")
FeaturePlot()

FeaturePlot(MERS.clustered, features = "DPP4")
FeaturePlot(MERS.clustered, features = "3UTR") 
FeaturePlot(MERS.clustered, features = "KRT4") #deuterosomal
FeaturePlot(MERS.clustered, features = "DEUP1") #deuterosomal - interesting

##################################################################################################################
#Renaming the found clusters to the respective cell types for downstream analysis differential expression analysis
new.cluster.ids <- c("Basal", "Club", "Club", "Basal", "Ciliated", "Deuterosomal")
names(new.cluster.ids) <- levels(MERS.clustered)
MERS.clustered <- RenameIdents(MERS.clustered, new.cluster.ids)
DimPlot(MERS.clustered, reduction = "umap", label.size = 6, repel = TRUE, label = TRUE, pt.size = 0.5) + NoLegend()

#Saving the seurat object to disk to access in other script for analysis
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.ferus")
SaveH5Seurat(MERS.clustered, "MERS.clustered", overwrite = TRUE)



