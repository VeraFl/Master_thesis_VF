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

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.dromedarius")

MERS.clustered <- LoadH5Seurat("MERS.combined.h5seurat")
MERS.clustered

##########################################################################################################################
#Analysis of the combined dataset
# Run the standard workflow for visualization and clustering
MERS.clustered <- ScaleData(MERS.clustered, vars.to.regress = c("Percent.Viral", "S.Score", "G2M.Score"), do.par = TRUE, num.cores = 16, verbose = FALSE)
MERS.clustered <- RunPCA(MERS.clustered, npcs = 30, verbose = FALSE)
MERS.clustered <- FindNeighbors(MERS.clustered, reduction = "pca", dims = 1:5)
MERS.clustered <- FindClusters(MERS.clustered, resolution = 0.2)
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

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Dromedary_camel_MERS")
write.table(Camel_MERS_Markers_clean,"Camel_Dromedary_MERS_Markers_clean.txt", row.names = FALSE)
#Plotting expression of different markers (making the clusters)

Camel_MERS_Markers_clean <- MERS.clustered_markers %>% select(6:7,)

Camel_MERS_Markers_clean_Cluster_0 <- Camel_MERS_Markers_clean[Camel_MERS_Markers_clean$cluster == 0, ]
Camel_MERS_Markers_clean_Cluster_1 <- Camel_MERS_Markers_clean[Camel_MERS_Markers_clean$cluster == 1, ]
Camel_MERS_Markers_clean_Cluster_2 <- Camel_MERS_Markers_clean[Camel_MERS_Markers_clean$cluster == 2, ]
Camel_MERS_Markers_clean_Cluster_3 <- Camel_MERS_Markers_clean[Camel_MERS_Markers_clean$cluster == 3, ]
Camel_MERS_Markers_clean_Cluster_4 <- Camel_MERS_Markers_clean[Camel_MERS_Markers_clean$cluster == 4, ]
Camel_MERS_Markers_clean_Cluster_5 <- Camel_MERS_Markers_clean[Camel_MERS_Markers_clean$cluster == 5, ]
Camel_MERS_Markers_clean_Cluster_6 <- Camel_MERS_Markers_clean[Camel_MERS_Markers_clean$cluster == 6, ]

Camel_MERS_Markers_clean_Cluster_0 <- Camel_MERS_Markers_clean_Cluster_0 %>% select(2,)
Camel_MERS_Markers_clean_Cluster_1 <- Camel_MERS_Markers_clean_Cluster_1 %>% select(2,)
Camel_MERS_Markers_clean_Cluster_2 <- Camel_MERS_Markers_clean_Cluster_2 %>% select(2,)
Camel_MERS_Markers_clean_Cluster_3 <- Camel_MERS_Markers_clean_Cluster_3 %>% select(2,)
Camel_MERS_Markers_clean_Cluster_4 <- Camel_MERS_Markers_clean_Cluster_4 %>% select(2,)
Camel_MERS_Markers_clean_Cluster_5 <- Camel_MERS_Markers_clean_Cluster_5 %>% select(2,)
Camel_MERS_Markers_clean_Cluster_6 <- Camel_MERS_Markers_clean_Cluster_6 %>% select(2,)


setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Dromedary_camel_MERS")
write.table(Camel_MERS_Markers_clean_Cluster_0,"Camel_MERS_Markers_clean_Cluster_0.txt", row.names = FALSE)
write.table(Camel_MERS_Markers_clean_Cluster_1,"Camel_MERS_Markers_clean_Cluster_1.txt", row.names = FALSE)
write.table(Camel_MERS_Markers_clean_Cluster_2,"Camel_MERS_Markers_clean_Cluster_2.txt", row.names = FALSE)
write.table(Camel_MERS_Markers_clean_Cluster_3,"Camel_MERS_Markers_clean_Cluster_3.txt", row.names = FALSE)
write.table(Camel_MERS_Markers_clean_Cluster_4,"Camel_MERS_Markers_clean_Cluster_4.txt", row.names = FALSE)
write.table(Camel_MERS_Markers_clean_Cluster_5,"Camel_MERS_Markers_clean_Cluster_5.txt", row.names = FALSE)
write.table(Camel_MERS_Markers_clean_Cluster_6,"Camel_MERS_Markers_clean_Cluster_6.txt", row.names = FALSE)


#############################################################################################################################################

#Cluster0
FeaturePlot(MERS.clustered, features = c("LOC105094590", "LOC105102682"))
FeaturePlot(MERS.clustered, features = "LOC105094590")
FeaturePlot(MERS.clustered, features = "LOC105102682")

#Cluster1
FeaturePlot(MERS.clustered, features = c("SCGB1A1", "IVL"))
FeaturePlot(MERS.clustered, features = "SCGB1A1") #secretory
FeaturePlot(MERS.clustered, features = "IVL") #secretory

#Cluster2
FeaturePlot(MERS.clustered, features = "SPP1")

#Cluster3
FeaturePlot(MERS.clustered, features = c("LOC105089543", "WFDC2"))
FeaturePlot(MERS.clustered, features = "LOC105089543")
FeaturePlot(MERS.clustered, features = "WFDC2")

#Cluster4
FeaturePlot(MERS.clustered, features = c("LOC105088388", "LOC105093359"))
FeaturePlot(MERS.clustered, features = "LOC105088388")
FeaturePlot(MERS.clustered, features = "LOC105093359")

#Cluster5
FeaturePlot(MERS.clustered, features = "TPPP3") #ciliated1

#Cluster6
FeaturePlot(MERS.clustered, features = c("TOP2A", "CENPF"))
FeaturePlot(MERS.clustered, features = "TOP2A")
FeaturePlot(MERS.clustered, features = "CENPF")
FeaturePlot()

FeaturePlot(MERS.clustered, features = "PLK4") #deuterosomal - interesting
FeaturePlot(MERS.clustered, features = "CCNO") #deuterosomal
FeaturePlot(MERS.clustered, features = "CEP78") #deuterosomal
FeaturePlot(MERS.clustered, features = "DEUP1") #deuterosomal - interesting

##################################################################################################################
#Renaming the found clusters to the respective cell types for downstream analysis differential expression analysis
# new.cluster.ids <- c("Secretory", "Ciliated", "Suprabasal", "Basal", "NA", "Deuterosomal")
# names(new.cluster.ids) <- levels(MERS.clustered)
# MERS.clustered <- RenameIdents(MERS.clustered, new.cluster.ids)
# DimPlot(MERS.clustered, reduction = "tsne", label.size = 6, repel = TRUE, label = TRUE, pt.size = 0.5) + NoLegend()


#Saving the seurat object to disk to access in other script for analysis
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/V.pacos")
SaveH5Seurat(MERS.clustered, "MERS.clustered", overwrite = TRUE)
