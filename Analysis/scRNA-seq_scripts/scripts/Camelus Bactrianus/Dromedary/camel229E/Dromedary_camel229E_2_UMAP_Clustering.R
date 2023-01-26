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

camel229E.clustered <- LoadH5Seurat("camel229E.combined.h5seurat")
camel229E.clustered

##########################################################################################################################
#Analysis of the combined dataset
# Run the standard workflow for visualization and clustering
camel229E.clustered <- ScaleData(camel229E.clustered, vars.to.regress = c("Percent.Viral", "S.Score", "G2M.Score"), do.par = TRUE, num.cores = 16, verbose = FALSE)
camel229E.clustered <- RunPCA(camel229E.clustered, npcs = 30, verbose = FALSE)
camel229E.clustered <- FindNeighbors(camel229E.clustered, reduction = "pca", dims = 1:5)
camel229E.clustered <- FindClusters(camel229E.clustered, resolution = 0.2)
camel229E.clustered <- RunUMAP(camel229E.clustered, reduction = "pca", dims = 1:5)
#camel229E.clustered <- RunTSNE(camel229E.clustered, reduction = "pca", dims = 1:5)


DimPlot(camel229E.clustered, reduction = "umap")
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


Camel_camel229E_Markers_clean <- camel229E.clustered_markers %>% select(7,)

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Dromedary_camel_camel229E")
write.table(Camel_camel229E_Markers_clean,"Camel_Dromedary_camel229E_Markers_clean.txt", row.names = FALSE)

Camel_camel229E_Markers_clean <- camel229E.clustered_markers %>% select(6:7,)

Camel_camel229E_Markers_clean_Cluster_0 <- Camel_camel229E_Markers_clean[Camel_camel229E_Markers_clean$cluster == 0, ]
Camel_camel229E_Markers_clean_Cluster_1 <- Camel_camel229E_Markers_clean[Camel_camel229E_Markers_clean$cluster == 1, ]
Camel_camel229E_Markers_clean_Cluster_2 <- Camel_camel229E_Markers_clean[Camel_camel229E_Markers_clean$cluster == 2, ]
Camel_camel229E_Markers_clean_Cluster_3 <- Camel_camel229E_Markers_clean[Camel_camel229E_Markers_clean$cluster == 3, ]
Camel_camel229E_Markers_clean_Cluster_4 <- Camel_camel229E_Markers_clean[Camel_camel229E_Markers_clean$cluster == 4, ]
Camel_camel229E_Markers_clean_Cluster_5 <- Camel_camel229E_Markers_clean[Camel_camel229E_Markers_clean$cluster == 5, ]

Camel_camel229E_Markers_clean_Cluster_0 <- Camel_camel229E_Markers_clean_Cluster_0 %>% select(2,)
Camel_camel229E_Markers_clean_Cluster_1 <- Camel_camel229E_Markers_clean_Cluster_1 %>% select(2,)
Camel_camel229E_Markers_clean_Cluster_2 <- Camel_camel229E_Markers_clean_Cluster_2 %>% select(2,)
Camel_camel229E_Markers_clean_Cluster_3 <- Camel_camel229E_Markers_clean_Cluster_3 %>% select(2,)
Camel_camel229E_Markers_clean_Cluster_4 <- Camel_camel229E_Markers_clean_Cluster_4 %>% select(2,)
Camel_camel229E_Markers_clean_Cluster_5 <- Camel_camel229E_Markers_clean_Cluster_5 %>% select(2,)


setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Dromedary_camel_camel229E")
write.table(Camel_camel229E_Markers_clean_Cluster_0,"Camel_camel229E_Markers_clean_Cluster_0.txt", row.names = FALSE)
write.table(Camel_camel229E_Markers_clean_Cluster_1,"Camel_camel229E_Markers_clean_Cluster_1.txt", row.names = FALSE)
write.table(Camel_camel229E_Markers_clean_Cluster_2,"Camel_camel229E_Markers_clean_Cluster_2.txt", row.names = FALSE)
write.table(Camel_camel229E_Markers_clean_Cluster_3,"Camel_camel229E_Markers_clean_Cluster_3.txt", row.names = FALSE)
write.table(Camel_camel229E_Markers_clean_Cluster_4,"Camel_camel229E_Markers_clean_Cluster_4.txt", row.names = FALSE)
write.table(Camel_camel229E_Markers_clean_Cluster_5,"Camel_camel229E_Markers_clean_Cluster_5.txt", row.names = FALSE)

#Cluster5 = ciliated

#############################################################################################################################################

#Cluster0
FeaturePlot(camel229E.clustered, features = c("LOC105094590", "LOC105102682"))
FeaturePlot(camel229E.clustered, features = "LOC105094590")
FeaturePlot(camel229E.clustered, features = "LOC105102682")

#Cluster1
FeaturePlot(camel229E.clustered, features = c("SCGB1A1", "IVL"))
FeaturePlot(camel229E.clustered, features = "SCGB1A1") #secretory
FeaturePlot(camel229E.clustered, features = "IVL") #secretory

#Cluster2
FeaturePlot(camel229E.clustered, features = "SPP1")

#Cluster3
FeaturePlot(camel229E.clustered, features = c("LOC105089543", "WFDC2"))
FeaturePlot(camel229E.clustered, features = "LOC105089543")
FeaturePlot(camel229E.clustered, features = "WFDC2")

#Cluster4
FeaturePlot(camel229E.clustered, features = c("LOC105088388", "LOC105093359"))
FeaturePlot(camel229E.clustered, features = "LOC105088388")
FeaturePlot(camel229E.clustered, features = "LOC105093359")

#Cluster5
FeaturePlot(camel229E.clustered, features = "TPPP3") #ciliated1

#Cluster6
FeaturePlot(camel229E.clustered, features = c("TOP2A", "CENPF"))
FeaturePlot(camel229E.clustered, features = "TOP2A")
FeaturePlot(camel229E.clustered, features = "CENPF")
FeaturePlot()

FeaturePlot(camel229E.clustered, features = "PLK4") #deuterosomal - interesting
FeaturePlot(camel229E.clustered, features = "CCNO") #deuterosomal
FeaturePlot(camel229E.clustered, features = "CEP78") #deuterosomal
FeaturePlot(camel229E.clustered, features = "DEUP1") #deuterosomal - interesting

##################################################################################################################
#Renaming the found clusters to the respective cell types for downstream analysis differential expression analysis
# new.cluster.ids <- c("Secretory", "Ciliated", "Suprabasal", "Basal", "NA", "Deuterosomal")
# names(new.cluster.ids) <- levels(camel229E.clustered)
# camel229E.clustered <- RenameIdents(camel229E.clustered, new.cluster.ids)
# DimPlot(camel229E.clustered, reduction = "tsne", label.size = 6, repel = TRUE, label = TRUE, pt.size = 0.5) + NoLegend()

