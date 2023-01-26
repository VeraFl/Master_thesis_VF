# UMAP and Clustering of Single-cell Seq data
# 20.06.2022

#load libraries
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

##########################################################################################################################
#Setting the color schemes for the script:

#Treatment colors (HCoV-229E, MERS-CoV, Mock)
cols_treat = c("#84A59D","#A5668B", "#96ADC8")
show_col(cols_treat)
#Status colors
cols_stat = c("#CD4631","#DEA47E")
show_col(cols_stat)
#Celltype colors
cols_cell = hcl.colors(7, "Temps")
show_col(cols_cell)

##########################################################################################################################
#Loading filtered and merged Seurat Object into Script
#Infected and Mock Samples merged

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens")

HCoV.clustered <- LoadH5Seurat("HCoV.combined.h5seurat")
HCoV.clustered

##########################################################################################################################
#Analysis of the combined dataset

# Run the standard workflow for visualization and clustering
HCoV.clustered <- RunPCA(HCoV.clustered, npcs = 30, verbose = FALSE)
HCoV.clustered <- FindNeighbors(HCoV.clustered, reduction = "pca", dims = 1:20)
HCoV.clustered <- FindClusters(HCoV.clustered, resolution = 0.3)
HCoV.clustered <- RunUMAP(HCoV.clustered, reduction = "pca", dims = 1:20)

png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_HCoV/umap_cluster.png",
    width=500, height=500)
DimPlot(HCoV.clustered, reduction = "umap", label.size = 6, repel = TRUE, label = TRUE, pt.size = 0.5, cols = cols_cell) + NoLegend()
dev.off()


png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_HCoV/umap_split_status_treat.png",
    width=800, height=300)
DimPlot(HCoV.clustered, reduction = "umap", group.by = "Status", cols = cols_stat)
dev.off()

png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_HCoV/umap_split_status_treat.png",
    width=800, height=300)
DimPlot(HCoV.clustered, reduction = "umap", pt.size = 0.5, group.by = "Treat", cols = cols_treat)
dev.off()

DimPlot(HCoV.clustered, reduction = "umap", split.by = "Treat", pt.size = 0.5, cols = cols_cell) + NoAxes() + NoLegend()

##################################################################################################################
# Cluster Annotation
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
HCoV.clustered_markers <- FindAllMarkers(HCoV.clustered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
HCoV.clustered_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

Find
#Saving markers gene lists
Human_HCoV_Markers_clean <- HCoV.clustered_markers %>% select(7,)

HCoV.clustered_markers_0 <- HCoV.clustered_markers[HCoV.clustered_markers$cluster == 0, ]

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Human_HCoV")
write.table(Human_HCoV_Markers_clean,"Human_HCoV_Markers_clean.txt", row.names = FALSE)

Human_HCoV_Markers_clean <- HCoV.clustered_markers %>% select(6:7,)

Human_HCoV_Markers_clean_Cluster_0 <- Human_HCoV_Markers_clean[Human_HCoV_Markers_clean$cluster == 0, ]
Human_HCoV_Markers_clean_Cluster_1 <- Human_HCoV_Markers_clean[Human_HCoV_Markers_clean$cluster == 1, ]
Human_HCoV_Markers_clean_Cluster_2 <- Human_HCoV_Markers_clean[Human_HCoV_Markers_clean$cluster == 2, ]
Human_HCoV_Markers_clean_Cluster_3 <- Human_HCoV_Markers_clean[Human_HCoV_Markers_clean$cluster == 3, ]
Human_HCoV_Markers_clean_Cluster_4 <- Human_HCoV_Markers_clean[Human_HCoV_Markers_clean$cluster == 4, ]
Human_HCoV_Markers_clean_Cluster_5 <- Human_HCoV_Markers_clean[Human_HCoV_Markers_clean$cluster == 5, ]


Human_HCoV_Markers_clean_Cluster_0 <- Human_HCoV_Markers_clean_Cluster_0 %>% select(2,)
Human_HCoV_Markers_clean_Cluster_1 <- Human_HCoV_Markers_clean_Cluster_1 %>% select(2,)
Human_HCoV_Markers_clean_Cluster_2 <- Human_HCoV_Markers_clean_Cluster_2 %>% select(2,)
Human_HCoV_Markers_clean_Cluster_3 <- Human_HCoV_Markers_clean_Cluster_3 %>% select(2,)
Human_HCoV_Markers_clean_Cluster_4 <- Human_HCoV_Markers_clean_Cluster_4 %>% select(2,)
Human_HCoV_Markers_clean_Cluster_5 <- Human_HCoV_Markers_clean_Cluster_5 %>% select(2,)


setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Human_HCoV")
write.table(Human_HCoV_Markers_clean_Cluster_0,"Human_HCoV_Markers_clean_Cluster_0.txt", row.names = FALSE)
write.table(Human_HCoV_Markers_clean_Cluster_1,"Human_HCoV_Markers_clean_Cluster_1.txt", row.names = FALSE)
write.table(Human_HCoV_Markers_clean_Cluster_2,"Human_HCoV_Markers_clean_Cluster_2.txt", row.names = FALSE)
write.table(Human_HCoV_Markers_clean_Cluster_3,"Human_HCoV_Markers_clean_Cluster_3.txt", row.names = FALSE)
write.table(Human_HCoV_Markers_clean_Cluster_4,"Human_HCoV_Markers_clean_Cluster_4.txt", row.names = FALSE)
write.table(Human_HCoV_Markers_clean_Cluster_5,"Human_HCoV_Markers_clean_Cluster_5.txt", row.names = FALSE)

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

#viral
FeaturePlot(HCoV.clustered, features = "HCoV") #viral

#All
FeaturePlot(HCoV.clustered, features = c("FN1", "KRT16", "FOXI1", "SCGB1A1", "SPDEF", "SCGB3A1", "PLK4", "FOXN4", "FOXJ1","DEUP1"))
?FeaturePlot

#BASAL - INF
FeaturePlot(HCoV.clustered, features = "KRT5") #basal1 -interesting
FeaturePlot(HCoV.clustered, features = "TP63") #basal2
FeaturePlot(HCoV.clustered, features = "FN1") #basal3 -interesting
FeaturePlot(HCoV.clustered, features = "DLK2") #basal4
FeaturePlot(HCoV.clustered, features = "KRT6A") #basal5

#SUPRABASAL
FeaturePlot(HCoV.clustered, features = "KRT4") #suprabasal1
FeaturePlot(HCoV.clustered, features = "KRT13") #suprabasal2 - interesting
FeaturePlot(HCoV.clustered, features = "KRT16") #suprabasal3 - interesting
FeaturePlot(HCoV.clustered, features = "KRT23") #suprabasal4

#TUFT
FeaturePlot(HCoV.clustered, features = "POUF2f3") #tuft1
FeaturePlot(HCoV.clustered, features = "TRPM5") #tuft2
FeaturePlot(HCoV.clustered, features = "LRMP") #tuft3

#IONOCYTES
FeaturePlot(HCoV.clustered, features = "FOXI1") #ionocytes1
FeaturePlot(HCoV.clustered, features = "CFTR") #ionocytes2

##PNEC
FeaturePlot(HCoV.clustered, features = c("PCSK1N")) #PNEC1
FeaturePlot(HCoV.clustered, features = c("SCGN")) #PNEC2
FeaturePlot(HCoV.clustered, features = c("NEB")) #PNEC3
FeaturePlot(HCoV.clustered, features = c("HOXB1")) #PNEC4
FeaturePlot(HCoV.clustered, features = c("ASCL1")) #PNEC5
FeaturePlot(HCoV.clustered, features = c("ASCL2")) #PNEC6
FeaturePlot(HCoV.clustered, features = c("FOXA2")) #PNEC7
FeaturePlot(HCoV.clustered, features = c("PSMD5")) #PNEC8
FeaturePlot(HCoV.clustered, features = c("NGF")) #PNEC9

#CLUB
FeaturePlot(HCoV.clustered, features = c("SCGB1A1")) #club -interesting
FeaturePlot(HCoV.clustered, features = c("BPIFA1")) #club
FeaturePlot(HCoV.clustered, features = c("KRT7")) #club
FeaturePlot(HCoV.clustered, features = c("MSMB")) #club
FeaturePlot(HCoV.clustered, features = c("KRT19")) #club

#GOBLET/SECRETORY
FeaturePlot(HCoV.clustered, features = "MUC5B") #secretory1 - interesting
FeaturePlot(HCoV.clustered, features = "SCGB3A1") #secretory2 -interesting
FeaturePlot(HCoV.clustered, features = "MUC5AC") #secretory3 (goblet)
FeaturePlot(HCoV.clustered, features = "FOXQ1") #goblet1
FeaturePlot(HCoV.clustered, features = "SPDEF") #goblet2

#DEUTEROSOMAL
FeaturePlot(HCoV.clustered, features = "PLK4") #deuterosomal - interesting
FeaturePlot(HCoV.clustered, features = "CCNO") #deuterosomal
FeaturePlot(HCoV.clustered, features = "CEP78") #deuterosomal
FeaturePlot(HCoV.clustered, features = "DEUP1") #deuterosomal - interesting


#CILIATED
FeaturePlot(HCoV.clustered, features = "FOXN4") #preciliated
FeaturePlot(HCoV.clustered, features = "FOXJ1") #ciliated1 - interesting
FeaturePlot(HCoV.clustered, features = "PIFO") #ciliated2
FeaturePlot(HCoV.clustered, features = "TPPP3") #ciliated3
FeaturePlot(HCoV.clustered, features = "SPEF2") #ciliated4 -interesting
FeaturePlot(HCoV.clustered, features = "DNAH5") #ciliated5
FeaturePlot(HCoV.clustered, features = "LRRC6") #ciliated6


###########################################################################################################################
#Cluster0
FeaturePlot(HCoV.clustered, features = c("orf1ab", "S", "orf3", "orf4a", "orf4b", "orf5", "E", "M", "orf8b")) #viral genes
FeaturePlot(HCoV.clustered, features = c("PCNA", "TOP2A", "MCM6", "MKI67")) #cell cycle state genes
FeaturePlot(HCoV.clustered, features = c("MUC5B", "MUC5AC")) #secretory
FeaturePlot(HCoV.clustered, features = "MUC5B") #secretory1
FeaturePlot(HCoV.clustered, features = "SCGB3A1") #secretory2
FeaturePlot(HCoV.clustered, features = "HCoV") #secretory3 (goblet)
#Cluster1
FeaturePlot(HCoV.clustered, features = c("FOXJ1", "TPPP3", "SNTN", "TUBA1A", "CETN2")) #ciliated
FeaturePlot(HCoV.clustered, features = "FOXJ1") #ciliated1
FeaturePlot(HCoV.clustered, features = "PIFO") #ciliated2
FeaturePlot(HCoV.clustered, features = "TPPP3") #ciliated3

#Cluster2
FeaturePlot(HCoV.clustered, features = c("KRT5", "TP63", "DLK2", "KRT6A")) #basal
FeaturePlot(HCoV.clustered, features = "KRT5") #basal1
FeaturePlot(HCoV.clustered, features = "TP63") #basal2
FeaturePlot(HCoV.clustered, features = "FN1") #basal3
#Cluster3
FeaturePlot(HCoV.clustered, features = "HYDIN")


FeaturePlot(HCoV.clustered, features = c("PCSK1N", "SCGN", "NEB", "HOXB1", "ASCL1", "ASCL2", "FOXA2")) #PNEC
#Cluster5
FeaturePlot(HCoV.clustered, features = c("ASCL3", "CFTR", "FOXI1", "DMRT2", "V-ATPase")) #ionocytes
FeaturePlot(HCoV.clustered, features = "IGFBP5") #ionocytes1
FeaturePlot(HCoV.clustered, features = "CFTR") #ionocytes2


##########################################################################################################################
#GO Annotation

BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("EnhancedVolcano")

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(EnhancedVolcano)
library(ggrepel)


Cluster_0 <- Human_HCoV_Markers_clean_Cluster_0 %>% select(2,)
genes_cluster_0 <- rownames(Cluster_0)
GO_results <- enrichGO(gene = genes_cluster_0, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
as.data.frame(GO_results)
fit <- plot(barplot(GO_results, showCategory = 20))

fit

HCoV.clustered_markers_0
EnhancedVolcano(HCoV.clustered_markers_0, x = "avg_log2FC", y = "p_val_adj", lab = HCoV.clustered_markers_0$gene)

##################################################################################################################
#Renaming the found clusters to the respective cell types for downstream analysis differential expression analysis
# new.cluster.ids <- c("Secretory", "Ciliated", "Suprabasal", "Basal", "NA", "Deuterosomal")
# names(new.cluster.ids) <- levels(HCoV.clustered)
# HCoV.clustered <- RenameIdents(HCoV.clustered, new.cluster.ids)
# DimPlot(HCoV.clustered, reduction = "tsne", label.size = 6, repel = TRUE, label = TRUE, pt.size = 0.5) + NoLegend()

##################################################################################################################

#Saving the seurat object to disk to access in other script for analysis
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects")
SaveH5Seurat(HCoV.clustered, "HCoV.clustered", overwrite = TRUE)



