# UMAP and Clustering of Single-cell Seq data
# 21.06.2022

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
library(glmGamPoi)

#Treatment colors(MERS, Mock)
cols_treat = c("#A5668B", "#96ADC8")
show_col(cols_treat)
#Status colors
cols_stat = c("#CD4631","#DEA47E")
show_col(cols_stat)
#Celltype colors
cols_cell = hcl.colors(8, "Temps")
show_col(cols_cell)

##########################################################################################################################
#Loading filtered and merged Seurat Object into Script
#Infected and Mock Samples merged

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/V.pacos")

MERS.clustered <- LoadH5Seurat("MERS.combined.h5seurat")

##########################################################################################################################
#Analysis of the combined dataset
DimHeatmap(MERS.clustered)
# Run the standard workflow for visualization and clustering
#MERS.clustered <- ScaleData(MERS.clustered, vars.to.regress = c("Percent.Viral", "S.Score", "G2M.Score"), do.par = TRUE, num.cores = 16, verbose = FALSE)
MERS.clustered <- RunPCA(MERS.clustered, npcs = 30, verbose = FALSE)
MERS.clustered <- FindNeighbors(MERS.clustered, reduction = "pca", dims = 1:30)
MERS.clustered <- FindClusters(MERS.clustered, resolution = 0.3)
MERS.clustered <- RunUMAP(MERS.clustered, reduction = "pca", dims = 1:30)
# MERS.clustered <- RunTSNE(MERS.clustered, reduction = "pca", dims = 1:5)

png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_MERS/UMAP_clusters.png",
    width=1000, height=500)
DimPlot(MERS.clustered, reduction = "umap", cols = cols_cell) + NoLegend()
dev.off()

png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_MERS/UMAP_by_status.png",
    width=1000, height=500)
DimPlot(MERS.clustered, reduction = "umap", group.by = "Status", cols = cols_stat)
dev.off()

png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_MERS/UMAP_by_treat.png",
    width=1000, height=500)
DimPlot(MERS.clustered, reduction = "umap", group.by = "Treat", cols = cols_treat)
dev.off()


##################################################################################################################
# Markers to find cell type clusters
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
MERS.clustered_markers <- FindAllMarkers(MERS.clustered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MERS.clustered_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_MERS")
write.table(Alpaca_MERS_Markers_clean,"Human_MERS_Markers_clean.txt", row.names = FALSE)

Alpaca_MERS_Markers_clean <- MERS.clustered_markers %>% select(6:7,)

Alpaca_MERS_Markers_clean_Cluster_0 <- Alpaca_MERS_Markers_clean[Alpaca_MERS_Markers_clean$cluster == 0, ]
Alpaca_MERS_Markers_clean_Cluster_1 <- Alpaca_MERS_Markers_clean[Alpaca_MERS_Markers_clean$cluster == 1, ]
Alpaca_MERS_Markers_clean_Cluster_2 <- Alpaca_MERS_Markers_clean[Alpaca_MERS_Markers_clean$cluster == 2, ]
Alpaca_MERS_Markers_clean_Cluster_3 <- Alpaca_MERS_Markers_clean[Alpaca_MERS_Markers_clean$cluster == 3, ]
Alpaca_MERS_Markers_clean_Cluster_4 <- Alpaca_MERS_Markers_clean[Alpaca_MERS_Markers_clean$cluster == 4, ]
Alpaca_MERS_Markers_clean_Cluster_5 <- Alpaca_MERS_Markers_clean[Alpaca_MERS_Markers_clean$cluster == 5, ]
Alpaca_MERS_Markers_clean_Cluster_6 <- Alpaca_MERS_Markers_clean[Alpaca_MERS_Markers_clean$cluster == 6, ]
Alpaca_MERS_Markers_clean_Cluster_7 <- Alpaca_MERS_Markers_clean[Alpaca_MERS_Markers_clean$cluster == 7, ]



Alpaca_MERS_Markers_clean_Cluster_0 <- Alpaca_MERS_Markers_clean_Cluster_0 %>% select(2,)
Alpaca_MERS_Markers_clean_Cluster_1 <- Alpaca_MERS_Markers_clean_Cluster_1 %>% select(2,)
Alpaca_MERS_Markers_clean_Cluster_2 <- Alpaca_MERS_Markers_clean_Cluster_2 %>% select(2,)
Alpaca_MERS_Markers_clean_Cluster_3 <- Alpaca_MERS_Markers_clean_Cluster_3 %>% select(2,)
Alpaca_MERS_Markers_clean_Cluster_4 <- Alpaca_MERS_Markers_clean_Cluster_4 %>% select(2,)
Alpaca_MERS_Markers_clean_Cluster_5 <- Alpaca_MERS_Markers_clean_Cluster_5 %>% select(2,)
Alpaca_MERS_Markers_clean_Cluster_6 <- Alpaca_MERS_Markers_clean_Cluster_6 %>% select(2,)
Alpaca_MERS_Markers_clean_Cluster_7 <- Alpaca_MERS_Markers_clean_Cluster_7 %>% select(2,)



setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_MERS")
write.table(Alpaca_MERS_Markers_clean_Cluster_0,"Alpaca_MERS_Cluster_0.txt", row.names = FALSE)
write.table(Alpaca_MERS_Markers_clean_Cluster_1,"Alpaca_MERS_Cluster_1.txt", row.names = FALSE)
write.table(Alpaca_MERS_Markers_clean_Cluster_2,"Alpaca_MERS_Cluster_2.txt", row.names = FALSE)
write.table(Alpaca_MERS_Markers_clean_Cluster_3,"Alpaca_MERS_Cluster_3.txt", row.names = FALSE)
write.table(Alpaca_MERS_Markers_clean_Cluster_4,"Alpaca_MERS_Cluster_4.txt", row.names = FALSE)
write.table(Alpaca_MERS_Markers_clean_Cluster_5,"Alpaca_MERS_Cluster_5.txt", row.names = FALSE)
write.table(Alpaca_MERS_Markers_clean_Cluster_6,"Alpaca_MERS_Cluster_6.txt", row.names = FALSE)
write.table(Alpaca_MERS_Markers_clean_Cluster_7,"Alpaca_MERS_Cluster_7.txt", row.names = FALSE)



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

FeaturePlot(MERS.clustered, features = "KRT16")
FeaturePlot(MERS.clustered, features = "KRT23")
FeaturePlot(MERS.clustered, features = "CCNO")
FeaturePlot(MERS.clustered, features = "DEUP1")
FeaturePlot(MERS.clustered, features = "DPP4")


#Cluster0
FeaturePlot(MERS.clustered, features = c("COL10A1", "CDH6")) #unidentifiable
FeaturePlot(MERS.clustered, features = "COL10A1") 
FeaturePlot(MERS.clustered, features = "CDH6") 


#Cluster1
FeaturePlot(MERS.clustered, features = c("LOC116662762", "LOC102523613")) #Club cells
FeaturePlot(MERS.clustered, features = "LOC1166627621")
FeaturePlot(MERS.clustered, features = "LOC102523613")


#Cluster2
FeaturePlot(MERS.clustered, features = c("TPPP3", "SEC14L3")) #Ciliated
FeaturePlot(MERS.clustered, features = "TPPP3")
FeaturePlot(MERS.clustered, features = "SEC14L3")

#Cluster3
FeaturePlot(MERS.clustered, features = c("CCDC17", "LOC102519002")) #Ciliated
FeaturePlot(MERS.clustered, features = "CCDC17")
FeaturePlot(MERS.clustered, features = "LOC102519002")

#Cluster4
FeaturePlot(MERS.clustered, features = c("TAF1D", "TTN"))
FeaturePlot(MERS.clustered, features = "TAF1D")
FeaturePlot(MERS.clustered, features = "TTN")


##################################################################################################################
#Renaming the found clusters to the respective cell types for downstream analysis differential expression analysis
new.cluster.ids <- c("Basal", "Club", "Ciliated", "Ciliated", "Club")
names(new.cluster.ids) <- levels(MERS.clustered)
MERS.clustered <- RenameIdents(MERS.clustered, new.cluster.ids)
DimPlot(MERS.clustered, reduction = "umap", label.size = 6, repel = TRUE, label = TRUE, pt.size = 0.5) + NoLegend()


#Saving the seurat object to disk to access in other script for analysis
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/V.pacos")
SaveH5Seurat(MERS.clustered, "MERS.clustered", overwrite = TRUE)



#################################################################################################################
#Mean Viral expression per Cluster
df <- MERS.clustered@meta.data %>% select(c("seurat_clusters","Percent.Viral", "Status", "Treat"))

#Average Viral Percentage in each Cluster
plot1 <- df %>%
  group_by(seurat_clusters) %>%
  summarise_at(vars(Percent.Viral), list(average_viral_percentage = mean))

#Plot1
p1 <- plot1 %>% ggplot(aes(x = seurat_clusters, y = average_viral_percentage))+
  geom_col()


#Mean cell number per Status and Treatment in each Cluster
plot2 <- df %>%
  group_by(seurat_clusters, Status, Treat) %>%
  count(Status, Treat)

#Plot2
plot2 %>% ggplot(aes(x = Treat, y = n, fill = Status))+
  geom_col()+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        legend.position="bottom")+
  facet_grid(cols = vars(seurat_clusters))


#Normalized bystander/infected/non-infected ratios
plot3 <- plot2 %>% 
  group_by(seurat_clusters) %>%
  summarize(Sum = sum(n),
            n = n,
            seurat_clusters = seurat_clusters,
            Treat = Treat,
            Status = Status) %>%
  mutate("Freq" = n/Sum)

#Plot3
plot3 %>% ggplot(aes(x = Treat, y = Freq, fill = Status))+
  geom_col()+
  ylim(0,1)+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        legend.position="bottom")+
  facet_grid(cols = vars(seurat_clusters))

plot4 <- plot3

plot4 %>% ggplot(aes(x = seurat_clusters, y = Freq, fill = Status))+
  geom_col(position="fill")+
  scale_fill_manual(values = cols)+
  ylim(0,1)+
  theme(legend.position="bottom")


#########################################################################################################
#Correlation DPP4 expression and Viral_counts
#Extract the info in which row te DPP4 expression data is
a <- as.data.frame(MERS.clustered@assays$RNA@counts@Dimnames[[1]])

#5269
#show the specific location in the matrix for control
DPP4 <- as.matrix(MERS.clustered@assays$RNA@counts[5269,])
#Count the DPP4 expression
MERS.clustered <- AddMetaData(object = MERS.clustered, metadata = DPP4, col.name = "Total.DPP4.Counts")

head(MERS.clustered@meta.data)

df1 <- MERS.clustered@meta.data %>% select(c("seurat_clusters","Total.DPP4.Counts", "Status", "Treat"))

#Average Viral Percentage in each Cluster
plot1_1 <- df1 %>%
  group_by(seurat_clusters) %>%
  summarise_at(vars(Total.DPP4.Counts), list(average_dpp4 = mean))


grid.arrange(p1, p1_1, ncol=2)

#Plot1
p1_1 <- plot1_1 %>% ggplot(aes(x = seurat_clusters, y = average_dpp4))+
  geom_col()


plot2_1 <- df1 %>%
  group_by(seurat_clusters, Status, Treat) %>%
  count(Status, Treat)


plot2_1 %>% ggplot(aes(x = Treat, y = n, fill = Status))+
  geom_col()+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        legend.position="bottom")+
  facet_grid(cols = vars(seurat_clusters))


