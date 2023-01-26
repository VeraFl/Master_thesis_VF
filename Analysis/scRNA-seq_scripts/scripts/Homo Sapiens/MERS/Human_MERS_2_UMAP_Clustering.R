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
library(glmGamPoi)


#Treatment colors(camel-229E, MERS, Mock)
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

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens")

MERS.clustered <- LoadH5Seurat("MERS.combined.h5seurat")
MERS.clustered

##########################################################################################################################
# head(MERS.clustered@meta.data)
# #How many cells are infected with that threshold?
# MERS_groups <- MERS.clustered@meta.data %>% select(c("Status", "Treat"))
# 
# MERS_groups <- MERS_groups %>%
#   group_by(Status, Treat) %>%
#   count(Status, Treat)
# 
# # MERS.clustered <- LoadH5Seurat("MERS.clustered.h5seurat")
# # MERS.clustered
##########################################################################################################################
#Analysis of the combined dataset
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# 
# MERS.clustered <- CellCycleScoring(MERS.clustered, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# Run the standard workflow for visualization and clustering
# MERS.clustered <- ScaleData(MERS.clustered, vars.to.regress = c("Percent.Viral", "S.Score", "G2M.Score"), do.par = TRUE, num.cores = 16, verbose = FALSE)
MERS.clustered <- RunPCA(MERS.clustered, npcs = 30, verbose = FALSE)

MERS.clustered <- JackStraw(MERS.clustered, num.replicate = 100)
# MERS.clustered <- ScoreJackStraw(MERS.clustered, dims = 1:20)
# JackStrawPlot(MERS.clustered, dims = 1:20)
# ElbowPlot(MERS.clustered)


MERS.clustered <- FindNeighbors(MERS.clustered, reduction = "pca", dims = 1:25)
MERS.clustered <- FindClusters(MERS.clustered, resolution = 0.3)
MERS.clustered <- RunUMAP(MERS.clustered, reduction = "pca", dims = 1:25)
# MERS.clustered <- RunTSNE(MERS.clustered, reduction = "pca", dims = 1:5)

png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_MERS/UMAP_clusters.png",
    width=500, height=600)
DimPlot(MERS.clustered, reduction = "umap",label.size = 6, repel = TRUE, label = TRUE, pt.size = 0.5, cols = cols_cell) + NoLegend()
dev.off()

# DimPlot(MERS.clustered, reduction = "tsne")
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_MERS/UMAP_status.png",
    width=500, height=600)
DimPlot(MERS.clustered, reduction = "umap", group.by = "Status", cols = cols_stat)
dev.off()

png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_MERS/UMAP_treat.png",
    width=500, height=600)
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


#Saving markers gene lists
Human_MERS_Markers_clean <- MERS.clustered_markers %>% select(7,)


setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Human_MERS")
write.table(Human_MERS_Markers_clean,"Human_MERS_Markers_clean.txt", row.names = FALSE)

Human_MERS_Markers_clean <- MERS.clustered_markers %>% select(6:7,)

Human_MERS_Markers_clean_Cluster_0 <- Human_MERS_Markers_clean[Human_MERS_Markers_clean$cluster == 0, ]
Human_MERS_Markers_clean_Cluster_1 <- Human_MERS_Markers_clean[Human_MERS_Markers_clean$cluster == 1, ]
Human_MERS_Markers_clean_Cluster_2 <- Human_MERS_Markers_clean[Human_MERS_Markers_clean$cluster == 2, ]
Human_MERS_Markers_clean_Cluster_3 <- Human_MERS_Markers_clean[Human_MERS_Markers_clean$cluster == 3, ]
Human_MERS_Markers_clean_Cluster_4 <- Human_MERS_Markers_clean[Human_MERS_Markers_clean$cluster == 4, ]
Human_MERS_Markers_clean_Cluster_5 <- Human_MERS_Markers_clean[Human_MERS_Markers_clean$cluster == 5, ]
Human_MERS_Markers_clean_Cluster_6 <- Human_MERS_Markers_clean[Human_MERS_Markers_clean$cluster == 6, ]
Human_MERS_Markers_clean_Cluster_7 <- Human_MERS_Markers_clean[Human_MERS_Markers_clean$cluster == 7, ]


Human_MERS_Markers_clean_Cluster_0 <- Human_MERS_Markers_clean_Cluster_0 %>% select(2,)
Human_MERS_Markers_clean_Cluster_1 <- Human_MERS_Markers_clean_Cluster_1 %>% select(2,)
Human_MERS_Markers_clean_Cluster_2 <- Human_MERS_Markers_clean_Cluster_2 %>% select(2,)
Human_MERS_Markers_clean_Cluster_3 <- Human_MERS_Markers_clean_Cluster_3 %>% select(2,)
Human_MERS_Markers_clean_Cluster_4 <- Human_MERS_Markers_clean_Cluster_4 %>% select(2,)
Human_MERS_Markers_clean_Cluster_5 <- Human_MERS_Markers_clean_Cluster_5 %>% select(2,)
Human_MERS_Markers_clean_Cluster_6 <- Human_MERS_Markers_clean_Cluster_6 %>% select(2,)
Human_MERS_Markers_clean_Cluster_7 <- Human_MERS_Markers_clean_Cluster_7 %>% select(2,)


setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Human_MERS")
write.table(Human_MERS_Markers_clean_Cluster_0,"Human_MERS_Markers_clean_Cluster_0.txt", row.names = FALSE)
write.table(Human_MERS_Markers_clean_Cluster_1,"Human_MERS_Markers_clean_Cluster_1.txt", row.names = FALSE)
write.table(Human_MERS_Markers_clean_Cluster_2,"Human_MERS_Markers_clean_Cluster_2.txt", row.names = FALSE)
write.table(Human_MERS_Markers_clean_Cluster_3,"Human_MERS_Markers_clean_Cluster_3.txt", row.names = FALSE)
write.table(Human_MERS_Markers_clean_Cluster_4,"Human_MERS_Markers_clean_Cluster_4.txt", row.names = FALSE)
write.table(Human_MERS_Markers_clean_Cluster_5,"Human_MERS_Markers_clean_Cluster_5.txt", row.names = FALSE)
write.table(Human_MERS_Markers_clean_Cluster_6,"Human_MERS_Markers_clean_Cluster_6.txt", row.names = FALSE)
write.table(Human_MERS_Markers_clean_Cluster_7,"Human_MERS_Markers_clean_Cluster_7.txt", row.names = FALSE)

#cluster1 = ciliated

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
FeaturePlot(MERS.clustered, features = "mers") #viral

#All

FeaturePlot(MERS.clustered, features = c("FN1", "KRT7", "MSMB", "TPPP3", "CFTR", "DNAAF1", "RRAD", "CCNO"))

DefaultAssay(MERS.clustered)
#integrated
?IntegrateData()
DefaultAssay(MERS.clustered) <- "integrated"
p1 <- FeaturePlot(MERS.clustered, slot = "scale.data", features = "FN1", split.by = )
p1 <- p1 + ggtitle("Scale Data")
p2 <- FeaturePlot(MERS.clustered, slot = "data", features = "FN1")
p2 <- p2 + ggtitle("Data - Back calculated normalized UMI counts")
plot_integrated <- p1 | p2
plot_integrated + ggtitle("Integrated")

DefaultAssay(MERS.clustered) <- "SCT"
p3 <- FeaturePlot(MERS.clustered, slot = "scale.data", features = "FN1")
p3 + ggtitle("Scale Data")
p4 <- FeaturePlot(MERS.clustered, slot = "data", features = "FN1")
p4 + ggtitle("Data")
p5 <- FeaturePlot(MERS.clustered, slot = "counts", features = "FN1")
p5 + ggtitle("Counts")
plot_sct <- p3 | p4 | p5
plot_sct + ggtitle("SCTransform")

DefaultAssay(MERS.clustered) <- "RNA"
p6 <- FeaturePlot(MERS.clustered, slot = "counts", features = "FN1")
p6 + ggitle("Counts")
p7 <- FeaturePlot(MERS.clustered, slot = "data", features = "FN1")
p7 + ggitle("Data")
plot_rna <- p6 | p7
plot_rna + ggtitle("RNA")

plot_integrated / plot_sct / plot_rna

str(MERS.clustered)

#BASAL - INF
FeaturePlot(MERS.clustered, features = "KRT5") #basal1 -interesting
FeaturePlot(MERS.clustered, features = "DAPL1")
FeaturePlot(MERS.clustered, features = "TP63") #basal2
FeaturePlot(MERS.clustered, features = "FN1") #basal3 -interesting
FeaturePlot(MERS.clustered, features = "DLK2") #basal4
FeaturePlot(MERS.clustered, features = "KRT6A") #basal5

#SUPRABASAL
FeaturePlot(MERS.clustered, features = "KRT4") #suprabasal1
FeaturePlot(MERS.clustered, features = "KRT13") #suprabasal2 - interesting
FeaturePlot(MERS.clustered, features = "KRT16") #suprabasal3 - interesting
FeaturePlot(MERS.clustered, features = "KRT23") #suprabasal4

#TUFT
FeaturePlot(MERS.clustered, features = "POUF2f3") #tuft1
FeaturePlot(MERS.clustered, features = "TRPM5") #tuft2
FeaturePlot(MERS.clustered, features = "LRMP") #tuft3

#IONOCYTES
FeaturePlot(MERS.clustered, features = "FOXI1") #ionocytes1
FeaturePlot(MERS.clustered, features = "CFTR") #ionocytes2

##PNEC
FeaturePlot(MERS.clustered, features = c("PCSK1N")) #PNEC1
FeaturePlot(MERS.clustered, features = c("SCGN")) #PNEC2
FeaturePlot(MERS.clustered, features = c("NEB")) #PNEC3
FeaturePlot(MERS.clustered, features = c("HOXB1")) #PNEC4
FeaturePlot(MERS.clustered, features = c("ASCL1")) #PNEC5
FeaturePlot(MERS.clustered, features = c("ASCL2")) #PNEC6
FeaturePlot(MERS.clustered, features = c("FOXA2")) #PNEC7
FeaturePlot(MERS.clustered, features = c("PSMD5")) #PNEC8
FeaturePlot(MERS.clustered, features = c("NGF")) #PNEC9

#CLUB
FeaturePlot(MERS.clustered, features = c("SCGB1A1")) #club -interesting
FeaturePlot(MERS.clustered, features = c("KRT15")) #club
FeaturePlot(MERS.clustered, features = c("BPIFA1")) #club
FeaturePlot(MERS.clustered, features = c("KRT7")) #club
FeaturePlot(MERS.clustered, features = c("MSMB")) #club
FeaturePlot(MERS.clustered, features = c("KRT19")) #club
FeaturePlot(MERS.clustered, features = c("CYP2F2")) #club
FeaturePlot(MERS.clustered, features = c("LYPD2")) #club
FeaturePlot(MERS.clustered, features = c("CBR2")) #club

#GOBLET/SECRETORY
FeaturePlot(MERS.clustered, features = "MUC5B") #secretory1 - interesting
FeaturePlot(MERS.clustered, features = "SCGB3A1") #secretory2 -interesting
FeaturePlot(MERS.clustered, features = "MUC5AC") #secretory3 (goblet)
FeaturePlot(MERS.clustered, features = "FOXQ1") #goblet1
FeaturePlot(MERS.clustered, features = "SPDEF") #goblet2
FeaturePlot(MERS.clustered, features = "GP2") #goblet2
FeaturePlot(MERS.clustered, features = "SPDEF") #goblet2

#DEUTEROSOMAL
FeaturePlot(MERS.clustered, features = "PLK4") #deuterosomal - interesting
FeaturePlot(MERS.clustered, features = "CCNO") #deuterosomal
FeaturePlot(MERS.clustered, features = "CEP78") #deuterosomal
FeaturePlot(MERS.clustered, features = "DEUP1") #deuterosomal - interesting


#CILIATED
FeaturePlot(MERS.clustered, features = "FOXN4") #preciliated
FeaturePlot(MERS.clustered, features = "FOXJ1") #ciliated1 - interesting
FeaturePlot(MERS.clustered, features = "PIFO") #ciliated2
FeaturePlot(MERS.clustered, features = "TPPP3") #ciliated3
FeaturePlot(MERS.clustered, features = "SPEF2") #ciliated4 -interesting
FeaturePlot(MERS.clustered, features = "DNAH5") #ciliated5
FeaturePlot(MERS.clustered, features = "LRRC6") #ciliated6
FeaturePlot(MERS.clustered, features = "CCDC153") #ciliated6
FeaturePlot(MERS.clustered, features = "CCDC113") #ciliated6
FeaturePlot(MERS.clustered, features = "MLF1") #ciliated6
FeaturePlot(MERS.clustered, features = "LZTFL1") #ciliated6

###########################################################################################################################

FeaturePlot(MERS.clustered, features = c("orf1ab", "S", "orf3", "orf4a", "orf4b", "orf5", "E", "M", "orf8b")) #viral genes
FeaturePlot(MERS.clustered, features = c("PCNA", "TOP2A", "MCM6", "MKI67")) #cell cycle state genes

#Cluster0
FeaturePlot(MERS.clustered, features = c("RIMS1", "FAM3D"))


#Cluster1
FeaturePlot(MERS.clustered, features = c("ADRB1", "ELK3 ")) 


#Cluster2
FeaturePlot(MERS.clustered, features = c("CRCT1", "TRERF1")) 


#Cluster3
FeaturePlot(MERS.clustered, features = c("TMPRSS4", "CCNL1")) 


#Cluster4
FeaturePlot(MERS.clustered, features = c("GJA1", "SPARC"))


#Cluster5
FeaturePlot(MERS.clustered, features = c("ODAD1", "PTRH1")) 


#Cluster6
FeaturePlot(MERS.clustered, features = c("BPIFA1", "SCGB3A1")) 


#Cluster7
FeaturePlot(MERS.clustered, features = c("ATP6V0B", "JADE1")) 

##################################################################################################################
#Renaming the found clusters to the respective cell types for downstream analysis differential expression analysis
levels(MERS.clustered)
new.cluster.ids <- c("Club", "Ciliated", "Club", "Basal", "Ciliated", "Deuterosomal")
names(new.cluster.ids) <- levels(MERS.clustered)
MERS.clustered <- RenameIdents(MERS.clustered, new.cluster.ids)
DimPlot(MERS.clustered, reduction = "umap", label.size = 6, repel = TRUE, label = TRUE, pt.size = 0.5) + NoLegend()

#Saving the seurat object to disk to access in other script for analysis
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens")
SaveH5Seurat(MERS.clustered, "MERS.clustered", overwrite = TRUE)

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens")
MERS <- LoadH5Seurat("MERS.h5seurat")



#################################################################################################################
#Mean Viral expression per Cluster
df <- MERS@meta.data %>% select(c("seurat_clusters","Percent.Viral", "Status", "Treat"))

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



