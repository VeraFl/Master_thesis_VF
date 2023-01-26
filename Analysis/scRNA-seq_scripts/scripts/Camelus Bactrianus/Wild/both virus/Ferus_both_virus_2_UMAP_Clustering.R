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

##########################################################################################################################
#Setting the color schemes for the script:

#Treatment colors(camel-229E, MERS, Mock)
cols_treat = c("#F39C6B","#A5668B", "#96ADC8")
show_col(cols_treat)

#Status colors infected vs uninfected
cols_stat = c("#BAA597", "#89043D","#BAA597")
show_col(cols_stat)

#Status colors infected, bystander and uninfected
cols_stat_2 = c("#519872", "#CD4631","#DEA47E")
show_col(cols_stat_2)

#Celltype colors for initial UMAP
cols_cells = hcl.colors(6, "Temps")
show_col(cols_cells)

cols = c("#807EBF", "#7EA8BE", "#089392", "#40B48B", "#9DCD84", "#EAE29C")
show_col(cols)


memory.limit(24000)
##########################################################################################################################
#Loading filtered and merged Seurat Object into Script
#Infected and Mock Samples merged (three samples)

#set the right working directory to the seurat objects folder
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.Ferus")

#Load the combined seurat object and name it "clustered"

#Ferus.clustered <- LoadH5Seurat("Ferus.combined.h5seurat")
Ferus.clustered <- LoadH5Seurat("Ferus.clustered.h5seurat")

?IntegrateData
#Check the right structure of the seurat object
str(Ferus.clustered)

##########################################################################################################################
###Perform an integrated analysis
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(Ferus.clustered) <- "integrated" #set the default assay as "integrated" - can be changed anytime

# Run the standard workflow for visualization and clustering
Ferus.clustered <- RunPCA(Ferus.clustered, npcs = 30, verbose = FALSE) #Principle component analysis
Ferus.clustered <- FindNeighbors(Ferus.clustered, reduction = "pca", dims = 1:30) #Finding neighbor within the defined dimensions
Ferus.clustered <- FindClusters(Ferus.clustered, resolution = 0.6) #Finding clusters with the set resolution
Ferus.clustered <- RunUMAP(Ferus.clustered, reduction = "pca", dims = 1:30) #Calculating the UMAP with the defined dimensions


# save the UMAP plots in different variations:
# Visualizing the initial unmodified UMAP
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_clustered_label.pdf",
    width=5, height=5)
DimPlot(Ferus.clustered, reduction = "umap", label.size = 6, repel = TRUE, label = TRUE, pt.size = 0.5) + NoLegend() +NoAxes()
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_clustered_legend.pdf",
    width=5, height=5)
DimPlot(Ferus.clustered, reduction = "umap", label.size = 6, repel = T, label = F, pt.size = 0.5)  +NoAxes()
dev.off()

# Visualizing the UMAP grouped by "Status"
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_group_stat_inf_uninf.pdf",
    width=10, height=5)
DimPlot(Ferus.clustered, reduction = "umap", split.by = "Treat", group.by = "Status", cols = cols_stat)
dev.off()

# Visualizing the UMAP split by "Treatment"
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_split_treat.pdf",
    width=15, height=5)
DimPlot(Ferus.clustered, reduction = "umap", split.by = "Treat", cols = cols)
dev.off()

# Visualizing the UMAP group by "Treatment"
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_group_treat.pdf",
    width=500, height=500)
DimPlot(Ferus.clustered, reduction = "umap", group.by = "Treat", cols = cols_treat)
dev.off()

# Visualizing the UMAP split by "Treat" and grouped by "Status" 
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_split_treat_group_stat.pdf",
    width=8, height=4)
DimPlot(Ferus.clustered, reduction = "umap", split.by = "Treat", group.by = "Status", cols = cols_stat)
dev.off()

##################################################################################################################
# Find celltype cluster markers for Annotation and Comparison
# find markers for every cluster compared to all remaining cells, report only the positive ones
# Threshold is a log fold change of 0.25
DefaultAssay(Ferus.clustered) <- "RNA"
Ferus.clustered_markers <- FindAllMarkers(Ferus.clustered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- Ferus.clustered_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

head(Ferus.clustered@meta.data)


row.names(Ferus.clustered_con_markers)

list <- unique(markers$gene)

DoHeatmap(Ferus.clustered, features = list)

DoHeatmap()
?FindAllMarkers
#dotplot for fresh unannotated clusters
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/Heatmap_initial_clusters.pdf",
    width=9, height=7)
DoHeatmap(Ferus.clustered, features = list)
dev.off()

##################################################################################################################
# Export the marker lists per cluster for cluster comparison to the other species (other R and python scripts)
# Select the cluster and the gene name columns from that list
Ferus_Markers_clean <- Ferus.clustered_markers %>% select(6:7,)

# save all the markers belonging to a cluster to the respective file
Ferus_Markers_clean_Cluster_0 <- Ferus_Markers_clean[Ferus_Markers_clean$cluster == 0, ]
Ferus_Markers_clean_Cluster_1 <- Ferus_Markers_clean[Ferus_Markers_clean$cluster == 1, ]
Ferus_Markers_clean_Cluster_2 <- Ferus_Markers_clean[Ferus_Markers_clean$cluster == 2, ]
Ferus_Markers_clean_Cluster_3 <- Ferus_Markers_clean[Ferus_Markers_clean$cluster == 3, ]
Ferus_Markers_clean_Cluster_4 <- Ferus_Markers_clean[Ferus_Markers_clean$cluster == 4, ]
Ferus_Markers_clean_Cluster_5 <- Ferus_Markers_clean[Ferus_Markers_clean$cluster == 5, ]
Ferus_Markers_clean_Cluster_6 <- Ferus_Markers_clean[Ferus_Markers_clean$cluster == 6, ]
Ferus_Markers_clean_Cluster_7 <- Ferus_Markers_clean[Ferus_Markers_clean$cluster == 7, ]
Ferus_Markers_clean_Cluster_8 <- Ferus_Markers_clean[Ferus_Markers_clean$cluster == 8, ]

# only choose the gene name column (for python script easier)
Ferus_Markers_clean_Cluster_0 <- Ferus_Markers_clean_Cluster_0 %>% select(2,)
Ferus_Markers_clean_Cluster_1 <- Ferus_Markers_clean_Cluster_1 %>% select(2,)
Ferus_Markers_clean_Cluster_2 <- Ferus_Markers_clean_Cluster_2 %>% select(2,)
Ferus_Markers_clean_Cluster_3 <- Ferus_Markers_clean_Cluster_3 %>% select(2,)
Ferus_Markers_clean_Cluster_4 <- Ferus_Markers_clean_Cluster_4 %>% select(2,)
Ferus_Markers_clean_Cluster_5 <- Ferus_Markers_clean_Cluster_5 %>% select(2,)
Ferus_Markers_clean_Cluster_6 <- Ferus_Markers_clean_Cluster_6 %>% select(2,)
Ferus_Markers_clean_Cluster_7 <- Ferus_Markers_clean_Cluster_7 %>% select(2,)
Ferus_Markers_clean_Cluster_8 <- Ferus_Markers_clean_Cluster_8 %>% select(2,)


# set the right working directory and save each gene marker list per cluster to one file

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Ferus_both")
write.table(Ferus_Markers_clean_Cluster_0,"Ferus_Cluster_0.txt", row.names = FALSE)
write.table(Ferus_Markers_clean_Cluster_1,"Ferus_Cluster_1.txt", row.names = FALSE)
write.table(Ferus_Markers_clean_Cluster_2,"Ferus_Cluster_2.txt", row.names = FALSE)
write.table(Ferus_Markers_clean_Cluster_3,"Ferus_Cluster_3.txt", row.names = FALSE)
write.table(Ferus_Markers_clean_Cluster_4,"Ferus_Cluster_4.txt", row.names = FALSE)
write.table(Ferus_Markers_clean_Cluster_5,"Ferus_Cluster_5.txt", row.names = FALSE)
write.table(Ferus_Markers_clean_Cluster_6,"Ferus_Cluster_6.txt", row.names = FALSE)
write.table(Ferus_Markers_clean_Cluster_7,"Ferus_Cluster_7.txt", row.names = FALSE)
write.table(Ferus_Markers_clean_Cluster_8,"Ferus_Cluster_8.txt", row.names = FALSE)

####################################################################################################################
#Plotting expression of different celltype specific markers 

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
DefaultAssay(Ferus.clustered) <- "RNA" 
# Visualization of the Feature plots of the main markers together
DimPlot(Ferus.clustered, reduction = "umap", label = TRUE, label.size = 5, repel = TRUE, pt.size = 0.5) + NoAxes() + NoLegend()
p2 <- FeaturePlot(Ferus.clustered, features = "FN1") + NoAxes() + NoLegend()
p3 <- FeaturePlot(Ferus.clustered, features = "MT1X") + NoAxes() + NoLegend()
p4 <- FeaturePlot(Ferus.clustered, features = "KRT7") + NoAxes() + NoLegend()
p5 <- FeaturePlot(Ferus.clustered, features = "MSMB") + NoAxes() + NoLegend()
p6 <- FeaturePlot(Ferus.clustered, features = "TPPP3") + NoAxes() + NoLegend()
p7 <- FeaturePlot(Ferus.clustered, features = "CCNO") + NoAxes() + NoLegend()
p8 <- FeaturePlot(Ferus.clustered, features = "CFTR") + NoAxes() + NoLegend()

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 2)

FeaturePlot(Ferus.clustered, features = "LOC102516515") + NoAxes() + NoLegend()
FeaturePlot(Ferus, features = "FOXJ1") + NoAxes() + NoLegend()
FeaturePlot(Alpaca, features = "FOXJ1") + NoAxes() + NoLegend()

FeaturePlot(Ferus.clustered, features = "TUB") + NoAxes() + NoLegend()

cols_cells = c("#807EBF", "#7EA8BE", "#089392", "#40B48B", "#9DCD84", "#EAE29C")

?FeaturePlot
str(Ferus.clustered)
DimPlot(Ferus.clustered, reduction = "umap", label = TRUE, label.size = 5, repel = TRUE, pt.size = 0.5) + NoAxes() + NoLegend()

DefaultAssay(Ferus.clustered) <- "RNA" 

p1 <- FeaturePlot(Ferus.clustered, slot = "data", features = "KRT5") + NoAxes() + NoLegend() #Basal
p2 <- FeaturePlot(Ferus.clustered, slot = "data", features = "SCGB1A1") + NoAxes() + NoLegend() #Club
p3 <- FeaturePlot(Ferus.clustered, slot = "data", features = "TPPP3") + NoAxes() + NoLegend() #Ciliated
p4 <- FeaturePlot(Ferus.clustered, slot = "data", features = "MSMB") + NoAxes() + NoLegend() #Secretory
p5 <- FeaturePlot(Ferus.clustered, slot = "data", features = "TOP2A") + NoAxes() + NoLegend() #Cluster 9



p7 <- VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "KRT5", pt.size = 0) + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25))) #Basal
p8 <- VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "SCGB1A1", pt.size = 0) + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)))#Club
p9 <- VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "TPPP3", pt.size = 0) + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)))#Ciliated
p10 <- VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "MSMB") + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)))#Secretory 
p11 <- VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "TOP2A", pt.size = 0) + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)))#Cluster 9


a <- as.data.frame(Ferus.clustered@assays$RNA@data@Dimnames[[1]])
VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "TPPP3", pt.size = 0) + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)))
VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "FOXJ1", pt.size = 0) + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)))
VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "TUBB", pt.size = 0) + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)))
VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "TUBB1", pt.size = 0) + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)))
VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "TUBB2A", pt.size = 0) + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)))
VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "TUBB3", pt.size = 0) + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)))
VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "TUBB4A", pt.size = 0) + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)))
VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "TUBB4B", pt.size = 0) + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)))
VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "TUBB6", pt.size = 0) + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)))
VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "TUBA4A", pt.size = 0) + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)))
VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "TUBA1A", pt.size = 0) + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)))
VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "TUBA8", pt.size = 0) + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)))
VlnPlot(Ferus.clustered, cols = cols_cells, slot = "data", features = "TUBAL3", pt.size = 0) + NoLegend() +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)))

?VlnPlot

head(Ferus_229E@meta.data)

Idents(Ferus.clustered) <- "Treat" 
Ferus_229E <- subset(x = Ferus.clustered, idents = c("Mock", "dcCoV-ACN4"))
Ferus_MERS <- subset(x = Ferus.clustered, idents = c("Mock", "MERS-CoV"))
Idents(Ferus_229E) <- "celltype" 
Idents(Ferus_MERS) <- "celltype" 
Idents(Ferus.clustered) <- "celltype" 

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/Vln_Plots_MERS_status_ciliated.pdf",
    width=17, height=15)
VlnPlot(Ferus_MERS, cols = cols_stat_2, slot = "data", features = c("DNAH7", "DNAH5", "LRRC6", "TPPP3", "FOXJ1", "TUBB", "TUBB2A",
                                                                         "LRRIQ1","TUBA4A","TUBA1A"), split.by = 'Status', pt.size = 0) +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)),
                                                                                                                                             legend.position = "top")
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/Vln_Plots_ACN4_status_ciliated.pdf",
    width=17, height=15)
VlnPlot(Ferus_229E, cols = cols_stat_2, slot = "data", features = c("DNAH7", "DNAH5", "LRRC6", "TPPP3", "FOXJ1", "TUBB", "TUBB2A",
                                                                    "LRRIQ1","TUBA4A","TUBA1A"), split.by = 'Status', pt.size = 0) +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)),
                                                                                                                                          legend.position = "top")
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/Vln_Plots_treat_ciliated.pdf",
    width=17, height=15)
VlnPlot(Ferus.clustered, cols = cols_treat, slot = "data", features = c("DNAH7", "DNAH5", "LRRC6", "TPPP3", "FOXJ1", "TUBB", "TUBB2A",
                                                                         "LRRIQ1","TUBA4A","TUBA1A"), split.plot = T, split.by = 'Treat', pt.size = 0) +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)),
                                                                                                                                               legend.position = "top")
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/Dot_Plots_RECEPTORS.pdf",
    width=6, height=4.2)
DotPlot(Ferus.clustered, features = c("DPP4", "ANPEP", "ACE2"), scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

dev.off()


cil_dotplot <- DotPlot(Ferus.clustered, features = c("DNAH5", "SPEF2", "PIFO", "LRRC6", "TPPP3", "FOXJ1", "TUBB", "TUBB1","TUBB2A",
                                      "TUBB4A","TUBB4B","TUBB6","TUBA4A","TUBA1A"), scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

sec_dotplot <- DotPlot(Ferus.clustered, features = c("FOXQ1", "SPDEF", "SCGB1A1", "MSMB", "MUC1", "MUC2", "MUC3A", 
                                                     "MUC4", "MUC5AC","MUC5B", "MUC6", "MUC7", "MUC15", "MUC19", "MUC20"), scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

deut_dotplot <- DotPlot(Ferus.clustered, features = c("PLK4", "CCNO", "CEP78", "DEUP1"), scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/camel_ferus_both_virus/Supplemenatary_fig_1.pdf",
    width=6, height=7)
cil_dotplot/sec_dotplot
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/Umaps_genes.pdf",
    width=3, height=15)
grid.arrange(p1, p2, p3, p4, p5, ncol = 1)
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/Expression_violin_genes.pdf",
    width=3, height=15)
grid.arrange(p7, p8, p9, p10, p11, ncol = 1)
dev.off()

DotPlot(Ferus.clustered, features = c("FN1","KRT5","KRT14","KRT17","MT1X","KRT6A","KRT16","KRT7","KRT4","KRT13","SCGB1A1","MSMB",
                                      "BPIFA1","CEACAM6", "SPDEF","MUC5B","FOXA2", "CYP26B1", "KRT8", "MUC5AC", "C5orf49", "FOXN4", "PLK4", "SYNE1", "SNTN", "CCNO", 
                                      "CFTR","MKI67","TOP2A"), scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


##################################################################################################################
#Renaming the found clusters to the respective cell types for downstream analysis differential expression analysis
levels(Ferus.clustered) #check the order of the unnamed clusters
new.cluster.ids <- c("Basal", "Secretory", "Secretory", "Basal", "Secretory", "Club", "Ciliated", "Club", "Basal", "Camel Cluster 1", "Camel Cluster 2","Basal") #enter new names in that order
names(new.cluster.ids) <- levels(Ferus.clustered) #assign the names to the new levels
levels(Ferus.clustered) #check the changed levels
Idents(Ferus.clustered)
Ferus.clustered <- RenameIdents(Ferus.clustered, new.cluster.ids) #change the identities to the celltype names
Ferus.clustered$celltype <- Idents(Ferus.clustered) #Make the previously defined Idents a new column in metadata (celltype)

Idents(Ferus.clustered) <- factor(x = Idents(Ferus.clustered), levels = c('Camel Cluster 1','Camel Cluster 2','Secretory', 'Ciliated', 'Club', 'Basal'))


levels(Ferus.clustered)
cols <- hcl.colors(15, "Temps")
cols_cells = c("#807EBF", "#7EA8BE", "#089392", "#40B48B", "#9DCD84", "#EAE29C")
show_col(cols_cells)

#Freshly assign colors to the celltype in the right order

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/umap_legend.pdf",
    width=7, height=5)
DimPlot(Ferus.clustered, reduction = "umap", label = FALSE, label.size = 6, repel = TRUE, pt.size = 0.5, cols = cols_cells) +NoAxes()
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/umap_label_2.pdf",
    width=5, height=5)
DimPlot(Ferus.clustered, reduction = "umap", label = TRUE, label.size = 6, repel = TRUE, pt.size = 0.5, cols = cols_cells) +NoAxes() +NoLegend()
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/umap_split_treat.pdf",
    width=8, height=4)
DimPlot(Ferus.clustered, reduction = "umap", split.by = "Treat", cols = cols_cells) + NoAxes() + NoLegend()
dev.off()

#self chosen markers
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/dot_plot.pdf",
    width=8, height=5)
DotPlot(Ferus.clustered, features = c("DDX5","LOC102518454","TUFT1","GOLGB1", "MKI67","TOP2A","CCNA2","MSMB", "MUC4", "PLK4","KRT4",
                                      "BPIFA1", "MUC5B", "KRT8","SYNE1","DEUP1", "SNTN","SCGB1A1", "KRT5","KRT17","KRT7"), scale.by = "radius", scale = T, cols = c("orange", "darkblue"))  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


DefaultAssay(Ferus.clustered) <- "RNA" 
head(Ferus.clustered@meta.data)

#Make dotplot for new annotated clusters
Ferus.clustered_markers <- FindAllMarkers(Ferus.clustered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- Ferus.clustered_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  print(n = 50)

Ferus.clustered <- Ferus
#Make dotplot with the new annotated clusters and their gene sets
#choose only the unique gene names
list_2 <- unique(markers$gene)

#Make plot
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/dot_plot_new_clusters_legend.pdf",
    width=9, height=6)
DotPlot(Ferus.clustered, features = list_2, scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                                                                                                                                  legend.position = "top")
dev.off()


#Save the gene lists of the newly annotated clusters to compare with human clusters
Camel_Markers_clean <- Ferus.clustered_markers %>% dplyr::select(6:7,)

Camel_Markers_clean_Cluster_0 <- Camel_Markers_clean[Camel_Markers_clean$cluster == "Secretory", ]
Camel_Markers_clean_Cluster_1 <- Camel_Markers_clean[Camel_Markers_clean$cluster == "Ciliated", ]
Camel_Markers_clean_Cluster_2 <- Camel_Markers_clean[Camel_Markers_clean$cluster == "Club", ]
Camel_Markers_clean_Cluster_3 <- Camel_Markers_clean[Camel_Markers_clean$cluster == "Basal", ]
Camel_Markers_clean_Cluster_4 <- Camel_Markers_clean[Camel_Markers_clean$cluster == "Camel Cluster 1", ]
Camel_Markers_clean_Cluster_5 <- Camel_Markers_clean[Camel_Markers_clean$cluster == "Camel Cluster 2", ]


Camel_Markers_clean_Cluster_0 <- Camel_Markers_clean_Cluster_0 %>% dplyr::select(2,)
Camel_Markers_clean_Cluster_1 <- Camel_Markers_clean_Cluster_1 %>% dplyr::select(2,)
Camel_Markers_clean_Cluster_2 <- Camel_Markers_clean_Cluster_2 %>% dplyr::select(2,)
Camel_Markers_clean_Cluster_3 <- Camel_Markers_clean_Cluster_3 %>% dplyr::select(2,)
Camel_Markers_clean_Cluster_4 <- Camel_Markers_clean_Cluster_4 %>% dplyr::select(2,)
Camel_Markers_clean_Cluster_5 <- Camel_Markers_clean_Cluster_5 %>% dplyr::select(2,)


setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Cluster_markers/Ferus_both")
write.csv(Ferus.clustered_markers, "Cluster_markers_Camel", row.names = F)

write.table(Camel_Markers_clean_Cluster_0,"Ferus_Cluster_0.txt", row.names = FALSE)
write.table(Camel_Markers_clean_Cluster_1,"Ferus_Cluster_1.txt", row.names = FALSE)
write.table(Camel_Markers_clean_Cluster_2,"Ferus_Cluster_2.txt", row.names = FALSE)
write.table(Camel_Markers_clean_Cluster_3,"Ferus_Cluster_3.txt", row.names = FALSE)
write.table(Camel_Markers_clean_Cluster_4,"Ferus_Cluster_4.txt", row.names = FALSE)
write.table(Camel_Markers_clean_Cluster_5,"Ferus_Cluster_5.txt", row.names = FALSE)




#Saving the seurat object to disk to access in other script for analysis
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.ferus")
SaveH5Seurat(Ferus.clustered, "Ferus.clustered", overwrite = TRUE)


