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
#Status colors
cols_stat = c("#DEA47E","#CD4631","#DEA47E")
show_col(cols_stat)

cols_stat_2 = c("#519872", "#CD4631","#DEA47E")
show_col(cols_stat_2)

#Celltype colors
cols_cell = hcl.colors(7, "Temps")
show_col(cols_cell)
cols_cell

cols = c("#8A6E7E", "#40679B", "#089392", "#40B48B", "#9DCD84", "#EAE29C")
show_col(cols)

memory.limit(24000)
##########################################################################################################################
#Loading filtered and merged Seurat Object into Script
#Infected and Mock Samples merged

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/V.pacos")

#Alpaca.clustered <- LoadH5Seurat("Alpaca.combined.h5seurat") #Import either seurat object
Alpaca.clustered <- LoadH5Seurat("Alpaca.clustered.h5seurat") #Import either seurat object

head(Alpaca.clustered@meta.data)

##########################################################################################################################
# Run the standard workflow for visualization and clustering

Alpaca.clustered <- RunPCA(Alpaca.clustered, npcs = 30, verbose = FALSE) #Run PCA as dimensionality reduction

#Visualize the dimensions
#VizDimLoadings(Alpaca.clustered, dims = 1:5, reduction = "pca")

#Using the chosen
Alpaca.clustered <- FindNeighbors(Alpaca.clustered, reduction = "pca", dims = 1:30)
Alpaca.clustered <- FindClusters(Alpaca.clustered, resolution = 0.6)
Alpaca.clustered <- RunUMAP(Alpaca.clustered, reduction = "pca", dims = 1:30)



show_col(cols)
levels(Alpaca.clustered)
#save the UMAP plots in different variations
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/umap_cluster_label.pdf",
    width=5, height=5)
DimPlot(Alpaca.clustered, reduction = "umap", label.size = 6, repel = TRUE, label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/umap_cluster_legend.pdf",
    width=5, height=5)
DimPlot(Alpaca.clustered, reduction = "umap", label.size = 6, repel = TRUE, label = TRUE, pt.size = 0.5)
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/umap_split_status_treat.pdf",
    width=8, height=4)
DimPlot(Alpaca.clustered, reduction = "umap", split.by = "Treat", pt.size = 0.5, group.by = "Status", cols = cols_stat)
dev.off()


DimPlot(Alpaca.clustered, reduction = "umap", split.by = "Treat", pt.size = 0.5, cols = cols_cell) + NoAxes() + NoLegend()
##################################################################################################################
# Markers to find cell type clusters

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
Alpaca.clustered_markers <- FindAllMarkers(Alpaca.clustered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- Alpaca.clustered_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

list <- unique(markers$gene)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/Heatmap_initial_clusters.pdf",
    width=9, height=7)
DoHeatmap(Alpaca.clustered, features = list)
dev.off()

#dotplot for fresh unannotated clusters
DotPlot(Ferus.clustered, features = list, scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



DotPlot(Alpaca.clustered, features = c("LOXL4", "MKI67","TOP2A","COL7A1", "TNC", "LOC102527556", "KRT19", "HES1",
                                       "KRT5", "TIMP3", "TP63", "LOC102527828", "CCN1",
                                       "ID1","KRT4", "LOC102526268", "MUC21","KRT7", "BPIFA1", "KRT23",
                                       "KRT8", "PLK4", "SYNE1","CFTR", "LRRIQ1", "TPPP3",
                                       "SCGB1A1","MSMB", "VMO1", "LOC116279189"), scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DotPlot(Alpaca.clustered, features = list, scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



# setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_both")
# write.table(Alpaca_Alpaca_Markers_clean,"Human_Alpaca_Markers_clean.txt", row.names = FALSE)

Alpaca_Markers_clean <- Alpaca.clustered_markers %>% dplyr::select(6:7,)

Alpaca_Markers_clean_Cluster_0 <- Alpaca_Markers_clean[Alpaca_Markers_clean$cluster == 0, ]
Alpaca_Markers_clean_Cluster_1 <- Alpaca_Markers_clean[Alpaca_Markers_clean$cluster == 1, ]
Alpaca_Markers_clean_Cluster_2 <- Alpaca_Markers_clean[Alpaca_Markers_clean$cluster == 2, ]
Alpaca_Markers_clean_Cluster_3 <- Alpaca_Markers_clean[Alpaca_Markers_clean$cluster == 3, ]
Alpaca_Markers_clean_Cluster_4 <- Alpaca_Markers_clean[Alpaca_Markers_clean$cluster == 4, ]
Alpaca_Markers_clean_Cluster_5 <- Alpaca_Markers_clean[Alpaca_Markers_clean$cluster == 5, ]
Alpaca_Markers_clean_Cluster_6 <- Alpaca_Markers_clean[Alpaca_Markers_clean$cluster == 6, ]
Alpaca_Markers_clean_Cluster_7 <- Alpaca_Markers_clean[Alpaca_Markers_clean$cluster == 7, ]
Alpaca_Markers_clean_Cluster_8 <- Alpaca_Markers_clean[Alpaca_Markers_clean$cluster == 8, ]
Alpaca_Markers_clean_Cluster_9 <- Alpaca_Markers_clean[Alpaca_Markers_clean$cluster == 9, ]
Alpaca_Markers_clean_Cluster_10 <- Alpaca_Markers_clean[Alpaca_Markers_clean$cluster == 10, ]



Alpaca_Markers_clean_Cluster_0 <- Alpaca_Markers_clean_Cluster_0 %>% dplyr::select(2,)
Alpaca_Markers_clean_Cluster_1 <- Alpaca_Markers_clean_Cluster_1 %>% dplyr::select(2,)
Alpaca_Markers_clean_Cluster_2 <- Alpaca_Markers_clean_Cluster_2 %>% dplyr::select(2,)
Alpaca_Markers_clean_Cluster_3 <- Alpaca_Markers_clean_Cluster_3 %>% dplyr::select(2,)
Alpaca_Markers_clean_Cluster_4 <- Alpaca_Markers_clean_Cluster_4 %>% dplyr::select(2,)
Alpaca_Markers_clean_Cluster_5 <- Alpaca_Markers_clean_Cluster_5 %>% dplyr::select(2,)
Alpaca_Markers_clean_Cluster_6 <- Alpaca_Markers_clean_Cluster_6 %>% dplyr::select(2,)
Alpaca_Markers_clean_Cluster_7 <- Alpaca_Markers_clean_Cluster_7 %>% dplyr::select(2,)
Alpaca_Markers_clean_Cluster_8 <- Alpaca_Markers_clean_Cluster_8 %>% dplyr::select(2,)
Alpaca_Markers_clean_Cluster_9 <- Alpaca_Markers_clean_Cluster_9 %>% dplyr::select(2,)
Alpaca_Markers_clean_Cluster_10 <- Alpaca_Markers_clean_Cluster_10 %>% dplyr::select(2,)



setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_both")
write.table(Alpaca_Markers_clean_Cluster_0,"Alpaca_Cluster_0.txt", row.names = FALSE)
write.table(Alpaca_Markers_clean_Cluster_1,"Alpaca_Cluster_1.txt", row.names = FALSE)
write.table(Alpaca_Markers_clean_Cluster_2,"Alpaca_Cluster_2.txt", row.names = FALSE)
write.table(Alpaca_Markers_clean_Cluster_3,"Alpaca_Cluster_3.txt", row.names = FALSE)
write.table(Alpaca_Markers_clean_Cluster_4,"Alpaca_Cluster_4.txt", row.names = FALSE)
write.table(Alpaca_Markers_clean_Cluster_5,"Alpaca_Cluster_5.txt", row.names = FALSE)
write.table(Alpaca_Markers_clean_Cluster_6,"Alpaca_Cluster_6.txt", row.names = FALSE)
write.table(Alpaca_Markers_clean_Cluster_7,"Alpaca_Cluster_7.txt", row.names = FALSE)
write.table(Alpaca_Markers_clean_Cluster_8,"Alpaca_Cluster_8.txt", row.names = FALSE)
write.table(Alpaca_Markers_clean_Cluster_9,"Alpaca_Cluster_9.txt", row.names = FALSE)
write.table(Alpaca_Markers_clean_Cluster_10,"Alpaca_Cluster_10.txt", row.names = FALSE)


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

DefaultAssay(Alpaca.clustered) <- "RNA"


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/umap_legend.pdf",
    width=5, height=5)
DimPlot(Alpaca.clustered, reduction = "umap", label = T, label.size = 5, repel = TRUE, pt.size = 0.5) + NoAxes()
dev.off()

Alpaca.clustered <- Alpaca
p1 <- FeaturePlot(Alpaca.clustered, slot = "data", features = "LOXL4") + NoAxes() + NoLegend() #Basal
p2 <- FeaturePlot(Alpaca.clustered, slot = "data", features = "MUC21") + NoAxes() + NoLegend() #early Club
p3 <- FeaturePlot(Alpaca.clustered, slot = "data", features = "BPIFA1") + NoAxes() + NoLegend() #Secretory
p4 <- FeaturePlot(Alpaca.clustered, slot = "data", features = "WDR66") + NoAxes() + NoLegend() #Ciliated and pre-ciliated


p5 <- VlnPlot(Alpaca.clustered, cols = cols, slot = "data", features = "LOXL4", pt.size = 0) + NoLegend()
p6 <- VlnPlot(Alpaca.clustered, cols = cols, slot = "data", features = "MUC21", pt.size = 0) + NoLegend()
p7 <- VlnPlot(Alpaca.clustered, cols = cols, slot = "data", features = "BPIFA1", pt.size = 0) + NoLegend()
p8 <- VlnPlot(Alpaca.clustered, cols = cols, slot = "data", features = "WDR66", pt.size = 0) + NoLegend()


FeaturePlot(Alpaca.clustered, slot = "data", features = "LOXL4") + NoAxes()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/Umaps_genes.pdf",
    width=3, height=12)
grid.arrange(p1, p2, p3, p4, ncol = 1)
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/Expression_violin_genes.pdf",
    width=3, height=15)
grid.arrange(p5, p6, p7, p8,  ncol = 1)
dev.off()

?FeaturePlot


FeaturePlot(Alpaca.clustered, slot = "data", features = "BPIFA1") + NoAxes()  #secretory
FeaturePlot(Alpaca.clustered, slot = "data", features = "KRT4") + NoAxes() + NoLegend() #Basal
FeaturePlot(Alpaca.clustered, slot = "data", features = "LOXL4") + NoAxes() + NoLegend() #Basal
FeaturePlot(Alpaca.clustered, slot = "data", features = "TP63") + NoAxes() + NoLegend() #Basal
FeaturePlot(Alpaca.clustered, slot = "data", features = "SCGB1A1") + NoAxes() + NoLegend() #Club
FeaturePlot(Alpaca.clustered, slot = "data", features = "LOC102527828") + NoAxes() + NoLegend() #KRT14
FeaturePlot(Alpaca.clustered, slot = "data", features = "LOC102526268") + NoAxes() + NoLegend() #KRT13
FeaturePlot(Alpaca.clustered, slot = "data", features = "KRT7") + NoAxes() + NoLegend() #KRT7
FeaturePlot(Alpaca.clustered, slot = "data", features = "MKI67") + NoAxes() + NoLegend() #Cycling
FeaturePlot(Alpaca.clustered, slot = "data", features = "LYN") + NoAxes() + NoLegend() #Secretory
FeaturePlot(Alpaca.clustered, slot = "data", features = "MUC21") + NoAxes() + NoLegend() #Club
FeaturePlot(Alpaca.clustered, slot = "data", features = "MSMB") + NoAxes() + NoLegend() #Secretory
FeaturePlot(Alpaca.clustered, slot = "data", features = "SYNE1") + NoAxes() + NoLegend() #Ciliated
FeaturePlot(Alpaca.clustered, slot = "data", features = "TPPP3") + NoAxes() + NoLegend() #Ciliated
FeaturePlot(Alpaca.clustered, slot = "data", features = "KRT17") + NoAxes()  #basal
FeaturePlot(Alpaca.clustered, slot = "data", features = "MT1X") + NoAxes()  #subrabasal
FeaturePlot(Alpaca.clustered, slot = "data", features = "KRT6A") + NoAxes()  #subrabasal
FeaturePlot(Alpaca.clustered, slot = "data", features = "KRT16") + NoAxes() #subrabasal
FeaturePlot(Alpaca.clustered, slot = "data", features = "CXCL8") + NoAxes()  #subrabasal
FeaturePlot(Alpaca.clustered, slot = "data", features = "KRT7") + NoAxes()  #Club
FeaturePlot(Alpaca.clustered, slot = "data", features = "KRT4") + NoAxes() #Club
FeaturePlot(Alpaca.clustered, slot = "data", features = "KRT13") + NoAxes()  #Club
FeaturePlot(Alpaca.clustered, slot = "data", features = "SCGB1A1") + NoAxes()  #Club
FeaturePlot(Alpaca.clustered, slot = "data", features = "MSMB") + NoAxes()  #secretory

FeaturePlot(Alpaca.clustered, slot = "data", features = "MKI67") + NoAxes()  #cycling
FeaturePlot(Alpaca.clustered, slot = "data", features = "TOP2A") + NoAxes()  #cycling
FeaturePlot(Alpaca.clustered, slot = "data", features = "MUC5AC") + NoAxes()  #secretory CEACAM6
FeaturePlot(Alpaca.clustered, slot = "data", features = "FAM111B") + NoAxes()  #secretory CEACAM6
FeaturePlot(Alpaca.clustered, slot = "data", features = "CDC20B") + NoAxes()  #secretory CEACAM6

head(Alpaca.clustered@meta.data)
Idents(Alpaca.clustered)
levels(Alpaca.clustered)


Idents(Alpaca.clustered) <- "Treat" 
Alpaca_MERS <- subset(x = Alpaca.clustered, idents = c("Mock", "MERS-CoV"))
Alpaca_ACN4 <- subset(x = Alpaca.clustered, idents = c("Mock", "dcCoV-ACN4"))
Idents(Alpaca_ACN4) <- "celltype" 
Idents(Alpaca_MERS) <- "celltype" 
 


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/Vln_Plots_MERS_Status_cilia.pdf",
    width=17, height=10)
VlnPlot(Alpaca_MERS, cols = cols_stat_2, slot = "data", features = c("DNAH7", "WDR66", "LRRC6","LRRIQ1", "TPPP3", "FOXJ1", "TUBB"), split.by = 'Status', pt.size = 0) +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)),
                                                                                                                                             legend.position = "top")
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/Vln_Plots_ACN4_Status_cilia.pdf",
    width=17, height=10)
VlnPlot(Alpaca_ACN4, cols = cols_stat_2, slot = "data", features = c("DNAH7", "WDR66", "LRRC6","LRRIQ1", "TPPP3", "FOXJ1", "TUBB"), split.by = 'Status', pt.size = 0) +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)),
                                                                                                                                                                                  legend.position = "top")
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/Vln_Plots_treat_cilia.pdf",
    width=17, height=10)
VlnPlot(Alpaca.clustered, cols = cols_treat, slot = "data", features = c("DNAH7","WDR66", "LRRC6","LRRIQ1", "TPPP3", "FOXJ1", "TUBB"), split.by = 'Treat', pt.size = 0) +theme(axis.title.y = element_text(hjust = 0.5, vjust = -0.5, margin = margin(r = 25)),
                                                                                                                                                                         legend.position = "top")
dev.off()

?DotPlot

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/Dot_Plots_RECEPTORS.pdf",
    width=6, height=4.2)
DotPlot(Alpaca.clustered, features = c("DPP4", "ANPEP", "ACE2", "camel229E-3UTR"), scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


FeaturePlot(Alpaca.clustered, slot = "data", features = "KRT18") + NoAxes()

DotPlot(Alpaca.clustered, features = c("DNAH5", "SPEF2", "PIFO", "LRRC6", "TPPP3", "FOXJ1", "TUBB", "TUBB1","TUBB2A",
                                                     "TUBB4A","TUBB4B","TUBB6","TUBA4A","TUBA1A"), scale.by = "radius", scale = T, cols = c("orange", "darkblue", "red")) +coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

sec_dotplot <- DotPlot(Alpaca.clustered, features = c("FOXQ1", "SPDEF", "SCGB1A1", "MSMB", "MUC1", "MUC2", "MUC3A", 
                                                     "MUC4", "MUC5AC","MUC5B", "MUC6", "MUC7", "MUC15", "MUC19", "MUC20"), scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

deut_dotplot <- DotPlot(Alpaca.clustered, features = c("PLK4", "CCNO", "CEP78", "DEUP1"), scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/Supplemenatary_fig_1.pdf",
    width=6, height=7)
cil_dotplot/sec_dotplot
dev.off()


DotPlot(Alpaca.clustered, features = c("LOXL4", "MKI67","TOP2A","COL7A1", "TNC", "LOC102527556", "KRT19", "HES1",
                                       "KRT5", "TIMP3", "TP63", "LOC102527828", "CCN1",
                                       "ID1","KRT4", "LOC102526268", "MUC21","KRT7", "BPIFA1", "KRT23",
                                       "KRT8", "PLK4", "SYNE1","CFTR", "LRRIQ1", "TPPP3",
                                       "SCGB1A1","MSMB", "VMO1", "LOC116279189"), scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


DefaultAssay(Alpaca.clustered) <- "SCT"
DoHeatmap(Alpaca.clustered, features = c("MUC4", "KRT4", "LOC102541537", "LRRIQ1", "LOC102545839","DYNLRB2", "DMBT1", "FN1","KRT5",
                                         "KRT7","KRT4","SCGB1A1","MSMB", "BPIFA1", "SPDEF", "KRT8", "MUC5AC", "FOXN4", "PLK4", "SYNE1", "CCNO", 
                                         "CFTR","MKI67","TOP2A"))


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/heatmap_cluster_first5.pdf",
    width=4, height=5)
DoHeatmap(Alpaca.clustered, features = Mark)
dev.off()
##################################################################################################################



##################################################################################################################
#Renaming the found clusters to the respective cell types for downstream analysis differential expression analysis
Idents(Alpaca.clustered) <- Alpaca.clustered$seurat_clusters
new.cluster.ids <- c("Club", "Ciliated", "Basal", "Secretory", "Club", "Basal", "Ciliated", "Secretory", "Club", "Llama Cluster 1", "Llama Cluster 2")
names(new.cluster.ids) <- levels(Alpaca.clustered)
levels(Alpaca.clustered)
Idents(Alpaca.clustered)
Alpaca.clustered <- RenameIdents(Alpaca.clustered, new.cluster.ids)
Alpaca.clustered$celltype <- Idents(Alpaca.clustered) #first make the celltype a column entry
Idents(Alpaca.clustered) <- factor(x = Idents(Alpaca.clustered), levels = c("Llama Cluster 1", "Llama Cluster 2", "Secretory", "Ciliated", "Club", "Basal"))

Alpaca.clustered@meta.data

Alpaca.clustered <- Alpaca
Alpaca.clustered_markers <- FindAllMarkers(Alpaca.clustered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Alpaca.clustered_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

levels(Alpaca.clustered)

cols = c("#8A6E7E", "#40679B", "#089392", "#40B48B", "#9DCD84", "#EAE29C")
show_col(cols)

Clus_1_vs_basal_markers <- FindMarkers(Alpaca.clustered, only.pos = TRUE, ident.1 = "Llama Cluster 1", ident.2 = "Basal", min.pct = 0.25, logfc.threshold = 0.25)

?FindMarkers

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/umap_cluster_named_label.pdf",
    width=5, height=5)
DimPlot(Alpaca.clustered, reduction = "umap", label.size = 6, repel = TRUE, label = TRUE, pt.size = 0.5, cols = cols) + NoLegend() + NoAxes()
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/umap_group_treat.pdf",
    width=8, height=4)
DimPlot(Alpaca.clustered, reduction = "umap", pt.size = 0.5, split.by = "Treat", cols = cols) + NoLegend()
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/umap_cluster_named_legend.pdf",
    width=8, height=7)
DimPlot(Alpaca.clustered, reduction = "umap", label.size = 6, repel = TRUE, label = F, pt.size = 0.5, cols = cols) + NoAxes()
dev.off()



DefaultAssay(Alpaca)
DefaultAssay(Alpaca.clustered) <- "integrated"

Alpaca.clustered <- Alpaca
Alpaca.clustered_markers <- FindAllMarkers(Alpaca.clustered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#Markers for new clusters and dotplot
markers <- Alpaca.clustered_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  print(n = 50)




#Make dotplot with the new annotated clusters and their gene sets
#choose only the unique gene names
list_2 <- unique(markers$gene)


#Make plot
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/dot_plot_new_clusters.pdf",
    width=4, height=9)
DotPlot(Alpaca.clustered, features = list_2, scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                                                                                                                                  legend.position = "top")
dev.off()



#Save the gene lists of the newly annotated clusters to compare with human clusters
Alpaca_Markers_clean <- Alpaca.clustered_markers %>% dplyr::select(6:7,)

Alpaca_Markers_clean_Cluster_0 <- Alpaca_Markers_clean[Alpaca_Markers_clean$cluster == "Secretory", ]
Alpaca_Markers_clean_Cluster_1 <- Alpaca_Markers_clean[Alpaca_Markers_clean$cluster == "Ciliated", ]
Alpaca_Markers_clean_Cluster_2 <- Alpaca_Markers_clean[Alpaca_Markers_clean$cluster == "Club", ]
Alpaca_Markers_clean_Cluster_3 <- Alpaca_Markers_clean[Alpaca_Markers_clean$cluster == "Basal", ]
Alpaca_Markers_clean_Cluster_4 <- Alpaca_Markers_clean[Alpaca_Markers_clean$cluster == "Llama Cluster 1", ]
Alpaca_Markers_clean_Cluster_5 <- Alpaca_Markers_clean[Alpaca_Markers_clean$cluster == "Llama Cluster 2", ]

# Alpaca_Markers_clean_Cluster_0 <- Alpaca_Markers_clean_Cluster_0 %>% dplyr::select(2,)
# Alpaca_Markers_clean_Cluster_1 <- Alpaca_Markers_clean_Cluster_1 %>% dplyr::select(2,)
# Alpaca_Markers_clean_Cluster_2 <- Alpaca_Markers_clean_Cluster_2 %>% dplyr::select(2,)
# Alpaca_Markers_clean_Cluster_3 <- Alpaca_Markers_clean_Cluster_3 %>% dplyr::select(2,)
# Alpaca_Markers_clean_Cluster_4 <- Alpaca_Markers_clean_Cluster_4 %>% dplyr::select(2,)
# Alpaca_Markers_clean_Cluster_5 <- Alpaca_Markers_clean_Cluster_5 %>% dplyr::select(2,)
# 

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Cluster_markers/Llama_both")
write.csv(Alpaca.clustered_markers, "Cluster_markers_Llama.csv", row.names = F)
write.table(Alpaca_Markers_clean_Cluster_0,"Alpaca_Cluster_0.txt", row.names = FALSE)
write.table(Alpaca_Markers_clean_Cluster_1,"Alpaca_Cluster_1.txt", row.names = FALSE)
write.table(Alpaca_Markers_clean_Cluster_2,"Alpaca_Cluster_2.txt", row.names = FALSE)
write.table(Alpaca_Markers_clean_Cluster_3,"Alpaca_Cluster_3.txt", row.names = FALSE)
write.table(Alpaca_Markers_clean_Cluster_4,"Alpaca_Cluster_4.txt", row.names = FALSE)
write.table(Alpaca_Markers_clean_Cluster_5,"Alpaca_Cluster_5.txt", row.names = FALSE)








#Saving the seurat object to disk to access in other script for analysis
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/V.pacos")
SaveH5Seurat(Alpaca.clustered, "Alpaca.clustered", overwrite = TRUE)



#########################################################################################################################################################

DefaultAssay(Alpaca.clustered) <- "RNA"
Idents(Alpaca.clustered) <- "celltype" 

Ferus_MERS <- subset(x = Alpaca.clustered, idents = c("Ciliated"))

levels(Ferus_MERS)
