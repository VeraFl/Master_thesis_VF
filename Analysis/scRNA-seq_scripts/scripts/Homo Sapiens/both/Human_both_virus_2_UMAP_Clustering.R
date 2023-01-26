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

#Treatment colors (HCoV-229E, MERS-CoV, Mock)
cols_treat = c("#84A59D","#A5668B", "#96ADC8")
show_col(cols_treat)
#Status colors
cols_stat = c("#DEA47E","#CD4631","#DEA47E")
show_col(cols_stat)

cols_stat_2 = c("#519872", "#CD4631","#DEA47E")
show_col(cols_stat_2)

cols_stat_2 = c("#519872", "#CD4631","#DEA47E")
show_col(cols_stat_2)

#Celltype colors
cols_cell = hcl.colors(7, "Temps")
show_col(cols_cell)
cols_cell


viridis(20)
#Secretory = "#089392"
#Ciliated = "#40B48B"
#Club = "#9DCD84"
#Basal = "#EAE29C"
#Suprabasal = "#EAB672"
#Deuterosomal = "#E6866A"
#Ionocytes = "#CF597E"

#colors from figure of cell development trajectory
cols_cell_2 = c("#2E9288", "#F9D57F", "#AA9AB4", "#C19BA6", "#B4AACD", "#DEA6A5")
memory.limit(24000)
##########################################################################################################################
#Loading filtered and merged Seurat Object into Script
#Infected and Mock Samples merged

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens")

Human.clustered <- LoadH5Seurat("Human.combined.h5seurat") #Import either seurat object
Human.clustered <- LoadH5Seurat("Human.clustered.h5seurat") #Import either seurat object
head(Human.clustered@meta.data)
str(Human.clustered@meta.data)
##########################################################################################################################
#Change Assay to integrated to make UMAP
DefaultAssay(Human.clustered) <- "integrated"
# Run the standard workflow for visualization and clustering
Human.clustered <- RunPCA(Human.clustered, npcs = 30, verbose = FALSE) #Run PCA as dimensionality reduction

#Using the chosen
Human.clustered <- FindNeighbors(Human.clustered, reduction = "pca", dims = 1:30)
Human.clustered <- FindClusters(Human.clustered, resolution = 0.6)
Human.clustered <- RunUMAP(Human.clustered, reduction = "pca", dims = 1:30)


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/umap_label.pdf",
    width=5, height=5)
DimPlot(Human.clustered, reduction = "umap", label = TRUE, label.size = 6) + NoLegend() +NoAxes()
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/umap_legend.pdf",
    width=5, height=5)
DimPlot(Human.clustered, reduction = "umap", label = F, label.size = 6,) +NoAxes()
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/umap_split_treat_status.pdf",
    width=8, height=4)
DimPlot(Human.clustered, reduction = "umap", split.by = "Treat", group.by = "Status", cols = cols_stat) + NoLegend()
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/umap_split_treat.pdf",
    width=8, height=4)
DimPlot(Human.clustered, reduction = "umap", split.by = "Treat") + NoLegend() +NoAxes()
dev.off()


##############################################################################################

##################################################################################################################

# Markers to find cell type clusters
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
Human.clustered_markers_wilcox <- FindAllMarkers(Human.clustered, only.pos = TRUE, test.use="wilcox", min.pct = 0.25, logfc.threshold = 0.25)
markers <- Human.clustered_markers_wilcox %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)%>%
  print(n = 60)

list <- unique(markers$gene)

#heatmap for fresh unannotated clusters
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/Heatmap_initial_clusters.pdf",
    width=10, height=7)
DoHeatmap(Human.clustered, features = list)
dev.off()

?DoHeatmap()

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_Human")
write.table(Human_Human_Markers_clean,"Human_Human_Markers_clean.txt", row.names = FALSE)

Human_Markers_clean <- Human.clustered_markers %>% select(6:7,)

Human_Markers_clean_Cluster_0 <- Human_Markers_clean[Human_Markers_clean$cluster == 0, ]
Human_Markers_clean_Cluster_1 <- Human_Markers_clean[Human_Markers_clean$cluster == 1, ]
Human_Markers_clean_Cluster_2 <- Human_Markers_clean[Human_Markers_clean$cluster == 2, ]
Human_Markers_clean_Cluster_3 <- Human_Markers_clean[Human_Markers_clean$cluster == 3, ]
Human_Markers_clean_Cluster_4 <- Human_Markers_clean[Human_Markers_clean$cluster == 4, ]
Human_Markers_clean_Cluster_5 <- Human_Markers_clean[Human_Markers_clean$cluster == 5, ]
Human_Markers_clean_Cluster_6 <- Human_Markers_clean[Human_Markers_clean$cluster == 6, ]
Human_Markers_clean_Cluster_7 <- Human_Markers_clean[Human_Markers_clean$cluster == 7, ]
Human_Markers_clean_Cluster_8 <- Human_Markers_clean[Human_Markers_clean$cluster == 8, ]


Human_Markers_clean_Cluster_0 <- Human_Markers_clean_Cluster_0 %>% select(2,)
Human_Markers_clean_Cluster_1 <- Human_Markers_clean_Cluster_1 %>% select(2,)
Human_Markers_clean_Cluster_2 <- Human_Markers_clean_Cluster_2 %>% select(2,)
Human_Markers_clean_Cluster_3 <- Human_Markers_clean_Cluster_3 %>% select(2,)
Human_Markers_clean_Cluster_4 <- Human_Markers_clean_Cluster_4 %>% select(2,)
Human_Markers_clean_Cluster_5 <- Human_Markers_clean_Cluster_5 %>% select(2,)
Human_Markers_clean_Cluster_6 <- Human_Markers_clean_Cluster_6 %>% select(2,)
Human_Markers_clean_Cluster_7 <- Human_Markers_clean_Cluster_7 %>% select(2,)
Human_Markers_clean_Cluster_8 <- Human_Markers_clean_Cluster_8 %>% select(2,)


setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Human_both")
write.table(Human_Markers_clean_Cluster_0,"Human_Cluster_0.txt", row.names = FALSE)
write.table(Human_Markers_clean_Cluster_1,"Human_Cluster_1.txt", row.names = FALSE)
write.table(Human_Markers_clean_Cluster_2,"Human_Cluster_2.txt", row.names = FALSE)
write.table(Human_Markers_clean_Cluster_3,"Human_Cluster_3.txt", row.names = FALSE)
write.table(Human_Markers_clean_Cluster_4,"Human_Cluster_4.txt", row.names = FALSE)
write.table(Human_Markers_clean_Cluster_5,"Human_Cluster_5.txt", row.names = FALSE)
write.table(Human_Markers_clean_Cluster_6,"Human_Cluster_6.txt", row.names = FALSE)
write.table(Human_Markers_clean_Cluster_7,"Human_Cluster_7.txt", row.names = FALSE)
write.table(Human_Markers_clean_Cluster_8,"Human_Cluster_8.txt", row.names = FALSE)



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
#basal = "#EAE29C"
#suprabasal = "#EAB672"
#Club = "#9DCD84"
#Deuterosomal = "#E6866A"
#Goblet/Secretory = "#089392"
#Ciliated = "#40B48B"
#Ionocytes = "#CF597E"


DefaultAssay(Human.clustered) <- "RNA"
# p1 <- DimPlot(Human.clustered, reduction = "umap", pt.size = 0.5, cols = cols_cell) + NoAxes() + NoLegend()
# p2 <- FeaturePlot(Human.clustered, features = c("FN1"), cols = c("lightgrey", "#EAE29C")) + NoAxes() + NoLegend()
# p3 <- FeaturePlot(Human.clustered, features = c("MT1X"), cols = c("lightgrey", "#EAB672")) + NoAxes() + NoLegend()
# p4 <- FeaturePlot(Human.clustered, features = c("KRT7"), cols = c("lightgrey", "#9DCD84")) + NoAxes() + NoLegend()
# p5 <- FeaturePlot(Human.clustered, features = c("MSMB"), cols = c("lightgrey", "#089392")) + NoAxes() + NoLegend()
# p6 <- FeaturePlot(Human.clustered, features = c("TPPP3"), cols = c("lightgrey", "#40B48B")) + NoAxes() + NoLegend()
# p7 <- FeaturePlot(Human.clustered, features = c("CCNO"), cols = c("lightgrey", "#E6866A")) + NoAxes() + NoLegend()
# p8 <- FeaturePlot(Human.clustered, features = c("CFTR"), cols = c("lightgrey", "#CF597E")) + NoAxes() + NoLegend()
# 
# grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 2)
?FeaturePlot
str(Human.clustered)
DimPlot(Human.clustered, reduction = "umap", label = TRUE, label.size = 5, repel = TRUE, pt.size = 0.5) + NoAxes() + NoLegend()





DefaultAssay(Human.clustered) <- "SCT"
FeaturePlot(Human.clustered, slot = "data", features = "FN1") + NoAxes()  #basal
FeaturePlot(Human.clustered, slot = "data", features = "KRT5") + NoAxes()  #basal
FeaturePlot(Human.clustered, slot = "data", features = "KRT17") + NoAxes()  #basal
FeaturePlot(Human.clustered, slot = "data", features = "MT1X") + NoAxes()  #subrabasal
FeaturePlot(Human.clustered, slot = "data", features = "KRT6A") + NoAxes()  #subrabasal
FeaturePlot(Human.clustered, slot = "data", features = "KRT16") + NoAxes() #subrabasal
FeaturePlot(Human.clustered, slot = "data", features = "CXCL8") + NoAxes()  #subrabasal
FeaturePlot(Human.clustered, slot = "data", features = "KRT7") + NoAxes()  #Club
FeaturePlot(Human.clustered, slot = "data", features = "KRT4") + NoAxes() #Club
FeaturePlot(Human.clustered, slot = "data", features = "KRT13") + NoAxes()  #Club
FeaturePlot(Human.clustered, slot = "data", features = "SCGB1A1") + NoAxes()  #Club
FeaturePlot(Human.clustered, slot = "data", features = "MSMB") + NoAxes()  #secretory
FeaturePlot(Human.clustered, slot = "data", features = "BPIFA1") + NoAxes()  #secretory
FeaturePlot(Human.clustered, slot = "data", features = "CEACAM6") + NoAxes()  #secretory
FeaturePlot(Human.clustered, slot = "data", features = "KRT8") + NoAxes()  #secretory
FeaturePlot(Human.clustered, slot = "data", features = "MUC5AC") + NoAxes()  #secretory CEACAM6
FeaturePlot(Human.clustered, slot = "data", features = "C5orf49") + NoAxes()  #early ciliated
FeaturePlot(Human.clustered, slot = "data", features = "FOXN4") + NoAxes()  #early ciliated
FeaturePlot(Human.clustered, slot = "data", features = "PLK4") + NoAxes()  #early ciliated
FeaturePlot(Human.clustered, slot = "data", features = "SYNE1") + NoAxes()  #later ciliated
FeaturePlot(Human.clustered, slot = "data", features = "SNTN") + NoAxes()  #mature ciliated
FeaturePlot(Human.clustered, slot = "data", features = "CCNO") + NoAxes()  #early ciliated/deuterosomal
FeaturePlot(Human.clustered, slot = "data", features = "CFTR") + NoAxes()  #Ionocytes

FeaturePlot(Human.clustered, slot = "data", features = "MKI67") + NoAxes()  #cycling
FeaturePlot(Human.clustered, slot = "data", features = "TOP2A") + NoAxes()  #cycling
FeaturePlot(Human.clustered, slot = "data", features = "MUC5AC") + NoAxes()  #secretory CEACAM6
FeaturePlot(Human.clustered, slot = "data", features = "FAM111B") + NoAxes()  #secretory CEACAM6
FeaturePlot(Human.clustered, slot = "data", features = "CDC20B") + NoAxes()  #secretory CEACAM6


FeaturePlot(Human.clustered, slot = "data", features = "KRT5") + NoAxes()  
FeaturePlot(Human.clustered, slot = "data", features = "TP63") + NoAxes()  
FeaturePlot(Human.clustered, slot = "data", features = "SCGB1A1") + NoAxes()  
FeaturePlot(Human.clustered, slot = "data", features = "MUC5AC") + NoAxes() 
FeaturePlot(Human.clustered, slot = "data", features = "FOXJ1") + NoAxes() 
FeaturePlot(Human.clustered, slot = "data", features = "DEUP1") + NoAxes()
FeaturePlot(Human.clustered, slot = "data", features = "CFTR") + NoAxes() 


FeaturePlot(Human.clustered, slot = "data", features = "RFX3")


Idents(Human.clustered) <- "celltype" 
Human.clustered_ciliated <- subset(x = Human.clustered, idents = c("Ciliated"))

p1 <- DotPlot(Human.clustered, features = c("FN1","KRT5","KRT14","KRT17","MT1X","KRT6A","KRT16","KRT7","KRT4",
                                      "KRT13","SCGB1A1","MSMB","BPIFA1","CEACAM6", "SPDEF","MUC5B","FOXA2", 
                                      "CYP26B1", "KRT8", "MUC5AC", "C5orf49", "FOXN4", "PLK4", "SYNE1", 
                                      "SNTN", "CCNO", "CFTR","MKI67","TOP2A"), scale.by = "radius", idents = "Ciliated", scale = F, cols = c("orange", "darkblue"), group.by = "Treat") + coord_flip() + theme(axis.text.x = element_text(angle = 90, size = 12)) +NoLegend()

p2 <- DotPlot(Human.clustered, features = c("FN1","KRT5","KRT14","KRT17","MT1X","KRT6A","KRT16","KRT7","KRT4",
                                            "KRT13","SCGB1A1","MSMB","BPIFA1","CEACAM6", "SPDEF","MUC5B","FOXA2", 
                                            "CYP26B1", "KRT8", "MUC5AC", "C5orf49", "FOXN4", "PLK4", "SYNE1", 
                                            "SNTN", "CCNO", "CFTR","MKI67","TOP2A"), scale.by = "radius", idents = "Secretory", scale = F, cols = c("orange", "darkblue"), group.by = "Treat") + coord_flip() + theme(axis.text.x = element_text(angle = 90, size = 12)) +NoLegend()

p3 <- DotPlot(Human.clustered, features = c("FN1","KRT5","KRT14","KRT17","MT1X","KRT6A","KRT16","KRT7","KRT4",
                                            "KRT13","SCGB1A1","MSMB","BPIFA1","CEACAM6", "SPDEF","MUC5B","FOXA2", 
                                            "CYP26B1", "KRT8", "MUC5AC", "C5orf49", "FOXN4", "PLK4", "SYNE1", 
                                            "SNTN", "CCNO", "CFTR","MKI67","TOP2A"), scale.by = "radius", idents = "Ionocytes", scale = F, cols = c("orange", "darkblue"), group.by = "Treat") + coord_flip() + title("Secretory") + theme(axis.text.x = element_text(angle = 90, size = 12),
                                                                                                                                                                                                                      axis.text.y = element_blank(),
                                                                                                                                                                                                                      axis.ticks.y = element_blank(),
                                                                                                                                                                                                                      axis.title.x = element_blank(),
                                                                                                                                                                                                                      title = element_text(size = 12)) +NoLegend()


p1 + p2 + p3



Idents(Human.clustered) <- "Treat" 
Human.clustered_mock <- subset(x = Human.clustered, idents = c("Mock"))

DotPlot(Human.clustered, features = c("FN1","KRT5","KRT14","KRT17","MT1X","KRT6A","KRT16","KRT7","KRT4",
                                               "KRT13","SCGB1A1","MSMB","BPIFA1","CEACAM6", "SPDEF","MUC5B","FOXA2", 
                                               "CYP26B1", "KRT8", "MUC5AC", "C5orf49", "FOXN4", "PLK4", "SYNE1", 
                                               "SNTN", "CCNO", "CFTR","MKI67","TOP2A"), scale.by = "radius", scale = T, cols=, split.by = "Treat", group.by = "celltype") + coord_flip() + theme(axis.text = element_text(angle = 25, size = 12))




?DotPlot()
p1 <- FeaturePlot(Human.clustered, slot = "data", features = "KRT5") + NoAxes()  + NoLegend()#basal
p2 <- FeaturePlot(Human.clustered, slot = "data", features = "MT1X") + NoAxes()  + NoLegend()#suprabasal
p3 <- FeaturePlot(Human.clustered, slot = "data", features = "SCGB1A1") + NoAxes()  + NoLegend()#Club
p4 <- FeaturePlot(Human.clustered, slot = "data", features = "CEACAM6") + NoAxes()  + NoLegend()#Secretory
p5 <- FeaturePlot(Human.clustered, slot = "data", features = "CCNO") + NoAxes()  + NoLegend()#early ciliated
p6 <- FeaturePlot(Human.clustered, slot = "data", features = "TPPP3") + NoAxes()  + NoLegend()#later ciliated
p7 <- FeaturePlot(Human.clustered, slot = "data", features = "CFTR") + NoAxes()  + NoLegend()#Ionocytes

p8 <- VlnPlot(Human.clustered, cols = cols_cell, slot = "data", features = "KRT5", pt.size = 0) + NoLegend()
p9 <- VlnPlot(Human.clustered, cols = cols_cell, slot = "data", features = "MT1X", pt.size = 0) + NoLegend()
p10 <- VlnPlot(Human.clustered, cols = cols_cell, slot = "data", features = "SCGB1A1", pt.size = 0) + NoLegend()
p11 <- VlnPlot(Human.clustered, cols = cols_cell, slot = "data", features = "CEACAM6", pt.size = 0) + NoLegend()
p12 <- VlnPlot(Human.clustered, cols = cols_cell, slot = "data", features = "CCNO", pt.size = 0) + NoLegend()
p13 <- VlnPlot(Human.clustered, cols = cols_cell, slot = "data", features = "TPPP3", pt.size = 0) + NoLegend()
p14 <- VlnPlot(Human.clustered, cols = cols_cell, slot = "data", features = "CFTR", pt.size = 0) + NoLegend()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/Umaps_genes.pdf",
    width=3, height=20)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, ncol = 1)
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/Expression_violin_genes.pdf",
    width=3, height=20)
grid.arrange(p8, p9, p10, p11, p12, p13, p14, ncol = 1)
dev.off()

?FeaturePlot

head(Human.clustered@meta.data)


#Make barplot about cluster proportion
Human_clustered <- Human.clustered@meta.data %>%  dplyr::select(c("Status", "Treat", "celltype"))


#Grouping by celltype and treatment allos wo count the cells in each celltype per treat
Human_celltype <- Human_clustered %>%
  group_by(celltype, Treat) %>%
  dplyr::count(celltype)

Plot_X <- Human_celltype %>%
  group_by(Treat) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            Treat = Treat) %>%
  mutate("Percentage" = n/Sum*100)


# Reordering Levels
Plot_X$Treat <- factor(Plot_X$Treat, levels=c('Mock','HCoV-229E','MERS-CoV'))

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/fraction_of_cells_per_celltype.pdf",
    width=6, height=5)
Plot_X %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = Treat, y = n, fill = celltype))+
  geom_bar(position = "fill", stat = "identity", width = 0.8, color = "#4F4A4B") +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_cell)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="right")
dev.off()

##################################################################################################################
#Renaming the found clusters to the respective cell types for downstream analysis differential expression analysis
new.cluster.ids <- c("Secretory", "Ciliated", "Secretory", "Basal", "Secretory", "Ciliated", "Secretory", "Club", "Suprabasal", "Ciliated", "Suprabasal", "Ciliated", "Suprabasal", "Deuterosomal", "Ionocytes")
names(new.cluster.ids) <- levels(Human.clustered)
Human.clustered <- RenameIdents(Human.clustered, new.cluster.ids)
Human.clustered$celltype <- Idents(Human.clustered)
Idents(Human.clustered)
levels(Human.clustered)
head(Human.clustered@meta.data)

Idents(Human.clustered) <- factor(x = Idents(Human.clustered), levels = c('Secretory', 'Ciliated', 'Club', 'Basal', 'Suprabasal', 'Deuterosomal', 'Ionocytes'))


cols_cell = hcl.colors(7, "Temps")
show_col(cols_cell)
cols_cell

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/umap_label.pdf",
    width=5, height=5)
DimPlot(Human.clustered, reduction = "umap", label = T, label.size = 6, repel = TRUE, pt.size = 0.5, cols = cols_cell) +NoAxes() +NoLegend()
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/umap_legend.pdf",
    width=6, height=5)
DimPlot(Human.clustered, reduction = "umap", label = F, label.size = 6, repel = TRUE, pt.size = 0.5, cols = cols_cell) +NoAxes()
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/umap_split_treat.pdf",
    width=8, height=4)
DimPlot(Human.clustered, reduction = "umap", split.by = "Treat", cols = cols_cell) + NoAxes() + NoLegend()
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/dot_plot.pdf",
    width=8, height=5)
DotPlot(Human.clustered, features = c("FN1","KRT5","KRT14","KRT17","MT1X","KRT6A","KRT7","KRT4","SCGB1A1","MSMB",
                                      "BPIFA1","CEACAM6", "MUC5B", "KRT8", "C5orf49", "FOXN4", "PLK4", "SYNE1", "SNTN", "CCNO", 
                                      "CFTR","MKI67","TOP2A"), scale.by = "radius", scale = T, cols = c("orange", "darkblue"))  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


?DotPlot

#Finding markers per assigned cell group
Human.celltype_markers <- FindAllMarkers(Human.clustered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- Human.celltype_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  print(n = 50)


#Make dotplot with the new annotated clusters and their gene sets
#choose only the unique gene names
list_2 <- unique(markers$gene)

#Make plot
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/dot_plot_new_clusters.pdf",
    width=6, height=9)
DotPlot(Human.clustered, features = list_2, scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


#Save the gene lists of the newly annotated clusters to compare with human clusters
Human_Markers_clean <- Human.celltype_markers %>% dplyr::select(6:7,)

Human_Markers_clean_Cluster_0 <- Human_Markers_clean[Human_Markers_clean$cluster == "Secretory", ]
Human_Markers_clean_Cluster_1 <- Human_Markers_clean[Human_Markers_clean$cluster == "Ciliated", ]
Human_Markers_clean_Cluster_2 <- Human_Markers_clean[Human_Markers_clean$cluster == "Club", ]
Human_Markers_clean_Cluster_3 <- Human_Markers_clean[Human_Markers_clean$cluster == "Basal", ]
Human_Markers_clean_Cluster_4 <- Human_Markers_clean[Human_Markers_clean$cluster == "Suprabasal", ]
Human_Markers_clean_Cluster_5 <- Human_Markers_clean[Human_Markers_clean$cluster == "Deuterosomal", ]
Human_Markers_clean_Cluster_6 <- Human_Markers_clean[Human_Markers_clean$cluster == "Ionocytes", ]

Human_Markers_clean_Cluster_0 <- Human_Markers_clean_Cluster_0 %>% dplyr::select(2,)
Human_Markers_clean_Cluster_1 <- Human_Markers_clean_Cluster_1 %>% dplyr::select(2,)
Human_Markers_clean_Cluster_2 <- Human_Markers_clean_Cluster_2 %>% dplyr::select(2,)
Human_Markers_clean_Cluster_3 <- Human_Markers_clean_Cluster_3 %>% dplyr::select(2,)
Human_Markers_clean_Cluster_4 <- Human_Markers_clean_Cluster_4 %>% dplyr::select(2,)
Human_Markers_clean_Cluster_5 <- Human_Markers_clean_Cluster_5 %>% dplyr::select(2,)
Human_Markers_clean_Cluster_6 <- Human_Markers_clean_Cluster_6 %>% dplyr::select(2,)

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Cluster_markers/Human_both")
write.csv(Human.celltype_markers, "Cluster_markers_Human.csv", row.names = F)
write.table(Human_Markers_clean_Cluster_0,"Human_Cluster_0.txt", row.names = FALSE)
write.table(Human_Markers_clean_Cluster_1,"Human_Cluster_1.txt", row.names = FALSE)
write.table(Human_Markers_clean_Cluster_2,"Human_Cluster_2.txt", row.names = FALSE)
write.table(Human_Markers_clean_Cluster_3,"Human_Cluster_3.txt", row.names = FALSE)
write.table(Human_Markers_clean_Cluster_4,"Human_Cluster_4.txt", row.names = FALSE)
write.table(Human_Markers_clean_Cluster_5,"Human_Cluster_5.txt", row.names = FALSE)
write.table(Human_Markers_clean_Cluster_6,"Human_Cluster_6.txt", row.names = FALSE)


#Saving the seurat object to disk to access in other script for analysis
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens")
SaveH5Seurat(Human.clustered, "Human.clustered", overwrite = TRUE)


