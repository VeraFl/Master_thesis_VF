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
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(EnhancedVolcano)
library(ggrepel)
library(DESeq2, quietly=T)
library(ComplexHeatmap, quietly = T)
library(circlize, quietly = T)
library(DEGreport, quietly=T)
library(tidyverse, quietly = T)
library(data.table, quietly = T)
library(clusterProfiler, quietly = T)
library(enrichplot, quietly = T)
library(ReactomePA, quietly = T)
library(ggVennDiagram,quietly = T)
library(PCAtools, quietly = T)
library(gprofiler2)
library(biomaRt)
library(gprofiler2)
library(ggVennDiagram)
library(ggupset)
library(MAST)
library(streamgraph)


#Setting the color schemes for the script:

#Treatment colors(camel-229E, MERS, Mock)
cols_treat = c("#F39C6B","#A5668B", "#96ADC8")
show_col(cols_treat)
#Status colors
cols_stat = c("#CD4631","#DEA47E")
show_col(cols_stat)

cols_stat_2 = c("#519872", "#CD4631","#DEA47E")
show_col(cols_stat_2)

#Celltype colors for initial UMAP
cols = hcl.colors(7, "Temps")


# The celltypes and their assigned color:
#Secretory = "#089392"
#Ciliated = "#40B48B"
#Club = "#9DCD84"
#Basal = "#EAE29C"
#Suprabasal = "#EAB672"
#Deuterosomal = "#E6866A"
#Ionocytes = "#CF597E"

memory.limit(24000)



##########################################################################################################################
#Loading filtered and merged Seurat Object into Script
#Infected and Mock Samples merged and cell clusters assigned

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.ferus")

Ferus <- LoadH5Seurat("Ferus.clustered.h5seurat")
Ferus <- Ferus.clustered
str(Ferus)
head(Ferus@meta.data)

###################################################################################################################################################
DefaultAssay(Ferus) <- "RNA"
Idents(Ferus) <- "Treat" 

Ferus_MERS <- subset(x = Ferus, idents = c("Mock", "MERS-CoV"))
Idents(Ferus_MERS) <- "celltype" 
levels(Ferus_MERS)



#MERS counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_mers_counts.pdf",
    width=5, height=5)
FeaturePlot(Ferus_MERS, features = "mers-3UTR", cols = c("orange", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T)
dev.off()

#DPP4 counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_mers_dpp4_counts.pdf",
    width=5, height=5)
FeaturePlot(Ferus_MERS, features = "DPP4", cols = c("orange", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T)
dev.off()

#ANPEP
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_mers_anpep_counts.pdf",
    width=5, height=5)
FeaturePlot(Ferus_MERS, features = "ANPEP", cols = c("orange", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T)
dev.off()

Idents(Ferus_MERS) <- "Status" 

#Reorder the status colors
cols_stat_2 = c("#DEA47E", "#519872", "#CD4631")
show_col(cols_stat_2)
cols_stat = c("#DEA47E","#DEA47E","#CD4631")
show_col(cols_stat)

#Infection
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_status.pdf",
    width=5, height=5)
DimPlot(Ferus_MERS, reduction = "umap", cols = cols_stat)
dev.off()





##################################################################################################################################################

#Threshold

##################################################################################################################################################
#########################################################################################################
#Correlation DPP4 expression and Viral_counts
#Extract the info in which row te DPP4 expression data is

a <- as.data.frame(Ferus_MERS@assays$RNA@data@Dimnames[[1]])

nrow(Ferus_MERS@assays$RNA@counts)
ncol(Ferus_MERS@assays$RNA@counts)

#5269
#show the specific location in the matrix for control
DPP4 <- as.matrix(Ferus_MERS@assays$RNA@counts[5269,])
Ferus_MERS <- AddMetaData(object = Ferus_MERS, metadata = DPP4, col.name = "Total.DPP4.Counts")


#Count the DPP4 expression

VlnPlot(Ferus_MERS, slot = "data", features = "DPP4")
#Threshold for DPP4 expression
Ferus_MERS@meta.data <- Ferus_MERS@meta.data %>% mutate(DPP4_Status = case_when(Total.DPP4.Counts > 0 ~ "DPP4_positive", Total.DPP4.Counts == 0 ~ "DPP4_negative"))


#26302
#show the specific location in the matrix for control
ANPEP <- as.matrix(Ferus_MERS@assays$RNA@counts[26302,])
Ferus_MERS <- AddMetaData(object = Ferus_MERS, metadata = ANPEP, col.name = "Total.ANPEP.Counts")
#Threshold for ANPEP expression
Ferus_MERS@meta.data <- Ferus_MERS@meta.data %>% mutate(ANPEP_Status = case_when(Total.ANPEP.Counts > 0 ~ "ANPEP_positive", Total.ANPEP.Counts == 0 ~ "ANPEP_negative"))

head(Ferus_MERS@meta.data)


df1 <- Ferus_MERS@meta.data %>% tidyseurat::select(c("celltype","Total.DPP4.Counts", "Total.ANPEP.Counts", "Status", "Treat", "DPP4_Status", "ANPEP_Status"))


#Plot 1: How many cells per celltype are Infected, Uninfected and Bystander
MERS_per_celltype <- df1 %>%
  group_by(celltype) %>%
  dplyr::count(Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            Status = Status) %>%
  mutate("Percentage" = n/Sum*100)


# Reordering Levels
MERS_per_celltype$Status <- factor(MERS_per_celltype$Status, levels=c('Infected','Bystander','Uninfected'))

cols_stat_2 = c("#CD4631","#519872","#DEA47E")
show_col(cols_stat_2)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_cells_status_per_celltype.pdf",
    width=6, height=5)
MERS_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = celltype, y = n, fill = Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 6, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_stat_2)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")
dev.off()


#Plot 2: How many infected cells have DPP4 expression

DPP4_per_celltype <- df1 %>%
  group_by(celltype, Status) %>%
  dplyr::count(DPP4_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            DPP4_Status = DPP4_Status,
            Status = Status) %>%
  mutate("Percentage" = n/Sum*100)

cols_1 = c("#F8A251","#9C4868")

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_inf_cells_dpp4_status_per_celltype.pdf",
    width=6, height=5)
DPP4_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  filter(Status == "Infected") %>%
  ggplot(aes(x = celltype, y = n, fill = DPP4_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 6, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_1)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")
dev.off()



#Plot 3: How many infected cells have ANPEP expression

ANPEP_per_celltype <- df1 %>%
  group_by(celltype) %>%
  count(ANPEP_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            ANPEP_Status = ANPEP_Status) %>%
  mutate("Percentage" = n/Sum*100)

cols_1 = c("orange","darkblue")

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_cells_dpp4_status_per_celltype.pdf",
    width=7, height=5)
ANPEP_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = celltype, y = n, fill = ANPEP_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_1)+
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


#Plot4: Fraction of cells which are DPP4 positive are infected?

DPP4_per_Status_per_celltype <- df1 %>%
  group_by(celltype, Status) %>%
  count(DPP4_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            DPP4_Status = DPP4_Status) %>%
  mutate("Percentage" = n/Sum*100)


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_cells_dpp4_status_per_celltype.pdf",
    width=15, height=5)
DPP4_per_Status_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = Status, y = n, fill = DPP4_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_stat_2)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="right")+
  facet_grid(cols = vars(celltype))
dev.off()










#Average Viral Percentage in each Cluster
dpp4 <- df1 %>%
  group_by(celltype) %>%
  summarise_at(vars(Total.DPP4.Counts), list(average_dpp4 = mean))

anpep <- df1 %>%
  group_by(celltype) %>%
  summarise_at(vars(Total.ANPEP.Counts), list(average_anpep = mean))


#Plot1
plot_dpp4 <- dpp4 %>% ggplot(aes(x = celltype, y = average_dpp4, fill = celltype))+
  geom_col() +
  ylab("Average DPP4 counts")+
  scale_fill_manual(values = cols)+
  ylim(c(0,0.3))+
  theme_classic()+
  NoLegend()

#Plot1
plot_anpep <- anpep %>% ggplot(aes(x = celltype, y = average_anpep, fill = celltype))+
  geom_col()+
  ylab("Average ANPEP counts") +
  scale_fill_manual(values = cols)+
  ylim(c(0,0.3))+
  theme_classic()+
  NoLegend()

grid.arrange(plot_dpp4, plot_anpep)




#Scatterplots for Correlation of MERS counts and DPP4

FeatureScatter(Ferus_MERS, feature1 = "Total.Norm.DPP4.Counts", feature2 = "Total.MERS.Counts")







##################################################################################################################################################
DefaultAssay(Ferus) <- "RNA"
Idents(Ferus) <- "Treat" 

Ferus_ACN4 <- subset(x = Ferus, idents = c("Mock", "dcCoV-ACN4"))
Idents(Ferus_ACN4) <- "celltype" 
levels(Ferus_ACN4)

#ACN4 counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_camel229E_counts.pdf",
    width=5, height=5)
FeaturePlot(Ferus_ACN4, features = "camel229E-3UTR", cols = c("orange", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T)
dev.off()

#DPP4 counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_camel229E_dpp4_counts.pdf",
    width=5, height=5)
FeaturePlot(Ferus_ACN4, features = "DPP4", cols = c("orange", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T)
dev.off()

#ANPEP
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_camel229E_anpep_counts.pdf",
    width=5, height=5)
FeaturePlot(Ferus_ACN4, features = "ANPEP", cols = c("orange", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T)
dev.off()
