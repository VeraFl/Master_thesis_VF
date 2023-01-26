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
cols = hcl.colors(7, "Temps")


memory.limit(24000)

##########################################################################################################################
#Loading filtered and merged Seurat Object into Script
#Infected and Mock Samples merged and cell clusters assigned

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.ferus")
Ferus <- LoadH5Seurat("Ferus.clustered.h5seurat")
#setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens")
#Human <- LoadH5Seurat("Human.clustered.h5seurat")
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/V.pacos")
Alpaca <- LoadH5Seurat("Alpaca.clustered.h5seurat")

head(Ferus@meta.data)
head(Human@meta.data)
head(Alpaca@meta.data)

###################################################################################################################################################

#Camel MERS , DPP4 and ANPEP counts in UMAP

###################################################################################################################################################
DefaultAssay(Ferus) <- "RNA"
Idents(Ferus) <- "Treat" 

Ferus_MERS <- subset(x = Ferus, idents = c("Mock", "MERS-CoV"))
Idents(Ferus_MERS) <- "celltype" 
levels(Ferus_MERS)

#MERS counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_mers_counts.pdf",
    width=5, height=5)
FeaturePlot(Ferus_MERS, features = "mers-3UTR", cols = c("#F8A251", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T)
dev.off()

#DPP4 counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_mers_dpp4_counts.pdf",
    width=5, height=5)
FeaturePlot(Ferus_MERS, features = "DPP4", cols = c("#F8A251", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T)
dev.off()

#ANPEP
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_mers_anpep_counts.pdf",
    width=5, height=5)
FeaturePlot(Ferus_MERS, features = "ANPEP", cols = c("#F8A251", "#06606E"), pt.size = 1, label.size = 5, label = T, repel = T)
dev.off()


##########################################################################################################################
#Receptor expression in cells - Quantitative
##########################################################################################################################


#Reorder the status colors

#Status colors infected vs uninfected
cols_stat = c("#BAA597", "#89043D","#BAA597")
show_col(cols_stat)

#Status colors infected, bystander and uninfected
cols_stat_2 = c("#519872", "#CD4631","#DEA47E")
show_col(cols_stat_2)

#define ANPEP and DPP4 status colors
cols_anpep = c("#F8A251","#06606E")
show_col(cols_anpep)

cols_dpp4 = c("#F8A251","#9C4868")
show_col(cols_dpp4)

Idents(Ferus_MERS) <- "Status" 

#Infection
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_Infection.pdf",
    width=7, height=4)
DimPlot(Ferus_MERS, reduction = "umap", split.by = "Treat", cols = cols_stat_2) +NoAxes()
dev.off()

##################################################################################################################################################

#Threshold

##################################################################################################################################################
#########################################################################################################
#Correlation DPP4 expression and Viral_counts
#Extract the info in which row te DPP4 expression data is

a <- as.data.frame(Ferus_MERS@assays$SCT@data@Dimnames[[1]])

nrow(Ferus_MERS@assays$RNA@counts)
ncol(Ferus_MERS@assays$RNA@counts)

#2825
#show the specific location in the matrix for control
DPP4 <- as.matrix(Ferus_MERS@assays$SCT@counts[2825,])
Ferus_MERS <- AddMetaData(object = Ferus_MERS, metadata = DPP4, col.name = "Total.DPP4.Counts")


#Threshold for DPP4 expression
Ferus_MERS@meta.data <- Ferus_MERS@meta.data %>% mutate(DPP4_Status = case_when(Total.DPP4.Counts > 0 ~ "DPP4_positive", Total.DPP4.Counts == 0 ~ "DPP4_negative"))


#14770
#show the specific location in the matrix for control
ANPEP <- as.matrix(Ferus_MERS@assays$SCT@counts[14770,])
Ferus_MERS <- AddMetaData(object = Ferus_MERS, metadata = ANPEP, col.name = "Total.ANPEP.Counts")

#Threshold for ANPEP expression
Ferus_MERS@meta.data <- Ferus_MERS@meta.data %>% mutate(ANPEP_Status = case_when(Total.ANPEP.Counts > 0 ~ "ANPEP_positive", Total.ANPEP.Counts == 0 ~ "ANPEP_negative"))
head(Ferus_MERS@meta.data)



Idents(Ferus_MERS) <- "DPP4_Status"

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_DPP4_status.pdf",
    width=5, height=5)
DimPlot(Ferus_MERS, reduction = "umap", label = F, label.size = 5, repel = TRUE, pt.size = 1, cols = c("#F8A251","#9C4868")) + NoAxes() 
dev.off()



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
MERS_per_celltype$celltype <- factor(MERS_per_celltype$celltype, levels = c('Camel Cluster 1','Camel Cluster 2','Secretory', 'Ciliated', 'Club', 'Basal'))


cols_stat_2 = c("#CD4631","#519872","#A06136")
show_col(cols_stat_2)

cols_stat_1 = c("#CD4631","#519872")
show_col(cols_stat_1)

cols_anpep = c("#F8A251","#06606E")
show_col(cols_anpep)

cols_dpp4 = c("#F8A251","#9C4868")
show_col(cols_dpp4)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_MERS_cells_status_per_celltype.pdf",
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
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")
dev.off()

MERS_per_celltype_inf_bys <- df1 %>%
  filter(Status != "Uninfected") %>%
  group_by(celltype) %>%
  dplyr::count(Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            Status = Status) %>%
  mutate("Percentage" = n/Sum*100)


# Reordering Levels
MERS_per_celltype_inf_bys$Status <- factor(MERS_per_celltype_inf_bys$Status, levels=c('Infected','Bystander'))
MERS_per_celltype_inf_bys$celltype <- factor(MERS_per_celltype_inf_bys$celltype, levels = c('Camel Cluster 1','Camel Cluster 2','Secretory', 'Ciliated', 'Club', 'Basal'))


#Plot 1.1: How many cells per celltype are Infected, Bystander
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_MERS_cells_inf_bys_per_celltype.pdf",
    width=7, height=6)
MERS_per_celltype_inf_bys %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = celltype, y = n, fill = Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 6, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_stat_1)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")
dev.off()


#Plot 2: How many cells have DPP4 expression

DPP4_per_celltype <- df1 %>%
  group_by(celltype) %>%
  dplyr::count(DPP4_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            DPP4_Status = DPP4_Status) %>%
  mutate("Percentage" = n/Sum*100)


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_MERS_cells_dpp4_status_per_celltype.pdf",
    width=6, height=5)
DPP4_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = celltype, y = n, fill = DPP4_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 6, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_dpp4)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")
dev.off()



#Plot 3: How many cells have ANPEP expression

ANPEP_per_celltype <- df1 %>%
  group_by(celltype) %>%
  dplyr::count(ANPEP_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            ANPEP_Status = ANPEP_Status) %>%
  mutate("Percentage" = n/Sum*100)


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_MERS_cells_anpep_status_per_celltype.pdf",
    width=7, height=5)
ANPEP_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = celltype, y = n, fill = ANPEP_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_anpep)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")
dev.off()


#Plot4: DPP4 downregulation in MERS-CoV exposed cells per Status

DPP4_per_DPP4_Status_per_celltype <- df1 %>%
  group_by(celltype, Status) %>%
  dplyr::count(DPP4_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            DPP4_Status = DPP4_Status) %>%
  mutate("Percentage" = n/Sum*100)


DPP4_per_DPP4_Status_per_celltype$Status <- factor(DPP4_per_DPP4_Status_per_celltype$Status, levels=c('Uninfected','Infected', 'Bystander'))


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_MERS_cells_dpp4_status_status_celltype.pdf",
    width=14, height=6)
DPP4_per_DPP4_Status_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  #filter(celltype == "Ciliated") %>%
  ggplot(aes(x = celltype, y = n, fill = DPP4_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_dpp4)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(colour= "black", angle = 45, hjust = 1, vjust = 1, size = 15),
        axis.text.y = element_text(colour= "black", size = 15),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")+
  facet_grid(cols = vars(Status))
dev.off()

#Plot4.1: DPP4 downregulation in MERS-CoV exposed cells per Treat

DPP4_Status_per_Treat_per_celltype <- df1 %>%
  group_by(celltype, Treat) %>%
  dplyr::count(DPP4_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            DPP4_Status = DPP4_Status,
            Treat = Treat) %>%
  mutate("Percentage" = n/Sum*100)


DPP4_Status_per_Treat_per_celltype$Treat <- factor(DPP4_Status_per_Treat_per_celltype$Treat, levels=c('Mock', 'MERS-CoV'))


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_MERS_cells_dpp4_status_treat_celltype.pdf",
    width=8, height=6)
DPP4_Status_per_Treat_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  #filter(celltype == "Ciliated") %>%
  ggplot(aes(x = celltype, y = n, fill = DPP4_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_dpp4)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(colour= "black", angle = 45, hjust = 1, vjust = 1, size = 15),
        axis.text.y = element_text(colour= "black", size = 15),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")+
  facet_grid(cols = vars(Treat))
dev.off()



#Plot5: ANPEP downregulation in MERS-CoV exposed cells

ANPEP_per_ANPEP_Status_per_celltype <- df1 %>%
  group_by(celltype, Status) %>%
  dplyr::count(ANPEP_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            ANPEP_Status = ANPEP_Status) %>%
  mutate("Percentage" = n/Sum*100)


ANPEP_per_ANPEP_Status_per_celltype$Status <- factor(ANPEP_per_ANPEP_Status_per_celltype$Status, levels=c('Uninfected','Infected', 'Bystander'))


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_MERS_cells_anpep_status_celltype.pdf",
    width=14, height=6)
ANPEP_per_ANPEP_Status_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  #filter(celltype == "Ciliated") %>%
  ggplot(aes(x = celltype, y = n, fill = ANPEP_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_anpep)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(colour= "black", angle = 45, hjust = 1, vjust = 1, size = 15),
        axis.text.y = element_text(colour= "black", size = 15),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")+
  facet_grid(cols = vars(Status))
dev.off()


###################################################################################################################################################

#Human MERS , DPP4 and ANPEP counts in UMAP

###################################################################################################################################################
Human <- Human.clustered
DefaultAssay(Human) <- "RNA"
Idents(Human) <- "Treat" 

Human_MERS <- subset(x = Human, idents = c("Mock", "MERS-CoV"))
Idents(Human_MERS) <- "celltype" 
levels(Human_MERS)

#MERS counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/UMAP_mers_counts.pdf",
    width=5, height=5)
FeaturePlot(Human_MERS, slot = "data", features = "mers-3UTR", cols = c("orange", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T)
dev.off()

#DPP4 counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/UMAP_mers_dpp4_counts.pdf",
    width=5, height=5)
FeaturePlot(Human_MERS, slot = "data", features = "DPP4", cols = c("orange", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T)
dev.off()

#ANPEP
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/UMAP_mers_anpep_counts.pdf",
    width=5, height=5)
FeaturePlot(Human_MERS, slot = "data", features = "ANPEP", cols = c("orange", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T)
dev.off()


##########################################################################################################################
#Receptor expression in cells - Quantitative
##########################################################################################################################


#Reorder the status colors
cols_stat_2 = c("#DEA47E", "#519872", "#CD4631")
show_col(cols_stat_2)
cols_stat = c("#DEA47E","#DEA47E","#CD4631")
show_col(cols_stat)

#define ANPEP and DPP4 status colors
cols_anpep = c("#F8A251","#06606E")
show_col(cols_anpep)

cols_dpp4 = c("#F8A251","#9C4868")
show_col(cols_dpp4)

Idents(Human_MERS) <- "Status" 

#Infection
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_Human_both_virus/UMAP_status.pdf",
    width=5, height=5)
DimPlot(Human_MERS, reduction = "umap", split.by = "Treat", cols = cols_stat_2)
dev.off()

##################################################################################################################################################

#Threshold

##################################################################################################################################################
#########################################################################################################
#Correlation DPP4 expression and Viral_counts
#Extract the info in which row te DPP4 expression data is

a <- as.data.frame(Human_MERS@assays$RNA@data@Dimnames[[1]])

nrow(Human_MERS@assays$RNA@counts)
ncol(Human_MERS@assays$RNA@counts)

#2779
#show the specific location in the matrix for control
DPP4 <- as.matrix(Human_MERS@assays$RNA@counts[2779,])
Human_MERS <- AddMetaData(object = Human_MERS, metadata = DPP4, col.name = "Total.DPP4.Counts")


#Threshold for DPP4 expression
Human_MERS@meta.data <- Human_MERS@meta.data %>% mutate(DPP4_Status = case_when(Total.DPP4.Counts > 0 ~ "DPP4_positive", Total.DPP4.Counts == 0 ~ "DPP4_negative"))


#14643
#show the specific location in the matrix for control
ANPEP <- as.matrix(Human_MERS@assays$RNA@counts[14643,])
Human_MERS <- AddMetaData(object = Human_MERS, metadata = ANPEP, col.name = "Total.ANPEP.Counts")
#Threshold for ANPEP expression
Human_MERS@meta.data <- Human_MERS@meta.data %>% mutate(ANPEP_Status = case_when(Total.ANPEP.Counts > 0 ~ "ANPEP_positive", Total.ANPEP.Counts == 0 ~ "ANPEP_negative"))

head(Human_MERS@meta.data)


df1 <- Human_MERS@meta.data %>% tidyseurat::select(c("celltype","Total.DPP4.Counts", "Total.ANPEP.Counts", "Status", "Treat", "DPP4_Status", "ANPEP_Status"))


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

cols_stat_2 = c("#CD4631","#519872","#A06136")
show_col(cols_stat_2)

cols_stat_1 = c("#CD4631","#519872")
show_col(cols_stat_1)

cols_anpep = c("#F8A251","#06606E")
show_col(cols_anpep)

cols_dpp4 = c("#F8A251","#9C4868")
show_col(cols_dpp4)


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/fraction_of_MERS_cells_status_per_celltype.pdf",
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
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")
dev.off()


MERS_per_celltype_inf_bys <- df1 %>%
  filter(Status != "Uninfected") %>%
  group_by(celltype) %>%
  dplyr::count(Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            Status = Status) %>%
  mutate("Percentage" = n/Sum*100)

# Reordering Levels
MERS_per_celltype_inf_bys$Status <- factor(MERS_per_celltype_inf_bys$Status, levels=c('Infected','Bystander'))
MERS_per_celltype_inf_bys$celltype <- factor(MERS_per_celltype_inf_bys$celltype, levels = c('Secretory', 'Ciliated', 'Club', 'Basal', 'Suprabasal', 'Deuterosomal', 'Ionocytes'))


#Plot 1.1: How many cells per celltype are Infected, Bystander
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/fraction_of_MERS_cells_inf_bys_per_celltype.pdf",
    width=7, height=6)
MERS_per_celltype_inf_bys %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = celltype, y = n, fill = Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 6, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_stat_1)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")
dev.off()

#Plot 2: How many cells have DPP4 expression

DPP4_per_celltype <- df1 %>%
  group_by(celltype) %>%
  dplyr::count(DPP4_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            DPP4_Status = DPP4_Status) %>%
  mutate("Percentage" = n/Sum*100)



pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/fraction_of_MERS_inf_cells_dpp4_status_per_celltype.pdf",
    width=6, height=5)
DPP4_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = celltype, y = n, fill = DPP4_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 6, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_dpp4)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")
dev.off()



#Plot 3: How many cells have ANPEP expression

ANPEP_per_celltype <- df1 %>%
  group_by(celltype) %>%
  dplyr::count(ANPEP_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            ANPEP_Status = ANPEP_Status) %>%
  mutate("Percentage" = n/Sum*100)


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/fraction_of_MERS_cells_anpep_status_per_celltype.pdf",
    width=7, height=5)
ANPEP_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = celltype, y = n, fill = ANPEP_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_anpep)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")
dev.off()


#Plot4: DPP4 downregulation in MERS-CoV exposed ciliated cells

DPP4_per_Status_per_celltype <- df1 %>%
  group_by(celltype, Status) %>%
  dplyr::count(DPP4_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            DPP4_Status = DPP4_Status) %>%
  mutate("Percentage" = n/Sum*100)


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/fraction_of_MERS_cells_dpp4_status_per_celltype.pdf",
    width=15, height=5)
DPP4_per_Status_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  #filter(celltype == "Ciliated") %>%
  ggplot(aes(x = celltype, y = n, fill = DPP4_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_dpp4)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="right")+
  facet_grid(cols = vars(Status))
dev.off()

#Plot4.1: DPP4 downregulation in MERS-CoV exposed cells

DPP4_per_Treat_per_celltype <- df1 %>%
  group_by(celltype, Treat) %>%
  dplyr::count(DPP4_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            DPP4_Status = DPP4_Status) %>%
  mutate("Percentage" = n/Sum*100)

DPP4_per_Treat_per_celltype$Treat <- factor(DPP4_per_Treat_per_celltype$Treat, levels=c('Mock','MERS-CoV'))


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/fraction_of_MERS_cells_dpp4_treat_per_celltype.pdf",
    width=15, height=5)
DPP4_per_Treat_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  #filter(celltype == "Ciliated") %>%
  ggplot(aes(x = celltype, y = n, fill = DPP4_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_dpp4)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")+
  facet_grid(cols = vars(Treat))
dev.off()

#Plot5: ANPEP downregulation in MERS-CoV exposed cells

ANPEP_per_ANPEP_Status_per_celltype <- df1 %>%
  group_by(celltype, Status) %>%
  dplyr::count(ANPEP_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            ANPEP_Status = ANPEP_Status) %>%
  mutate("Percentage" = n/Sum*100)


ANPEP_per_ANPEP_Status_per_celltype$Status <- factor(ANPEP_per_ANPEP_Status_per_celltype$Status, levels=c('Uninfected','Infected', 'Bystander'))


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_MERS_cells_anpep_status_celltype.pdf",
    width=14, height=6)
ANPEP_per_ANPEP_Status_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  #filter(celltype == "Ciliated") %>%
  ggplot(aes(x = celltype, y = n, fill = ANPEP_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_anpep)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")+
  facet_grid(cols = vars(Status))
dev.off()

#Plot5.1: ANPEP downregulation in MERS-CoV exposed cells

ANPEP_per_ANPEP_Treat_per_celltype <- df1 %>%
  group_by(celltype, Treat) %>%
  dplyr::count(ANPEP_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            ANPEP_Status = ANPEP_Status) %>%
  mutate("Percentage" = n/Sum*100)


ANPEP_per_ANPEP_Treat_per_celltype$Treat <- factor(ANPEP_per_ANPEP_Treat_per_celltype$Treat, levels=c('Mock','MERS-CoV'))


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_MERS_cells_anpep_treat_celltype.pdf",
    width=14, height=6)
ANPEP_per_ANPEP_Treat_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  #filter(celltype == "Ciliated") %>%
  ggplot(aes(x = celltype, y = n, fill = ANPEP_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_anpep)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")+
  facet_grid(cols = vars(Treat))
dev.off()

###################################################################################################################################################

#Llama MERS , DPP4 and ANPEP counts in UMAP

###################################################################################################################################################
Alpaca <- Alpaca.clustered
DefaultAssay(Alpaca) <- "RNA"
Idents(Alpaca) <- "Treat" 

Alpaca_MERS <- subset(x = Alpaca, idents = c("Mock", "MERS-CoV"))
Idents(Alpaca_MERS) <- "celltype" 
levels(Alpaca_MERS)


#MERS counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/UMAP_mers_counts.pdf",
    width=7, height=7)
FeaturePlot(Alpaca_MERS, slot = "data", features = "mers-3UTR", cols = c("#F8A251", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T) +NoAxes()
dev.off()

#DPP4 counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/UMAP_mers_dpp4_counts.pdf",
    width=5, height=5)
FeaturePlot(Alpaca_MERS, slot = "data", features = "DPP4", cols = c("#F8A251", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T) +NoAxes()
dev.off()

#ANPEP
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/UMAP_mers_anpep_counts.pdf",
    width=5, height=5)
FeaturePlot(Alpaca_MERS, slot = "data", features = "ANPEP", cols = c("#F8A251", "#06606E"), pt.size = 1, label.size = 5, label = T, repel = T) +NoAxes()
dev.off()


##########################################################################################################################
#Receptor expression in cells - Quantitative
##########################################################################################################################


#Reorder the status colors
cols_stat_2 = c("#DEA47E", "#519872", "#CD4631")
show_col(cols_stat_2)
cols_stat = c("#DEA47E","#DEA47E","#CD4631")
show_col(cols_stat)


#define ANPEP and DPP4 status colors
cols_anpep = c("#F8A251","#06606E")
show_col(cols_anpep)

cols_dpp4 = c("#F8A251","#9C4868")
show_col(cols_dpp4)


Idents(Alpaca_MERS) <- "Status" 

#Infection
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/UMAP_status.pdf",
    width=5, height=5)
DimPlot(Alpaca_MERS, reduction = "umap", split.by = "Treat", cols = cols_stat) +NoAxes()
dev.off()

##################################################################################################################################################

#Threshold

##################################################################################################################################################
#########################################################################################################
#Correlation DPP4 expression and Viral_counts
#Extract the info in which row te DPP4 expression data is

a <- as.data.frame(Alpaca_MERS@assays$SCT@data@Dimnames[[1]])

nrow(Alpaca_MERS@assays$RNA@counts)
ncol(Alpaca_MERS@assays$RNA@counts)

#2475
#show the specific location in the matrix for control
DPP4 <- as.matrix(Alpaca_MERS@assays$SCT@counts[2475,])
Alpaca_MERS <- AddMetaData(object = Alpaca_MERS, metadata = DPP4, col.name = "Total.DPP4.Counts")


#Threshold for DPP4 expression
Alpaca_MERS@meta.data <- Alpaca_MERS@meta.data %>% mutate(DPP4_Status = case_when(Total.DPP4.Counts > 0 ~ "DPP4_positive", Total.DPP4.Counts == 0 ~ "DPP4_negative"))


#11227
#show the specific location in the matrix for control
ANPEP <- as.matrix(Alpaca_MERS@assays$SCT@counts[11227,])
Alpaca_MERS <- AddMetaData(object = Alpaca_MERS, metadata = ANPEP, col.name = "Total.ANPEP.Counts")
#Threshold for ANPEP expression
Alpaca_MERS@meta.data <- Alpaca_MERS@meta.data %>% mutate(ANPEP_Status = case_when(Total.ANPEP.Counts > 0 ~ "ANPEP_positive", Total.ANPEP.Counts == 0 ~ "ANPEP_negative"))

head(Alpaca_MERS@meta.data)

Idents(Alpaca_MERS) <- "DPP4_Status"

#DPP4_Status
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/UMAP_DPP4_status.pdf",
    width=8, height=7)
DimPlot(Alpaca_MERS, reduction = "umap", label = F, label.size = 5, repel = TRUE, pt.size = 1, cols = c("#F8A251","#9C4868")) + NoAxes()
dev.off()


df1 <- Alpaca_MERS@meta.data %>% tidyseurat::select(c("celltype","Total.DPP4.Counts", "Total.ANPEP.Counts", "Status", "Treat", "DPP4_Status", "ANPEP_Status"))


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

cols_stat_2 = c("#CD4631","#519872","#A06136")
show_col(cols_stat_2)

cols_stat_1 = c("#CD4631","#519872")
show_col(cols_stat_1)

cols_anpep = c("#F8A251","#06606E")
show_col(cols_anpep)

cols_dpp4 = c("#F8A251","#9C4868")
show_col(cols_dpp4)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/fraction_of_MERS_cells_status_per_celltype.pdf",
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
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")
dev.off()

MERS_per_celltype_inf_bys <- df1 %>%
  filter(Status != "Uninfected") %>%
  group_by(celltype) %>%
  dplyr::count(Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            Status = Status) %>%
  mutate("Percentage" = n/Sum*100)

# Reordering Levels
MERS_per_celltype_inf_bys$Status <- factor(MERS_per_celltype_inf_bys$Status, levels=c('Infected','Bystander'))
MERS_per_celltype_inf_bys$celltype <- factor(MERS_per_celltype_inf_bys$celltype, levels = c("Llama Cluster 1", "Llama Cluster 2", "Secretory", "Ciliated", "Club", "Basal"))


#Plot 1.1: How many cells per celltype are Infected, Bystander
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/fraction_of_MERS_cells_inf_bys_per_celltype.pdf",
    width=7, height=6)
MERS_per_celltype_inf_bys %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = celltype, y = n, fill = Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 6, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_stat_1)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")
dev.off()


#Plot 2: How many cells have DPP4 expression

DPP4_per_celltype <- df1 %>%
  group_by(celltype) %>%
  dplyr::count(DPP4_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            DPP4_Status = DPP4_Status) %>%
  mutate("Percentage" = n/Sum*100)


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/fraction_of_MERS_cells_dpp4_status_per_celltype.pdf",
    width=6, height=5)
DPP4_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = celltype, y = n, fill = DPP4_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 6, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_dpp4)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")
dev.off()



#Plot 3: How many cells have ANPEP expression

ANPEP_per_celltype <- df1 %>%
  group_by(celltype) %>%
  dplyr::count(ANPEP_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            ANPEP_Status = ANPEP_Status) %>%
  mutate("Percentage" = n/Sum*100)


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/fraction_of_MERS_cells_anpep_status_per_celltype.pdf",
    width=7, height=5)
ANPEP_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = celltype, y = n, fill = ANPEP_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_anpep)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")
dev.off()


#Plot4: DPP4 downregulation in MERS-CoV exposed ciliated cells

DPP4_per_Status_per_celltype <- df1 %>%
  group_by(celltype, Status) %>%
  dplyr::count(DPP4_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            DPP4_Status = DPP4_Status) %>%
  mutate("Percentage" = n/Sum*100)

DPP4_per_Status_per_celltype$Status <- factor(DPP4_per_Status_per_celltype$Status, levels=c('Uninfected','Infected','Bystander'))
DPP4_per_Status_per_celltype$celltype <- factor(DPP4_per_Status_per_celltype$celltype, levels=c('Basal', 'Secretory','Club', 'Ciliated', 'Llama Cluster 1', 'Llama Cluster 2'))


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/fraction_of_MERS_cells_dpp4_status_per_celltype.pdf",
    width=14, height=6)
DPP4_per_Status_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  #filter(celltype == "Ciliated") %>%
  ggplot(aes(x = celltype, y = n, fill = DPP4_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_dpp4)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(colour= "black", angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(colour= "black", size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")+
  facet_grid(cols = vars(Status))
dev.off()


#Plot5: ANPEP downregulation in MERS-CoV exposed cells

ANPEP_per_ANPEP_Status_per_celltype <- df1 %>%
  group_by(celltype, Status) %>%
  dplyr::count(ANPEP_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            ANPEP_Status = ANPEP_Status) %>%
  mutate("Percentage" = n/Sum*100)


ANPEP_per_ANPEP_Status_per_celltype$Status <- factor(ANPEP_per_ANPEP_Status_per_celltype$Status, levels=c('Uninfected','Infected', 'Bystander'))
ANPEP_per_ANPEP_Status_per_celltype$celltype <- factor(ANPEP_per_ANPEP_Status_per_celltype$celltype, levels=c('Basal', 'Secretory','Club', 'Ciliated', 'Llama Cluster 1', 'Llama Cluster 2'))


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/fraction_of_MERS_cells_anpep_status_celltype.pdf",
    width=14, height=6)
ANPEP_per_ANPEP_Status_per_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  #filter(celltype == "Ciliated") %>%
  ggplot(aes(x = celltype, y = n, fill = ANPEP_Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_anpep)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")+
  facet_grid(cols = vars(Status))
dev.off()

