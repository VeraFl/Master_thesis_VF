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
library(gprofiler2)
library(biomaRt)
library(gprofiler2)
library(ggVennDiagram)
library(ggupset)
library(MAST)

#Setting the color schemes for the script:

#Treatment colors(camel-229E, 229E, Mock)
cols_treat = c("#F39C6B","#A5668B", "#96ADC8")
show_col(cols_treat)
#Status colors
cols_stat = c("#CD4631","#DEA47E")
show_col(cols_stat)

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

?log1p
###################################################################################################################################################

#Camel 229E , DPP4 and ANPEP counts in UMAP

###################################################################################################################################################
Ferus <- Ferus.clustered
DefaultAssay(Ferus) <- "RNA"
Idents(Ferus) <- "Treat" 

Ferus_229E <- subset(x = Ferus, idents = c("Mock", "dcCoV-ACN4"))
Idents(Ferus_229E) <- "celltype" 
levels(Ferus_229E)

#229E counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_229E_counts.pdf",
    width=5, height=5)
FeaturePlot(Ferus_229E, features = "camel229E-3UTR", cols = c("#F8A251", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T)
dev.off()

#DPP4 counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_229E_dpp4_counts.pdf",
    width=5, height=5)
FeaturePlot(Ferus_229E, features = "DPP4", cols = c("#F8A251", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T)
dev.off()

#ANPEP
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_229E_anpep_counts.pdf",
    width=5, height=5)
FeaturePlot(Ferus_229E, features = "ANPEP", cols = c("#F8A251", "#06606E"), pt.size = 1, label.size = 5, label = T, repel = T)
dev.off()


##########################################################################################################################
#Receptor expression in cells - Quantitative
##########################################################################################################################


#Reorder the status colors
cols_stat_2 = c("#DEA47E", "#519872", "#CD4631")
show_col(cols_stat_2)
cols_stat = c("#DEA47E","#DEA47E","#CD4631")
show_col(cols_stat)

#define receptor color
cols_anpep = c("#F8A251","#06606E")
show_col(cols_anpep)

cols_dpp4 = c("#F8A251","#9C4868")
show_col(cols_dpp4)

Idents(Ferus_229E) <- "Status" 

#Infection
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_status.pdf",
    width=5, height=5)
DimPlot(Ferus_229E, reduction = "umap", split.by = "Treat", cols = cols_stat_2)
dev.off()

##################################################################################################################################################

#Threshold

##################################################################################################################################################
#########################################################################################################
#Correlation DPP4 expression and Viral_counts
#Extract the info in which row te DPP4 expression data is

a <- as.data.frame(Ferus@assays$SCT@data@Dimnames[[1]])

nrow(Ferus_229E@assays$RNA@counts)
ncol(Ferus_229E@assays$RNA@counts)

#2825
#show the specific location in the matrix for control
DPP4 <- as.matrix(Ferus_229E@assays$SCT@counts[2825,])
Ferus_229E <- AddMetaData(object = Ferus_229E, metadata = DPP4, col.name = "Total.DPP4.Counts")


#Threshold for DPP4 expression
Ferus_229E@meta.data <- Ferus_229E@meta.data %>% mutate(DPP4_Status = case_when(Total.DPP4.Counts > 0 ~ "DPP4_positive", Total.DPP4.Counts == 0 ~ "DPP4_negative"))


#14770
#show the specific location in the matrix for control
ANPEP <- as.matrix(Ferus_229E@assays$SCT@counts[14770,])
Ferus_229E <- AddMetaData(object = Ferus_229E, metadata = ANPEP, col.name = "Total.ANPEP.Counts")


#Threshold for ANPEP expression
Ferus_229E@meta.data <- Ferus_229E@meta.data %>% mutate(ANPEP_Status = case_when(Total.ANPEP.Counts > 0 ~ "ANPEP_positive", Total.ANPEP.Counts == 0 ~ "ANPEP_negative"))

head(Ferus_229E@meta.data)

Idents(Ferus_229E) <- "ANPEP_Status"

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_ANPEP_status.pdf",
    width=5, height=5)
DimPlot(Ferus_229E, reduction = "umap", label = F, label.size = 5, repel = TRUE, pt.size = 1, cols = c("#F8A251","#0E7C7B")) + NoAxes()
dev.off()


df1 <- Ferus_229E@meta.data %>% tidyseurat::select(c("celltype","Total.DPP4.Counts", "Total.ANPEP.Counts", "Status", "Treat", "DPP4_Status", "ANPEP_Status"))


#Plot 1: How many cells per celltype are Infected, Uninfected and Bystander
camel229E_per_celltype <- df1 %>%
  group_by(celltype) %>%
  dplyr::count(Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            Status = Status) %>%
  mutate("Percentage" = n/Sum*100)


# Reordering Levels
camel229E_per_celltype$Status <- factor(camel229E_per_celltype$Status, levels=c('Infected','Bystander','Uninfected'))

cols_stat_2 = c("#CD4631","#519872","#A06136")
show_col(cols_stat_2)

cols_stat_1 = c("#CD4631","#519872")
show_col(cols_stat_1)

cols_anpep = c("#F8A251","#06606E")
show_col(cols_anpep)

cols_dpp4 = c("#F8A251","#9C4868")
show_col(cols_dpp4)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_229E_cells_status_per_celltype.pdf",
    width=6, height=5)
camel229E_per_celltype %>%
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

camel229E_per_celltype_inf_bys <- df1 %>%
  filter(Status != "Uninfected") %>%
  group_by(celltype) %>%
  dplyr::count(Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            Status = Status) %>%
  mutate("Percentage" = n/Sum*100)


# Reordering Levels
camel229E_per_celltype_inf_bys$Status <- factor(camel229E_per_celltype_inf_bys$Status, levels=c('Infected','Bystander'))
camel229E_per_celltype_inf_bys$celltype <- factor(camel229E_per_celltype_inf_bys$celltype, levels = c('Camel Cluster 1','Camel Cluster 2','Secretory', 'Ciliated', 'Club', 'Basal'))


#Plot 1.1: How many cells per celltype are Infected, Bystander
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_camel229E_cells_inf_bys_per_celltype.pdf",
    width=7, height=6)
camel229E_per_celltype_inf_bys %>%
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


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_229E_inf_cells_dpp4_status_per_celltype.pdf",
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



#Plot 3: How many infected cells have ANPEP expression

ANPEP_per_celltype <- df1 %>%
  group_by(celltype) %>%
  dplyr::count(ANPEP_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            ANPEP_Status = ANPEP_Status) %>%
  mutate("Percentage" = n/Sum*100)


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_229E_inf_cells_anpep_status_per_celltype.pdf",
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


#Plot4: DPP4 downregulation in 229E-CoV exposed cells per Status

DPP4_per_Status_per_celltype <- df1 %>%
  group_by(celltype, Status) %>%
  dplyr::count(DPP4_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            DPP4_Status = DPP4_Status) %>%
  mutate("Percentage" = n/Sum*100)


DPP4_per_Status_per_celltype$Status <- factor(DPP4_per_Status_per_celltype$Status, levels=c('Uninfected','Infected', 'Bystander'))


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/Plot_229E_DPP4_camel.pdf",
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
        axis.text.x = element_text(colour= "black", angle = 45, hjust = 1, vjust = 1, size = 15),
        axis.text.y = element_text(colour= "black", size = 15),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="right")+
  facet_grid(cols = vars(Status))
dev.off()



#Plot5: ANPEP downregulation in 229E-CoV exposed cells per Status

ANPEP_per_Status_per_celltype <- df1 %>%
  group_by(celltype, Status) %>%
  dplyr::count(ANPEP_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            ANPEP_Status = ANPEP_Status) %>%
  mutate("Percentage" = n/Sum*100)


ANPEP_per_Status_per_celltype$Status <- factor(ANPEP_per_Status_per_celltype$Status, levels=c('Uninfected','Infected', 'Bystander'))


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/Plot_229E_ANPEP_camel.pdf",
    width=14, height=6)
ANPEP_per_Status_per_celltype %>%
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
        legend.position="right")+
  facet_grid(cols = vars(Status))
dev.off()




###################################################################################################################################################

#Human 229E , DPP4 and ANPEP counts in UMAP

###################################################################################################################################################
Human <- Human.clustered
DefaultAssay(Human) <- "RNA"
Idents(Human) <- "Treat" 

Human_229E <- subset(x = Human, idents = c("Mock", "HCoV-229E"))
Idents(Human_229E) <- "celltype" 
levels(Human_229E)


#229E counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/UMAP_229E_counts.pdf",
    width=5, height=5)
FeaturePlot(Human_229E, slot = "data", features = "hcov-3UTR", cols = c("#F8A251", "#06606E"), pt.size = 1, label.size = 5, label = T, repel = T)
dev.off()

#DPP4 counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/UMAP_229E_dpp4_counts.pdf",
    width=5, height=5)
FeaturePlot(Human_229E, slot = "data", features = "DPP4", cols = c("#F8A251", "#06606E"), pt.size = 1, label.size = 5, label = T, repel = T)
dev.off()

#ANPEP
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/UMAP_229E_anpep_counts.pdf",
    width=5, height=5)
FeaturePlot(Human_229E, slot = "data", features = "ANPEP", cols = c("#F8A251", "#06606E"), pt.size = 1, label.size = 5, label = T, repel = T)
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


Idents(Human_229E) <- "Status" 

#Infection
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_Human_both_virus/UMAP_status.pdf",
    width=5, height=5)
DimPlot(Human_229E, reduction = "umap", split.by = "Treat", cols = cols_stat_2)
dev.off()

##################################################################################################################################################

#Threshold

##################################################################################################################################################
#########################################################################################################
#Correlation DPP4 expression and Viral_counts
#Extract the info in which row te DPP4 expression data is

a <- as.data.frame(Human_229E@assays$RNA@data@Dimnames[[1]])

nrow(Human_229E@assays$RNA@counts)
ncol(Human_229E@assays$RNA@counts)

#2779
#show the specific location in the matrix for control
DPP4 <- as.matrix(Human_229E@assays$RNA@counts[2779,])
Human_229E <- AddMetaData(object = Human_229E, metadata = DPP4, col.name = "Total.DPP4.Counts")


#Threshold for DPP4 expression
Human_229E@meta.data <- Human_229E@meta.data %>% mutate(DPP4_Status = case_when(Total.DPP4.Counts > 0 ~ "DPP4_positive", Total.DPP4.Counts == 0 ~ "DPP4_negative"))


#14643
#show the specific location in the matrix for control
ANPEP <- as.matrix(Human_229E@assays$RNA@counts[14643,])
Human_229E <- AddMetaData(object = Human_229E, metadata = ANPEP, col.name = "Total.ANPEP.Counts")
#Threshold for ANPEP expression
Human_229E@meta.data <- Human_229E@meta.data %>% mutate(ANPEP_Status = case_when(Total.ANPEP.Counts > 0 ~ "ANPEP_positive", Total.ANPEP.Counts == 0 ~ "ANPEP_negative"))

head(Human_229E@meta.data)

Idents(Alpaca_MERS) <- "ANPEP_Status"

#DPP4_Status
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/UMAP_ANPEP_status.pdf",
    width=8, height=7)
DimPlot(Human_MERS, reduction = "umap", label = F, label.size = 5, repel = TRUE, pt.size = 1, cols = c("#F8A251","#06606E")) + NoAxes()
dev.off()


df1 <- Human_229E@meta.data %>% tidyseurat::select(c("celltype","Total.DPP4.Counts", "Total.ANPEP.Counts", "Status", "Treat", "DPP4_Status", "ANPEP_Status"))


#Plot 1: How many cells per celltype are Infected, Uninfected and Bystander
HCoV229E_per_celltype <- df1 %>%
  group_by(celltype) %>%
  dplyr::count(Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            Status = Status) %>%
  mutate("Percentage" = n/Sum*100)


# Reordering Levels
HCoV229E_per_celltype$Status <- factor(HCoV229E_per_celltype$Status, levels=c('Infected','Bystander','Uninfected'))

cols_stat_2 = c("#CD4631","#519872","#A06136")
show_col(cols_stat_2)

cols_stat_1 = c("#CD4631","#519872")
show_col(cols_stat_1)

cols_anpep = c("#F8A251","#06606E")
show_col(cols_anpep)

cols_dpp4 = c("#F8A251","#9C4868")
show_col(cols_dpp4)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/fraction_of_229E_cells_status_per_celltype.pdf",
    width=6, height=5)
HCoV229E_per_celltype %>%
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



HCoV229E_per_celltype_inf_bys <- df1 %>%
  filter(Status != "Uninfected") %>%
  group_by(celltype) %>%
  dplyr::count(Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            Status = Status) %>%
  mutate("Percentage" = n/Sum*100)

# Reordering Levels
HCoV229E_per_celltype_inf_bys$Status <- factor(HCoV229E_per_celltype_inf_bys$Status, levels=c('Infected','Bystander'))
HCoV229E_per_celltype_inf_bys$celltype <- factor(HCoV229E_per_celltype_inf_bys$celltype, levels = c('Secretory', 'Ciliated', 'Club', 'Basal', 'Suprabasal', 'Deuterosomal', 'Ionocytes'))


#Plot 1.1: How many cells per celltype are Infected, Bystander
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/fraction_of_HCoV229E_cells_inf_bys_per_celltype.pdf",
    width=7, height=6)
HCoV229E_per_celltype_inf_bys %>%
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




pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/fraction_of_229E_inf_cells_dpp4_status_per_celltype.pdf",
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


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/fraction_of_229E_cells_anpep_status_per_celltype.pdf",
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


#Plot4: ANPEP downregulation in 229E-CoV exposed cells per Status

ANPEP_per_Status_per_celltype <- df1 %>%
  group_by(celltype, Status) %>%
  dplyr::count(ANPEP_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            ANPEP_Status = ANPEP_Status) %>%
  mutate("Percentage" = n/Sum*100)

ANPEP_per_Status_per_celltype$Status <- factor(ANPEP_per_Status_per_celltype$Status, levels=c('Infected','Bystander','Uninfected'))


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/fraction_of_229E_cells_anpep_status_per_celltype.pdf",
    width=15, height=5)
ANPEP_per_Status_per_celltype %>%
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



#Plot4.1: ANPEP downregulation in 229E-CoV exposed cells per Treat

ANPEP_per_Treat_per_celltype <- df1 %>%
  group_by(celltype, Treat) %>%
  dplyr::count(ANPEP_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            ANPEP_Status = ANPEP_Status) %>%
  mutate("Percentage" = n/Sum*100)

ANPEP_per_Treat_per_celltype$Treat <- factor(ANPEP_per_Treat_per_celltype$Treat, levels=c('Mock','HCoV-229E'))

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/fraction_of_229E_cells_anpep_treat_per_celltype.pdf",
    width=15, height=5)
ANPEP_per_Treat_per_celltype %>%
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



#Plot5: DPP4 downregulation in MERS-CoV exposed cells per Status

DPP4_per_DPP4_Status_per_celltype <- df1 %>%
  group_by(celltype, Status) %>%
  dplyr::count(DPP4_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            DPP4_Status = DPP4_Status) %>%
  mutate("Percentage" = n/Sum*100)

DPP4_per_DPP4_Status_per_celltype$Status <- factor(DPP4_per_DPP4_Status_per_celltype$Status, levels=c('Uninfected','Infected', 'Bystander'))

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_MERS_cells_dpp4_status_celltype.pdf",
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
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")+
  facet_grid(cols = vars(Status))
dev.off()



#Plot5.1: DPP4 downregulation in MERS-CoV exposed cells per Treat

DPP4_per_DPP4_Treat_per_celltype <- df1 %>%
  group_by(celltype, Treat) %>%
  dplyr::count(DPP4_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            DPP4_Status = DPP4_Status) %>%
  mutate("Percentage" = n/Sum*100)


DPP4_per_DPP4_Treat_per_celltype$Treat <- factor(DPP4_per_DPP4_Treat_per_celltype$Treat, levels=c('Mock','HCoV-229E'))


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_MERS_cells_dpp4_treat_celltype.pdf",
    width=14, height=6)
DPP4_per_DPP4_Treat_per_celltype %>%
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

###################################################################################################################################################

#Llama 229E , DPP4 and ANPEP counts in UMAP

###################################################################################################################################################
Alpaca <- Alpaca.clustered
DefaultAssay(Alpaca) <- "RNA"
Idents(Alpaca) <- "Treat" 

Alpaca_229E <- subset(x = Alpaca, idents = c("Mock", "dcCoV-ACN4"))
Idents(Alpaca_229E) <- "celltype" 
levels(Alpaca_229E)


#229E counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/UMAP_229E_counts.pdf",
    width=7, height=7)
FeaturePlot(Alpaca_229E, slot = "data", features = "camel229E-3UTR", cols = c("#F8A251", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T) +NoAxes()
dev.off()

#DPP4 counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/UMAP_229E_dpp4_counts.pdf",
    width=5, height=5)
FeaturePlot(Alpaca_229E, slot = "data", features = "DPP4", cols = c("#F8A251", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T) +NoAxes()
dev.off()

#ANPEP
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/UMAP_229E_anpep_counts.pdf",
    width=5, height=5)
FeaturePlot(Alpaca_229E, slot = "data", features = "ANPEP", cols = c("#F8A251", "#06606E"), pt.size = 1, label.size = 5, label = T, repel = T) +NoAxes()
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

Idents(Alpaca_229E) <- "Status" 

#Infection
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_Alpaca_both_virus/UMAP_status.pdf",
    width=5, height=5)
DimPlot(Alpaca_229E, reduction = "umap", split.by = "Treat", cols = cols_stat)
dev.off()

##################################################################################################################################################

#Threshold

##################################################################################################################################################
#########################################################################################################
#Correlation DPP4 expression and Viral_counts
#Extract the info in which row te DPP4 expression data is

a <- as.data.frame(Alpaca_229E@assays$SCT@data@Dimnames[[1]])

nrow(Alpaca_229E@assays$RNA@counts)
ncol(Alpaca_229E@assays$RNA@counts)

#2475
#show the specific location in the matrix for control
DPP4 <- as.matrix(Alpaca_229E@assays$SCT@counts[2475,])
Alpaca_229E <- AddMetaData(object = Alpaca_229E, metadata = DPP4, col.name = "Total.DPP4.Counts")


#Threshold for DPP4 expression
Alpaca_229E@meta.data <- Alpaca_229E@meta.data %>% mutate(DPP4_Status = case_when(Total.DPP4.Counts > 0 ~ "DPP4_positive", Total.DPP4.Counts == 0 ~ "DPP4_negative"))


#11227
#show the specific location in the matrix for control
ANPEP <- as.matrix(Alpaca_229E@assays$SCT@counts[11227,])
Alpaca_229E <- AddMetaData(object = Alpaca_229E, metadata = ANPEP, col.name = "Total.ANPEP.Counts")
#Threshold for ANPEP expression
Alpaca_229E@meta.data <- Alpaca_229E@meta.data %>% mutate(ANPEP_Status = case_when(Total.ANPEP.Counts > 0 ~ "ANPEP_positive", Total.ANPEP.Counts == 0 ~ "ANPEP_negative"))

head(Alpaca_229E@meta.data)


Idents(Alpaca_229E) <- "ANPEP_Status"

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/UMAP_ANPEP_status.pdf",
    width=8, height=7)
DimPlot(Alpaca_229E, reduction = "umap", label = F, label.size = 5, repel = TRUE, pt.size = 1, cols = c("#F8A251","#06606E")) + NoAxes()
dev.off()




df1 <- Alpaca_229E@meta.data %>% tidyseurat::select(c("celltype","Total.DPP4.Counts", "Total.ANPEP.Counts", "Status", "Treat", "DPP4_Status", "ANPEP_Status"))


#Plot 1: How many cells per celltype are Infected, Uninfected and Bystander
dcCoV229E_per_celltype <- df1 %>%
  group_by(celltype) %>%
  dplyr::count(Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            Status = Status) %>%
  mutate("Percentage" = n/Sum*100)


# Reordering Levels
dcCoV229E_per_celltype$Status <- factor(dcCoV229E_per_celltype$Status, levels=c('Infected','Bystander','Uninfected'))

cols_stat_2 = c("#CD4631","#519872","#A06136")
show_col(cols_stat_2)

cols_stat_1 = c("#CD4631","#519872")
show_col(cols_stat_1)

cols_anpep = c("#F8A251","#06606E")
show_col(cols_anpep)

cols_dpp4 = c("#F8A251","#9C4868")
show_col(cols_dpp4)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_Alpaca_both_virus/fraction_of_229E_cells_status_per_celltype.pdf",
    width=6, height=5)
dcCoV229E_per_celltype %>%
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



dcCoV229E_per_celltype_inf_bys <- df1 %>%
  filter(Status != "Uninfected") %>%
  group_by(celltype) %>%
  dplyr::count(Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            Status = Status) %>%
  mutate("Percentage" = n/Sum*100)

# Reordering Levels
dcCoV229E_per_celltype_inf_bys$Status <- factor(dcCoV229E_per_celltype_inf_bys$Status, levels=c('Infected','Bystander'))
dcCoV229E_per_celltype_inf_bys$celltype <- factor(dcCoV229E_per_celltype_inf_bys$celltype, levels = c("Llama Cluster 1", "Llama Cluster 2", "Secretory", "Ciliated", "Club", "Basal"))


#Plot 1.1: How many cells per celltype are Infected, Bystander
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/fraction_of_dcCoV229E_cells_inf_bys_per_celltype.pdf",
    width=7, height=6)
dcCoV229E_per_celltype_inf_bys %>%
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



pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_Alpaca_both_virus/fraction_of_229E_cells_dpp4_status_per_celltype.pdf",
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



pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_Alpaca_both_virus/fraction_of_229E_cells_anpep_status_per_celltype.pdf",
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




#Plot4: ANPEP downregulation in 229E-CoV exposed cells per Status

ANPEP_per_Status_per_celltype <- df1 %>%
  group_by(celltype, Status) %>%
  dplyr::count(ANPEP_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            ANPEP_Status = ANPEP_Status) %>%
  mutate("Percentage" = n/Sum*100)

ANPEP_per_Status_per_celltype$Status <- factor(ANPEP_per_Status_per_celltype$Status, levels=c('Uninfected', 'Infected','Bystander'))
ANPEP_per_Status_per_celltype$celltype <- factor(ANPEP_per_Status_per_celltype$celltype, levels=c('Basal', 'Secretory','Club', 'Ciliated', 'Llama Cluster 1', 'Llama Cluster 2'))


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/Plot_229E_ANPEP_llama.pdf",
    width=14, height=6)
ANPEP_per_Status_per_celltype %>%
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
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 15),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")+
  facet_grid(cols = vars(Status))
dev.off()



#Plot4.1: ANPEP downregulation in 229E-CoV exposed cells per Treat

ANPEP_per_Treat_per_celltype <- df1 %>%
  group_by(celltype, Treat) %>%
  dplyr::count(ANPEP_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            ANPEP_Status = ANPEP_Status) %>%
  mutate("Percentage" = n/Sum*100)

ANPEP_per_Status_per_celltype$Status <- factor(ANPEP_Status_per_celltype$Status, levels=c('Uninfected','Infected', 'Bystander'))


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_virus/fraction_of_229E_cells_anpep_status_per_celltype.pdf",
    width=15, height=5)
ANPEP_per_Treat_per_celltype %>%
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
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 15),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="top")+
  facet_grid(cols = vars(Treat))
dev.off()



#Plot5: DPP4 downregulation in MERS-CoV exposed cells per Status

DPP4_Status_per_celltype <- df1 %>%
  group_by(celltype, Status) %>%
  dplyr::count(DPP4_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            DPP4_Status = DPP4_Status) %>%
  mutate("Percentage" = n/Sum*100)

DPP4_Status_per_celltype$Status <- factor(DPP4_Status_per_celltype$Status, levels=c('Uninfected','Infected', 'Bystander'))
DPP4_Status_per_celltype$celltype <- factor(DPP4_Status_per_celltype$celltype, levels=c('Basal', 'Secretory','Club', 'Ciliated', 'Llama Cluster 1', 'Llama Cluster 2'))


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/Plot_229E_DPP4_llama.pdf",
    width=14, height=6)
DPP4_Status_per_celltype %>%
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
  facet_grid(cols = vars(Status))
dev.off()



#Plot5.1: DPP4 downregulation in MERS-CoV exposed cells per Treat

DPP4_per_DPP4_Treat_per_celltype <- df1 %>%
  group_by(celltype, Treat) %>%
  dplyr::count(DPP4_Status) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            DPP4_Status = DPP4_Status) %>%
  mutate("Percentage" = n/Sum*100)


DPP4_per_DPP4_Treat_per_celltype$Treat <- factor(DPP4_per_DPP4_Treat_per_celltype$Treat, levels=c('Mock','HCoV-229E'))


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_MERS_cells_dpp4_treat_celltype.pdf",
    width=14, height=6)
DPP4_per_DPP4_Treat_per_celltype %>%
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

