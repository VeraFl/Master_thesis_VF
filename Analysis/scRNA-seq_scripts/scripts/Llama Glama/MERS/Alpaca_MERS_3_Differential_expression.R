# Differential expression of Single-cell Seq data
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
cols_cell = hcl.colors(7, "Temps")
show_col(cols_cell)


##########################################################################################################################
#Loading filtered and merged Seurat Object into Script
#Infected and Mock Samples merged and cell clusters assigned

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/V.pacos")

MERS <- LoadH5Seurat("MERS.combined.h5seurat")
Assays(MERS)

##################################################################################################################
# "Bulk" differential expression in the Status groups "Infected", "Bystander", and "Non-infected" without the cell groups
Idents(MERS) <- "Status" #Parameters "Status" needs to be marked as a Idents in the object
MERS_inf_markers <- FindAllMarkers(MERS, only.pos = TRUE, min.pct = 0.25, grouping.var = "Status", logfc.threshold = 0.25)
MERS_inf_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)



# Find differentially expressed features between Infected and Non-infected cells
Infection_vs_Noninfection.markers <- FindMarkers(MERS, ident.1 = "Infected", ident.2 = "Non-infected")
# view results
head(Infection_vs_Noninfection.markers)

# Find differentially expressed features between Infected and Bystander cells
Infection_vs_bystander.markers <- FindMarkers(MERS, ident.1 = "Infected", ident.2 = "Bystander")
# view results
head(Infection_vs_bystander.markers)

# Find differentially expressed features between Non-Infected and Bystander cells
NonInfection_vs_bystander.markers <- FindMarkers(MERS, ident.1 = "Non-infected", ident.2 = "Bystander")
# view results
head(NonInfection_vs_bystander.markers)

# Find differentially expressed features between Non-Infected and the other two
NonInfection.markers <- FindMarkers(MERS, ident.1 = "Non-infected")
# view results
head(NonInfection.markers)
#############################################################################################################################################
DoHeatmap(MERS.combined, features = c("APOE", "KRT5"), group.by = "Status")

#############################################################################################################################################
#Difference between Status (Infected, Non-infected, Bystander) within celltypes
MERS@meta.data
MERS$celltype.stim <- paste(Idents(MERS), MERS$Status, sep = "_")
MERS$celltype <- Idents(MERS)
Idents(MERS) <- "celltype.stim"
ciliated.response <- FindMarkers(MERS, ident.1 = "Ciliated_Non-infected", ident.2 = "Ciliated_Bystander", verbose = FALSE)
head(ciliated.response, n = 15)

#Of markers differently expressed in non-infected ciliated cells and bystander ciliated cells we can see different expression patterns in all the cell types:
plots <- VlnPlot(MERS, features = c("MYC", "ISG15", "ADSS2", "STC1", "LOC102523265"), split.by = "Status", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)

VlnPlot(MERS, features = c("ISG15","LOC102523265"), split.by = "Status", group.by = "celltype", 
        pt.size = 0, combine = FALSE)



#How many of each celltype are infected?
Idents(MERS) <- "Status"
Infected <- subset(MERS, idents = "Infected")
df <- MERS@meta.data


plot_count <- df %>%
  group_by(celltype, Status) %>%
  count(Status)



#Plot2
plot_count %>% ggplot(aes(x = Status, y = n, fill = Status))+
  geom_col()+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        legend.position="bottom")+
  facet_grid(cols = vars(celltype))


