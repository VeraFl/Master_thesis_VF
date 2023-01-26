# Differential expression of Single-cell Seq data
# 20.06.2022

#load libraries
library(dplyr)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(patchwork)
library(readr)
library(gdata)
library(Matrix)
library(writexl)

##########################################################################################################################
#Loading filtered and merged Seurat Object into Script
#Infected and Mock Samples merged and cell clusters assigned

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens")

MERS <- LoadH5Seurat("MERS.h5seurat")
MERS

MERS@meta.data

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

DoHeatmap(MERS, features = c("VIM", "APOE", "CALD1"), group.by = "Status")

#Of markers differently expressed in non-infected ciliated cells and bystander ciliated cells we can see different expression patterns in all the cell types:
plots <- VlnPlot(MERS, features = c("TENM2", "PGF", "C1S", "MMP10", "COL3A1", "ARHGAP22"), split.by = "Status", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)

##################################################################################################################
#Difference between Status (Infected, Non-infected, Bystander) within celltypes
MERS$celltype.stim <- paste(Idents(MERS), MERS$Status, sep = "_")
MERS$celltype <- Idents(MERS)
Idents(MERS) <- "celltype.stim"
ciliated.response <- FindMarkers(MERS, ident.1 = "Ciliated_Non-infected", ident.2 = "Ciliated_Bystander", verbose = FALSE)
head(ciliated.response, n = 15)

#Of markers differently expressed in non-infected ciliated cells and bystander ciliated cells we can see different expression patterns in all the cell types:
plots <- VlnPlot(MERS, features = c("CXCL17", "SLPI", "CXCL8"), split.by = "Status", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)


#Of the low number of infected cells, which cell types are these cells?
Idents(MERS) <- "Status"
Infected <- subset(MERS, idents = "Infected")
Infected@meta.data


