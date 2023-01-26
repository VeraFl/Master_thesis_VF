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

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.ferus")

MERS <- LoadH5Seurat("MERS.clustered.h5seurat")
MERS

##################################################################################################################
# differential expression in the Status groups "Infected", "Bystander", and "Non-infected"
Idents(MERS) <- "Status" #Parameters "Status" needs to be marked as a Idents in the object
MERS_inf_markers <- FindAllMarkers(MERS, only.pos = TRUE, min.pct = 0.25, grouping.var = "Status", logfc.threshold = 0.25)
MERS_inf_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results")
write.xlsx(MERS.combined_inf_markers,"Camel_ferus_MERS_Infection_markers.xlsx")

#Ask what the levels are of the object
levels(MERS)


# Find differentially expressed features between Infected and Non-infected cells
Infection_vs_Noninfection.markers <- FindMarkers(MERS, ident.1 = "Infected", ident.2 = "Non-infected")
# view results
head(Infection_vs_Noninfection.markers)

# Find differentially expressed features between Infected and Bystander cells
Infection_vs_bystander.markers <- FindMarkers(MERS, ident.1 = "Infected", ident.2 = "Bystander")
# view results
head(Infection_vs_bystander.markers)

# Find differentially expressed features between Non-Infected and Bystander cells
NonInfection_vs_bystander.markers <- FindMarkers(MERS.combined, ident.1 = "Non-infected", ident.2 = "Bystander")
# view results
head(NonInfection_vs_bystander.markers)

# Find differentially expressed features between Non-Infected and the other two
NonInfection.markers <- FindMarkers(MERS.combined, ident.1 = "Non-infected")
# view results
head(NonInfection.markers)
#############################################################################################################################################

DoHeatmap(MERS.combined, features = c("RSAD2", "IFIT3"), group.by = "Status")

