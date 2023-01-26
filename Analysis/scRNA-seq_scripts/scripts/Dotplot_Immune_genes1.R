#Dotplot of immune genes over both species and all conditions and all celltypes
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
library(EnhancedVolcano)
library(ggrepel)
library(gprofiler2)
library(biomaRt)
library(gprofiler2)
library(ggVennDiagram)
library(ggupset)
library(MAST)

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.ferus")
Ferus <- LoadH5Seurat("Ferus.clustered.h5seurat")
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens")
Human <- LoadH5Seurat("Human.clustered.h5seurat")
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/V.pacos")
Alpaca <- LoadH5Seurat("Alpaca.clustered.h5seurat")




features1 <- c("IFNB1", "IFNL1", "IFNL2", "IFNL3", "IFNAR1", "IFNLR1", "IFI27", "IFITM3",
               "IFI6", "IFIT1", "MX1", "ISG15", "CXCL9", "CXCL10", "CXCL11","CXCL16", "CCL2",
               "IL1A", "IL1B", "IL1RN", "IL6", "IL10", "TNF")

features2 <- c("DPP4", "ANPEP", "ACN2")

Ferus@meta.data


DefaultAssay(Ferus) <- "SCT"
DefaultAssay(Alpaca) <- "SCT"
#Increase memory limit
memory.limit(24000)
# Find differentially expressed features between Infected and uninfected for one Treatment group
# MERS-CoV
str(Ferus)
Ferus@meta.data
Alpaca@meta.data
Idents(Ferus) <- "celltype"
Idents(Alpaca) <- "celltype"

#Make a fused new identity of each treatment and status combination
Ferus$cell.treat.status <- paste(Idents(Ferus), Ferus$Treat, Ferus$Status, sep = "_")
Alpaca$cell.treat.status <- paste(Idents(Alpaca), Alpaca$Treat, Alpaca$Status, sep = "_")


#make it the new identity
Idents(Ferus) <- "cell.treat.status"
Idents(Alpaca) <- "cell.treat.status"
#check if levels are correct
levels(Ferus)
levels(Alpaca)

DotPlot(Ferus, features = features1, scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "top")

DotPlot(Alpaca, features = features1, scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "top")



DotPlot(Ferus, features = features2, scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "top")

DotPlot(Alpaca, features = features2, scale.by = "radius", scale = T, cols = c("orange", "darkblue")) +coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "top")

?DotPlot
