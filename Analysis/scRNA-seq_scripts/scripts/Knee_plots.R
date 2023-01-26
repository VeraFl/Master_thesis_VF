#Barcode Inflection plot

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
library(DropletUtils)

#Treatment colors
cols_treat = c("#84A59D","#A5668B", "#96ADC8")
show_col(cols_treat)
#Status colors
cols_stat = c("#7C5869","#6E969F")
show_col(cols_stat)
#Celltype colors
cols_cell = hcl.colors(7, "Temps")
show_col(cols_cell)

###################################################################################################################################

# Import data which was downloaded from UBELIX after cellranger count run
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/")
# Create each individual Seurat object (if multiple samples)
for (file in c("mock_raw_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("data/H. Sapiens/both_virus/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   #min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

#"mock_filtered_feature_bc_matrix", "hcov_inf_filtered_feature_bc_matrix", "mers_inf_filtered_feature_bc_matrix"

# # Check the metadata in the new Seurat objects
# head(hcov_inf_filtered_feature_bc_matrix@meta.data)
# head(mers_inf_filtered_feature_bc_matrix@meta.data)
# head(mock_filtered_feature_bc_matrix@meta.data)


mock_raw_feature_bc_matrix <- CalculateBarcodeInflections(mock_raw_feature_bc_matrix,threshold.high = 15000, threshold.low = 10)
p1 <- BarcodeInflectionsPlot(mock_raw_feature_bc_matrix) + NoLegend()
p1 + 
  scale_x_log10() +
  scale_y_log10() +
  theme(axis.text = element_text(size = 12))




