if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")


install.packages("remotes")
library(remotes)
remotes::install_github("satijalab/seurat-wrappers")

devtools::install_github('cole-trapnell-lab/monocle3')

library(dplyr)
library(SeuratData)
library(SeuratDisk)
library(SeuratWrappers)
library(Seurat)
library(patchwork)
library(readr)
library(gdata)
library(Matrix)
library(ggplot2)
library(magrittr)
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
library(cowplot)
library(monocle3)
library(devtools)


#Import integrated dataset
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens/")
integrated <- LoadH5Seurat("Human.clustered.h5seurat")
integrated <- Ferus

cds <- as.cell_data_set(integrated)
## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

## Step 4: Cluster the cells
cds <- cluster_cells(cds)

## Step 5: Learn a graph
cds <- learn_graph(cds)

## Step 6: Order cells
cds <- order_cells(cds)

plot_cells(cds)

cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)
p1
p2

as.Seurat(integrated)

integrated.sub <- base::subset(SeuratWrappers::as.Seurat(cds, assay = "integrated"), monocle3_partitions == 1)
cds <- as.cell_data_set(cds)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

?as.Seurat()
?subset
subset


max.avp <- which.max(unlist(FetchData(cds, "AVP")))
max.avp <- colnames(integrated.sub)[max.avp]
cds <- order_cells(cds, root_cells = max.avp)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)


# Set the assay back as 'integrated'
integrated.sub <- as.Seurat(cds, assay = "integrated")
FeaturePlot(integrated.sub, "monocle3_pseudotime")

