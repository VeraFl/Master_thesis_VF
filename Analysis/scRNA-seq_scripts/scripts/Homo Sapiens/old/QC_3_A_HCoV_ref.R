# Analysis of Single-cell Seq data
# 17.03.2022
# Count data of Mock Samples (3_A) with concatanated custom Reference Genome "Homo_sapiens.GRCh38.105.HCoV_229E.gtf"
# Count matrix produced by "Cell_Ranger_Count_3_A.sh"

#load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(readr)
library(gdata)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(tidyseurat)
###################################################################################################################################
#Assign a theme

# Use colourblind-friendly colours
friendly_cols <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51")

# Set theme
my_theme <-
  list(
    scale_fill_manual(values = friendly_cols),
    scale_color_manual(values = friendly_cols),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 12),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_blank(),
        #axis.title.x = element_text(margin = margin(t = 1, r = 5, b = 5, l = 10)),
        #axis.title.y = element_text(margin = margin(t = 1, r = 5, b = 5, l = 10))
      )
  )

###################################################################################################################################
getwd()
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/")
# import data which was downloaded from ubelix after cellranger count run
# filtered_feature_bc_matrix


for (file in c("mock_filtered_feature_bc_matrix")){
  HCoV_mock_data <- Read10X(data.dir = paste0("data/H. sapiens/HCoV/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}



#Importing Cell cycle marker genes to regress out later
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Giving a cell cycle score to every cell
mock_filtered_feature_bc_matrix <- CellCycleScoring(mock_filtered_feature_bc_matrix, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(mock_filtered_feature_bc_matrix[[]])

# #Making the seurat object to a tidyseurat
# HCoV_mock <- tidy(mock_filtered_feature_bc_matrix)
# HCoV_mock

###################################################################################################################################
# HCoV_mock %>%
#   tidyseurat::ggplot(aes(nFeature_RNA, fill = Phase)) +
#   geom_histogram() +
#   my_theme


# QC

#Calculating some QC statistics
#goes through every row and adds all the counts of the columns (molecules per cell)
counts_per_cell <- Matrix::colSums(HCoV_mock_data) 
#goes through every column and adds all the counts of the rows (molecules per gene)
counts_per_gene <- Matrix::rowSums(HCoV_mock_data) 
#count gene itself only if it has non-zero reads mapped (not value).
genes_per_cell <- Matrix::colSums(HCoV_mock_data>0) # only count gene if it has non-zero reads mapped (not value).
cells_per_gene <- Matrix::rowSums(HCoV_mock_data>0) # only count cells where the gene is expressed

# #Counts per cell before filtering and before excluding cell cycle markers
# HCoV_mock %>%
#   tidyseurat::ggplot(aes(log10(counts_per_cell+1),main='counts per cell', fill = Phase)) +
#   geom_histogram() +
#   my_theme
# 
# #Genes per cell before filtering and before excluding cell cycle markers
# HCoV_mock %>%
#   tidyseurat::ggplot(aes(log10(genes_per_cell+1), main='genes per cell', fill = Phase)) +
#   geom_histogram() +
#   my_theme
# 
# 
# HCoV_mock %>%
#   tidyseurat::ggplot() +
#   geom_point(aes(x = counts_per_cell, y = genes_per_cell), group = "Phase") +
#   my_theme


#Genes per cell ordered
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered) - HCoV - MOCK')


#################################################################################################################
# Add column for Mitochondrial Content in percentage to the meta data
mock_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(mock_filtered_feature_bc_matrix, pattern = "^MT-")
head(mock_filtered_feature_bc_matrix@meta.data)


# Make violin plot to show number of features, number of counts and mitochondiral content
VlnPlot(mock_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

?PercentageFeatureSet
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

?FeatureScatter

# We filter cells that have unique feature counts over
# We filter cells that have >5% mitochondrial counts

mock_filtered_feature_bc_matrix <- subset(mock_filtered_feature_bc_matrix, subset = nFeature_RNA > 1000 & percent.mt < 30)

# Make violin plot to show number of features, number of counts and mitochondiral content after filtering
VlnPlot(mock_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


#Normalization

Sample_3_A <- NormalizeData(Sample_3_A, normalization.method = "LogNormalize", scale.factor = 10000)


# Identification of highly variable features (feature selection)

Sample_3_A <- FindVariableFeatures(Sample_3_A, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Sample_3_A), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Sample_3_A)
plot2 <- LabelPoints(plot = plot1, points = top10)
plot1 + plot2

#Scaling the data
all.genes <- rownames(Sample_3_A)
Sample_3_A <- ScaleData(Sample_3_A, features = all.genes)

###################################################################################################################################
# Dimensionality Reduction - PCA Run
Sample_3_A <- RunPCA(Sample_3_A, features = VariableFeatures(object = Sample_3_A))

# Examine and visualize PCA results a few different ways
print(Sample_3_A[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Sample_3_A, dims = 1:2, reduction = "pca")
DimPlot(Sample_3_A, reduction = "pca")

#Heatmap
DimHeatmap(Sample_3_A, dims = 1, cells = 500, balanced = TRUE)


###################################################################################################################################
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Regressing out cell-cycle markers (same procedure for viral ones later)
Sample_3_A <- CellCycleScoring(Sample_3_A, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(Sample_3_A[[]])

# Visualize the distribution of cell cycle markers across
RidgePlot(Sample_3_A, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

#Running a PCA on cell cycle genes reveals, unsurprisingly, that cells seperate entirely by phase
Sample_3_A <- RunPCA(Sample_3_A, features = c(s.genes, g2m.genes))
DimPlot(Sample_3_A)

Sample_3_A <- ScaleData(Sample_3_A, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Sample_3_A))

Sample_3_A <- RunPCA(Sample_3_A, features = c(s.genes, g2m.genes))
DimPlot(Sample_3_A)

#Finding Neighbors and Clusters
Sample_3_A <- FindNeighbors(Sample_3_A, dims = 1:15)
Sample_3_A <- FindClusters(Sample_3_A, resolution = 0.2)

Sample_3_A <- RunPCA(Sample_3_A, features = c(s.genes, g2m.genes))
DimPlot(Sample_3_A)
###################################################################################################################################
# Dimensionality Reduction - UMAP Run

#Finding Neighbors and Clusters
Sample_3_A <- FindNeighbors(Sample_3_A, dims = 1:15)
Sample_3_A <- FindClusters(Sample_3_A, resolution = 0.5)

#Run UMAP
Sample_3_A <- RunUMAP(Sample_3_A, dims = 1:15)
DimPlot(Sample_3_A, reduction = "umap")

?RunUMAP

#Finding Markers to respective Clusters made above to differentiate between clusters
cluster2.markers <- FindMarkers(Sample_3_A, ident.1  = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

cluster5.markers <- FindMarkers(Sample_3_A, ident.1 = 2, ident.2 = 0, min.pct = 0.25)
head(cluster5.markers, n = 5)



# Find all Markers
# find markers for every cluster compared to all remaining cells, report only the positive ones
Sample_3_A.markers <- FindAllMarkers(Sample_3_A, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Sample_3_A.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#Plotting expression of different markers (making the clusters)

VlnPlot(Sample_3_A, features = c("SCGB3A1", "GPRC5A", "SLPI"))
VlnPlot(Sample_3_A, features = c("SCGB3A1", "GPRC5A"), slot = "counts", log = TRUE)

#KRT5, TP63, DLK2, KRT6A = Basal
#MUC5B, MUC5AC = Secretory cells
#FOXJ1, TPPP3, SNTN, TUBA1A, CETN2 = ciliated
#LRMP, RGS13, HOXC5, HMX2, ANXA4 = Tuft
#PCSK1N, SCGN, NEB, HOXB1, ASCL1/2, FOXA2 = PNEC
#ASCL3, CFTR, FOXI1, DMRT2, V-ATPase = ionocytes

#Cluster0
FeaturePlot(Sample_3_A, features = c("MUC5B", "MUC5AC")) #secretory
FeaturePlot(Sample_3_A, features = "MUC5B") #secretory1
FeaturePlot(Sample_3_A, features = "SCGB3A1") #secretory2
FeaturePlot(Sample_3_A, features = "MUC5AC") #secretory3 (goblet)
FeaturePlot(Sample_3_A, features = "SCGB1A1") #secretory4 (club)

#Cluster1
FeaturePlot(Sample_3_A, features = c("FOXJ1", "TPPP3", "SNTN", "TUBA1A", "CETN2")) #ciliated
FeaturePlot(Sample_3_A, features = "FOXJ1") #ciliated1
FeaturePlot(Sample_3_A, features = "PIFO") #ciliated2
FeaturePlot(Sample_3_A, features = "TPPP3") #ciliated3

#Cluster2
FeaturePlot(Sample_3_A, features = c("KRT5", "TP63", "DLK2", "KRT6A")) #basal
FeaturePlot(Sample_3_A, features = "KRT5") #basal1
FeaturePlot(Sample_3_A, features = "TP63") #basal2
FeaturePlot(Sample_3_A, features = "FN1") #basal3
#Cluster3
FeaturePlot(Sample_3_A, features = "BPIFB1")



FeaturePlot(Sample_3_A, features = c("PCSK1N", "SCGN", "NEB", "HOXB1", "ASCL1", "ASCL2", "FOXA2")) #PNEC
#Cluster5
FeaturePlot(Sample_3_A, features = c("ASCL3", "CFTR", "FOXI1", "DMRT2", "V-ATPase")) #ionocytes
FeaturePlot(Sample_3_A, features = "IGFBP5") #ionocytes1
FeaturePlot(Sample_3_A, features = "CFTR") #ionocytes2

FeaturePlot(Sample_3_A, features = "BPIFB1")

RidgePlot(Sample_3_A, features = "SCGB3A1")


grid.arrange(p1 + p2, ncol = 2)

#all the best markers together
FeaturePlot(Sample_3_A, features = c("KRT5", "MUC5B", "FOXJ1","AcTub", "PCSK1N", "ASCL3", "FOXI1", "TP63","KRT5","CFTR"))


#Assign new clusters
new.cluster.ids <- c("Secretory", "Ciliated", "Basal", "Tuft", "PNEC","Ionocytes")
names(new.cluster.ids) <- levels(Sample_3_A)
Sample_3_A <- RenameIdents(Sample_3_A, new.cluster.ids)
DimPlot(Sample_3_A, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

levels(Sample_3_A)


