# Analysis of Single-cell Seq data
# 10.05.2022
# Count data of Infected Samples (3_A) with concatanated custom Reference Genome "Homo_sapiens.GRCh38.105.MERS_CoV.gtf"
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
# #Assign a theme
# 
# # Use colourblind-friendly colours
# friendly_cols <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51")
# 
# # Set theme
# my_theme <-
#   list(
#     scale_fill_manual(values = friendly_cols),
#     scale_color_manual(values = friendly_cols),
#     theme_bw() +
#       theme(
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         panel.grid.major = element_line(size = 0.2),
#         panel.grid.minor = element_line(size = 0.1),
#         text = element_text(size = 12),
#         legend.position = "bottom",
#         aspect.ratio = 1,
#         strip.background = element_blank(),
#         #axis.title.x = element_text(margin = margin(t = 1, r = 5, b = 5, l = 10)),
#         #axis.title.y = element_text(margin = margin(t = 1, r = 5, b = 5, l = 10))
#       )
#   )

###################################################################################################################################

# Import data which was downloaded from UBELIX after cellranger count run
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/")
# Create each individual Seurat object (if multiple samples)
for (file in c("mock_filtered_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("data/H. sapiens/MERS/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

###################################################################################################################################
# QC
counts_per_gene <- Matrix::rowSums(seurat_data)
genes_per_cell <- Matrix::colSums(seurat_data)
nrow(seurat_data)
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered) - MERS - MOCK')

# Add column for Mitochondrial Content in percentage to the meta data
mock_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(mock_filtered_feature_bc_matrix, pattern = "^MT-")
head(mock_filtered_feature_bc_matrix@meta.data)

?AddMetaData

# Make violin plot to show number of features, number of counts and mitochondiral content
VlnPlot(mock_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

?PercentageFeatureSet
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# We filter cells that have unique feature counts over
# We filter cells that have >5% mitochondrial counts

mock_filtered_feature_bc_matrix <- subset(mock_filtered_feature_bc_matrix, subset = nFeature_RNA > 1000 & percent.mt < 30)

# Make violin plot to show number of features, number of counts and mitochondiral content after filtering
VlnPlot(mock_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2



#Normalization

mock_filtered_feature_bc_matrix <- NormalizeData(mock_filtered_feature_bc_matrix, normalization.method = "LogNormalize", scale.factor = 10000)


# Identification of highly variable features (feature selection)

mock_filtered_feature_bc_matrix <- FindVariableFeatures(mock_filtered_feature_bc_matrix, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mock_filtered_feature_bc_matrix), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mock_filtered_feature_bc_matrix)
plot2 <- LabelPoints(plot = plot1, points = top10)
plot1 + plot2

#Scaling the data
all.genes <- rownames(mock_filtered_feature_bc_matrix)
mock_filtered_feature_bc_matrix <- ScaleData(mock_filtered_feature_bc_matrix, features = all.genes)

##################################################################################################################
#PCA
mock_filtered_feature_bc_matrix <- RunPCA(mock_filtered_feature_bc_matrix, features = VariableFeatures(object = mock_filtered_feature_bc_matrix))

# Examine and visualize PCA results a few different ways
print(mock_filtered_feature_bc_matrix[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mock_filtered_feature_bc_matrix, dims = 1:2, reduction = "pca")
DimPlot(mock_filtered_feature_bc_matrix, reduction = "pca")

#Heatmap
DimHeatmap(mock_filtered_feature_bc_matrix, dims = 1, cells = 500, balanced = TRUE)

##################################################################################################################
#Excluding cell cycle state markers
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Regressing out cell-cycle markers (same procedure for viral ones later)
mock_filtered_feature_bc_matrix <- CellCycleScoring(mock_filtered_feature_bc_matrix, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(mock_filtered_feature_bc_matrix[[]])

# Visualize the distribution of cell cycle markers across
RidgePlot(mock_filtered_feature_bc_matrix, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

#Running a PCA on cell cycle genes reveals, it looks ok
mock_filtered_feature_bc_matrix <- RunPCA(mock_filtered_feature_bc_matrix, features = c(s.genes, g2m.genes))
DimPlot(mock_filtered_feature_bc_matrix)

# Scaling the data
mock_filtered_feature_bc_matrix <- ScaleData(mock_filtered_feature_bc_matrix, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(mock_filtered_feature_bc_matrix))


##################################################################################################################
#UMAP
mock_filtered_feature_bc_matrix <- FindNeighbors(mock_filtered_feature_bc_matrix, dims = 1:20)
mock_filtered_feature_bc_matrix <- FindClusters(mock_filtered_feature_bc_matrix, resolution = 0.2)

mock_filtered_feature_bc_matrix <- RunUMAP(mock_filtered_feature_bc_matrix, dims = 1:20)

DimPlot(mock_filtered_feature_bc_matrix, reduction = "umap")

##################################################################################################################
# Markers
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
mock_filtered_feature_bc_matrix_markers <- FindAllMarkers(mock_filtered_feature_bc_matrix, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mock_filtered_feature_bc_matrix_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#Plotting expression of different markers (making the clusters)

#BASAL: KRT5, TP63, DLK2, KRT6A
#Suprabasal: KRT4, KRT13, KRT16, KRT23
#Club: SCGB1A1, KRT7, KRT19
#Deuterosomal: PLK4, CCNO, CEP78, DEUP1
#Secretory: MUC5B, MUC5AC
#Goblet: FOXQ1, SPDEF, MUC5AC
#Ciliated: FOXJ1, TPPP3, SNTN, TUBA1A, CETN2, SPEF2, PIFO, LRRC6
#Tuft: POUF2f3, TRPM5, LRMP, RGS13, HOXC5, HMX2, ANXA4
#PNEC: PSMD5, NGF, PCSK1N, SCGN, NEB, HOXB1, ASCL1/2, FOXA2
#Ionocytes: ASCL3, CFTR, FOXI1, DMRT2, V-ATPase

FeaturePlot(mock_filtered_feature_bc_matrix, features = "APN")

#BASAL
FeaturePlot(mock_filtered_feature_bc_matrix, features = "KRT5") #basal1
FeaturePlot(mock_filtered_feature_bc_matrix, features = "TP63") #basal2
FeaturePlot(mock_filtered_feature_bc_matrix, features = "FN1") #basal3
FeaturePlot(mock_filtered_feature_bc_matrix, features = "DLK2") #basal4
FeaturePlot(mock_filtered_feature_bc_matrix, features = "KRT6A") #basal5

#SUPRABASAL
FeaturePlot(mock_filtered_feature_bc_matrix, features = "KRT4") #suprabasal1
FeaturePlot(mock_filtered_feature_bc_matrix, features = "KRT13") #suprabasal2
FeaturePlot(mock_filtered_feature_bc_matrix, features = "KRT16") #suprabasal3
FeaturePlot(mock_filtered_feature_bc_matrix, features = "KRT23") #suprabasal4

#TUFT
FeaturePlot(mock_filtered_feature_bc_matrix, features = "POUF2f3") #tuft1
FeaturePlot(mock_filtered_feature_bc_matrix, features = "TRPM5") #tuft2
FeaturePlot(mock_filtered_feature_bc_matrix, features = "LRMP") #tuft3

#IONOCYTES
FeaturePlot(mock_filtered_feature_bc_matrix, features = "FOXI1") #ionocytes1
FeaturePlot(mock_filtered_feature_bc_matrix, features = "CFTR") #ionocytes2

##PNEC
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("PCSK1N")) #PNEC1
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("SCGN")) #PNEC2
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("NEB")) #PNEC3
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("HOXB1")) #PNEC4
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("ASCL1")) #PNEC5
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("ASCL2")) #PNEC6
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("FOXA2")) #PNEC7
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("PSMD5")) #PNEC8
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("NGF")) #PNEC9

#CLUB
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("PCSK1N")) #PNEC1
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("SCGN")) #PNEC2
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("NEB")) #PNEC3

#GOBLET/SECRETORY
FeaturePlot(mock_filtered_feature_bc_matrix, features = "MUC5B") #secretory1
FeaturePlot(mock_filtered_feature_bc_matrix, features = "SCGB3A1") #secretory2
FeaturePlot(mock_filtered_feature_bc_matrix, features = "MUC5AC") #secretory3 (goblet)
FeaturePlot(mock_filtered_feature_bc_matrix, features = "FOXQ1") #goblet1
FeaturePlot(mock_filtered_feature_bc_matrix, features = "SPDEF") #goblet2

#DEUTEROSOMAL
FeaturePlot(mock_filtered_feature_bc_matrix, features = "PLK4") #ciliated1
FeaturePlot(mock_filtered_feature_bc_matrix, features = "CCNO") #ciliated2
FeaturePlot(mock_filtered_feature_bc_matrix, features = "CEP78") #ciliated3
FeaturePlot(mock_filtered_feature_bc_matrix, features = "DEUP1") #ciliated1


#CILIATED
FeaturePlot(mock_filtered_feature_bc_matrix, features = "FOXJ1") #ciliated1
FeaturePlot(mock_filtered_feature_bc_matrix, features = "PIFO") #ciliated2
FeaturePlot(mock_filtered_feature_bc_matrix, features = "TPPP3") #ciliated3
FeaturePlot(mock_filtered_feature_bc_matrix, features = "SPEF2") #ciliated1
FeaturePlot(mock_filtered_feature_bc_matrix, features = "DNAH5") #ciliated2
FeaturePlot(mock_filtered_feature_bc_matrix, features = "LRRC6") #ciliated3


#############################################################################################################################################

#Cluster0
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("orf1ab", "S", "orf3", "orf4a", "orf4b", "orf5", "E", "M", "orf8b")) #viral genes
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("PCNA", "TOP2A", "MCM6", "MKI67")) #cell cycle state genes
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("MUC5B", "MUC5AC")) #secretory
FeaturePlot(mock_filtered_feature_bc_matrix, features = "MUC5B") #secretory1
FeaturePlot(mock_filtered_feature_bc_matrix, features = "SCGB3A1") #secretory2
FeaturePlot(mock_filtered_feature_bc_matrix, features = "MUC5AC") #secretory3 (goblet)
#Cluster1
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("FOXJ1", "TPPP3", "SNTN", "TUBA1A", "CETN2")) #ciliated
FeaturePlot(mock_filtered_feature_bc_matrix, features = "FOXJ1") #ciliated1
FeaturePlot(mock_filtered_feature_bc_matrix, features = "PIFO") #ciliated2
FeaturePlot(mock_filtered_feature_bc_matrix, features = "TPPP3") #ciliated3

#Cluster2
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("KRT5", "TP63", "DLK2", "KRT6A")) #basal
FeaturePlot(mock_filtered_feature_bc_matrix, features = "KRT5") #basal1
FeaturePlot(mock_filtered_feature_bc_matrix, features = "TP63") #basal2
FeaturePlot(mock_filtered_feature_bc_matrix, features = "FN1") #basal3
#Cluster3
FeaturePlot(mock_filtered_feature_bc_matrix, features = "BPIFB1")


FeaturePlot(mock_filtered_feature_bc_matrix, features = c("PCSK1N", "SCGN", "NEB", "HOXB1", "ASCL1", "ASCL2", "FOXA2")) #PNEC
#Cluster5
FeaturePlot(mock_filtered_feature_bc_matrix, features = c("ASCL3", "CFTR", "FOXI1", "DMRT2", "V-ATPase")) #ionocytes
FeaturePlot(mock_filtered_feature_bc_matrix, features = "IGFBP5") #ionocytes1
FeaturePlot(mock_filtered_feature_bc_matrix, features = "CFTR") #ionocytes2
