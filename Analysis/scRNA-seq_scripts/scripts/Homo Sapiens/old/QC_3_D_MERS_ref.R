# Analysis of Single-cell Seq data
# 10.05.2022
# Count data of Infected Samples (3_D) with concatanated custom Reference Genome "Homo_sapiens.GRCh38.105.MERS_CoV.gtf"
# Count matrix produced by "Cell_Ranger_Count_3_D.sh"

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

# Import data which was downloaded from UBELIX after cellranger count run
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/")
# Create each individual Seurat object (if multiple samples)
for (file in c("inf_filtered_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("data/H. sapiens/MERS/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

# for (file in c("mock_raw_feature_bc_matrix", "inf_raw_feature_bc_matrix")){
#   seurat_data <- Read10X(data.dir = paste0("data/H. sapiens/MERS/", file))
#   seurat_obj <- CreateSeuratObject(counts = seurat_data, 
#                                    min.features = 100, 
#                                    project = file)
#   assign(file, seurat_obj)
# }

inf_raw_data <- Read10X(data.dir = "data/H. sapiens/MERS/inf_raw_feature_bc_matrix")
inf_filtered_data <- Read10X(data.dir = "data/H. sapiens/MERS/inf_filtered_feature_bc_matrix")

inf_raw <- CreateSeuratObject(counts = inf_raw_data, project = "inf", min.cells = 3, min.features = 200)
inf_filtered <- CreateSeuratObject(counts = inf_filtered_data, project = "inf", min.cells = 3, min.features = 200)


#calculating counts per gene and checking the viral genes if we have counts for them.
#looking at raw and filtered matrices to see a difference
raw_counts_per_gene <- Matrix::rowSums(inf_raw_data)
filtered_counts_per_gene <- Matrix::rowSums(inf_filtered_data)
virus_genes_raw <- tail(raw_counts_per_gene, 11) %>% data.frame()
virus_genes_filtered <- tail(filtered_counts_per_gene, 11) %>% data.frame()

#Calculating some QC statistics
#goes through every row and adds all the counts of the columns (molecules per cell)
counts_per_cell <- Matrix::colSums(seurat_data) 
nrow(seurat_data)
#goes through every column and adds all the counts of the rows (molecules per gene)
counts_per_gene <- Matrix::rowSums(seurat_data) 
#count gene itself only if it has non-zero reads mapped (not value).
genes_per_cell <- Matrix::colSums(seurat_data>0) # only count gene if it has non-zero reads mapped (not value).
cells_per_gene <- Matrix::rowSums(seurat_data>0) # only count cells where the gene is expressed




###################################################################################################################################
# QC

plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered) - MERS - INF')

# Add column for Mitochondrial Content in percentage to the meta data
inf_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(inf_filtered_feature_bc_matrix, pattern = "^MT-")
head(inf_filtered_feature_bc_matrix@meta.data)

?AddMetaData

# Make violin plot to show number of features, number of counts and mitochondiral content
VlnPlot(inf_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

?PercentageFeatureSet
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# We filter cells that have unique feature counts over
# We filter cells that have >5% mitochondrial counts

inf_filtered_feature_bc_matrix <- subset(inf_filtered_feature_bc_matrix, subset = nFeature_RNA > 1000 & percent.mt < 30)

# Make violin plot to show number of features, number of counts and mitochondiral content after filtering
VlnPlot(inf_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


#Normalization

inf_filtered_feature_bc_matrix <- NormalizeData(inf_filtered_feature_bc_matrix, normalization.method = "LogNormalize", scale.factor = 10000)


# Identification of highly variable features (feature selection)

inf_filtered_feature_bc_matrix <- FindVariableFeatures(inf_filtered_feature_bc_matrix, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(inf_filtered_feature_bc_matrix), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(inf_filtered_feature_bc_matrix)
plot2 <- LabelPoints(plot = plot1, points = top10)
plot1 + plot2

#Scaling the data
all.genes <- rownames(inf_filtered_feature_bc_matrix)
inf_filtered_feature_bc_matrix <- ScaleData(inf_filtered_feature_bc_matrix, features = all.genes)

##################################################################################################################
#PCA
inf_filtered_feature_bc_matrix <- RunPCA(inf_filtered_feature_bc_matrix, features = VariableFeatures(object = inf_filtered_feature_bc_matrix))

# Examine and visualize PCA results a few different ways
print(inf_filtered_feature_bc_matrix[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(inf_filtered_feature_bc_matrix, dims = 1:2, reduction = "pca")
DimPlot(inf_filtered_feature_bc_matrix, reduction = "pca")

#Heatmap
DimHeatmap(inf_filtered_feature_bc_matrix, dims = 1, cells = 500, balanced = TRUE)

##################################################################################################################
#Excluding cell cycle state markers
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Regressing out cell-cycle markers (same procedure for viral ones later)
inf_filtered_feature_bc_matrix <- CellCycleScoring(inf_filtered_feature_bc_matrix, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(inf_filtered_feature_bc_matrix[[]])

# Visualize the distribution of cell cycle markers across
RidgePlot(inf_filtered_feature_bc_matrix, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

#Running a PCA on cell cycle genes reveals, it looks ok
inf_filtered_feature_bc_matrix <- RunPCA(inf_filtered_feature_bc_matrix, features = c(s.genes, g2m.genes))
DimPlot(inf_filtered_feature_bc_matrix)

inf_filtered_feature_bc_matrix <- ScaleData(inf_filtered_feature_bc_matrix, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(inf_filtered_feature_bc_matrix))


inf_filtered_feature_bc_matrix <- FindNeighbors(inf_filtered_feature_bc_matrix, dims = 1:20)
inf_filtered_feature_bc_matrix <- FindClusters(inf_filtered_feature_bc_matrix, resolution = 0.4)



inf_filtered_feature_bc_matrix <- RunUMAP(inf_filtered_feature_bc_matrix, dims = 1:20)

DimPlot(inf_filtered_feature_bc_matrix, reduction = "umap")



# find markers for every cluster compared to all remaining cells, report only the positive
# ones
inf_filtered_feature_bc_matrix <- FindAllMarkers(inf_filtered_feature_bc_matrix, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
inf_filtered_feature_bc_matrix %>%
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

#BASAL
FeaturePlot(inf_filtered_feature_bc_matrix, features = "KRT5") #basal1
FeaturePlot(inf_filtered_feature_bc_matrix, features = "TP63") #basal2
FeaturePlot(inf_filtered_feature_bc_matrix, features = "FN1") #basal3
FeaturePlot(inf_filtered_feature_bc_matrix, features = "DLK2") #basal4
FeaturePlot(inf_filtered_feature_bc_matrix, features = "KRT6A") #basal5

#SUPRABASAL
FeaturePlot(inf_filtered_feature_bc_matrix, features = "KRT4") #suprabasal1
FeaturePlot(inf_filtered_feature_bc_matrix, features = "KRT13") #suprabasal2
FeaturePlot(inf_filtered_feature_bc_matrix, features = "KRT16") #suprabasal3
FeaturePlot(inf_filtered_feature_bc_matrix, features = "KRT23") #suprabasal4

#TUFT
FeaturePlot(inf_filtered_feature_bc_matrix, features = "POUF2f3") #tuft1
FeaturePlot(inf_filtered_feature_bc_matrix, features = "TRPM5") #tuft2
FeaturePlot(inf_filtered_feature_bc_matrix, features = "LRMP") #tuft3

#IONOCYTES
FeaturePlot(inf_filtered_feature_bc_matrix, features = "FOXI1") #ionocytes1
FeaturePlot(inf_filtered_feature_bc_matrix, features = "CFTR") #ionocytes2

##PNEC
FeaturePlot(inf_filtered_feature_bc_matrix, features = c("PCSK1N")) #PNEC1
FeaturePlot(inf_filtered_feature_bc_matrix, features = c("SCGN")) #PNEC2
FeaturePlot(inf_filtered_feature_bc_matrix, features = c("NEB")) #PNEC3
FeaturePlot(inf_filtered_feature_bc_matrix, features = c("HOXB1")) #PNEC4
FeaturePlot(inf_filtered_feature_bc_matrix, features = c("ASCL1")) #PNEC5
FeaturePlot(inf_filtered_feature_bc_matrix, features = c("ASCL2")) #PNEC6
FeaturePlot(inf_filtered_feature_bc_matrix, features = c("FOXA2")) #PNEC7
FeaturePlot(inf_filtered_feature_bc_matrix, features = c("PSMD5")) #PNEC8
FeaturePlot(inf_filtered_feature_bc_matrix, features = c("NGF")) #PNEC9

#CLUB
FeaturePlot(inf_filtered_feature_bc_matrix, features = c("PCSK1N")) #PNEC1
FeaturePlot(inf_filtered_feature_bc_matrix, features = c("SCGN")) #PNEC2
FeaturePlot(inf_filtered_feature_bc_matrix, features = c("NEB")) #PNEC3

#GOBLET/SECRETORY
FeaturePlot(inf_filtered_feature_bc_matrix, features = "MUC5B") #secretory1
FeaturePlot(inf_filtered_feature_bc_matrix, features = "SCGB3A1") #secretory2
FeaturePlot(inf_filtered_feature_bc_matrix, features = "MUC5AC") #secretory3 (goblet)
FeaturePlot(inf_filtered_feature_bc_matrix, features = "FOXQ1") #goblet1
FeaturePlot(inf_filtered_feature_bc_matrix, features = "SPDEF") #goblet2

#DEUTEROSOMAL
FeaturePlot(inf_filtered_feature_bc_matrix, features = "PLK4") #ciliated1
FeaturePlot(inf_filtered_feature_bc_matrix, features = "CCNO") #ciliated2
FeaturePlot(inf_filtered_feature_bc_matrix, features = "CEP78") #ciliated3
FeaturePlot(inf_filtered_feature_bc_matrix, features = "DEUP1") #ciliated1


#CILIATED
FeaturePlot(inf_filtered_feature_bc_matrix, features = "FOXJ1") #ciliated1
FeaturePlot(inf_filtered_feature_bc_matrix, features = "PIFO") #ciliated2
FeaturePlot(inf_filtered_feature_bc_matrix, features = "TPPP3") #ciliated3
FeaturePlot(inf_filtered_feature_bc_matrix, features = "SPEF2") #ciliated1
FeaturePlot(inf_filtered_feature_bc_matrix, features = "DNAH5") #ciliated2
FeaturePlot(inf_filtered_feature_bc_matrix, features = "LRRC6") #ciliated3







#Cluster0
FeaturePlot(inf_filtered_feature_bc_matrix, features = c("orf1ab", "S", "orf3", "orf4a", "orf4b", "orf5", "E", "M", "orf8b")) #viral genes
FeaturePlot(inf_filtered_feature_bc_matrix, features = c("PCNA", "TOP2A", "MCM6", "MKI67")) #cell cycle state genes
FeaturePlot(inf_filtered_feature_bc_matrix, features = c("MUC5B", "MUC5AC")) #secretory
FeaturePlot(inf_filtered_feature_bc_matrix, features = "MUC5B") #secretory1
FeaturePlot(inf_filtered_feature_bc_matrix, features = "SCGB3A1") #secretory2
FeaturePlot(inf_filtered_feature_bc_matrix, features = "MUC5AC") #secretory3 (goblet)
#Cluster1
FeaturePlot(inf_filtered_feature_bc_matrix, features = c("FOXJ1", "TPPP3", "SNTN", "TUBA1A", "CETN2")) #ciliated
FeaturePlot(inf_filtered_feature_bc_matrix, features = "FOXJ1") #ciliated1
FeaturePlot(MERS.combined, features = "PIFO") #ciliated2
FeaturePlot(MERS.combined, features = "TPPP3") #ciliated3

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
