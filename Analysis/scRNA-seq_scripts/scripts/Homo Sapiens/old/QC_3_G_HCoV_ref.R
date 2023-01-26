# Analysis of Single-cell Seq data
# 06.05.2022
# Count data of Infected Samples (3_G) with concatanated custom Reference Genome "Homo_sapiens.GRCh38.105.HCoV_229E.gtf"
# Count matrix produced by "Cell_Ranger_Count_3_G.sh"

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


for (file in c("inf_filtered_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("data/H. sapiens/HCoV/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

#Importing the infected samples independently
inf_raw_data <- Read10X(data.dir = "data/H. sapiens/HCoV/inf_raw_feature_bc_matrix")
inf_filtered_data <- Read10X(data.dir = "data/H. sapiens/HCoV/inf_filtered_feature_bc_matrix")

inf_raw <- CreateSeuratObject(counts = inf_raw_data, project = "inf", min.cells = 3, min.features = 200)
inf_filtered <- CreateSeuratObject(counts = inf_filtered_data, project = "inf", min.cells = 3, min.features = 200)


#calculating counts per gene and checking the viral genes if we have counts for them.
#looking at raw and filtered matrices to see a difference
raw_counts_per_gene <- Matrix::rowSums(inf_raw_data)
filtered_counts_per_gene <- Matrix::rowSums(inf_filtered_data)
virus_genes_raw <- tail(raw_counts_per_gene, 11) %>% data.frame()
virus_genes_filtered <- tail(filtered_counts_per_gene, 11) %>% data.frame()


#Cell cycle marker genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Giving a cell cycle score to every cell
HCoV_inf <- CellCycleScoring(HCoV_inf, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(HCoV_inf[[]])

HCoV_inf <- HCoV_inf

#Making the seurat object to a tidyseurat
# HCoV_inf <- tidy(HCoV_inf)
# HCoV_inf

###################################################################################################################################
# HCoV_mock %>%
#   tidyseurat::ggplot(aes(nFeature_RNA, fill = Phase)) +
#   geom_histogram() +
#   my_theme

###################################################################################################################################
# QC

#Calculating some QC statistics
#goes through every row and adds all the counts of the columns (molecules per cell)
counts_per_cell <- Matrix::colSums(seurat_data) 
#goes through every column and adds all the counts of the rows (molecules per gene)
counts_per_gene <- Matrix::rowSums(seurat_data)
tail(counts_per_gene, 20)
#count gene itself only if it has non-zero reads mapped (not value).
genes_per_cell <- Matrix::colSums(seurat_data>0) # only count gene if it has non-zero reads mapped (not value).
cells_per_gene <- Matrix::rowSums(seurat_data>0) # only count cells where the gene is expressed

#calculating counts per gene and checking the viral genes if we have counts for them.
#looking at raw and filteres matrices to see a difference

tail(counts_per_gene, 20)


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
# 

#Genes per cell ordered
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered) - HCoV - INF')


#################################################################################################################
# Add column for Mitochondrial Content in percentage to the meta data
inf_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(inf_filtered_feature_bc_matrix, pattern = "^MT-")
head(inf_filtered_feature_bc_matrix@meta.data)


# Make violin plot to show number of features, number of counts and mitochondiral content
VlnPlot(inf_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

?PercentageFeatureSet
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

?FeatureScatter

# We filter cells that have unique feature counts over
# We filter cells that have >5% mitochondrial counts

inf_filtered_feature_bc_matrix <- subset(inf_filtered_feature_bc_matrix, subset = nFeature_RNA > 1000 & percent.mt < 30)

# Make violin plot to show number of features, number of counts and mitochondiral content after filtering
VlnPlot(inf_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


#Normalization

HCoV_inf <- NormalizeData(HCoV_inf, normalization.method = "LogNormalize", scale.factor = 10000)


# Identification of highly variable features (feature selection)

HCoV_inf <- FindVariableFeatures(HCoV_inf, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(HCoV_inf), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(HCoV_inf)
plot2 <- LabelPoints(plot = plot1, points = top10)
plot1 + plot2

#Scaling the data
all.genes <- rownames(HCoV_inf)
HCoV_inf <- ScaleData(HCoV_inf, features = all.genes)


###################################################################################################################################
# Dimensionality Reduction - PCA Run
HCoV_inf <- RunPCA(HCoV_inf, features = VariableFeatures(object = HCoV_inf))

# Examine and visualize PCA results a few different ways
print(HCoV_inf[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(HCoV_inf, dims = 1:2, reduction = "pca")
DimPlot(HCoV_inf, reduction = "pca")

#Heatmap
DimHeatmap(HCoV_inf, dims = 1, cells = 500, balanced = TRUE)


###################################################################################################################################
#Regressing out the cell cycle markers. The cell cycle markers we defined at the beginning are not plotted and regressed out.

# Visualize the distribution of cell cycle markers across
RidgePlot(HCoV_inf, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2) # no expression of these markers

#Running a PCA on cell cycle genes reveals, unsurprisingly, that cells seperate entirely by phase
HCoV_inf <- RunPCA(HCoV_inf, features = c(s.genes, g2m.genes))
DimPlot(HCoV_inf)

#Regressing out the cell cycle state markers should improve the bias towards the cell cycle markers (not that strong in that sample anyways)
HCoV_inf <- ScaleData(HCoV_inf, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(HCoV_inf))

HCoV_inf <- RunPCA(HCoV_inf, features = c(s.genes, g2m.genes))
DimPlot(HCoV_inf)

ElbowPlot(HCoV_inf)

###################################################################################################################################
# Dimensionality Reduction - UMAP Run

#Finding Neighbors and Clusters
HCoV_inf <- FindNeighbors(HCoV_inf, dims = 1:5)
HCoV_inf <- FindClusters(HCoV_inf, resolution = 0.2)

#Run UMAP
HCoV_inf <- RunUMAP(HCoV_inf, dims = 1:5)
DimPlot(HCoV_inf, reduction = "umap")

?RunUMAP

# #Finding Markers to respective Clusters made above to differentiate between clusters
# cluster2.markers <- FindMarkers(HCoV_inf, ident.1  = 2, min.pct = 0.25)
# head(cluster2.markers, n = 5)
# 
# cluster5.markers <- FindMarkers(Sample_3_A, ident.1 = 2, ident.2 = 0, min.pct = 0.25)
# head(cluster5.markers, n = 5)
# 


# Find all Markers
# find markers for every cluster compared to all remaining cells, report only the positive ones
HCoV_inf.markers <- FindAllMarkers(HCoV_inf, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
HCoV_inf.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#Plotting expression of different markers (making the clusters)

VlnPlot(HCoV_inf, features = c("SCGB3A1", "GPRC5A", "SLPI"))
VlnPlot(HCoV_inf, features = c("SCGB3A1", "GPRC5A"), slot = "counts", log = TRUE)


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
FeaturePlot(HCoV_inf, features = "KRT5") #basal1
FeaturePlot(HCoV_inf, features = "TP63") #basal2
FeaturePlot(HCoV_inf, features = "FN1") #basal3
FeaturePlot(HCoV_inf, features = "DLK2") #basal4
FeaturePlot(HCoV_inf, features = "KRT6A") #basal5

#SUPRABASAL
FeaturePlot(HCoV_inf, features = "KRT4") #suprabasal1
FeaturePlot(HCoV_inf, features = "KRT13") #suprabasal2
FeaturePlot(HCoV_inf, features = "KRT16") #suprabasal3
FeaturePlot(HCoV_inf, features = "KRT23") #suprabasal4

#TUFT
FeaturePlot(HCoV_inf, features = "POUF2f3") #tuft1
FeaturePlot(HCoV_inf, features = "TRPM5") #tuft2
FeaturePlot(HCoV_inf, features = "LRMP") #tuft3

#IONOCYTES
FeaturePlot(HCoV_inf, features = "FOXI1") #ionocytes1
FeaturePlot(HCoV_inf, features = "CFTR") #ionocytes2

##PNEC
FeaturePlot(HCoV_inf, features = c("PCSK1N")) #PNEC1
FeaturePlot(HCoV_inf, features = c("SCGN")) #PNEC2
FeaturePlot(HCoV_inf, features = c("NEB")) #PNEC3
FeaturePlot(HCoV_inf, features = c("HOXB1")) #PNEC4
FeaturePlot(HCoV_inf, features = c("ASCL1")) #PNEC5
FeaturePlot(HCoV_inf, features = c("ASCL2")) #PNEC6
FeaturePlot(HCoV_inf, features = c("FOXA2")) #PNEC7
FeaturePlot(HCoV_inf, features = c("PSMD5")) #PNEC8
FeaturePlot(HCoV_inf, features = c("NGF")) #PNEC9

#CLUB
FeaturePlot(HCoV_inf, features = c("PCSK1N")) #PNEC1
FeaturePlot(HCoV_inf, features = c("SCGN")) #PNEC2
FeaturePlot(HCoV_inf, features = c("NEB")) #PNEC3

#GOBLET/SECRETORY
FeaturePlot(HCoV_inf, features = "MUC5B") #secretory1
FeaturePlot(HCoV_inf, features = "SCGB3A1") #secretory2
FeaturePlot(HCoV_inf, features = "MUC5AC") #secretory3 (goblet)
FeaturePlot(HCoV_inf, features = "FOXQ1") #goblet1
FeaturePlot(HCoV_inf, features = "SPDEF") #goblet2

#DEUTEROSOMAL
FeaturePlot(HCoV_inf, features = "PLK4") #ciliated1
FeaturePlot(HCoV_inf, features = "CCNO") #ciliated2
FeaturePlot(HCoV_inf, features = "CEP78") #ciliated3
FeaturePlot(HCoV_inf, features = "DEUP1") #ciliated1


#CILIATED
FeaturePlot(HCoV_inf, features = "FOXJ1") #ciliated1
FeaturePlot(HCoV_inf, features = "PIFO") #ciliated2
FeaturePlot(HCoV_inf, features = "TPPP3") #ciliated3
FeaturePlot(HCoV_inf, features = "SPEF2") #ciliated1
FeaturePlot(HCoV_inf, features = "DNAH5") #ciliated2
FeaturePlot(HCoV_inf, features = "LRRC6") #ciliated3




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


