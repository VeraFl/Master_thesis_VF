# Analysis of Single-cell Seq data
# 22.06.2022


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

###################################################################################################################################

# Import data which was downloaded from UBELIX after cellranger count run
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/")
# Create each individual Seurat object (if multiple samples)
for (file in c("mock_filtered_feature_bc_matrix", "inf_filtered_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("data/C. Dromedarius/camel229E_exon/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}


# Check the metadata in the new Seurat objects
head(inf_filtered_feature_bc_matrix@meta.data)
head(mock_filtered_feature_bc_matrix@meta.data)

#####################################################################################
# QC
#Calculating some QC statistics
#goes through every row and adds all the counts of the columns (molecules per cell)
counts_per_cell_inf <- Matrix::colSums(inf_filtered_feature_bc_matrix@assays$RNA@counts)
counts_per_cell_mock <- Matrix::colSums(mock_filtered_feature_bc_matrix@assays$RNA@counts)
#goes through every column and adds all the counts of the rows (molecules per gene)
counts_per_gene_inf <- Matrix::rowSums(inf_filtered_feature_bc_matrix@assays$RNA@counts)
counts_per_gene_mock <- Matrix::rowSums(mock_filtered_feature_bc_matrix@assays$RNA@counts)
#count gene itself only if it has non-zero reads mapped (not value).
genes_per_cell_inf <- Matrix::colSums(inf_filtered_feature_bc_matrix@assays$RNA@counts>0) # only count gene if it has non-zero reads mapped (not value).
cells_per_gene_inf <- Matrix::rowSums(inf_filtered_feature_bc_matrix@assays$RNA@counts>0) # only count cells where the gene is expressed
genes_per_cell_mock <- Matrix::colSums(mock_filtered_feature_bc_matrix@assays$RNA@counts>0) # only count gene if it has non-zero reads mapped (not value).
cells_per_gene_mock <- Matrix::rowSums(mock_filtered_feature_bc_matrix@assays$RNA@counts>0) # only count cells where the gene is expressed

#plotting how many genes per cell we have in the samples
plot(sort(genes_per_cell_inf), ylim = c(1, 10000), xlim = c(0, 4500), xlab='cell', log='y', main='genes per cell (ordered) - camel229E infected')
plot(sort(genes_per_cell_mock), ylim = c(1, 10000), xlim = c(0, 4500), xlab='cell', log='y', main='genes per cell (ordered) - camel229E mock')


#####################################################################################
#show the number of rows of the matrix to locate the viral genes
nrow(mock_filtered_feature_bc_matrix@assays$RNA@counts)
#checking for the viral gene "camel229E" at the end of the matrix
tail(counts_per_gene_inf, 15)
tail(counts_per_gene_mock, 15)
#show the specific location in the matrix for control
inf_filtered_feature_bc_matrix@assays$RNA@counts[29364:29375,]

#take the right rows in the matrix
Infected_viral_counts <- as.matrix(inf_filtered_feature_bc_matrix@assays$RNA@counts[29364:29375,])
Mock_viral_counts <- as.matrix(mock_filtered_feature_bc_matrix@assays$RNA@counts[29364:29375,])
#Count the viral counts
infected.total.viral.counts <- Matrix::colSums(Infected_viral_counts)
mock.total.viral.counts <- Matrix::colSums(Mock_viral_counts)
#Count the total gene counts
infected.total.counts <- Matrix::colSums(inf_filtered_feature_bc_matrix@assays$RNA@counts)
mock.total.counts <- Matrix::colSums(mock_filtered_feature_bc_matrix@assays$RNA@counts)

infected.percent.viral <- infected.total.viral.counts/infected.total.counts
mock.percent.viral <- mock.total.viral.counts/mock.total.counts


#Adding the calculated variables to the metadata of the corresponding sample
inf_filtered_feature_bc_matrix <- AddMetaData(object = inf_filtered_feature_bc_matrix, metadata = infected.total.viral.counts, col.name = "Total.Viral.Counts")
inf_filtered_feature_bc_matrix <- AddMetaData(object = inf_filtered_feature_bc_matrix, metadata = infected.percent.viral, col.name = "Percent.Viral")
inf_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(inf_filtered_feature_bc_matrix, pattern = "^MT-")
inf_filtered_feature_bc_matrix <- AddMetaData(inf_filtered_feature_bc_matrix, "INFECTED", col.name = "Treat")

head(inf_filtered_feature_bc_matrix@meta.data)

mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.total.viral.counts, col.name = "Total.Viral.Counts")
mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.percent.viral, col.name = "Percent.Viral")
mock_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(mock_filtered_feature_bc_matrix, pattern = "^MT-")
mock_filtered_feature_bc_matrix <- AddMetaData(mock_filtered_feature_bc_matrix, "MOCK", col.name = "Treat")
head(mock_filtered_feature_bc_matrix@meta.data)

#Extracting the Viral-Percentage column for every cell and transpose it to plot it below
Infected_viral_percentage <- inf_filtered_feature_bc_matrix@meta.data %>% select(5) %>% t()
Mock_viral_percentage <- mock_filtered_feature_bc_matrix@meta.data %>% select(5) %>% t()

#plotting how many virus positive cells there are in each treatment
plot(sort(Infected_viral_percentage), xlim = c(0, 4500), xlab='cell', main='Viral.Percentage per cell (ordered) - camel229E infected')
plot(sort(Mock_viral_percentage), xlim = c(0, 4500), xlab='cell', main='Viral.Percentage per cell (ordered) - camel229E mock')

#Threshold above 0.0002 Percent Viral is infected

inf_filtered_feature_bc_matrix@meta.data <- inf_filtered_feature_bc_matrix@meta.data %>% mutate(Status = case_when(Percent.Viral > 0.0002 ~ "Infected",
                                                                                                                   Percent.Viral <= 0.0002 ~ "Bystander"))
mock_filtered_feature_bc_matrix <- AddMetaData(mock_filtered_feature_bc_matrix, "Non-infected", col.name = "Status")


##################################################################################################################
#QC
# Make violin plot to show number of features, number of counts and mitochondrial content
#infected
VlnPlot(inf_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(inf_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), ncol = 2)
VlnPlot(inf_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), ncol = 2)
VlnPlot(inf_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), y.max = 0.0005, ncol = 2)
#mock
VlnPlot(mock_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(mock_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), ncol = 2)
VlnPlot(mock_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), ncol = 2)
VlnPlot(mock_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), y.max = 0.0005, ncol = 2)


#FeatureScatter is typically used to visualize feature-feature relationships, but can be used
#for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
#infected
plot1 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "Percent.Viral", feature2 = "percent.mt")
plot1 + plot2 + plot3
#mock
plot1 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# We filter cells that have unique feature counts over
# We filter cells that have >5% mitochondrial counts
#infected
inf_filtered_feature_bc_matrix <- subset(inf_filtered_feature_bc_matrix, subset = nFeature_RNA > 1000 & percent.mt < 30)
#mock
mock_filtered_feature_bc_matrix <- subset(mock_filtered_feature_bc_matrix, subset = nFeature_RNA > 1000 & percent.mt < 30)

###################################################################################################################
# #cell cycle state genes from the Seurat package
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# 
# inf_filtered_feature_bc_matrix <- CellCycleScoring(inf_filtered_feature_bc_matrix, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# mock_filtered_feature_bc_matrix <- CellCycleScoring(mock_filtered_feature_bc_matrix, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#################################################################################################################################
#Integrate the two sample matrices to one seurat object
camel229E.list <- c(mock_filtered_feature_bc_matrix, inf_filtered_feature_bc_matrix)
# normalize and identify variable features for each dataset independently
camel229E.list <- lapply(X = camel229E.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = camel229E.list)


camel229E.anchors <- FindIntegrationAnchors(object.list = camel229E.list, anchor.features = features)

# this command creates an 'integrated' data assay
camel229E.combined <- IntegrateData(anchorset = camel229E.anchors)


###Perform an integrated analysis
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(camel229E.combined) <- "integrated"

#Saving the seurat object to disk to access in other script for analysis
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.dromedarius")
SaveH5Seurat(camel229E.combined, "camel229E.combined", overwrite = TRUE)


