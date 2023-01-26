# Analysis of Single-cell Seq data
# 04.08.2022


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
library(cowplot)

#Treatment colors(camel-229E, Mock)
cols_treat = c("#F39C6B", "#96ADC8")
show_col(cols_treat)
#Status colors
cols_stat = c("#CD4631","#DEA47E")
show_col(cols_stat)
#Celltype colors
cols_cell = hcl.colors(7, "Temps")
show_col(cols_cell)


###################################################################################################################################

# Import data which was downloaded from UBELIX after cellranger count run
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/")
# Create each individual Seurat object (if multiple samples)
for (file in c("mock_filtered_feature_bc_matrix", "inf_filtered_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("data/C. Ferus/camel229E_exon_mito/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}


# Check the metadata in the new Seurat objects
head(inf_filtered_feature_bc_matrix@meta.data)
head(mock_filtered_feature_bc_matrix@meta.data)

#####################################################################################
#QC
#when there are more than one (all the camel229E genes separately):
#calculating some info about the expression of the camel229E genes and gene content in general
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
plot(sort(genes_per_cell_inf), ylim = c(1, 10000), xlim = c(0, 5000), xlab='cell', log='y', main='genes per cell (ordered) - camel229E infected')
plot(sort(genes_per_cell_mock), ylim = c(1, 10000), xlim = c(0, 5000), xlab='cell', log='y', main='genes per cell (ordered) - camel229E mock')

#####################################################################################
# Adding important info to the meta data column
#show the number of rows of the matrix to locate the camel229E genes
nrow(mock_filtered_feature_bc_matrix@assays$RNA@counts)
nrow(inf_filtered_feature_bc_matrix@assays$RNA@counts)
tail(counts_per_gene_inf, 15)
tail(counts_per_gene_mock, 15)
#show the specific location in the matrix for control
inf_filtered_feature_bc_matrix@assays$RNA@counts[30848:30857,]

#take the right rows in the matrix
Infected_camel229E_counts <- as.matrix(inf_filtered_feature_bc_matrix@assays$RNA@counts[30848:30857,])
Mock_camel229E_counts <- as.matrix(mock_filtered_feature_bc_matrix@assays$RNA@counts[30848:30857,])
#Count the camel229E counts
infected.total.camel229E.counts <- Matrix::colSums(Infected_camel229E_counts)
mock.total.camel229E.counts <- Matrix::colSums(Mock_camel229E_counts)
#Count the total gene counts
infected.total.counts <- Matrix::colSums(inf_filtered_feature_bc_matrix@assays$RNA@counts)
mock.total.counts <- Matrix::colSums(mock_filtered_feature_bc_matrix@assays$RNA@counts)

infected.percent.camel229E <- infected.total.camel229E.counts/infected.total.counts
mock.percent.camel229E <- mock.total.camel229E.counts/mock.total.counts


#Adding the calculated variables to the metadata of the corresponding sample
inf_filtered_feature_bc_matrix <- AddMetaData(object = inf_filtered_feature_bc_matrix, metadata = infected.total.camel229E.counts, col.name = "Total.camel229E.Counts")
inf_filtered_feature_bc_matrix <- AddMetaData(object = inf_filtered_feature_bc_matrix, metadata = infected.percent.camel229E, col.name = "Percent.camel229E")
inf_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(inf_filtered_feature_bc_matrix, pattern = "^KEG")
inf_filtered_feature_bc_matrix <- AddMetaData(inf_filtered_feature_bc_matrix, "dcCoV-ACN4", col.name = "Treat")

head(inf_filtered_feature_bc_matrix@meta.data)

mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.total.camel229E.counts, col.name = "Total.camel229E.Counts")
mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.percent.camel229E, col.name = "Percent.camel229E")
mock_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(mock_filtered_feature_bc_matrix, pattern = "^KEG")
mock_filtered_feature_bc_matrix <- AddMetaData(mock_filtered_feature_bc_matrix, "Mock", col.name = "Treat")
head(mock_filtered_feature_bc_matrix@meta.data)

#Adding a Meta data column with information about the Status of each cell. If the cell is Infected or Bystander
#in the camel229E infected Treatment group or non-infected in the Mock Group

#Extracting the camel229E-Percentage column for every cell and transpose it to plot it below
Infected_camel229E_percentage <- inf_filtered_feature_bc_matrix@meta.data %>% select(5) %>% t()
Mock_camel229E_percentage <- mock_filtered_feature_bc_matrix@meta.data %>% select(5) %>% t()

#plotting how many virus positive cells there are in each treatment
plot(sort(Infected_camel229E_percentage), ylim = c(0,0.001), xlim = c(0, 4500), xlab='cell', main='camel229E.Percentage vs cell (ordered) - camel229E infected')
plot(sort(Mock_camel229E_percentage), xlim = c(0, 4500), xlab='cell', main='camel229E.Percentage vs cell (ordered) - camel229E Mock')

#############################################################################################################################
#Calculation and visualization of infection threshold for the Human samples
#infected camel229E for camel229E
camel229E_infected_percentage <- inf_filtered_feature_bc_matrix@meta.data %>% select(5)

#calculate minimal number after 0 to set instead of 0, so they are included into plot
halfmin_camel229E <- min(camel229E_infected_percentage$Percent.camel229E[camel229E_infected_percentage$Percent.camel229E>0])
halfmin_camel229E

camel229E_infected_percentage$Percent.camel229E_log10 <- log10(camel229E_infected_percentage$Percent.camel229E+halfmin_camel229E)
range(camel229E_infected_percentage$Percent.camel229E_log10)

des.all <- density(camel229E_infected_percentage$Percent.camel229E_log10)
min.all <- des.all$x[which(diff(sign(diff(des.all$y)))==2)+1]
min.all
#-3.172746

#mock for camel229E
mock_percentage <- mock_filtered_feature_bc_matrix@meta.data %>% select(5)
mock_percentage$Percent.camel229E_log10 <- log10(mock_percentage$Percent.camel229E+halfmin_camel229E)
range(mock_percentage$Percent.camel229E_log10)


#set colors for infected and non-infected
camel229E_infected_percentage$group = ifelse(camel229E_infected_percentage$Percent.camel229E_log10 <= min.all[2], 'Non-infected', ifelse(camel229E_infected_percentage$Percent.camel229E_log10, 'Infected'))
mock_percentage$group = ifelse(mock_percentage$Percent.camel229E_log10 <= min.all[2], 'Non-infected', ifelse(mock_percentage$Percent.camel229E_log10, 'Infected'))

p1 <- camel229E_infected_percentage %>%
  ggplot(aes(x = Percent.camel229E_log10, fill = group))+
  geom_histogram(bins = 500)+
  scale_fill_manual(values = cols_stat) +
  xlim(-5,1)+
  coord_cartesian(ylim=c(0, 200))+
  ylab("Cell Number")+
  xlab("log10(% of mRNA of camel229E-CoV)")+
  ggtitle("camel229E-CoV Infected")+
  theme_classic()+
  geom_vline(xintercept = min.all[2], linetype = "dashed", size = 0.9)+
  theme(axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15)))


p2 <- mock_percentage %>%
  ggplot(aes(x = Percent.camel229E_log10, fill = group)) +
  geom_histogram(bins = 500) +
  scale_fill_manual(values = cols_stat) +
  xlim(-5,1) +
  coord_cartesian(ylim=c(0, 200))+
  ylab("Cell Number")+
  xlab("")+
  ggtitle("Mock")+
  theme_classic() +
  geom_vline(xintercept = min.all[2], linetype = "dashed", size = 0.9)+
  theme(axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15)))+
  NoLegend()

png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_camel229E/threshold_camel229E_infection.png",
    width=500, height=600)
p2 / p1
dev.off()

#the infection threshold calculated back from log10 is
10^(-3.172746)
#camel229E Infection threshold: 0.0006718217


#Threshold above 0.0002 Percent camel229E is infected

inf_filtered_feature_bc_matrix@meta.data <- inf_filtered_feature_bc_matrix@meta.data %>% mutate(Status = case_when(Percent.camel229E > 0.0006718217 ~ "Infected",
                                                                                                                   Percent.camel229E <= 0.0006718217 ~ "Non-infected"))

mock_filtered_feature_bc_matrix <- AddMetaData(mock_filtered_feature_bc_matrix, "Non-infected", col.name = "Status")

##################################################################################################################
# QC

# # Make violin plot to show number of features, number of counts and mitochondrial content
# #infected
# VlnPlot(inf_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# VlnPlot(inf_filtered_feature_bc_matrix, features = c("Percent.camel229E", "Total.camel229E.Counts"), ncol = 2)
# VlnPlot(inf_filtered_feature_bc_matrix, features = c("Percent.camel229E", "Total.camel229E.Counts"), ncol = 2)
# VlnPlot(inf_filtered_feature_bc_matrix, features = c("Percent.camel229E", "Total.camel229E.Counts"), y.max = 0.8, ncol = 2)
# #mock
# VlnPlot(mock_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# VlnPlot(mock_filtered_feature_bc_matrix, features = c("Percent.camel229E", "Total.camel229E.Counts"), ncol = 2)
# VlnPlot(mock_filtered_feature_bc_matrix, features = c("Percent.camel229E", "Total.camel229E.Counts"), ncol = 2)
# VlnPlot(mock_filtered_feature_bc_matrix, features = c("Percent.camel229E", "Total.camel229E.Counts"), y.max = 0.8, ncol = 2)
# 
# 
# #FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# #for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# #infected
# plot1 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot3 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "Total.camel229E.Counts")
# plot3
# #mock
# plot1 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

# Add number of genes per UMI for each cell to metadata
inf_filtered_feature_bc_matrix$log10GenesPerUMI <- log10(inf_filtered_feature_bc_matrix$nFeature_RNA) / log10(inf_filtered_feature_bc_matrix$nCount_RNA)
mock_filtered_feature_bc_matrix$log10GenesPerUMI <- log10(mock_filtered_feature_bc_matrix$nFeature_RNA) / log10(mock_filtered_feature_bc_matrix$nCount_RNA)


metadata <- bind_rows(inf_filtered_feature_bc_matrix@meta.data, mock_filtered_feature_bc_matrix@meta.data)
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)



metadata_count <- metadata %>%
  group_by(Status, Treat) %>%
  count(Status)


#Plots for QC
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_camel229E/cell_numbers_per_virus_before_filtering.png",
    width=600, height=500)
metadata_count %>% 
  ggplot(aes(x=Treat, y=n, fill=Status)) + 
  geom_col() +
  theme_classic() +
  scale_fill_manual(values = cols_stat)+
  xlab("Condition")+
  ylab("Number of cells")+
  ylim(0,5000)+
  coord_flip()+
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold"))
dev.off()


#Number of UMI (transcripts) per cell
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_camel229E/UMI_per_cell_density_log_before_filtering.png",
    width=600, height=500)
metadata %>% 
  ggplot(aes(color=Treat, x=nUMI, fill= Treat)) + 
  geom_density(alpha = 0.6) + 
  scale_color_manual(values= cols_treat) +
  scale_fill_manual(values = cols_treat) +
  scale_x_log10(labels = comma_format(big.mark = "'",
                                      decimal.mark = ",")) + 
  theme_classic() +
  ylab("Cell density") +
  xlab("Number of intracellular mRNA transcripts")+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold"))
dev.off()

#Number of Genes per cell
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_camel229E/genes_per_cell_density_log_before_filtering.png",
    width=600, height=500)
metadata %>% 
  ggplot(aes(color=Treat, x=nGene, fill= Treat)) + 
  geom_density(alpha = 0.6) +
  scale_color_manual(values= cols_treat) +
  scale_fill_manual(values = cols_treat) +
  theme_classic() +
  scale_x_log10(labels = comma_format(big.mark = "'",
                                      decimal.mark = ",")) + 
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold"))+
  ylab("Cell density")+
  xlab("Number of detected genes")
dev.off()

#Complexity
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_camel229E/complexity_log_before_filtering.png",
    width=700, height=500)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = Treat, fill=Treat)) +
  geom_density(alpha = 0.6) +
  scale_color_manual(values= cols_treat) +
  scale_fill_manual(values = cols_treat) +
  xlim(0.7,1)+
  theme_classic() +
  ylab("Cell density") +
  xlab("Complexity") +
  geom_vline(xintercept = 0.8, linetype = "dashed", size = 0.9) +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold"))
dev.off()


# Visualize the distribution of mitochondrial gene expression detected per cell
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_camel229E/Mitochondrial_content_density_before_filtering.png",
    width=700, height=500)
metadata %>% 
  ggplot(aes(color=Treat, x=percent.mt, fill=Treat)) + 
  geom_density(alpha = 0.6) + 
  scale_color_manual(values= cols_treat) +
  scale_fill_manual(values = cols_treat) +
  xlim(0,50)+
  theme_classic() +
  geom_vline(xintercept = 30)
dev.off()


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_camel229E/correlation_gene_UMI_log_before_filtering.png",
    width=1000, height=500)
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color = Treat)) + 
  geom_point() + 
  stat_smooth(method=lm, color = "black") +
  scale_color_manual(values= cols_treat) +
  xlab("Number of transcripts per cell") +
  ylab("Number of genes per cell")+
  scale_x_log10(labels = comma_format(big.mark = "'",
                                      decimal.mark = ",")) + 
  scale_y_log10(labels = comma_format(big.mark = "'",
                                      decimal.mark = ",")) + 
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold"))+
  facet_wrap(~Treat)
dev.off()

# We filter cells that have unique feature counts over
# We filter cells that have >5% mitochondrial counts
#infected
inf_filtered_feature_bc_matrix <- subset(inf_filtered_feature_bc_matrix, subset = nFeature_RNA > 1000 & percent.mt < 30)
#mock
mock_filtered_feature_bc_matrix <- subset(mock_filtered_feature_bc_matrix, subset = nFeature_RNA > 1000 & percent.mt < 30)


###################################################################################################################
# cell cycle state genes from the Seurat package
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes
g2m.genes <- c(g2m.genes, "PIMREG", "JPT1")
s.genes <- c(s.genes, "CENPU")
memory.limit(24000)
###################################################################################################################
#Integrate the two sample matrices to one seurat object
camel229E.list <- c(mock_filtered_feature_bc_matrix, inf_filtered_feature_bc_matrix)
# normalize and identify variable features for each dataset independently
camel229E.list <- lapply(X = camel229E.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  x <- SCTransform(x, vars.to.regress = c("Percent.camel229E", "S.Score", "G2M.Score", "percent.mt"), verbose = FALSE)
})


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = camel229E.list, nfeatures = 3000)
camel229E.list <- PrepSCTIntegration(object.list = camel229E.list, anchor.features = features)
camel229E.anchors <- FindIntegrationAnchors(object.list = camel229E.list, anchor.features = features, normalization.method = "SCT")

# this command creates an 'integrated' data assay
camel229E.combined <- IntegrateData(anchorset = camel229E.anchors, normalization.method = "SCT")


###Perform an integrated analysis
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(camel229E.combined) <- "integrated"

###################################################################################################################

#Saving the seurat object to disk to access in other script for analysis
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.ferus")
SaveH5Seurat(camel229E.combined, "camel229E.combined", overwrite = TRUE)

#########################################################################################################
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.ferus")
camel229E.combined <- LoadH5Seurat("camel229E.combined.h5seurat")
camel229E.combined


# Add number of genes per UMI for each cell to metadata
camel229E.combined$log10GenesPerUMI <- log10(camel229E.combined$nFeature_RNA) / log10(camel229E.combined$nCount_RNA)


metadata <- camel229E.combined@meta.data
# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

metadata_count <- metadata %>%
  group_by(Status, Treat) %>%
  count(Status)


#Plots for QC
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_camel229E/cell_numbers_per_virus_after_filtering.png",
    width=600, height=500)
metadata_count %>% 
  ggplot(aes(x=Treat, y=n, fill=Status)) + 
  geom_col() +
  theme_classic() +
  scale_fill_manual(values = cols_stat)+
  xlab("Condition")+
  ylab("Number of cells")+
  ylim(0,5000)+
  coord_flip()+
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), col = "black", size = 7, hjust = 0.5, check_overlap=TRUE)+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(r = 20)),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 20, face = "bold"))
dev.off()


#Number of UMI (transcripts) per cell
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_camel229E/UMI_per_cell_density_log_after_filtering.png",
    width=600, height=500)
metadata %>% 
  ggplot(aes(color=Treat, x=nUMI, fill= Treat)) + 
  geom_density(alpha = 0.6) + 
  scale_color_manual(values= cols_treat) +
  scale_fill_manual(values = cols_treat) +
  scale_x_log10(labels = comma_format(big.mark = "'",
                                      decimal.mark = ",")) + 
  theme_classic() +
  ylab("Cell density") +
  xlab("Number of intracellular mRNA transcripts")+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(r = 20)),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 20, face = "bold"))
dev.off()

#Number of Genes per cell
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_camel229E/genes_per_cell_density_log_after_filtering.png",
    width=600, height=500)
metadata %>% 
  ggplot(aes(color=Treat, x=nGene, fill= Treat)) + 
  geom_density(alpha = 0.6) +
  scale_color_manual(values= cols_treat) +
  scale_fill_manual(values = cols_treat) +
  theme_classic() +
  scale_x_log10(labels = comma_format(big.mark = "'",
                                      decimal.mark = ",")) + 
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(r = 20)),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 20, face = "bold"))+
  ylab("Cell density")+
  xlab("Number of detected genes")
dev.off()

#Complexity
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_camel229E/complexity_log_after_filtering.png",
    width=700, height=500)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = Treat, fill=Treat)) +
  geom_density(alpha = 0.6) +
  scale_color_manual(values= cols_treat) +
  scale_fill_manual(values = cols_treat) +
  xlim(0.7,1)+
  theme_classic() +
  ylab("Cell density") +
  xlab("Complexity (log10 Genes per UMI)") +
  geom_vline(xintercept = 0.8, linetype = "dashed", size = 0.9) +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(r = 20)),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 20, face = "bold"))
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per cell
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_camel229E/Mitochondrial_content_density_after_filtering.png",
    width=700, height=500)
metadata %>% 
  ggplot(aes(color=Treat, x=percent.mt, fill=Treat)) + 
  geom_density(alpha = 0.6) + 
  scale_color_manual(values= cols_treat) +
  scale_fill_manual(values = cols_treat) +
  xlim(0,50)+
  theme_classic() +
  ylab("Cell density") +
  xlab("Percentage of mitochondrial transcripts") +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(r = 20)),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 20, face = "bold"))+
  geom_vline(xintercept = 30, linetype = "dashed", colour = "red")
dev.off()


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_camel229E/correlation_gene_UMI_log_after_filtering.png",
    width=1000, height=500)
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color = Treat)) + 
  geom_point() + 
  stat_smooth(method=lm, color = "black") +
  scale_color_manual(values= cols_treat) +
  xlab("Number of transcripts per cell") +
  ylab("Number of genes per cell")+
  scale_x_log10(labels = comma_format(big.mark = "'",
                                      decimal.mark = ",")) + 
  scale_y_log10(labels = comma_format(big.mark = "'",
                                      decimal.mark = ",")) + 
  theme_classic() +
  facet_wrap(~Treat)+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(r = 20)),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 20, face = "bold"),
        strip.text = element_text(size = 17))
dev.off()

