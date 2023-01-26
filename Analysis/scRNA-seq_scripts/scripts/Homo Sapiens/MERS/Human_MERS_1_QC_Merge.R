# QC and Merging of Single-cell Seq data
# 13.05.2022


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


#Treatment colors(camel-229E, MERS, Mock)
cols_treat = c("#A5668B", "#96ADC8")
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
  seurat_data <- Read10X(data.dir = paste0("data/H. sapiens/MERS/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.features = 100,
                                   project = file)
  assign(file, seurat_obj)
}


# Check the metadata in the new Seurat objects
head(inf_filtered_feature_bc_matrix@meta.data)
head(mock_filtered_feature_bc_matrix@meta.data)

#####################################################################################
#when there are more than one (all the viral genes seperately):

#calculating some info about the expression of the viral genes and gene content in general
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
plot(sort(genes_per_cell_inf), ylim = c(1, 10000), xlim = c(0, 4500), xlab='cell', log='y', main='genes per cell (ordered) - MERS infected')
plot(sort(genes_per_cell_mock), ylim = c(1, 10000), xlim = c(0, 4500), xlab='cell', log='y', main='genes per cell (ordered) - MERS mock')

hist(colSums(counts(pbmc.sce)), breaks = 100, xlab="#features per cell", main = "Features per cell distribution")

#####################################################################################
# Adding important info to the meta data column
#show the number of rows of the matrix to locate the viral genes
nrow(mock_filtered_feature_bc_matrix@assays$RNA@counts)
nrow(inf_filtered_feature_bc_matrix@assays$RNA@counts)
tail(counts_per_gene_inf, 15)
tail(counts_per_gene_mock, 15)
#show the specific location in the matrix for control
inf_filtered_feature_bc_matrix@assays$RNA@counts[19989:20000,]
inf_filtered_feature_bc_matrix@assays$RNA@counts[19989:20000,]

#take the right rows in the matrix
Infected_viral_counts <- as.matrix(inf_filtered_feature_bc_matrix@assays$RNA@counts[19989:20000,])
Mock_viral_counts <- as.matrix(mock_filtered_feature_bc_matrix@assays$RNA@counts[19989:20000,])
#Count the viral counts
infected.total.viral.counts <- Matrix::colSums(Infected_viral_counts)
mock.total.viral.counts <- Matrix::colSums(Mock_viral_counts)
#Count the total gene counts
infected.total.counts <- Matrix::colSums(inf_filtered_feature_bc_matrix@assays$RNA@counts)
mock.total.counts <- Matrix::colSums(mock_filtered_feature_bc_matrix@assays$RNA@counts)

infected.percent.viral <- infected.total.viral.counts/infected.total.counts
mock.percent.viral <- mock.total.viral.counts/mock.total.counts


#Adding the calculated variables to the metadata of the corresponding sample
inf_filtered_feature_bc_matrix <- AddMetaData(object = inf_filtered_feature_bc_matrix, metadata = infected.total.viral.counts, col.name = "Total.MERS.Counts")
inf_filtered_feature_bc_matrix <- AddMetaData(object = inf_filtered_feature_bc_matrix, metadata = infected.percent.viral, col.name = "Percent.MERS")
inf_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(inf_filtered_feature_bc_matrix, pattern = "^MT-")
inf_filtered_feature_bc_matrix <- AddMetaData(inf_filtered_feature_bc_matrix, "MERS-CoV", col.name = "Treat")

head(inf_filtered_feature_bc_matrix@meta.data)

mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.total.viral.counts, col.name = "Total.MERS.Counts")
mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.percent.viral, col.name = "Percent.MERS")
mock_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(mock_filtered_feature_bc_matrix, pattern = "^MT-")
mock_filtered_feature_bc_matrix <- AddMetaData(mock_filtered_feature_bc_matrix, "Mock", col.name = "Treat")
head(mock_filtered_feature_bc_matrix@meta.data)

#Adding a Meta data column with information about the Status of each cell. If the cell is Infected or Bystander
#in the MERS infected Treatment group or non-infected in the Mock Group

#Extracting the Viral-Percentage column for every cell and transpose it to plot it below
Infected_viral_percentage <- inf_filtered_feature_bc_matrix@meta.data %>% select(5) %>% t()
Mock_viral_percentage <- mock_filtered_feature_bc_matrix@meta.data %>% select(5) %>% t()

#plotting how many virus positive cells there are in each treatment
plot(sort(Infected_viral_percentage), ylim = c(0, 0.0002), xlim = c(0, 4500), xlab='cell', main='Viral.Percentage vs cell (ordered) - MERS infected')
plot(sort(Mock_viral_percentage), xlim = c(0, 4500), xlab='cell', main='Viral.Percentage vs cell (ordered) - MERS Mock')

#############################################################################################################################
#Calculation and visualization of infection threshold for the Human samples
#infected mers for mers
mers_infected_percentage <- inf_filtered_feature_bc_matrix@meta.data %>% select(5)

#calculate minimal number after 0 to set instead of 0, so they are included into plot
halfmin_mers <- min(mers_infected_percentage$Percent.MERS[mers_infected_percentage$Percent.MERS>0])
halfmin_mers

mers_infected_percentage$Percent.MERS_log10 <- log10(mers_infected_percentage$Percent.MERS+halfmin_mers)
range(mers_infected_percentage$Percent.MERS_log10)

des.all <- density(mers_infected_percentage$Percent.MERS_log10)
min.all <- des.all$x[which(diff(sign(diff(des.all$y)))==2)+1]
min.all
#-2.9801647

#mock for mers
mock_percentage <- mock_filtered_feature_bc_matrix@meta.data %>% select(5)
mock_percentage$Percent.MERS_log10 <- log10(mock_percentage$Percent.MERS+halfmin_mers)
range(mock_percentage$Percent.MERS_log10)


#set colors for infected and non-infected
mers_infected_percentage$group = ifelse(mers_infected_percentage$Percent.MERS_log10 <= min.all[2], 'Non-infected', ifelse(mers_infected_percentage$Percent.MERS_log10, 'Infected'))
mock_percentage$group = ifelse(mock_percentage$Percent.MERS_log10 <= min.all[2], 'Non-infected', ifelse(mock_percentage$Percent.MERS_log10 > min.all[2], 'Infected'))

p1 <- mers_infected_percentage %>%
  ggplot(aes(x = Percent.MERS_log10, fill = group))+
  geom_histogram(bins = 500)+
  scale_fill_manual(values = cols_stat) +
  xlim(-5,1)+
  coord_cartesian(ylim=c(0, 200))+
  ylab("Cell Number")+
  xlab("log10(% of mRNA of MERS-CoV)")+
  ggtitle("MERS-CoV Infected")+
  theme_classic()+
  geom_vline(xintercept = min.all[2], linetype = "dashed", size = 0.9)+
  theme(axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15)))


p2 <- mock_percentage %>%
  ggplot(aes(x = Percent.MERS_log10, fill = group)) +
  geom_histogram(bins = 500) +
  scale_fill_manual(values = "#DEA47E") +
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

png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_MERS/threshold_mers_infection.png",
    width=500, height=600)
p2 / p1
dev.off()

#the infection threshold calculated back from log10 is
10^(-2.9801647)
#mers Infection threshold: 0.001046732




#Threshold above 0.002 Percent Viral is infected

inf_filtered_feature_bc_matrix@meta.data <- inf_filtered_feature_bc_matrix@meta.data %>% mutate(Status = case_when(Percent.MERS > 0.001046732 ~ "Infected",
                                                                                                                   Percent.MERS <= 0.001046732 ~ "Non-infected"))

mock_filtered_feature_bc_matrix <- AddMetaData(mock_filtered_feature_bc_matrix, "Non-infected", col.name = "Status")


##################################################################################################################################
# QC
# Make violin plot to show number of features, number of counts and mitochondrial content
#infected
# VlnPlot(inf_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# VlnPlot(inf_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), ncol = 2)
# VlnPlot(inf_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), ncol = 2)
# VlnPlot(inf_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), y.max = 0.0005, ncol = 2)
# #mock
# VlnPlot(mock_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# VlnPlot(mock_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), ncol = 2)
# VlnPlot(mock_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), ncol = 2)
# VlnPlot(mock_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), y.max = 0.0005, ncol = 2)


#FeatureScatter is typically used to visualize feature-feature relationships, but can be used
#for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# #infected
# plot1 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot3 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "Percent.Viral", feature2 = "percent.mt")
# plot3
# plot1 + plot2 + plot3
# #mock
# plot1 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2


# We filter cells that have unique feature counts over
# We filter cells that have >5% mitochondrial counts


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
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_MERS/cell_numbers_per_virus_before_filtering.png",
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
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_MERS/UMI_per_cell_density_log_before_filtering.png",
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
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_MERS/genes_per_cell_density_log_before_filtering.png",
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
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_MERS/complexity_log_before_filtering.png",
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
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_MERS/Mitochondrial_content_density_before_filtering.png",
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
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_MERS/correlation_gene_UMI_log_before_filtering.png",
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

memory.limit(24000)
#################################################################################################################################
#Integrate the two sample matrices to one seurat object
MERS.list <- c(mock_filtered_feature_bc_matrix, inf_filtered_feature_bc_matrix)
# normalize and identify variable features for each dataset independently
MERS.list <- lapply(X = MERS.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = c("Percent.MERS", "S.Score", "G2M.Score", "percent.mt"), verbose = FALSE)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = MERS.list, nfeatures = 3000)
MERS.list <- PrepSCTIntegration(object.list = MERS.list, anchor.features = features)
MERS.anchors <- FindIntegrationAnchors(object.list = MERS.list, anchor.features = features, normalization.method = "SCT")

# this command creates an 'integrated' data assay
MERS.combined <- IntegrateData(anchorset = MERS.anchors, normalization.method = "SCT")


###Perform an integrated analysis
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(MERS.combined) <- "integrated"

head(MERS.combined@meta.data)
###################################################################################################################

#Saving the seurat object to disk to access in other script for analysis
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens")
SaveH5Seurat(MERS.combined, "MERS.combined", overwrite = TRUE)

#########################################################################################################
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens")
MERS.combined <- LoadH5Seurat("MERS.combined.h5seurat")
MERS.combined


# Add number of genes per UMI for each cell to metadata
MERS.combined$log10GenesPerUMI <- log10(MERS.combined$nFeature_RNA) / log10(MERS.combined$nCount_RNA)

metadata <- MERS.combined@meta.data

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

metadata_count <- metadata %>%
  group_by(Status, Treat) %>%
  count(Status)


#Plots for QC
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_MERS/cell_numbers_per_virus_after_filtering.png",
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
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_MERS/UMI_per_cell_density_log_after_filtering.png",
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
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_MERS/genes_per_cell_density_log_after_filtering.png",
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
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_MERS/complexity_log_after_filtering.png",
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
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_MERS/Mitochondrial_content_density_after_filtering.png",
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
png(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_MERS/correlation_gene_UMI_log_after_filtering.png",
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



