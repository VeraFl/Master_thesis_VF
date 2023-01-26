# Analysis of Single-cell Seq data
# 21.06.2022


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

##########################################################################################################################
#Setting the color schemes for the script:
#Treatment colors (HCoV-229E, MERS-CoV, Mock)
cols_treat = c("#84A59D","#A5668B", "#96ADC8")
show_col(cols_treat)
#Status colors

cols_stat = c("#CD4631","#DEA47E")
show_col(cols_stat)

cols_stat_1 = c("#CD4631","#DEA47E","#DEA47E")
show_col(cols_stat_1)

cols_stat_2 = c("#519872", "#CD4631","#DEA47E")
show_col(cols_stat_2)

#Celltype colors
cols_cell = hcl.colors(7, "Temps")
show_col(cols_cell)

###################################################################################################################################
# Import data which was downloaded from UBELIX after cellranger count run
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/")
# Create each individual Seurat object (if multiple samples)
for (file in c("mock_filtered_feature_bc_matrix", "hcov_inf_filtered_feature_bc_matrix", "mers_inf_filtered_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("data/H. Sapiens/both_virus/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}


# Check the metadata in the new Seurat objects
head(hcov_inf_filtered_feature_bc_matrix@meta.data)
head(mers_inf_filtered_feature_bc_matrix@meta.data)
head(mock_filtered_feature_bc_matrix@meta.data)

#####################################################################################
# QC
#Calculating some QC statistics
#goes through every row and adds all the counts of the columns (molecules per cell)
counts_per_cell_inf_mers <- Matrix::colSums(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts)
counts_per_cell_inf_hcov <- Matrix::colSums(hcov_inf_filtered_feature_bc_matrix@assays$RNA@counts)
counts_per_cell_mock <- Matrix::colSums(mock_filtered_feature_bc_matrix@assays$RNA@counts)

#goes through every column and adds all the counts of the rows (molecules per gene)
counts_per_gene_inf_mers <- Matrix::rowSums(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts)
counts_per_gene_inf_hcov <- Matrix::rowSums(hcov_inf_filtered_feature_bc_matrix@assays$RNA@counts)
counts_per_gene_mock <- Matrix::rowSums(mock_filtered_feature_bc_matrix@assays$RNA@counts)

#count gene itself only if it has non-zero reads mapped (not value).
genes_per_cell_inf_mers <- Matrix::colSums(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts>0)
genes_per_cell_inf_hcov <- Matrix::colSums(hcov_inf_filtered_feature_bc_matrix@assays$RNA@counts>0)

# only count gene if it has non-zero reads mapped (not value).
cells_per_gene_inf_mers <- Matrix::rowSums(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts>0)
cells_per_gene_inf_hcov <- Matrix::rowSums(hcov_inf_filtered_feature_bc_matrix@assays$RNA@counts>0)

# only count cells where the gene is expressed
genes_per_cell_mock <- Matrix::colSums(mock_filtered_feature_bc_matrix@assays$RNA@counts>0) # only count gene if it has non-zero reads mapped (not value).
cells_per_gene_mock <- Matrix::rowSums(mock_filtered_feature_bc_matrix@assays$RNA@counts>0) # only count cells where the gene is expressed

#plotting how many genes per cell we have in the samples
plot(sort(genes_per_cell_inf_mers), ylim = c(1, 10000), xlim = c(0, 5000), xlab='cell', log='y', main='genes per cell (ordered) - Human MERS infected')
plot(sort(genes_per_cell_inf_hcov), ylim = c(1, 10000), xlim = c(0, 5000), xlab='cell', log='y', main='genes per cell (ordered) - Human HCoV infected')
plot(sort(genes_per_cell_inf_mers), ylim = c(1, 10000), xlim = c(0, 5000), xlab='cell', log='y', main='genes per cell (ordered) - Human MERS infected')
plot(sort(genes_per_cell_mock), ylim = c(1, 10000), xlim = c(0, 5000), xlab='cell', log='y', main='genes per cell (ordered) - Human Mock')

#####################################################################################
# Adding important info to the meta data column
#show the number of rows of the matrix to locate the viral genes
nrow(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts)
nrow(hcov_inf_filtered_feature_bc_matrix@assays$RNA@counts)
nrow(mock_filtered_feature_bc_matrix@assays$RNA@counts)
tail(counts_per_gene_inf_mers, 25)
tail(counts_per_gene_inf_hcov, 25)
tail(counts_per_gene_mock, 25)


#show the specific location in the matrix for control
mers_inf_filtered_feature_bc_matrix@assays$RNA@counts[19989:20000,]
hcov_inf_filtered_feature_bc_matrix@assays$RNA@counts[20001:20009,]
mock_filtered_feature_bc_matrix@assays$RNA@counts[19989:20000,] #mers mock
mock_filtered_feature_bc_matrix@assays$RNA@counts[20001:20009,] #hcov mock

#take the right rows in the matrix
mers_infected_mers_counts <- as.matrix(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts[19989:20000,])
mers_infected_hcov_counts <- as.matrix(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts[20001:20009,])

hcov_infected_mers_counts <- as.matrix(hcov_inf_filtered_feature_bc_matrix@assays$RNA@counts[19989:20000,])
hcov_infected_hcov_counts <- as.matrix(hcov_inf_filtered_feature_bc_matrix@assays$RNA@counts[20001:20009,])

mock_mers_counts <- as.matrix(mock_filtered_feature_bc_matrix@assays$RNA@counts[19989:20000,])
mock_hcov_counts <- as.matrix(mock_filtered_feature_bc_matrix@assays$RNA@counts[20001:20009,])

#Count the viral counts
mers.infected.total.mers.counts <- Matrix::colSums(mers_infected_mers_counts)
mers.infected.total.hcov.counts <- Matrix::colSums(mers_infected_hcov_counts)

hcov.infected.total.mers.counts <- Matrix::colSums(hcov_infected_mers_counts)
hcov.infected.total.hcov.counts <- Matrix::colSums(hcov_infected_hcov_counts)

mock.total.mers.counts <- Matrix::colSums(mock_mers_counts)
mock.total.hcov.counts <- Matrix::colSums(mock_hcov_counts)

#Count the total gene counts
mers.infected.total.counts <- Matrix::colSums(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts)
hcov.infected.total.counts <- Matrix::colSums(hcov_inf_filtered_feature_bc_matrix@assays$RNA@counts)
mock.total.counts <- Matrix::colSums(mock_filtered_feature_bc_matrix@assays$RNA@counts)

#Calculate Percent viral for all combinations
mers.infected.percent.mers <- mers.infected.total.mers.counts/mers.infected.total.counts
mers.infected.percent.hcov <- mers.infected.total.hcov.counts/mers.infected.total.counts
hcov.infected.percent.mers <- hcov.infected.total.mers.counts/hcov.infected.total.counts
hcov.infected.percent.hcov <- hcov.infected.total.hcov.counts/hcov.infected.total.counts
mock.percent.mers <- mock.total.mers.counts/mock.total.counts
mock.percent.hcov <- mock.total.hcov.counts/mock.total.counts

#Adding the calculated variables to the metadata of the corresponding sample
mers_inf_filtered_feature_bc_matrix <- AddMetaData(object = mers_inf_filtered_feature_bc_matrix, metadata = mers.infected.total.mers.counts, col.name = "Total.MERS.Counts")
mers_inf_filtered_feature_bc_matrix <- AddMetaData(object = mers_inf_filtered_feature_bc_matrix, metadata = mers.infected.total.hcov.counts, col.name = "Total.hcov.Counts")
mers_inf_filtered_feature_bc_matrix <- AddMetaData(object = mers_inf_filtered_feature_bc_matrix, metadata = mers.infected.percent.mers, col.name = "Percent.MERS")
mers_inf_filtered_feature_bc_matrix <- AddMetaData(object = mers_inf_filtered_feature_bc_matrix, metadata = mers.infected.percent.hcov, col.name = "Percent.hcov")
mers_inf_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(mers_inf_filtered_feature_bc_matrix, pattern = "^MT-")
mers_inf_filtered_feature_bc_matrix <- AddMetaData(mers_inf_filtered_feature_bc_matrix, "MERS-CoV", col.name = "Treat")

head(mers_inf_filtered_feature_bc_matrix@meta.data)

hcov_inf_filtered_feature_bc_matrix <- AddMetaData(object = hcov_inf_filtered_feature_bc_matrix, metadata = hcov.infected.total.mers.counts, col.name = "Total.MERS.Counts")
hcov_inf_filtered_feature_bc_matrix <- AddMetaData(object = hcov_inf_filtered_feature_bc_matrix, metadata = hcov.infected.total.hcov.counts, col.name = "Total.hcov.Counts")
hcov_inf_filtered_feature_bc_matrix <- AddMetaData(object = hcov_inf_filtered_feature_bc_matrix, metadata = hcov.infected.percent.mers, col.name = "Percent.MERS")
hcov_inf_filtered_feature_bc_matrix <- AddMetaData(object = hcov_inf_filtered_feature_bc_matrix, metadata = hcov.infected.percent.hcov, col.name = "Percent.hcov")
hcov_inf_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(hcov_inf_filtered_feature_bc_matrix, pattern = "^MT-")
hcov_inf_filtered_feature_bc_matrix <- AddMetaData(hcov_inf_filtered_feature_bc_matrix, "HCoV-229E", col.name = "Treat")

head(hcov_inf_filtered_feature_bc_matrix@meta.data)

mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.total.mers.counts, col.name = "Total.MERS.Counts")
mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.total.hcov.counts, col.name = "Total.hcov.Counts")
mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.percent.mers, col.name = "Percent.MERS")
mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.percent.hcov, col.name = "Percent.hcov")
mock_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(mock_filtered_feature_bc_matrix, pattern = "^MT-")
mock_filtered_feature_bc_matrix <- AddMetaData(mock_filtered_feature_bc_matrix, "Mock", col.name = "Treat")
head(mock_filtered_feature_bc_matrix@meta.data)



median(mock_filtered_feature_bc_matrix@meta.data$nFeature_RNA)
median(hcov_inf_filtered_feature_bc_matrix@meta.data$nFeature_RNA)
median(mers_inf_filtered_feature_bc_matrix@meta.data$nFeature_RNA)

mean(mock_filtered_feature_bc_matrix@meta.data$nFeature_RNA)
mean(hcov_inf_filtered_feature_bc_matrix@meta.data$nFeature_RNA)
mean(mers_inf_filtered_feature_bc_matrix@meta.data$nFeature_RNA)

#Adding a Meta data column with information about the Status of each cell. If the cell is Infected or Bystander
#in the MERS infected Treatment group or non-infected in the Mock Group


# #Extracting the Viral-Percentage column for every cell and transpose it to plot it below
# mers_infected_mers_percentage <- mers_inf_filtered_feature_bc_matrix@meta.data %>% select(6) %>% t()
# mers_infected_hcov_percentage <- mers_inf_filtered_feature_bc_matrix@meta.data %>% select(7) %>% t()
# hcov_infected_mers_percentage <- hcov_inf_filtered_feature_bc_matrix@meta.data %>% select(6) %>% t()
# hcov_infected_hcov_percentage <- hcov_inf_filtered_feature_bc_matrix@meta.data %>% select(7) %>% t()
# Mock_mers_percentage <- mock_filtered_feature_bc_matrix@meta.data %>% select(6) %>% t()
# Mock_hcov_percentage <- mock_filtered_feature_bc_matrix@meta.data %>% select(7) %>% t()
# 
# #plotting how many virus positive cells there are in each treatment
# plot(sort(mers_infected_mers_percentage), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Human MERS infected MERS counts')
# plot(sort(mers_infected_hcov_percentage), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Human MERS infected hcov counts')
# plot(sort(hcov_infected_mers_percentage), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Human hcov infected MERS counts')
# plot(sort(hcov_infected_hcov_percentage), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Human hcov infected hcov counts')
# plot(sort(Mock_mers_percentage), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Human Mock MERS counts')
# plot(sort(Mock_hcov_percentage), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Human Mock hcov counts')


#############################################################################################################################
#Calculation and visualization of infection threshold for the Human samples
#infected mers for mers
mers_infected_mers_percentage <- mers_inf_filtered_feature_bc_matrix@meta.data %>% select(6)

#calculate minimal number after 0 to set instead of 0, so they are included into plot
halfmin_mers <- min(mers_infected_mers_percentage$Percent.MERS[mers_infected_mers_percentage$Percent.MERS>0])
halfmin_mers
#1.479837e-05
mers_infected_mers_percentage$Percent.MERS_log10 <- log10(mers_infected_mers_percentage$Percent.MERS+halfmin_mers)
range(mers_infected_mers_percentage$Percent.MERS_log10)
#-4.8297861 -0.3359357
des.all <- density(mers_infected_mers_percentage$Percent.MERS_log10)
min.all <- des.all$x[which(diff(sign(diff(des.all$y)))==2)+1]
min.all
#-4.5565287 -2.9599637 -1.7258091 -0.9128342 -0.5308340
#chosen minimum: -2.9599637

#infected hcov for mers (fake mock)
hcov_infected_mers_percentage <- hcov_inf_filtered_feature_bc_matrix@meta.data %>% select(6)
hcov_infected_mers_percentage$Percent.MERS_log10 <- log10(hcov_infected_mers_percentage$Percent.MERS+halfmin_mers)
range(hcov_infected_mers_percentage$Percent.MERS_log10)
#-4.829786 -3.766074

#mock for mers
mock_percentage <- mock_filtered_feature_bc_matrix@meta.data %>% select(6)
mock_percentage$Percent.MERS_log10 <- log10(mock_percentage$Percent.MERS+halfmin_mers)
range(mock_percentage$Percent.MERS_log10)
#-4.829786 -3.844454

#set colors for infected and non-infected
mers_infected_mers_percentage$group = ifelse(mers_infected_mers_percentage$Percent.MERS_log10 <= min.all[2], 'Uninfected', ifelse(mers_infected_mers_percentage$Percent.MERS_log10, 'Infected'))
hcov_infected_mers_percentage$group = ifelse(hcov_infected_mers_percentage$Percent.MERS_log10 <= min.all[2], 'Uninfected', ifelse(hcov_infected_mers_percentage$Percent.MERS_log10, 'Infected'))
mock_percentage$group = ifelse(mock_percentage$Percent.MERS_log10 <= min.all[2], 'Uninfected', ifelse(mock_percentage$Percent.MERS_log10 > min.all[2], 'Infected'))

p1 <- mers_infected_mers_percentage %>% 
  ggplot(aes(x = Percent.MERS_log10, fill = group))+
  geom_histogram(bins = 500)+
  scale_fill_manual(values = cols_stat) +
  xlim(-5,1)+
  coord_cartesian(ylim=c(0, 60))+
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

p2 <- hcov_infected_mers_percentage %>% 
  ggplot(aes(x = Percent.MERS_log10, fill = group)) +
  geom_histogram(bins = 500) +
  scale_fill_manual(values = "#DEA47E") +
  xlim(-5,1) +
  coord_cartesian(ylim=c(0, 60))+
  ylab("Cell Number")+
  xlab("")+
  ggtitle("HCoV-229E Infected")+
  theme_classic() +
  geom_vline(xintercept = min.all[2], linetype = "dashed", size = 0.9)+
  theme(axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15)))+
  NoLegend()

p3 <- mock_percentage %>% 
  ggplot(aes(x = Percent.MERS_log10, fill = group)) +
  geom_histogram(bins = 500) +
  scale_fill_manual(values = "#DEA47E") +
  xlim(-5,1) +
  coord_cartesian(ylim=c(0, 60))+
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

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/threshold_mers_infection.pdf",
    width=5, height=8)
p3 / p2 / p1
dev.off()

#the infection threshold calculated back from log10 is
10^(-2.9599637)

###################################################################################################################
#infected hcov for hcov
hcov_infected_hcov_percentage <- hcov_inf_filtered_feature_bc_matrix@meta.data %>% select(7)

#calculate minimal number after 0 to set instead of 0, so they are included into plot
halfmin_hcov <- min(hcov_infected_hcov_percentage$Percent.hcov[hcov_infected_hcov_percentage$Percent.hcov>0])
halfmin_hcov
#1.563159e-05
log10(1.563159e-05)

hcov_infected_hcov_percentage$Percent.hcov_log10 <- log10(hcov_infected_hcov_percentage$Percent.hcov+halfmin_hcov)
range(hcov_infected_hcov_percentage$Percent.hcov_log10)
#-4.8059967 -0.9366344

des.all <- density(hcov_infected_hcov_percentage$Percent.hcov_log10)
min.all <- des.all$x[which(diff(sign(diff(des.all$y)))==2)+1]
min.all
#-4.603279 -4.118814 -3.973474 -3.432488 -3.222553 -3.044916 -2.600822 -1.736859 -1.470404 -1.446180 -1.413883 -1.381585 -1.316990
#-3.432488

#infected mers for hcov (fake mock)
mers_infected_hcov_percentage <- mers_inf_filtered_feature_bc_matrix@meta.data %>% select(7)
mers_infected_hcov_percentage$Percent.hcov_log10 <- log10(mers_infected_hcov_percentage$Percent.hcov+halfmin_hcov)
range(mers_infected_hcov_percentage$Percent.hcov_log10)

#mock for hcov
mock_percentage <- mock_filtered_feature_bc_matrix@meta.data %>% select(7)
mock_percentage$Percent.hcov_log10 <- log10(mock_percentage$Percent.hcov+halfmin_hcov)
range(mock_percentage$Percent.hcov_log10)


#set colors for infected and non-infected
hcov_infected_hcov_percentage$group = ifelse(hcov_infected_hcov_percentage$Percent.hcov_log10 <= min.all[4], 'Uninfected', ifelse(hcov_infected_hcov_percentage$Percent.hcov_log10, 'Infected'))
mers_infected_hcov_percentage$group = ifelse(mers_infected_hcov_percentage$Percent.hcov_log10 <= min.all[4], 'Uninfected', ifelse(mers_infected_hcov_percentage$Percent.hcov_log10, 'Infected'))
mock_percentage$group = ifelse(mock_percentage$Percent.hcov_log10 <= min.all[4], 'Uninfected', ifelse(mock_percentage$Percent.hcov_log10, 'Infected'))

p1 <- hcov_infected_hcov_percentage %>% 
  ggplot(aes(x = Percent.hcov_log10, fill = group))+
  geom_histogram(bins = 500)+
  scale_fill_manual(values = cols_stat) +
  xlim(-5,1)+
  coord_cartesian(ylim=c(0, 60))+
  ylab("Cell Number")+
  xlab("log10(% of mRNA of hcov-CoV)")+
  ggtitle("HCoV-229E Infected")+
  theme_classic()+
  geom_vline(xintercept = min.all[4], linetype = "dashed", size = 0.9)+
  theme(axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15)))

p2 <- mers_infected_hcov_percentage %>% 
  ggplot(aes(x = Percent.hcov_log10, fill = group)) +
  geom_histogram(bins = 500) +
  scale_fill_manual(values = cols_stat) +
  xlim(-5,1) +
  coord_cartesian(ylim=c(0, 60))+
  ylab("Cell Number")+
  xlab("")+
  ggtitle("MERS-CoV Infected")+
  theme_classic() +
  geom_vline(xintercept = min.all[4], linetype = "dashed", size = 0.9)+
  theme(axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15)))+
  NoLegend()

p3 <- mock_percentage %>% 
  ggplot(aes(x = Percent.hcov_log10, fill = group)) +
  geom_histogram(bins = 500) +
  scale_fill_manual(values = cols_stat) +
  xlim(-5,1) +
  coord_cartesian(ylim=c(0, 60))+
  ylab("Cell Number")+
  xlab("")+
  ggtitle("Mock")+
  theme_classic() +
  geom_vline(xintercept = min.all[4], linetype = "dashed", size = 0.9)+
  theme(axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15)))+
  NoLegend()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/threshold_hcov_infection.pdf",
    width=5, height=8)
p3 / p2 / p1
dev.off()

10^(-3.432488)
#Infection threshold: 0.0003694129


#Threshold above 0.002 Percent Viral is infected

mers_inf_filtered_feature_bc_matrix@meta.data <- mers_inf_filtered_feature_bc_matrix@meta.data %>% mutate(Status = case_when(Percent.MERS > 0.00109657 ~ "Infected",
                                                                                                                             Percent.MERS <= 0.00109657 ~ "Bystander"))
hcov_inf_filtered_feature_bc_matrix@meta.data <- hcov_inf_filtered_feature_bc_matrix@meta.data %>% mutate(Status = case_when(Percent.hcov > 0.0003694129 ~ "Infected",
                                                                                                                                       Percent.hcov <= 0.0003694129 ~ "Bystander"))

mock_filtered_feature_bc_matrix <- AddMetaData(mock_filtered_feature_bc_matrix, "Uninfected", col.name = "Status")

##################################################################################################################
# # QC
# # Make violin plot to show number of features, number of counts and mitochondrial content
# #mers infected
# VlnPlot(mers_inf_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# VlnPlot(mers_inf_filtered_feature_bc_matrix, features = c("Percent.MERS", "Total.MERS.Counts"), y.max = 0.8, ncol = 2)
# VlnPlot(mers_inf_filtered_feature_bc_matrix, features = c("Percent.hcov", "Total.hcov.Counts"), y.max = 0.8, ncol = 2)
# 
# #hcov infected
# VlnPlot(hcov_inf_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# VlnPlot(hcov_inf_filtered_feature_bc_matrix, features = c("Percent.MERS", "Total.MERS.Counts"), ncol = 2)
# VlnPlot(hcov_inf_filtered_feature_bc_matrix, features = c("Percent.hcov", "Total.hcov.Counts"), ncol = 2)
# #mock
# VlnPlot(mock_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# VlnPlot(mock_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), ncol = 2)
# VlnPlot(mock_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), ncol = 2)
# VlnPlot(mock_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), y.max = 0.8, ncol = 2)
# 
# ?VlnPlot()
# #FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# #for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# #infected
# plot1 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot3 <- FeatureScatter(inf_filtered_feature_bc_matrix, feature1 = "Percent.Viral", feature2 = "percent.mt")
# plot1 + plot2 + plot3
# #mock
# plot1 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2


# Add number of genes per UMI for each cell to metadata
mers_inf_filtered_feature_bc_matrix$log10GenesPerUMI <- log10(mers_inf_filtered_feature_bc_matrix$nFeature_RNA) / log10(mers_inf_filtered_feature_bc_matrix$nCount_RNA)
hcov_inf_filtered_feature_bc_matrix$log10GenesPerUMI <- log10(hcov_inf_filtered_feature_bc_matrix$nFeature_RNA) / log10(hcov_inf_filtered_feature_bc_matrix$nCount_RNA)
mock_filtered_feature_bc_matrix$log10GenesPerUMI <- log10(mock_filtered_feature_bc_matrix$nFeature_RNA) / log10(mock_filtered_feature_bc_matrix$nCount_RNA)


metadata <- bind_rows(mers_inf_filtered_feature_bc_matrix@meta.data, hcov_inf_filtered_feature_bc_matrix@meta.data, mock_filtered_feature_bc_matrix@meta.data)
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


metadata_count <- metadata %>%
  group_by(Status, Treat) %>%
  count(Status)

levels(metadata_count$Status)

metadata_count$Status <- factor(x = metadata_count$Status, levels = c('Infected','Bystander', 'Uninfected'))



#Plots for QC
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/cell_numbers_per_virus_before_filtering.pdf",
    width=6, height=5)
metadata_count %>% 
  ggplot(aes(x=Treat, y=n, fill=Status)) + 
  geom_col() +
  theme_classic() +
  scale_fill_manual(values = cols_stat_1)+
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

?geom_density

#Number of UMI (transcripts) per cell
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/UMI_per_cell_density_log_before_filtering.pdf",
    width=9, height=6)
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
        axis.text.x = element_text(colour = "black", size = 20),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold"))
dev.off()

?geom_density
#Number of Genes per cell
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/genes_per_cell_density_log_before_filtering.pdf",
    width=9, height=6)
metadata %>% 
  ggplot(aes(color=Treat, x=nGene, fill= Treat)) + 
  geom_density(alpha = 0.6) +
  scale_color_manual(values= cols_treat) +
  scale_fill_manual(values = cols_treat) +
  theme_classic() +
  scale_x_log10(labels = comma_format(big.mark = "'",
                                      decimal.mark = ",")) + 
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.text.x = element_text(colour = "black", size = 20),
        axis.text.y = element_text(colour = "black", size = 20),
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/complexity_log_before_filtering.pdf",
    width=9, height=6)
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
        axis.text.x = element_text(colour = "black", size = 20),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold"))
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per cell
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/Mitochondrial_content_density_before_filtering.pdf",
    width=9, height=6)
metadata %>% 
  ggplot(aes(color=Treat, x=percent.mt, fill=Treat)) + 
  geom_density(alpha = 0.6) + 
  scale_color_manual(values= cols_treat) +
  scale_fill_manual(values = cols_treat) +
  xlim(0,50)+
  theme_classic() +
  geom_vline(xintercept = 30)+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.text.x = element_text(colour = "black", size = 20),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold")) +
  ylab("Cell density") +
  xlab("Percentage of mitochondrial transcripts")
dev.off()


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/correlation_gene_UMI_log_before_filtering.pdf",
    width=10, height=5)
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



# We filter cells that have unique feature counts over 1000
# We filter cells that have <30% mitochondrial counts
#infected
mers_inf_filtered_feature_bc_matrix <- subset(mers_inf_filtered_feature_bc_matrix, subset = nFeature_RNA > 1000 & percent.mt < 30)
hcov_inf_filtered_feature_bc_matrix <- subset(hcov_inf_filtered_feature_bc_matrix, subset = nFeature_RNA > 1000 & percent.mt < 30)
#mock
mock_filtered_feature_bc_matrix <- subset(mock_filtered_feature_bc_matrix, subset = nFeature_RNA > 1000 & percent.mt < 30)

##################################################################################################################

# #Plot to show Violin Plots of the infected cells
# MERS <- mers_inf_filtered_feature_bc_matrix
# hcov<- hcov_inf_filtered_feature_bc_matrix
# Mock <- mock_filtered_feature_bc_matrix
# 
# Idents(MERS) <- "Treat"
# Idents(hcov) <- "Treat"
# Idents(Mock) <- "Treat"
# 
# p1 <- VlnPlot(MERS, features = "Percent.MERS", y.max = 0.8, ncol = 1) + NoLegend() + FontSize(x.text = 0, x.title = 0)
# p2 <- VlnPlot(hcov, features = "Percent.hcov", y.max = 0.8, ncol = 1) + NoLegend() + FontSize(x.text = 0, x.title = 0)
# p3 <- VlnPlot(Mock, features = c("Percent.MERS"), y.max = 0.8, ncol = 1) + NoLegend() + FontSize(x.text = 0, x.title = 0)
# p4 <- VlnPlot(Mock, features = c("Percent.hcov"), y.max = 0.8, ncol = 1) + NoLegend() + FontSize(x.text = 0, x.title = 0)
# 
# grid.arrange(p1, p2, p3, p4, ncol = 2)
# 

###################################################################################################################
#cell cycle state genes from the Seurat package
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

g2m.genes <- c(g2m.genes, "PIMREG", "JPT1") #renamed genes have to be added to the list
s.genes <- c(s.genes, "CENPU") #renamed genes have to be added to the list

memory.limit(24000)
#################################################################################################################################
#Integrate the two sample matrices to one seurat object
Human.list <- c(mock_filtered_feature_bc_matrix, mers_inf_filtered_feature_bc_matrix, hcov_inf_filtered_feature_bc_matrix)
# normalize and identify variable features for each dataset independently
Human.list <- lapply(X = Human.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  x <- SCTransform(x, vars.to.regress = c("Percent.MERS", "Percent.hcov", "S.Score", "G2M.Score", "percent.mt"), verbose = FALSE)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Human.list, nfeatures = 3000)
Human.list <- PrepSCTIntegration(object.list = Human.list, anchor.features = features)
Human.anchors <- FindIntegrationAnchors(object.list = Human.list, anchor.features = features, normalization.method = "SCT")

# this command creates an 'integrated' data assay
Human.combined <- IntegrateData(anchorset = Human.anchors, normalization.method = "SCT")

###Perform an integrated analysis
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(Human.combined) <- "integrated"

head(Human.combined@meta.data)
###################################################################################################################

#Saving the seurat object to disk to access in other script for analysis
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens")
SaveH5Seurat(Human.combined, "Human.combined", overwrite = TRUE)


###################################################################################################################
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens")
Human.combined <- LoadH5Seurat("Human.combined.h5seurat")
Human.combined


# Add number of genes per UMI for each cell to metadata
Human.combined$log10GenesPerUMI <- log10(Human.combined$nFeature_RNA) / log10(Human.combined$nCount_RNA)


metadata <- Human.combined@meta.data
# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


metadata_count <- metadata %>%
  group_by(Status, Treat) %>%
  count(Status)

levels(metadata_count$Status)

metadata_count$Status <- factor(x = metadata_count$Status, levels = c('Infected','Bystander', 'Uninfected'))


#Plots for QC
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/cell_numbers_per_virus_after_filtering.pdf",
    width=9, height=6)
metadata_count %>% 
  ggplot(aes(x=Treat, y=n, fill=Status)) + 
  geom_col() +
  theme_classic() +
  scale_fill_manual(values = cols_stat_1)+
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/UMI_per_cell_density_log_after_filtering.pdf",
    width=9, height=6)
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/genes_per_cell_density_log_after_filtering.pdf",
    width=9, height=6)
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/complexity_log_after_filtering.pdf",
    width=9, height=6)
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
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(r = 20)),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 20, face = "bold"))
dev.off()


# Visualize the distribution of mitochondrial gene expression detected per cell
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/Mitochondrial_content_density_after_filtering.pdf",
    width=9, height=6)
metadata %>% 
  ggplot(aes(color=Treat, x=percent.mt, fill=Treat)) + 
  geom_density(alpha = 0.6) + 
  theme_classic() +
  scale_color_manual(values= cols_treat) +
  scale_fill_manual(values = cols_treat) +
  xlim(0,50)+
  geom_vline(xintercept = 30, linetype ="dashed", size = 0.9)+
  ylab("Cell density") +
  xlab("Percentage of mitochondrial transcripts") +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(r = 20)),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 20, face = "bold"))
dev.off()

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/correlation_gene_UMI_log_filtering.pdf",
    width=15, height=6)
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
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(r = 20)),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 20, face = "bold"),
        strip.text = element_text(size = 15))+
  facet_wrap(~Treat)
dev.off()
