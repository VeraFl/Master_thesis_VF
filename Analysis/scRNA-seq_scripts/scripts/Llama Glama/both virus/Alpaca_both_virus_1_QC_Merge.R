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
library(glmGamPoi)

?Devices
###################################################################################################################################
#Color palettes
#Treatment colors(camel-229E, MERS, Mock)
cols_treat = c("#F39C6B","#A5668B", "#96ADC8")
show_col(cols_treat)
#Status colors
cols_stat = c("#CD4631","#DEA47E")
show_col(cols_stat)

cols_stat_2 = c("#CD4631", "#DEA47E", "#DEA47E")

#Celltype colors
cols_cell = hcl.colors(7, "Temps")
show_col(cols_cell)


###################################################################################################################################

# Import data which was downloaded from UBELIX after cellranger count run
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/")
# Create each individual Seurat object (if multiple samples)
for (file in c("mock_filtered_feature_bc_matrix", "camel229E_inf_filtered_feature_bc_matrix", "mers_inf_filtered_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("data/L. Glama/both_virus/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}


# Check the metadata in the new Seurat objects
head(camel229E_inf_filtered_feature_bc_matrix@meta.data)
head(mers_inf_filtered_feature_bc_matrix@meta.data)
head(mock_filtered_feature_bc_matrix@meta.data)


#####################################################################################
# QC
#Calculating some QC statistics
#goes through every row and adds all the counts of the columns (molecules per cell)
counts_per_cell_inf_mers <- Matrix::colSums(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts)
counts_per_cell_inf_camel229E <- Matrix::colSums(camel229E_inf_filtered_feature_bc_matrix@assays$RNA@counts)
counts_per_cell_mock <- Matrix::colSums(mock_filtered_feature_bc_matrix@assays$RNA@counts)

#goes through every column and adds all the counts of the rows (molecules per gene)
counts_per_gene_inf_mers <- Matrix::rowSums(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts)
counts_per_gene_inf_camel229E <- Matrix::rowSums(camel229E_inf_filtered_feature_bc_matrix@assays$RNA@counts)
counts_per_gene_mock <- Matrix::rowSums(mock_filtered_feature_bc_matrix@assays$RNA@counts)

#count gene itself only if it has non-zero reads mapped (not value).
genes_per_cell_inf_mers <- Matrix::colSums(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts>0)
genes_per_cell_inf_camel229E <- Matrix::colSums(camel229E_inf_filtered_feature_bc_matrix@assays$RNA@counts>0)

# only count gene if it has non-zero reads mapped (not value).
cells_per_gene_inf_mers <- Matrix::rowSums(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts>0)
cells_per_gene_inf_camel229E <- Matrix::rowSums(camel229E_inf_filtered_feature_bc_matrix@assays$RNA@counts>0)

# only count cells where the gene is expressed
genes_per_cell_mock <- Matrix::colSums(mock_filtered_feature_bc_matrix@assays$RNA@counts>0) # only count gene if it has non-zero reads mapped (not value).
cells_per_gene_mock <- Matrix::rowSums(mock_filtered_feature_bc_matrix@assays$RNA@counts>0) # only count cells where the gene is expressed
genes_per_cell_mock <- Matrix::colSums(mock_filtered_feature_bc_matrix@assays$RNA@counts>0) # only count gene if it has non-zero reads mapped (not value).
cells_per_gene_mock <- Matrix::rowSums(mock_filtered_feature_bc_matrix@assays$RNA@counts>0) # only count cells where the gene is expressed

#plotting how many genes per cell we have in the samples
plot(sort(genes_per_cell_inf_mers), ylim = c(1, 10000), xlim = c(0, 5000), xlab='cell', log='y', main='genes per cell (ordered) - Alpaca MERS infected')
plot(sort(genes_per_cell_inf_camel229E), ylim = c(1, 10000), xlim = c(0, 5000), xlab='cell', log='y', main='genes per cell (ordered) - Alpaca camel229E infected')

plot(sort(genes_per_cell_mock), ylim = c(1, 10000), xlim = c(0, 5000), xlab='cell', log='y', main='genes per cell (ordered) - Alpaca  Mock')


#####################################################################################
# Adding important info to the meta data column
#show the number of rows of the matrix to locate the viral genes
nrow(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts)
nrow(camel229E_inf_filtered_feature_bc_matrix@assays$RNA@counts)
nrow(mock_filtered_feature_bc_matrix@assays$RNA@counts)

# display the last 25 rows in the counts per gene (from above) to see the total viral counts in all the cells
tail(counts_per_gene_inf_mers, 25)
tail(counts_per_gene_inf_camel229E, 25)
tail(counts_per_gene_mock, 25)


#show the specific location in the matrix for control
mers_inf_filtered_feature_bc_matrix@assays$RNA@counts[28366:28377,]
camel229E_inf_filtered_feature_bc_matrix@assays$RNA@counts[28378:28387,]
mock_filtered_feature_bc_matrix@assays$RNA@counts[28366:28377,] #mers mock
mock_filtered_feature_bc_matrix@assays$RNA@counts[28378:28387,] #camel229E mock

#take the right rows in the matrix
mers_infected_mers_counts <- as.matrix(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts[28366:28377,])
mers_infected_camel229E_counts <- as.matrix(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts[28378:28387,])

camel229E_infected_mers_counts <- as.matrix(camel229E_inf_filtered_feature_bc_matrix@assays$RNA@counts[28366:28377,])
camel229E_infected_camel229E_counts <- as.matrix(camel229E_inf_filtered_feature_bc_matrix@assays$RNA@counts[28378:28387,])

mock_mers_counts <- as.matrix(mock_filtered_feature_bc_matrix@assays$RNA@counts[28366:28377,])
mock_camel229E_counts <- as.matrix(mock_filtered_feature_bc_matrix@assays$RNA@counts[28378:28387,])

#Count the viral counts
mers.infected.total.mers.counts <- Matrix::colSums(mers_infected_mers_counts)
mers.infected.total.camel229E.counts <- Matrix::colSums(mers_infected_camel229E_counts)

camel229E.infected.total.mers.counts <- Matrix::colSums(camel229E_infected_mers_counts)
camel229E.infected.total.camel229E.counts <- Matrix::colSums(camel229E_infected_camel229E_counts)

mock.total.mers.counts <- Matrix::colSums(mock_mers_counts)
mock.total.camel229E.counts <- Matrix::colSums(mock_camel229E_counts)

#Count the total gene counts
mers.infected.total.counts <- Matrix::colSums(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts)
camel229E.infected.total.counts <- Matrix::colSums(camel229E_inf_filtered_feature_bc_matrix@assays$RNA@counts)
mock.total.counts <- Matrix::colSums(mock_filtered_feature_bc_matrix@assays$RNA@counts)

#Calculate Percent viral for all combinations
mers.infected.percent.mers <- mers.infected.total.mers.counts/mers.infected.total.counts
mers.infected.percent.camel229E <- mers.infected.total.camel229E.counts/mers.infected.total.counts
camel229E.infected.percent.mers <- camel229E.infected.total.mers.counts/camel229E.infected.total.counts
camel229E.infected.percent.camel229E <- camel229E.infected.total.camel229E.counts/camel229E.infected.total.counts
mock.percent.mers <- mock.total.mers.counts/mock.total.counts
mock.percent.camel229E <- mock.total.camel229E.counts/mock.total.counts

#Adding the calculated variables to the metadata of the corresponding sample
mers_inf_filtered_feature_bc_matrix <- AddMetaData(object = mers_inf_filtered_feature_bc_matrix, metadata = mers.infected.total.mers.counts, col.name = "Total.MERS.Counts")
mers_inf_filtered_feature_bc_matrix <- AddMetaData(object = mers_inf_filtered_feature_bc_matrix, metadata = mers.infected.total.camel229E.counts, col.name = "Total.camel229E.Counts")
mers_inf_filtered_feature_bc_matrix <- AddMetaData(object = mers_inf_filtered_feature_bc_matrix, metadata = mers.infected.percent.mers, col.name = "Percent.MERS")
mers_inf_filtered_feature_bc_matrix <- AddMetaData(object = mers_inf_filtered_feature_bc_matrix, metadata = mers.infected.percent.camel229E, col.name = "Percent.camel229E")
mers_inf_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(mers_inf_filtered_feature_bc_matrix, pattern = "^KEF")
mers_inf_filtered_feature_bc_matrix <- AddMetaData(mers_inf_filtered_feature_bc_matrix, "MERS-CoV", col.name = "Treat")

head(mers_inf_filtered_feature_bc_matrix@meta.data)

camel229E_inf_filtered_feature_bc_matrix <- AddMetaData(object = camel229E_inf_filtered_feature_bc_matrix, metadata = camel229E.infected.total.mers.counts, col.name = "Total.MERS.Counts")
camel229E_inf_filtered_feature_bc_matrix <- AddMetaData(object = camel229E_inf_filtered_feature_bc_matrix, metadata = camel229E.infected.total.camel229E.counts, col.name = "Total.camel229E.Counts")
camel229E_inf_filtered_feature_bc_matrix <- AddMetaData(object = camel229E_inf_filtered_feature_bc_matrix, metadata = camel229E.infected.percent.mers, col.name = "Percent.MERS")
camel229E_inf_filtered_feature_bc_matrix <- AddMetaData(object = camel229E_inf_filtered_feature_bc_matrix, metadata = camel229E.infected.percent.camel229E, col.name = "Percent.camel229E")
camel229E_inf_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(camel229E_inf_filtered_feature_bc_matrix, pattern = "^KEF")
camel229E_inf_filtered_feature_bc_matrix <- AddMetaData(camel229E_inf_filtered_feature_bc_matrix, "dcCoV-ACN4", col.name = "Treat")

head(camel229E_inf_filtered_feature_bc_matrix@meta.data)

mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.total.mers.counts, col.name = "Total.MERS.Counts")
mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.total.camel229E.counts, col.name = "Total.camel229E.Counts")
mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.percent.mers, col.name = "Percent.MERS")
mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.percent.camel229E, col.name = "Percent.camel229E")
mock_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(mock_filtered_feature_bc_matrix, pattern = "^KEF")
mock_filtered_feature_bc_matrix <- AddMetaData(mock_filtered_feature_bc_matrix, "Mock", col.name = "Treat")
head(mock_filtered_feature_bc_matrix@meta.data)


median(mock_filtered_feature_bc_matrix@meta.data$nFeature_RNA)
median(camel229E_inf_filtered_feature_bc_matrix@meta.data$nFeature_RNA)
median(mers_inf_filtered_feature_bc_matrix@meta.data$nFeature_RNA)

mean(mock_filtered_feature_bc_matrix@meta.data$nFeature_RNA)
mean(camel229E_inf_filtered_feature_bc_matrix@meta.data$nFeature_RNA)
mean(mers_inf_filtered_feature_bc_matrix@meta.data$nFeature_RNA)

#Adding a Meta data column with information about the Status of each cell. If the cell is Infected or Bystander
#in the MERS infected Treatment group or non-infected in the Mock Group


# #Extracting the Viral-Percentage column for every cell and transpose it to plot it below
# mers_infected_mers_percentage <- mers_inf_filtered_feature_bc_matrix@meta.data %>% select(6) %>% t()
# mers_infected_camel229E_percentage <- mers_inf_filtered_feature_bc_matrix@meta.data %>% select(7) %>% t()
# camel229E_infected_mers_percentage <- camel229E_inf_filtered_feature_bc_matrix@meta.data %>% select(6) %>% t()
# camel229E_infected_camel229E_percentage <- camel229E_inf_filtered_feature_bc_matrix@meta.data %>% select(7) %>% t()
# Mock_mers_percentage <- mock_filtered_feature_bc_matrix@meta.data %>% select(6) %>% t()
# Mock_camel229E_percentage <- mock_filtered_feature_bc_matrix@meta.data %>% select(7) %>% t()
# 
# #plotting how many virus positive cells there are in each treatment
# plot(sort(mers_infected_mers_percentage), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Alpaca MERS infected MERS counts')
# plot(sort(mers_infected_camel229E_percentage), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Alpaca MERS infected camel229E counts')
# plot(sort(camel229E_infected_mers_percentage), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Alpaca camel229E infected MERS counts')
# plot(sort(camel229E_infected_camel229E_percentage), ylim = c(0,0.002), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Alpaca camel229E infected camel229E counts')
# plot(sort(Mock_mers_percentage), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Alpaca Mock MERS counts')
# plot(sort(Mock_camel229E_percentage), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Alpaca Mock camel229E counts')
# 

#############################################################################################################################
#Calculation and visualization of infection threshold for the Human samples
#infected mers for mers
mers_infected_mers_percentage <- mers_inf_filtered_feature_bc_matrix@meta.data %>% select(6)

#calculate minimal number after 0 to set instead of 0, so they are included into plot
halfmin_mers <- min(mers_infected_mers_percentage$Percent.MERS[mers_infected_mers_percentage$Percent.MERS>0])
halfmin_mers
#7.638252e-05

mers_infected_mers_percentage$Percent.MERS_log10 <- log10(mers_infected_mers_percentage$Percent.MERS+halfmin_mers)
range(mers_infected_mers_percentage$Percent.MERS_log10)
#-4.11700600 -0.07863714

des.all <- density(mers_infected_mers_percentage$Percent.MERS_log10)
min.all <- des.all$x[which(diff(sign(diff(des.all$y)))==2)+1]
min.all
#-3.9112885 -2.4898061 -2.0676689 -1.6713768 -1.2406246 -0.7581821 -0.4394254
#-2.4898061 [2]

#infected camel229E for mers (fake mock)
camel229E_infected_mers_percentage <- camel229E_inf_filtered_feature_bc_matrix@meta.data %>% select(6)
camel229E_infected_mers_percentage$Percent.MERS_log10 <- log10(camel229E_infected_mers_percentage$Percent.MERS+halfmin_mers)
range(camel229E_infected_mers_percentage$Percent.MERS_log10)

#mock for mers
mock_percentage <- mock_filtered_feature_bc_matrix@meta.data %>% select(6)
mock_percentage$Percent.MERS_log10 <- log10(mock_percentage$Percent.MERS+halfmin_mers)
range(mock_percentage$Percent.MERS_log10)


#set colors for infected and non-infected
mers_infected_mers_percentage$group = ifelse(mers_infected_mers_percentage$Percent.MERS_log10 <= min.all[2], 'Uninfected', ifelse(mers_infected_mers_percentage$Percent.MERS_log10, 'Infected'))
camel229E_infected_mers_percentage$group = ifelse(camel229E_infected_mers_percentage$Percent.MERS_log10 <= min.all[2], 'Uninfected', ifelse(camel229E_infected_mers_percentage$Percent.MERS_log10, 'Infected'))
mock_percentage$group = ifelse(mock_percentage$Percent.MERS_log10 <= min.all[2], 'Uninfected', ifelse(mock_percentage$Percent.MERS_log10 > min.all[2], 'Infected'))

p1 <- mers_infected_mers_percentage %>%
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

p2 <- camel229E_infected_mers_percentage %>%
  ggplot(aes(x = Percent.MERS_log10, fill = group)) +
  geom_histogram(bins = 500) +
  scale_fill_manual(values = "#DEA47E") +
  xlim(-5,1) +
  coord_cartesian(ylim=c(0, 200))+
  ylab("Cell Number")+
  xlab("")+
  ggtitle("Camel229E Infected")+
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

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/threshold_mers_infection.pdf",
    width=5, height=8)
p3 / p2 / p1
dev.off()

#the infection threshold calculated back from log10 is
10^(-2.4898061)

#mers Infection threshold: 0.003237382

###################################################################################################################
#infected camel229E for camel229E
camel229E_infected_camel229E_percentage <- camel229E_inf_filtered_feature_bc_matrix@meta.data %>% select(7)

#calculate minimal number after 0 to set instead of 0, so they are included into plot
halfmin_camel229E <- min(camel229E_infected_camel229E_percentage$Percent.camel229E[camel229E_infected_camel229E_percentage$Percent.camel229E>0])
halfmin_camel229E
#6.799483e-05

#log10 transform and add halfmin to all values
camel229E_infected_camel229E_percentage$Percent.camel229E_log10 <- log10(camel229E_infected_camel229E_percentage$Percent.camel229E+halfmin_camel229E)
range(camel229E_infected_camel229E_percentage$Percent.camel229E_log10)
#-4.1675241 -0.1844866

#calculate density
des.all <- density(camel229E_infected_camel229E_percentage$Percent.camel229E_log10)
min.all <- des.all$x[which(diff(sign(diff(des.all$y)))==2)+1]
min.all
#-3.9220216 -3.4691894 -2.6506083 -2.3893590 -1.9974850 -1.7014024 -1.3704866 -1.0656958 -0.8218631
#-2.6506083 [3]

#infected mers for camel229E (fake mock)
mers_infected_camel229E_percentage <- mers_inf_filtered_feature_bc_matrix@meta.data %>% select(7)
mers_infected_camel229E_percentage$Percent.camel229E_log10 <- log10(mers_infected_camel229E_percentage$Percent.camel229E+halfmin_camel229E)
range(mers_infected_camel229E_percentage$Percent.camel229E_log10)

#mock for camel229E
mock_percentage <- mock_filtered_feature_bc_matrix@meta.data %>% select(7)
mock_percentage$Percent.camel229E_log10 <- log10(mock_percentage$Percent.camel229E+halfmin_camel229E)
range(mock_percentage$Percent.camel229E_log10)


#set colors for infected and non-infected
camel229E_infected_camel229E_percentage$group = ifelse(camel229E_infected_camel229E_percentage$Percent.camel229E_log10 <= min.all[3], 'Uninfected', ifelse(camel229E_infected_camel229E_percentage$Percent.camel229E_log10, 'Infected'))
mers_infected_camel229E_percentage$group = ifelse(mers_infected_camel229E_percentage$Percent.camel229E_log10 <= min.all[3], 'Uninfected', ifelse(mers_infected_camel229E_percentage$Percent.camel229E_log10, 'Infected'))
mock_percentage$group = ifelse(mock_percentage$Percent.camel229E_log10 <= min.all[3], 'Uninfected', ifelse(mock_percentage$Percent.camel229E_log10, 'Infected'))

p1 <- camel229E_infected_camel229E_percentage %>%
  ggplot(aes(x = Percent.camel229E_log10, fill = group))+
  geom_histogram(bins = 500)+
  scale_fill_manual(values = cols_stat) +
  xlim(-5,1)+
  coord_cartesian(ylim=c(0, 200))+
  ylab("Cell Number")+
  xlab("log10(% of mRNA of Camel229E-CoV)")+
  ggtitle("Camel-229E Infected")+
  theme_classic()+
  geom_vline(xintercept = min.all[3], linetype = "dashed", size = 0.9)+
  theme(axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15)))

p2 <- mers_infected_camel229E_percentage %>%
  ggplot(aes(x = Percent.camel229E_log10, fill = group)) +
  geom_histogram(bins = 500) +
  scale_fill_manual(values = cols_stat) +
  xlim(-5,1) +
  coord_cartesian(ylim=c(0, 200))+
  ylab("Cell Number")+
  xlab("")+
  ggtitle("MERS-CoV Infected")+
  theme_classic() +
  geom_vline(xintercept = min.all[3], linetype = "dashed", size = 0.9)+
  theme(axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15)))+
  NoLegend()

p3 <- mock_percentage %>%
  ggplot(aes(x = Percent.camel229E_log10, fill = group)) +
  geom_histogram(bins = 500) +
  scale_fill_manual(values = cols_stat) +
  xlim(-5,1) +
  coord_cartesian(ylim=c(0, 200))+
  ylab("Cell Number")+
  xlab("")+
  ggtitle("Mock")+
  theme_classic() +
  geom_vline(xintercept = min.all[3], linetype = "dashed", size = 0.9)+
  theme(axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 15, margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15)))+
  NoLegend()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/threshold_camel229E_infection.pdf",
    width=5, height=8)
p3 / p2 / p1
dev.off()

10^(-2.6506083)

#camel229E Infection threshold: 0.002235588


#Threshold above 0.002 Percent Viral is infected

mers_inf_filtered_feature_bc_matrix@meta.data <- mers_inf_filtered_feature_bc_matrix@meta.data %>% mutate(Status = case_when(Percent.MERS > 0.003237382 ~ "Infected",
                                                                                                                             Percent.MERS <= 0.003237382 ~ "Bystander"))
camel229E_inf_filtered_feature_bc_matrix@meta.data <- camel229E_inf_filtered_feature_bc_matrix@meta.data %>% mutate(Status = case_when(Percent.camel229E > 0.002235588 ~ "Infected",
                                                                                                                                       Percent.camel229E <= 0.002235588 ~ "Bystander"))

mock_filtered_feature_bc_matrix <- AddMetaData(mock_filtered_feature_bc_matrix, "Uninfected", col.name = "Status")

##################################################################################################################
# # QC
# # Make violin plot to show number of features, number of counts and mitochondrial content
# #mers infected
# VlnPlot(mers_inf_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# VlnPlot(mers_inf_filtered_feature_bc_matrix, features = c("Percent.MERS", "Total.MERS.Counts"), y.max = 0.8, ncol = 2)
# VlnPlot(mers_inf_filtered_feature_bc_matrix, features = c("Percent.camel229E", "Total.camel229E.Counts"), y.max = 0.8, ncol = 2)
# 
# #camel229E infected
# VlnPlot(camel229E_inf_filtered_feature_bc_matrix, features = c("nFeature_RNA"), ncol = 3) + NoLegend() 
# p2 <- VlnPlot(camel229E_inf_filtered_feature_bc_matrix, features = c("nCount_RNA"), ncol = 3) + NoLegend()
# p3 <- VlnPlot(camel229E_inf_filtered_feature_bc_matrix, features = c("percent.mt"), ncol = 3) + NoLegend()
# p1 + p2 + p3
# 
# VlnPlot(camel229E_inf_filtered_feature_bc_matrix, features = c("Percent.MERS", "Total.MERS.Counts"), ncol = 2)
# VlnPlot(camel229E_inf_filtered_feature_bc_matrix, features = c("Percent.camel229E", "Total.camel229E.Counts"), ncol = 2)
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
# plot1 <- FeatureScatter(mers_inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
# plot2 <- FeatureScatter(mers_inf_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
# plot3 <- FeatureScatter(mers_inf_filtered_feature_bc_matrix, feature1 = "Percent.MERS", feature2 = "percent.mt") + NoLegend()
# plot1 / plot2 /plot3
# #mock
# plot1 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(mock_filtered_feature_bc_matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# 


# Add number of genes per UMI for each cell to metadata
mers_inf_filtered_feature_bc_matrix$log10GenesPerUMI <- log10(mers_inf_filtered_feature_bc_matrix$nFeature_RNA) / log10(mers_inf_filtered_feature_bc_matrix$nCount_RNA)
camel229E_inf_filtered_feature_bc_matrix$log10GenesPerUMI <- log10(camel229E_inf_filtered_feature_bc_matrix$nFeature_RNA) / log10(camel229E_inf_filtered_feature_bc_matrix$nCount_RNA)
mock_filtered_feature_bc_matrix$log10GenesPerUMI <- log10(mock_filtered_feature_bc_matrix$nFeature_RNA) / log10(mock_filtered_feature_bc_matrix$nCount_RNA)



metadata <- bind_rows(mers_inf_filtered_feature_bc_matrix@meta.data, camel229E_inf_filtered_feature_bc_matrix@meta.data, mock_filtered_feature_bc_matrix@meta.data)
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


metadata_count <- metadata %>%
  group_by(Status, Treat) %>%
  count(Status)

metadata_count$Status <- as.factor(metadata_count$Status)
str(metadata_count)

levels(metadata_count$Status)

metadata_count$Status <- factor(x = metadata_count$Status, levels = c('Infected','Bystander', 'Uninfected'))

str(metadata_count$Status)
#Plots for QC
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/cell_numbers_per_virus_before_filtering.pdf",
    width=9, height=6)
metadata_count %>%
  ggplot(aes(x=Treat, y=n, fill=Status)) +
  geom_col(position = "stack") +
  theme_classic() +
  scale_fill_manual(values = cols_stat_2)+
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

?geom_col

#Number of UMI (transcripts) per cell
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/UMI_per_cell_density_log_before_filtering.pdf",
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
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(r = 20)),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 20, face = "bold"))
dev.off()

#Number of Genes per cell
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/genes_per_cell_density_log_before_filtering.pdf",
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/complexity_log_before_filtering.pdf",
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
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(t = 20)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 20, margin = margin(r = 20)),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 20, face = "bold"))
dev.off()


# Visualize the distribution of mitochondrial gene expression detected per cell
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/Mitochondrial_content_density_before_filtering.pdf",
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/correlation_gene_UMI_log_before_filtering.pdf",
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

# We filter cells that have unique feature counts over 1000
# We filter cells that have <30% mitochondrial counts
#infected
mers_inf_filtered_feature_bc_matrix <- subset(mers_inf_filtered_feature_bc_matrix, subset = nFeature_RNA > 1000 & percent.mt < 30)
camel229E_inf_filtered_feature_bc_matrix <- subset(camel229E_inf_filtered_feature_bc_matrix, subset = nFeature_RNA > 1000 & percent.mt < 30)
#mock
mock_filtered_feature_bc_matrix <- subset(mock_filtered_feature_bc_matrix, subset = nFeature_RNA > 1000 & percent.mt < 30)

##################################################################################################################

# #Plot to show Violin Plots of the infected cells
# MERS <- mers_inf_filtered_feature_bc_matrix
# Camel229E<- camel229E_inf_filtered_feature_bc_matrix
# Mock <- mock_filtered_feature_bc_matrix
# 
# Idents(MERS) <- "Treat"
# Idents(Camel229E) <- "Treat"
# Idents(Mock) <- "Treat"
# 
# 
# p1 <- VlnPlot(MERS, features = "Percent.MERS", y.max = 0.8, ncol = 1) + NoLegend() + FontSize(x.text = 0, x.title = 0)
# p2 <- VlnPlot(Camel229E, features = "Percent.camel229E", y.max = 0.8, ncol = 1) + NoLegend() + FontSize(x.text = 0, x.title = 0)
# p3 <- VlnPlot(Mock, features = "Percent.MERS", y.max = 0.8, ncol = 1) + NoLegend() + FontSize(x.text = 0, x.title = 0)
# 
# 
# p1 <- VlnPlot(MERS, features = "Percent.MERS", y.max = 0.8, ncol = 1) + NoLegend() + FontSize(x.text = 0, x.title = 0)
# p2 <- VlnPlot(Camel229E, features = "Percent.camel229E", y.max = 0.8, ncol = 1) + NoLegend() + FontSize(x.text = 0, x.title = 0)
# p3 <- VlnPlot(Mock, features = "Percent.MERS", y.max = 0.8, ncol = 1) + NoLegend() + FontSize(x.text = 0, x.title = 0)
# 
# grid.arrange(p1, p2, p3, ncol = 3)


###################################################################################################################
#cell cycle state genes from the Seurat package
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

g2m.genes <- c(g2m.genes, "PIMREG", "JPT1") #renamed genes have to be added to the list
s.genes <- c(s.genes, "CENPU") #renamed genes have to be added to the list


# a <- as.data.frame(mock_filtered_feature_bc_matrix@assays$RNA@counts@Dimnames[[1]])
# mito <- a[grep("^KEG", a$`mock_filtered_feature_bc_matrix@assays$RNA@counts@Dimnames[[1]]`), ]
# inf_filtered_feature_bc_matrix <- CellCycleScoring(inf_filtered_feature_bc_matrix, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# mock_filtered_feature_bc_matrix <- CellCycleScoring(mock_filtered_feature_bc_matrix, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

memory.limit(24000)
#################################################################################################################################
#Integrate the two sample matrices to one seurat object
Alpaca.list <- c(mock_filtered_feature_bc_matrix, mers_inf_filtered_feature_bc_matrix, camel229E_inf_filtered_feature_bc_matrix)
# normalize and identify variable features for each dataset independently
Alpaca.list <- lapply(X = Alpaca.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = c("Percent.MERS", "Percent.camel229E", "S.Score", "G2M.Score", "percent.mt"), verbose = FALSE)
})


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Alpaca.list, nfeatures = 3000)
Alpaca.list <- PrepSCTIntegration(object.list = Alpaca.list, anchor.features = features)
Alpaca.anchors <- FindIntegrationAnchors(object.list = Alpaca.list, anchor.features = features, normalization.method = "SCT")

# this command creates an 'integrated' data assay
Alpaca.combined <- IntegrateData(anchorset = Alpaca.anchors, normalization.method = "SCT")

###Perform an integrated analysis
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(Alpaca.combined) <- "integrated"

#head(Alpaca.combined@meta.data)
###################################################################################################################

#Saving the seurat object to disk to access in other script for analysis
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/V.pacos")
SaveH5Seurat(Alpaca.combined, "Alpaca.combined", overwrite = TRUE)

#########################################################################################################
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/V.pacos")
Alpaca.combined <- LoadH5Seurat("Alpaca.combined.h5seurat")
head(Alpaca.combined@meta.data)




# Add number of genes per UMI for each cell to metadata
Alpaca.combined$log10GenesPerUMI <- log10(Alpaca.combined$nFeature_RNA) / log10(Alpaca.combined$nCount_RNA)


metadata <- Alpaca.combined@meta.data
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/cell_numbers_per_virus_after_filtering.pdf",
    width=9, height=6)
metadata_count %>% 
  ggplot(aes(x=Treat, y=n, fill=Status)) + 
  geom_col() +
  theme_classic() +
  scale_fill_manual(values = cols_stat_2)+
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/UMI_per_cell_density_log_after_filtering.pdf",
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/genes_per_cell_density_log_after_filtering.pdf",
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/complexity_log_after_filtering.pdf",
    width=10, height=6)
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
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, face="bold", size = 17, margin = margin(r = 20)),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 20, face = "bold"))
dev.off()


# Visualize the distribution of mitochondrial gene expression detected per cell
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/Mitochondrial_content_density_after_filtering.pdf",
    width=10, height=6)
metadata %>% 
  ggplot(aes(color=Treat, x=percent.mt, fill=Treat)) + 
  geom_density(alpha = 0.6) + 
  theme_classic() +
  scale_color_manual(values= cols_treat) +
  scale_fill_manual(values = cols_treat) +
  xlim(0,50)+
  geom_vline(xintercept = 30)+
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/correlation_gene_UMI_log_after_filtering.pdf",
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
