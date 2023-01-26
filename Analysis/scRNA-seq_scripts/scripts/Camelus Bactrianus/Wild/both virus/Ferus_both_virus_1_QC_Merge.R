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
library(cowplot)

#Treatment colors(camel-229E, MERS, Mock)
cols_treat = c("#F39C6B","#A5668B", "#96ADC8")
show_col(cols_treat) #show the colors
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
  seurat_data <- Read10X(data.dir = paste0("data/C. Ferus/both_virus/", file))
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
# Checking how many genes per cell each samples has. The genes per cell values are calculated for each cell and displayed ordered.

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

#plotting how many genes per cell we have in the samples
plot(sort(genes_per_cell_inf_mers), ylim = c(1, 10000), xlim = c(0, 5000), xlab='cell', log='y', main='genes per cell (ordered) - Ferus MERS infected')
plot(sort(genes_per_cell_inf_camel229E), ylim = c(1, 10000), xlim = c(0, 5000), xlab='cell', log='y', main='genes per cell (ordered) - Ferus camel229E infected')

plot(sort(genes_per_cell_mock), ylim = c(1, 10000), xlim = c(0, 5000), xlab='cell', log='y', main='genes per cell (ordered) - Ferus Mock')

#####################################################################################
# Adding important info (Viral counts for both viruses, viral percentage for each virus, mitochondrial content and grouping variables) to the meta data column
# show the number of rows of the matrix to know how many genes we have
nrow(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts)
nrow(camel229E_inf_filtered_feature_bc_matrix@assays$RNA@counts)
nrow(mock_filtered_feature_bc_matrix@assays$RNA@counts)

# display the last 25 rows in the counts per gene (from above) to see the total viral counts in all the cells
tail(counts_per_gene_inf_mers, 20)
tail(counts_per_gene_inf_camel229E, 25)
tail(counts_per_gene_mock, 25)


#show the specific location in the matrix for control
mers_inf_filtered_feature_bc_matrix@assays$RNA@counts[30848:30859,]
camel229E_inf_filtered_feature_bc_matrix@assays$RNA@counts[30860:30869,]
mock_filtered_feature_bc_matrix@assays$RNA@counts[30848:30859,] #mers mock
mock_filtered_feature_bc_matrix@assays$RNA@counts[30860:30869,] #camel229E mock

#take the right rows in the matrix
mers_infected_mers_counts <- as.matrix(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts[30848:30859,])
mers_infected_camel229E_counts <- as.matrix(mers_inf_filtered_feature_bc_matrix@assays$RNA@counts[30860:30869,])

camel229E_infected_mers_counts <- as.matrix(camel229E_inf_filtered_feature_bc_matrix@assays$RNA@counts[30848:30859,])
camel229E_infected_camel229E_counts <- as.matrix(camel229E_inf_filtered_feature_bc_matrix@assays$RNA@counts[30860:30869,])

mock_mers_counts <- as.matrix(mock_filtered_feature_bc_matrix@assays$RNA@counts[30848:30859,])
mock_camel229E_counts <- as.matrix(mock_filtered_feature_bc_matrix@assays$RNA@counts[30860:30869,])

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
#Adding a Meta data column with information about the Status of each cell. If the cell is Infected or Bystander
#in the MERS infected Treatment group or non-infected in the Mock Group
mers_inf_filtered_feature_bc_matrix <- AddMetaData(object = mers_inf_filtered_feature_bc_matrix, metadata = mers.infected.total.mers.counts, col.name = "Total.MERS.Counts")
mers_inf_filtered_feature_bc_matrix <- AddMetaData(object = mers_inf_filtered_feature_bc_matrix, metadata = mers.infected.total.camel229E.counts, col.name = "Total.camel229E.Counts")
mers_inf_filtered_feature_bc_matrix <- AddMetaData(object = mers_inf_filtered_feature_bc_matrix, metadata = mers.infected.percent.mers, col.name = "Percent.MERS")
mers_inf_filtered_feature_bc_matrix <- AddMetaData(object = mers_inf_filtered_feature_bc_matrix, metadata = mers.infected.percent.camel229E, col.name = "Percent.camel229E")
mers_inf_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(mers_inf_filtered_feature_bc_matrix, pattern = "^KEG")
mers_inf_filtered_feature_bc_matrix <- AddMetaData(mers_inf_filtered_feature_bc_matrix, "MERS-CoV", col.name = "Treat")

head(mers_inf_filtered_feature_bc_matrix@meta.data)

camel229E_inf_filtered_feature_bc_matrix <- AddMetaData(object = camel229E_inf_filtered_feature_bc_matrix, metadata = camel229E.infected.total.mers.counts, col.name = "Total.MERS.Counts")
camel229E_inf_filtered_feature_bc_matrix <- AddMetaData(object = camel229E_inf_filtered_feature_bc_matrix, metadata = camel229E.infected.total.camel229E.counts, col.name = "Total.camel229E.Counts")
camel229E_inf_filtered_feature_bc_matrix <- AddMetaData(object = camel229E_inf_filtered_feature_bc_matrix, metadata = camel229E.infected.percent.mers, col.name = "Percent.MERS")
camel229E_inf_filtered_feature_bc_matrix <- AddMetaData(object = camel229E_inf_filtered_feature_bc_matrix, metadata = camel229E.infected.percent.camel229E, col.name = "Percent.camel229E")
camel229E_inf_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(camel229E_inf_filtered_feature_bc_matrix, pattern = "^KEG")
camel229E_inf_filtered_feature_bc_matrix <- AddMetaData(camel229E_inf_filtered_feature_bc_matrix, "dcCoV-ACN4", col.name = "Treat")

head(camel229E_inf_filtered_feature_bc_matrix@meta.data)

mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.total.mers.counts, col.name = "Total.MERS.Counts")
mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.total.camel229E.counts, col.name = "Total.camel229E.Counts")
mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.percent.mers, col.name = "Percent.MERS")
mock_filtered_feature_bc_matrix <- AddMetaData(object = mock_filtered_feature_bc_matrix, metadata = mock.percent.camel229E, col.name = "Percent.camel229E")
mock_filtered_feature_bc_matrix[["percent.mt"]] <- PercentageFeatureSet(mock_filtered_feature_bc_matrix, pattern = "^KEG")
mock_filtered_feature_bc_matrix <- AddMetaData(mock_filtered_feature_bc_matrix, "Mock", col.name = "Treat")
head(mock_filtered_feature_bc_matrix@meta.data)


############################################################################################################################################################
# Visualize "infection" in the three species samples

# #Extracting the Viral-Percentage column for every cell and transpose it to plot it below
# mers_infected_mers_percentage <- mers_inf_filtered_feature_bc_matrix@meta.data %>% select(6) %>% t()
# mers_infected_camel229E_percentage <- mers_inf_filtered_feature_bc_matrix@meta.data %>% select(7) %>% t()
# camel229E_infected_mers_percentage <- camel229E_inf_filtered_feature_bc_matrix@meta.data %>% select(6) %>% t()
# camel229E_infected_camel229E_percentage <- camel229E_inf_filtered_feature_bc_matrix@meta.data %>% select(7) %>% t()
# Mock_mers_percentage <- mock_filtered_feature_bc_matrix@meta.data %>% select(6) %>% t()
# Mock_camel229E_percentage <- mock_filtered_feature_bc_matrix@meta.data %>% select(7) %>% t()
# 
# #plotting how many virus positive cells there are in each treatment
# plot(sort(mers_infected_mers_percentage), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Ferus MERS infected MERS counts')
# plot(sort(mers_infected_camel229E_percentage), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Ferus MERS infected camel229E counts')
# plot(sort(camel229E_infected_mers_percentage), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Ferus camel229E infected MERS counts')
# plot(sort(camel229E_infected_camel229E_percentage), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Ferus camel229E infected camel229E counts')
# plot(sort(Mock_mers_percentage), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Ferus Mock MERS counts')
# plot(sort(Mock_camel229E_percentage), xlim = c(0, 5000), xlab='cell', main='Viral.Percentage vs cell (ordered) - Ferus Mock camel229E counts')
# 
# head(mers_inf_filtered_feature_bc_matrix@meta.data)

#############################################################################################################################
#Calculation and visualization of infection threshold for the camel samples
#infected mers for mers
mers_infected_mers_percentage <- mers_inf_filtered_feature_bc_matrix@meta.data %>% select(6)

#calculate minimal number after 0 to set instead of 0, so they are included into plot
halfmin_mers <- min(mers_infected_mers_percentage$Percent.MERS[mers_infected_mers_percentage$Percent.MERS>0])
halfmin_mers

#Add the minimal MERS percentage to all the values and log10 transform all the values
mers_infected_mers_percentage$Percent.MERS_log10 <- log10(mers_infected_mers_percentage$Percent.MERS+halfmin_mers)
#give out the range of log10 transformed values
range(mers_infected_mers_percentage$Percent.MERS_log10)

#calculate the density distribution of all the transformed values
des.all <- density(mers_infected_mers_percentage$Percent.MERS_log10)
#take out all the local minimas
min.all <- des.all$x[which(diff(sign(diff(des.all$y)))==2)+1]
min.all # choose the right minimum for this virus and this sample
#-2.3159777

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


#Plot the densities with colored set to "infected" or "not-infected" and a line at the threshold
p1 <- mers_infected_mers_percentage %>%
  ggplot(aes(x = Percent.MERS_log10, fill = group))+
  geom_histogram(bins = 500)+
  scale_fill_manual(values = c("#CD4631","#DEA47E")) +
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

p2 <- camel229E_infected_mers_percentage %>%
  ggplot(aes(x = Percent.MERS_log10, fill = group)) +
  geom_histogram(bins = 500) +
  scale_fill_manual(values = c("#DEA47E", "#CD4631")) +
  xlim(-5,1) +
  coord_cartesian(ylim=c(0, 60))+
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
  scale_fill_manual(values = c("#DEA47E", "#CD4631")) +
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

#save the plots
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/threshold_mers_infection.pdf",
    width=5, height=8)
p3 / p2 / p1
dev.off()

#the infection threshold calculated back from log10 is
10^(-2.3159777)

###################################################################################################################
#infected camel229E for camel229E
camel229E_infected_camel229E_percentage <- camel229E_inf_filtered_feature_bc_matrix@meta.data %>% select(7)

#calculate minimal number after 0 to set instead of 0, so they are included into plot
halfmin_camel229E <- min(camel229E_infected_camel229E_percentage$Percent.camel229E[camel229E_infected_camel229E_percentage$Percent.camel229E>0])
halfmin_camel229E

camel229E_infected_camel229E_percentage$Percent.camel229E_log10 <- log10(camel229E_infected_camel229E_percentage$Percent.camel229E+halfmin_camel229E)
range(camel229E_infected_camel229E_percentage$Percent.camel229E_log10)

des.all <- density(camel229E_infected_camel229E_percentage$Percent.camel229E_log10)
min.all <- des.all$x[which(diff(sign(diff(des.all$y)))==2)+1]
min.all
#-2.243342

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
  scale_fill_manual(values = c("#CD4631","#DEA47E")) +
  xlim(-5,1)+
  coord_cartesian(ylim=c(0, 60))+
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
  scale_fill_manual(values = c("#DEA47E", "#CD4631")) +
  xlim(-5,1) +
  coord_cartesian(ylim=c(0, 60))+
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
  scale_fill_manual(values = c("#DEA47E","#CD4631")) +
  xlim(-5,1) +
  coord_cartesian(ylim=c(0, 60))+
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

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/threshold_camel229E_infection.pdf",
    width=5, height=8)
p3 / p2 / p1
dev.off()

10^(-2.243342)
#Infection threshold: 0.006483851

#Threshold above 0.004 Percent Viral is infected

mers_inf_filtered_feature_bc_matrix@meta.data <- mers_inf_filtered_feature_bc_matrix@meta.data %>% mutate(Status = case_when(Percent.MERS > 0.004830836 ~ "Infected",
                                                                                                                             Percent.MERS <= 0.004830836 ~ "Bystander"))
camel229E_inf_filtered_feature_bc_matrix@meta.data <- camel229E_inf_filtered_feature_bc_matrix@meta.data %>% mutate(Status = case_when(Percent.camel229E > 0.005710288 ~ "Infected",
                                                                                                                                       Percent.camel229E <= 0.005710288 ~ "Bystander"))

mock_filtered_feature_bc_matrix <- AddMetaData(mock_filtered_feature_bc_matrix, "Uninfected", col.name = "Status")

##################################################################################################################
# QC
# Make violin plot to show number of features, number of counts and mitochondrial content
#mers infected
VlnPlot(mers_inf_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(mers_inf_filtered_feature_bc_matrix, features = c("Percent.MERS", "Total.MERS.Counts"), y.max = 0.8, ncol = 2)
VlnPlot(mers_inf_filtered_feature_bc_matrix, features = c("Percent.camel229E", "Total.camel229E.Counts"), y.max = 0.8, ncol = 2)

#camel229E infected
VlnPlot(camel229E_inf_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(camel229E_inf_filtered_feature_bc_matrix, features = c("Percent.MERS", "Total.MERS.Counts"), ncol = 2)
VlnPlot(camel229E_inf_filtered_feature_bc_matrix, features = c("Percent.camel229E", "Total.camel229E.Counts"), ncol = 2)
#mock
VlnPlot(mock_filtered_feature_bc_matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(mock_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), ncol = 2)
VlnPlot(mock_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), ncol = 2)
VlnPlot(mock_filtered_feature_bc_matrix, features = c("Percent.Viral", "Total.Viral.Counts"), y.max = 0.8, ncol = 2)

?VlnPlot()
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


# Add number of genes per UMI for each cell to metadata
mers_inf_filtered_feature_bc_matrix$log10GenesPerUMI <- log10(mers_inf_filtered_feature_bc_matrix$nFeature_RNA) / log10(mers_inf_filtered_feature_bc_matrix$nCount_RNA)
camel229E_inf_filtered_feature_bc_matrix$log10GenesPerUMI <- log10(camel229E_inf_filtered_feature_bc_matrix$nFeature_RNA) / log10(camel229E_inf_filtered_feature_bc_matrix$nCount_RNA)
mock_filtered_feature_bc_matrix$log10GenesPerUMI <- log10(mock_filtered_feature_bc_matrix$nFeature_RNA) / log10(mock_filtered_feature_bc_matrix$nCount_RNA)

mers_inf_filtered_feature_bc_matrix@meta.data

metadata <- bind_rows(mers_inf_filtered_feature_bc_matrix@meta.data, camel229E_inf_filtered_feature_bc_matrix@meta.data, mock_filtered_feature_bc_matrix@meta.data)
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

#cellnumber for each status and treatment
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/cell_numbers_per_virus_before_filtering.pdf",
    width=6, height=5)
metadata_count %>% 
  ggplot(aes(x=Treat, y=n, fill=Status)) + 
  geom_col() +
  theme_classic() +
  scale_fill_manual(values = cols_stat_2)+
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMI_per_cell_density_log_before_filtering.pdf",
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

#Number of Genes per cell
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/genes_per_cell_density_log_before_filtering.pdf",
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/complexity_log_before_filtering.pdf",
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/Mitochondrial_content_density_before_filtering.pdf",
    width=9, height=6)
metadata %>% 
  ggplot(aes(color=Treat, x=percent.mt, fill=Treat)) + 
  geom_density(alpha = 0.6) + 
  scale_color_manual(values= cols_treat) +
  scale_fill_manual(values = cols_treat) +
  xlim(0,50)+
  theme_classic() +
  geom_vline(xintercept = 30) +
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/correlation_gene_UMI_log_before_filtering.pdf",
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
# p1 <- VlnPlot(MERS, features = "Percent.Viral", y.max = 0.8, ncol = 1) + NoLegend() + FontSize(x.text = 0, x.title = 0)
# p2 <- VlnPlot(Camel229E, features = "Percent.Viral", y.max = 0.8, ncol = 1) + NoLegend() + FontSize(x.text = 0, x.title = 0)
# p3 <-VlnPlot(Mock, features = "Percent.Viral", y.max = 0.8, ncol = 1) + NoLegend() + FontSize(x.text = 0, x.title = 0)
# 
# grid.arrange(p1, p2, p3, ncol = 3)


###################################################################################################################
#cell cycle state genes from the Seurat package
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

g2m.genes <- c(g2m.genes, "PIMREG", "JPT1") #renamed genes have to be added to the list
s.genes <- c(s.genes, "CENPU") #renamed genes have to be added to the list

memory.limit(24000) #seurat size cannot be allocated anymore, memory limit at 24000

#################################################################################################################################
#Integrate the two sample matrices to one seurat object
Ferus.list <- c(mock_filtered_feature_bc_matrix, mers_inf_filtered_feature_bc_matrix, camel229E_inf_filtered_feature_bc_matrix)
# normalize and identify variable features for each dataset independently
Ferus.list <- lapply(X = Ferus.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  x <- SCTransform(x, vars.to.regress = c("Percent.MERS", "Percent.camel229E", "S.Score", "G2M.Score", "percent.mt"), verbose = FALSE)
})



# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Ferus.list, nfeatures = 3000)
Ferus.list <- PrepSCTIntegration(object.list = Ferus.list, anchor.features = features)
Ferus.anchors <- FindIntegrationAnchors(object.list = Ferus.list, anchor.features = features, normalization.method = "SCT")

# this command creates an 'integrated' data assay with the SCT normalized data
Ferus.combined <- IntegrateData(anchorset = Ferus.anchors, normalization.method = "SCT")

##Perform an integrated analysis
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(Ferus.combined) <- "integrated"

head(Ferus.combined@meta.data)
###################################################################################################################

#Saving the seurat object to disk to access in other script for analysis
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.ferus")
SaveH5Seurat(Ferus.combined, "Ferus.combined", overwrite = TRUE)


#########################################################################################################
#Plot the QC plots after filtering

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.ferus")
Ferus.combined <- LoadH5Seurat("Ferus.combined.h5seurat")
Ferus.combined


# Add number of genes per UMI for each cell to metadata
Ferus.combined$log10GenesPerUMI <- log10(Ferus.combined$nFeature_RNA) / log10(Ferus.combined$nCount_RNA)


metadata <- Ferus.combined@meta.data
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
#cell number for each status and treatment
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/cell_numbers_per_virus_after_filtering.pdf",
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMI_per_cell_density_log_after_filtering.pdf",
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/genes_per_cell_density_log_after_filtering.pdf",
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/complexity_log_after_filtering.pdf",
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/Mitochondrial_content_density_after_filtering.pdf",
    width=9, height=6)
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
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/correlation_gene_UMI_log_after_filtering.pdf",
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



