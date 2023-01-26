# Differential expression of Single-cell Seq data
# 12.08.2022

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
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(EnhancedVolcano)
library(ggrepel)


#Treatment colors (HCoV-229E, MERS-CoV, Mock)
cols_treat = c("#84A59D","#A5668B", "#96ADC8")
show_col(cols_treat)
#Status colors
cols_stat = c("#CD4631","#DEA47E")
show_col(cols_stat)


##########################################################################################################################
#Loading filtered and merged Seurat Object into Script
#Either load just combined or also clustered data set

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens")

Human <- LoadH5Seurat("Human.clustered.h5seurat") # for bulk resolution analysis
#Human <- LoadH5Seurat("Human.combined.h5seurat") # for single cell resolution analysis
head(Human@meta.data)
summary(Human@meta.data)


##################################################################################################################
#Plot relationship between stimulated and unstimulated before to see outliers
Idents(Human) <- "Treat" #Parameters "Treat" needs to be marked as a Idents in the object
levels(Human)
# avg.Human <- log1p(AverageExpression(Human, verbose = FALSE)$RNA)
# avg.Human <- as.data.frame(avg.Human)
# 
# 
# #Plot genes which have a different average expression in Mock as in treated
# p1 <- ggplot(avg.Human, aes(Mock, `MERS-CoV`)) + geom_point() + ggtitle("CD4 Naive T Cells")
# p2 <- ggplot(avg.Human, aes(Mock, `HCoV-229E`)) + geom_point() + ggtitle("CD4 Naive T Cells")
# p1 | p2


##################################################################################################################################################

# Difference between Treatment (Mock, HCoV-229E, MERS-CoV) for all cells

##################################################################################################################################################

# "Bulk" differential expression in the Treatment groups Mock, HcoV-229E, MERS-CoV
DefaultAssay(Human) <- "RNA"
# Find differentially expressed features between Mock and treated with HCoV-229E
Mock_vs_HCoV.markers_wilcox <- FindMarkers(Human, slot = "data", ident.2 = "Mock", ident.1 = "HCoV-229E", test.use = "wilcox")
# view results
#head(Mock_vs_HCoV.markers)

# Find differentially expressed features between Mock and treated with MERS-CoV
Mock_vs_MERS.markers_wilcox <- FindMarkers(Human, slot = "data", ident.2 = "Mock", ident.1 = "MERS-CoV", test.use = "wilcox")
# view results
#head(Mock_vs_MERS.markers)

#Volcano#####################################################################################################################
x1 <- FindMarkers(Human, ident.2 = "Mock", ident.1 = "HCoV-229E", test.use = "wilcox")
p1 <- EnhancedVolcano(x1,
                      title = 'Mock versus HCoV-229E treated',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x1),
                      col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 0.5,
                      pCutoff = 0.05,
                      colAlpha = 1)

x2 <- FindMarkers(Human, ident.2 = "Mock", ident.1 = "MERS-CoV", test.use = "wilcox")
p2 <- EnhancedVolcano(x2,
                      title = 'Mock versus MERS-CoV treated',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x2),
                      col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 0.5,
                      pCutoff = 0.05,
                      colAlpha = 1)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/DEG_Mock_v_Virus.pdf",
    width=15, height=7)
p1|p2
dev.off()




#save the lists
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/DEG/Human_both")
write.csv(Mock_vs_HCoV.markers_wilcox,"Human_Mock_vs_HCoV_markers_wilcox.csv")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/DEG/Human_both")
write.csv(Mock_vs_MERS.markers_wilcox,"Human_Mock_vs_MERS_markers_wilcox.csv")


#############################################################################################################################################
# Visualizations
Marker_Plots <- function(a, b) {
  a_1 <- rownames(a)
  a_2 <- a_1[1:35]
  # p_1 <- DoHeatmap(b, features = a_2, group.by = "Treat")
  # print(p_1)
  return(a_2)
}

Marker_list_Human <- Marker_Plots(Human_inf.markers, Human)
Marker_list_Mock_vs_MERS <- Marker_Plots(Mock_vs_MERS.markers, Human)
Marker_list_Mock_vs_HCoV <- Marker_Plots(Mock_vs_HCoV.markers, Human)



##########################################################################################################################
#GO Annotation

#GO_results_Human <- enrichGO(gene = Marker_list_Human, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO_results_Mock_vs_MERS <- enrichGO(gene = Marker_list_Mock_vs_MERS, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO_results_Mock_vs_HCoV <- enrichGO(gene = Marker_list_Mock_vs_HCoV, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
#GO_results_MERS_vs_HCoV <- enrichGO(gene = Marker_list_MERS_vs_HCoV, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
#GO_results_Mock <- enrichGO(gene = Marker_list_Mock, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")

#as.data.frame(GO_results_Human)
as.data.frame(GO_results_Mock_vs_MERS)
as.data.frame(GO_results_Mock_vs_HCoV)
#as.data.frame(GO_results_MERS_vs_HCoV)
#as.data.frame(GO_results_Mock)

#plot(barplot(GO_results_Human, showCategory = 20))
plot(barplot(GO_results_Mock_vs_MERS, showCategory = 20))
plot(barplot(GO_results_Mock_vs_HCoV, showCategory = 20))
#plot(barplot(GO_results_MERS_vs_HCoV, showCategory = 20))
#plot(barplot(GO_results_Mock, showCategory = 20))


#Function for later
GO_Annotation <- function(a) {
  a_1 <- enrichGO(gene = a, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
  a_1 <- as.data.frame(a_1)
  return(a_1)
}


EnhancedVolcano(Human_inf.markers, x = "avg_log2FC", y = "p_val_adj", lab = rownames(Human_inf.markers))
EnhancedVolcano(Mock_vs_MERS.markers, x = "avg_log2FC", y = "p_val_adj", lab = rownames(Mock_vs_MERS.markers))
EnhancedVolcano(Mock_vs_HCoV.markers, x = "avg_log2FC", y = "p_val_adj", lab = rownames(Mock_vs_HCoV.markers))
EnhancedVolcano(MERS_vs_HCoV.markers, x = "avg_log2FC", y = "p_val_adj", lab = rownames(MERS_vs_HCoV.markers))
EnhancedVolcano(Mock.markers, x = "avg_log2FC", y = "p_val_adj", lab = rownames(Mock.markers))



#############################################################################################################################################
#Violin plots for interesting markers found above

VlnPlot(Human, features = c("GPRC5A", "PLK2", "CEACAM6", "S100A16", "VCL"), group.by = "Treat")
?VlnPlot()






##################################################################################################################################################

# Difference between Status (Infected, Uninfected, Bystander) for each Virus (HcoV-229E, MERS-CoV)

##################################################################################################################################################
#Increase memory limit
memory.limit(24000)
# Find differentially expressed features between Infected and uninfected for one Treatment group
# MERS-CoV
Human@meta.data
Idents(Human) <- "Treat"

#Make a fused new identity of each treatment and status combination
Human$treat.status <- paste(Idents(Human), Human$Status, sep = "_")

#make it the new identity
Idents(Human) <- "treat.status"
#check if levels are correct
levels(Human)


#Find markers for all combinations of Treatment/Status combinations that make sense

Uninf_vs_MERS_bys.markers <- FindMarkers(Human, ident.1 = "Mock_Uninfected", ident.2 = "MERS-CoV_Bystander", min.pct = 0.5, logfc.threshold = 0.5)
# view results
head(Uninf_vs_MERS_bys.markers)

Uninf_vs_MERS_inf.markers <- FindMarkers(Human, ident.1 = "Mock_Uninfected", ident.2 = "MERS-CoV_Infected", min.pct = 0.5, logfc.threshold = 0.5)
# view results
head(Uninf_vs_MERS_inf.markers)

Uninf_vs_HCOV_bys.markers <- FindMarkers(Human, ident.1 = "Mock_Uninfected", ident.2 = "HCoV-229E_Bystander", min.pct = 0.5, logfc.threshold = 0.5)
# view results
head(Uninf_vs_HCOV_bys.markers)

#not enough cells in the HCOV infected group
Uninf_vs_HCOV_inf.markers <- FindMarkers(Human, ident.1 = "Mock_Uninfected", ident.2 = "HCoV-229E_Infected", min.pct = 0.5, logfc.threshold = 0.5)
# view results
head(Uninf_vs_HCOV_inf.markers)

##################################################################################################################################################

Marker_Plots <- function(a, b) {
  a_1 <- rownames(a)
  a_2 <- a_1[1:35]
  p_1 <- DoHeatmap(b, features = a_2, group.by = "treat.status")
  print(p_1)
  return(a_2)
}

Marker_list_Uninf_vs_MERS_bys <- Marker_Plots(Uninf_vs_MERS_bys.markers, Human)
Marker_list_Uninf_vs_MERS_inf <- Marker_Plots(Uninf_vs_MERS_inf.markers, Human)
Marker_list_Uninf_vs_HCOV_bys <- Marker_Plots(Uninf_vs_HCOV_bys.markers, Human)



##########################################################################################################################
#GO Annotation

#GO_results_Human <- enrichGO(gene = Marker_list_Human, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO_results_Uninf_vs_MERS_bys <- enrichGO(gene = Marker_list_Uninf_vs_MERS_bys, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO_results_Uninf_vs_MERS_inf <- enrichGO(gene = Marker_list_Uninf_vs_MERS_inf, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO_results_Uninf_vs_HCOV_bys <- enrichGO(gene = Marker_list_Uninf_vs_HCOV_bys, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")


#as.data.frame(GO_results_Human)
as.data.frame(GO_results_Uninf_vs_MERS_bys)
as.data.frame(GO_results_Uninf_vs_MERS_inf)
as.data.frame(GO_results_Uninf_vs_HCOV_bys)


plot(barplot(GO_results_Uninf_vs_MERS_bys, showCategory = 20))
plot(barplot(GO_results_Uninf_vs_MERS_inf, showCategory = 20))
plot(barplot(GO_results_Uninf_vs_HCOV_bys, showCategory = 20))



# #Function for later
# GO_Annotation <- function(a) {
#   a_1 <- enrichGO(gene = a, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
#   a_1 <- as.data.frame(a_1)
#   return(a_1)
# }

#Volcano plots
EnhancedVolcano(Uninf_vs_MERS_bys.markers, x = "avg_log2FC", y = "p_val_adj", lab = rownames(Uninf_vs_MERS_bys.markers))
EnhancedVolcano(Uninf_vs_MERS_inf.markers, x = "avg_log2FC", y = "p_val_adj", lab = rownames(Uninf_vs_MERS_inf.markers))
EnhancedVolcano(Uninf_vs_HCOV_bys.markers, x = "avg_log2FC", y = "p_val_adj", lab = rownames(Uninf_vs_HCOV_bys.markers))



########################################################################################################################################################

# Difference between Status (Infected, uninfected, Bystander) for each celltype and virus

########################################################################################################################################################
DefaultAssay(Human) <- "integrated"
#Making new metadata groups
Human@meta.data
Idents(Human) <- "celltype"
Human$celltype.status.treat <- paste(Idents(Human), Human$Status, Human$Treat, sep = "_")
#Human$celltype.treat <- paste(Idents(Human), Human$Treat, sep = "_")


#Analyzing status groups over celltypes
Idents(Human) <- "celltype.status.treat"
levels(Human)

##############
#Secretory
#MERS-COV
#Uninfected vs Infected
# secretory.response_Uninf_vs_Inf_MERS <- FindMarkers(Human, ident.1 = "Secretory_Uninfected_Mock", ident.2 = "Secretory_Infected_MERS-CoV", verbose = FALSE, min.pct = 0.5, logfc.threshold = 0.5)
# head(secretory.response_Uninf_vs_Inf_MERS, n = 15)
#Uninfected vs Bystander
secretory.response_Uninf_vs_Bys_MERS <- FindMarkers(Human, ident.1 = "Secretory_Uninfected_Mock", ident.2 = "Secretory_Bystander_MERS-CoV", verbose = FALSE, min.pct = 0.5, logfc.threshold = 1)
head(secretory.response_Uninf_vs_Bys_MERS, n = 15)

?FindMarkers()
#HCoV-229E
#Uninfected vs Infected
# secretory.response_Uninf_vs_Inf_HCOV <- FindMarkers(Human, ident.1 = "Secretory_Uninfected_Mock", ident.2 = "Secretory_Infected_HCoV-229E", verbose = FALSE, min.pct = 0.5, logfc.threshold = 0.5)
# head(secretory.response_Uninf_vs_Inf_HCOV, n = 15)
#Uninfected vs Bystander
secretory.response_Uninf_vs_Bys_HCOV <- FindMarkers(Human, ident.1 = "Secretory_Uninfected_Mock", ident.2 = "Secretory_Bystander_HCoV-229E", verbose = FALSE, min.pct = 0.5, logfc.threshold = 1)
head(secretory.response_Uninf_vs_Bys_HCOV, n = 15)


##############
#Ciliated
#MERS-COV
#Uninfected vs Infected
# ciliated.response_Uninf_vs_Inf_MERS <- FindMarkers(Human, ident.1 = "Ciliated_Uninfected_Mock", ident.2 = "Ciliated_Infected_MERS-CoV", verbose = FALSE, min.pct = 0.5, logfc.threshold = 0.5)
# head(ciliated.response_Uninf_vs_Inf_MERS, n = 15)
#Uninfected vs Bystander
ciliated.response_Uninf_vs_Bys_MERS <- FindMarkers(Human, ident.1 = "Ciliated_Uninfected_Mock", ident.2 = "Ciliated_Bystander_MERS-CoV", verbose = FALSE, min.pct = 0.5, logfc.threshold = 0.5)
head(ciliated.response_Uninf_vs_Bys_MERS, n = 15)

#HCoV-229E
#Uninfected vs Infected
# ciliated.response_Uninf_vs_Inf_HCOV <- FindMarkers(Human, ident.1 = "Ciliated_Uninfected_Mock", ident.2 = "Ciliated_Infected_HCoV-229E", verbose = FALSE, min.pct = 0.5, logfc.threshold = 0.5)
# head(ciliated.response_Uninf_vs_Inf_HCOV, n = 15)
#Uninfected vs Bystander
ciliated.response_Uninf_vs_Bys_HCOV <- FindMarkers(Human, ident.1 = "Ciliated_Uninfected_Mock", ident.2 = "Ciliated_Bystander_HCoV-229E", verbose = FALSE, min.pct = 0.5, logfc.threshold = 0.5)
head(ciliated.response_Uninf_vs_Bys_HCOV, n = 15)


##############
#Club
#MERS-COV
#Uninfected vs Infected
# club.response_Uninf_vs_Inf_MERS <- FindMarkers(Human, ident.1 = "Club_Uninfected_Mock", ident.2 = "Club_Infected_MERS-CoV", verbose = FALSE)
# head(club.response_Uninf_vs_Inf_MERS, n = 15)
#Uninfected vs Bystander
club.response_Uninf_vs_Bys_MERS <- FindMarkers(Human, ident.1 = "Club_Uninfected_Mock", ident.2 = "Club_Bystander_MERS-CoV", verbose = FALSE)
head(club.response_Uninf_vs_Bys_MERS, n = 15)

#HCoV-229E
#Uninfected vs Infected
# club.response_Uninf_vs_Inf_HCOV <- FindMarkers(Human, ident.1 = "Club_Uninfected_Mock", ident.2 = "Club_Infected_HCoV-229E", verbose = FALSE)
# head(club.response_Uninf_vs_Inf_HCOV, n = 15)
#Uninfected vs Bystander
club.response_Uninf_vs_Bys_HCOV <- FindMarkers(Human, ident.1 = "Club_Uninfected_Mock", ident.2 = "Club_Bystander_HCoV-229E", verbose = FALSE)
head(club.response_Uninf_vs_Bys_HCOV, n = 15)



##############
#Deuterosomal
#MERS-COV
#Uninfected vs Infected
# deuterosomal.response_Uninf_vs_Inf_MERS <- FindMarkers(Human, ident.1 = "Deuterosomal_Uninfected_Mock", ident.2 = "Deuterosomal_Infected_MERS-CoV", verbose = FALSE)
# head(deuterosomal.response_Uninf_vs_Inf_MERS, n = 15)
#Uninfected vs Bystander
deuterosomal.response_Uninf_vs_Bys_MERS <- FindMarkers(Human, ident.1 = "Deuterosomal_Uninfected_Mock", ident.2 = "Deuterosomal_Bystander_MERS-CoV", verbose = FALSE)
head(deuterosomal.response_Uninf_vs_Bys_MERS, n = 15)

#HCoV-229E
#Uninfected vs Infected
# deuterosomal.response_Uninf_vs_Inf_HCOV <- FindMarkers(Human, ident.1 = "Deuterosomal_Uninfected_Mock", ident.2 = "Deuterosomal_Infected_HCoV-229E", verbose = FALSE)
# head(deuterosomal.response_Uninf_vs_Inf_HCOV, n = 15)
#Uninfected vs Bystander
deuterosomal.response_Uninf_vs_Bys_HCOV <- FindMarkers(Human, ident.1 = "Deuterosomal_Uninfected_Mock", ident.2 = "Deuterosomal_Bystander_HCoV-229E", verbose = FALSE)
head(deuterosomal.response_Uninf_vs_Bys_HCOV, n = 15)


##############
#Basal
#MERS-COV
#Uninfected vs Infected
# basal.response_Uninf_vs_Inf_MERS <- FindMarkers(Human, ident.1 = "Basal_Uninfected_Mock", ident.2 = "Basal_Infected_MERS-CoV", verbose = FALSE)
# head(basal.response_Uninf_vs_Inf_MERS, n = 15)
#Uninfected vs Bystander
basal.response_Uninf_vs_Bys_MERS <- FindMarkers(Human, ident.1 = "Basal_Uninfected_Mock", ident.2 = "Basal_Bystander_MERS-CoV", verbose = FALSE)
head(basal.response_Uninf_vs_Bys_MERS, n = 15)

#HCoV-229E
#Uninfected vs Infected
# basal.response_Uninf_vs_Inf_HCOV <- FindMarkers(Human, ident.1 = "Basal_Uninfected_Mock", ident.2 = "Basal_Infected_HCoV-229E", verbose = FALSE)
# head(basal.response_Uninf_vs_Inf_HCOV, n = 15)
#Uninfected vs Bystander
basal.response_Uninf_vs_Bys_HCOV <- FindMarkers(Human, ident.1 = "Basal_Uninfected_Mock", ident.2 = "Basal_Bystander_HCoV-229E", verbose = FALSE)
head(basal.response_Uninf_vs_Bys_HCOV, n = 15)



Marker_Plots <- function(a, b) {
  a_1 <- rownames(a)
  a_2 <- a_1[1:10]
  p_1 <- DoHeatmap(b, features = a_2, group.by = "celltype.status.treat")
  print(p_1)
  return(a_2)
}


#Secretory
Marker_list_secretory.response_Uninf_vs_Bys_MERS <- Marker_Plots(secretory.response_Uninf_vs_Bys_MERS, Human)
Marker_list_secretory.response_Uninf_vs_Bys_HCOV <- Marker_Plots(secretory.response_Uninf_vs_Bys_HCOV, Human)
#Ciliated
Marker_list_ciliated.response_Uninf_vs_Bys_MERS <- Marker_Plots(ciliated.response_Uninf_vs_Bys_MERS, Human)
Marker_list_ciliated.response_Uninf_vs_Bys_HCOV <- Marker_Plots(ciliated.response_Uninf_vs_Bys_HCOV, Human)
#Club
Marker_list_club.response_Uninf_vs_Bys_MERS <- Marker_Plots(club.response_Uninf_vs_Bys_MERS, Human)
Marker_list_club.response_Uninf_vs_Bys_HCOV <- Marker_Plots(club.response_Uninf_vs_Bys_HCOV, Human)
#Deuterosomal
Marker_list_deuterosomal.response_Uninf_vs_Bys_MERS <- Marker_Plots(deuterosomal.response_Uninf_vs_Bys_MERS, Human)
Marker_list_deuterosomal.response_Uninf_vs_Bys_HCOV <- Marker_Plots(deuterosomal.response_Uninf_vs_Bys_HCOV, Human)
#Basal
Marker_list_basal.response_Uninf_vs_Bys_MERS <- Marker_Plots(basal.response_Uninf_vs_Bys_MERS, Human)
Marker_list_basal.response_Uninf_vs_Bys_HCOV <- Marker_Plots(basal.response_Uninf_vs_Bys_HCOV, Human)



##########################################################################################################################
#GO Annotation

#GO_results_Human <- enrichGO(gene = Marker_list_Human, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO_results_secretory.response_Uninf_vs_Bys_MERS <- enrichGO(gene = Marker_list_secretory.response_Uninf_vs_Bys_MERS, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO_results_secretory.response_Uninf_vs_Bys_HCOV <- enrichGO(gene = Marker_list_secretory.response_Uninf_vs_Bys_HCOV, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
#GO_results_MERS_vs_HCoV <- enrichGO(gene = Marker_list_MERS_vs_HCoV, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
#GO_results_Mock <- enrichGO(gene = Marker_list_Mock, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")

#as.data.frame(GO_results_Human)
as.data.frame(GO_results_secretory.response_Uninf_vs_Bys_MERS)
as.data.frame(GO_results_secretory.response_Uninf_vs_Bys_HCOV)
#as.data.frame(GO_results_MERS_vs_HCoV)
#as.data.frame(GO_results_Mock)

#plot(barplot(GO_results_Human, showCategory = 20))
plot(barplot(GO_results_secretory.response_Uninf_vs_Bys_MERS, showCategory = 20))
plot(barplot(GO_results_secretory.response_Uninf_vs_Bys_HCOV, showCategory = 20))
#plot(barplot(GO_results_MERS_vs_HCoV, showCategory = 20))
#plot(barplot(GO_results_Mock, showCategory = 20))


# #Function for later
# GO_Annotation <- function(a) {
#   a_1 <- enrichGO(gene = a, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
#   a_1 <- as.data.frame(a_1)
#   return(a_1)
# }


EnhancedVolcano(Marker_list_secretory.response_Uninf_vs_Bys_MERS, x = "avg_log2FC", y = "p_val_adj", lab = rownames(Marker_list_secretory.response_Uninf_vs_Bys_MERS))
EnhancedVolcano(Marker_list_secretory.response_Uninf_vs_Bys_HCOV, x = "avg_log2FC", y = "p_val_adj", lab = rownames(Marker_list_secretory.response_Uninf_vs_Bys_HCOV))




#Of markers differently expressed in non-infected ciliated cells and bystander ciliated cells we can see different expression patterns in all the cell types:
plots <- VlnPlot(Human, features = c("MYC", "ISG15", "ADSS2", "STC1", "LOC102523265"), split.by = "Status", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)

VlnPlot(Human, features = c("ISG15","LOC102523265"), split.by = "Status", group.by = "celltype", 
        pt.size = 0, combine = FALSE)



#How many of each celltype are infected?
Idents(Human) <- "Status"
Infected <- subset(Human, idents = "Infected")
df <- Human@meta.data


plot_count <- df %>%
  group_by(celltype, Status) %>%
  count(Status)



#Plot2
plot_count %>% ggplot(aes(x = Status, y = n, fill = Status))+
  geom_col()+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        legend.position="bottom")+
  facet_grid(cols = vars(celltype))


