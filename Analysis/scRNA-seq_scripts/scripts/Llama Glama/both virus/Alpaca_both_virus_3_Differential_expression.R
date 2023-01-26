# Differential expression of Single-cell Seq data
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
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(AnnotationHub)
library(AnnotationForge)
library(EnhancedVolcano)
library(ggrepel)
library(stringr)
library(ggvenn)

#Treatment colors(camel-229E, MERS, Mock)
cols_treat = c("#F39C6B","#A5668B", "#96ADC8")
show_col(cols_treat)
#Status colors
cols_stat = c("#CD4631","#DEA47E")
show_col(cols_stat)
#Celltype colors
cols = hcl.colors(7, "Temps")
memory.limit(24000)

##########################################################################################################################
#Loading filtered and merged Seurat Object into Script
#Either load just combined or also clustered data set

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/V.pacos")

Alpaca <- LoadH5Seurat("Alpaca.clustered.h5seurat") # for bulk resolution analysis
#Alpaca <- LoadH5Seurat("Alpaca.combined.h5seurat") # for single cell resolution analysis
head(Alpaca@meta.data)

##################################################################################################################
#Plot relationship between stimulated and unstimulated before to see outliers
Idents(Alpaca) <- "Treat" #Parameters "Treat" needs to be marked as a Idents in the object
# avg.Alpaca <- log1p(AverageExpression(Alpaca, verbose = FALSE)$RNA) #take the data out of RNA assay and not integrated
# avg.Alpaca <- as.data.frame(avg.Alpaca)
# colnames(avg.Alpaca) <- c("Mock", "MERS", "camel229E")
levels(Alpaca)
# p1 <- ggplot(avg.Alpaca, aes(Mock, MERS)) + geom_point() + ggtitle("Mock vs MERS")
# p2 <- ggplot(avg.Alpaca, aes(Mock, camel229E)) + geom_point() + ggtitle("Mock vs camel229E")
# p1 | p2
# 

#DPP4 counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/UMAP_both_dpp4_counts.pdf",
    width=5, height=5)
FeaturePlot(Alpaca, features = "DPP4", cols = c("#F8A251", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T) + NoAxes()
dev.off()

#ANPEP
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/UMAP_both_anpep_counts.pdf",
    width=5, height=5)
FeaturePlot(Alpaca, features = "ANPEP", cols = c("#F8A251", "#034B56"), pt.size = 1, label.size = 5, label = T, repel = T) + NoAxes()
dev.off()

##################################################################################################################################################

# Difference between Treatment (Mock, dcCoV-ACN4, MERS-CoV) for all cells

##################################################################################################################################################
# "Bulk" differential expression in the Treat groups "Mock", "dcCoV-ACN4" "MERS-CoV" 
DefaultAssay(Alpaca) <- "RNA"

# Find differentially expressed features between Mock and treated with camel229E
Mock_vs_ACN4.markers_wilcox <- FindMarkers(Alpaca, slot = "data", ident.2 = "Mock", ident.1 = "dcCoV-ACN4", test.use = "wilcox")
# view results
#head(Mock_vs_ACN4.markers)

# Find differentially expressed features between Mock and treated with MERS-CoV
Mock_vs_MERS.markers_wilcox <- FindMarkers(Alpaca, slot = "data", ident.2 = "Mock", ident.1 = "MERS-CoV", test.use = "wilcox")
# view results
#head(Mock_vs_MERS.markers)


#save the lists
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_both")
write.csv(Mock_vs_ACN4.markers_wilcox,"Llama_Mock_vs_ACN4_markers_wilcox.csv")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_both")
write.csv(Mock_vs_MERS.markers_wilcox,"Llama_Mock_vs_MERS_markers_wilcox.csv")

####################################################################################################################################################
#Volcano Plots

x1 <- FindMarkers(Alpaca, ident.2 = "Mock", ident.1 = "dcCoV-ACN4", test.use = "wilcox")
p1 <- EnhancedVolcano(x1,
                      title = 'Mock versus dcCoV-ACN4 treated',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x1),
                      #col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 0.5,
                      pCutoff = 0.05,
                      colAlpha = 1)

x2 <- FindMarkers(Alpaca, ident.2 = "Mock", ident.1 = "MERS-CoV", test.use = "wilcox")
p2 <- EnhancedVolcano(x2,
                      title = 'Mock versus MERS-CoV treated',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x2),
                      #col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 0.5,
                      pCutoff = 0.05,
                      colAlpha = 1)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/DEG_Mock_v_Virus.pdf",
    width=15, height=7)
p1|p2
dev.off()

?FindMarkers
##################################################################################################################################################

# Difference between Status (Infected, Uninfected, Bystander) for each Virus (dcCoV-ACN4, MERS-CoV)

##################################################################################################################################################
DefaultAssay(Alpaca) <- "RNA"
#Increase memory limit
memory.limit(24000)
# Find differentially expressed features between Infected and uninfected for one Treatment group
# MERS-CoV
Alpaca@meta.data
Idents(Alpaca) <- "Treat"

#Make a fused new identity of each treatment and status combination
Alpaca$treat.status <- paste(Idents(Alpaca), Alpaca$Status, sep = "_")

#make it the new identity
Idents(Alpaca) <- "treat.status"
#check if levels are correct
levels(Alpaca)



# Find differentially expressed features between Mock and MERS-CoV Bystanders
x3 <- FindMarkers(Alpaca, ident.2 = "Mock_Uninfected", ident.1 = "MERS-CoV_Bystander", test.use = "wilcox")
p3 <- EnhancedVolcano(x3,
                      title = 'Uninfected versus MERS-CoV Bystanders',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x3),
                      #col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 0.5,
                      pCutoff = 0.05,
                      colAlpha = 1)


# Find differentially expressed features between Mock and MERS-CoV Infected
x4 <- FindMarkers(Alpaca, ident.2 = "Mock_Uninfected", ident.1 = "MERS-CoV_Infected", test.use = "wilcox")
p4 <- EnhancedVolcano(x4,
                      title = 'Uninfected versus MERS-CoV infected',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x4),
                      #col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 0.5,
                      pCutoff = 0.05,
                      colAlpha = 1)

# Find differentially expressed features between Mock and ACN4 Bystanders
x5 <- FindMarkers(Alpaca, ident.2 = "Mock_Uninfected", ident.1 = "dcCoV-ACN4_Bystander", test.use = "wilcox")
p5 <- EnhancedVolcano(x5,
                      title = 'Uninfected versus ACN4 Bystanders',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x5),
                      #col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 0.5,
                      pCutoff = 0.05,
                      colAlpha = 1)

# Find differentially expressed features between Mock and ACN4 Infected
x6 <- FindMarkers(Alpaca, ident.2 = "Mock_Uninfected", ident.1 = "dcCoV-ACN4_Infected", test.use = "wilcox")
p6 <- EnhancedVolcano(x6,
                      title = 'Uninfected versus ACN4 infected',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x6),
                      #col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 0.5,
                      pCutoff = 0.05,
                      colAlpha = 1)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/DEG_mock_v_MERS_by_inf.pdf",
    width=15, height=7)
p3|p4
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/DEG_mock_v_ACN4_by_inf.pdf",
    width=15, height=7)
p5|p6
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/DEG_MERS_inf_v_ACN4_inf.pdf",
    width=15, height=7)
p4|p6
dev.off()


?FindMarkers()

# Find differentially expressed features between Mock and MERS-CoV Bystanders
Uninf_vs_MERS_bys.markers_wilcox_df <- FindMarkers(Alpaca, ident.2 = "Mock_Uninfected", ident.1 = "MERS-CoV_Bystander", test.use = "wilcox", logfc.threshold = 0.5)

# Find differentially expressed features between Mock and MERS-CoV Infected
Uninf_vs_MERS_inf.markers_wilcox_df <- FindMarkers(Alpaca, ident.2 = "Mock_Uninfected", ident.1 = "MERS-CoV_Infected", test.use = "wilcox", logfc.threshold = 0.5)

# Find differentially expressed features between Mock and ACN4 Bystanders
Uninf_vs_ACN4_bys.markers_wilcox_df <- FindMarkers(Alpaca, ident.2 = "Mock_Uninfected", ident.1 = "dcCoV-ACN4_Bystander", test.use = "wilcox", logfc.threshold = 0.5)

# Find differentially expressed features between Mock and ACN4 Infected
Uninf_vs_ACN4_inf.markers_wilcox_df <- FindMarkers(Alpaca, ident.2 = "Mock_Uninfected", ident.1 = "dcCoV-ACN4_Infected", test.use = "wilcox", logfc.threshold = 0.5)



#save the lists
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_both")
write.csv(Uninf_vs_MERS_bys.markers_wilcox_df,"Llama_Uninf_vs_MERS_bys_markers_wilcox.csv")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_both")
write.csv(Uninf_vs_MERS_inf.markers_wilcox_df,"Llama_Uninf_vs_MERS_inf_markers_wilcox.csv")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_both")
write.csv(Uninf_vs_ACN4_bys.markers_wilcox_df,"Llama_Uninf_vs_ACN4_bys_markers_wilcox.csv")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_both")
write.csv(Uninf_vs_ACN4_inf.markers_wilcox_df,"Llama_Uninf_vs_ACN4_inf_markers_wilcox.csv")


#############################################################################################################################################

#Translate gene lists obtained by DE into the human equivalent to get GO terms with human database

##############################################################################################################################################################################################
#Retrieving the gene lists as vectors insteadt of the dataframe

#up-and-down regulated mers infected
Uninf_vs_MERS_inf.markers_down <- Uninf_vs_MERS_inf.markers_wilcox_df %>%
  filter(Uninf_vs_MERS_inf.markers_wilcox_df$avg_log2FC < 0) %>%
  rownames()

Uninf_vs_MERS_inf.markers_up <- Uninf_vs_MERS_inf.markers_wilcox_df %>%
  filter(Uninf_vs_MERS_inf.markers_wilcox_df$avg_log2FC > 0) %>%
  rownames()

#up-and-down regulated acn4 infected
Uninf_vs_ACN4_inf.markers_down <- Uninf_vs_ACN4_inf.markers_wilcox_df %>%
  filter(Uninf_vs_ACN4_inf.markers_wilcox_df$avg_log2FC < 0) %>%
  rownames()

Uninf_vs_ACN4_inf.markers_up <- Uninf_vs_ACN4_inf.markers_wilcox_df %>%
  filter(Uninf_vs_ACN4_inf.markers_wilcox_df$avg_log2FC > 0) %>%
  rownames()


#venn diagrams

All_genes <- list(`Llama MERS-CoV downregulated` = Uninf_vs_MERS_inf.markers_down,
                  `Llama MERS-CoV upregulated` = Uninf_vs_MERS_inf.markers_up,
                  `Llama dcCoV-ACN4 upregulated` = Uninf_vs_ACN4_inf.markers_up,
                  `Llama dcCoV-ACN4 downregulated` = Uninf_vs_ACN4_inf.markers_down)


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_llama_up_down.pdf",
    width=10, height=9)
ggvenn(All_genes, 
       fill_color = c("#457B9D", "#002C8A","#371F6F","#AB89DD"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)
dev.off()

#together
Uninf_vs_MERS_inf.markers_wilcox <- Uninf_vs_MERS_inf.markers_wilcox_df %>%
  rownames()
Uninf_vs_ACN4_inf.markers_wilcox <- Uninf_vs_ACN4_inf.markers_wilcox_df %>%
  rownames()

features = c("LOC102543804","FOXJ1","LOC116283662","BBOF1","CASC1","TXNIP","AK7","PLAT",
             "CUNH22orf23","SAXO2","RAMP1","DNAAF1","LRRIQ1","WDR66","DYNLRB2","NME9",
             "CDHR4","SEC14L3","BASP1","CCDC173","CETN2","BPIFA1","TPPP3","LOC102528737")

MERS_down_ACN4_down <- intersect(Uninf_vs_MERS_inf.markers_down, Uninf_vs_ACN4_inf.markers_down)

Alpaca_MERS <- subset(x = Alpaca, idents = c("Mock_Uninfected", "MERS-CoV_Infected", "dcCoV-ACN4_Infected"))

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/dotplot_shared_genes_both_virus.pdf",
    width=5, height=7)
DotPlot(Alpaca_MERS, features = c(MERS_down_ACN4_down), scale.by = "radius", scale = T, cols = c("orange", "darkblue"))+coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

DotPlot(Alpaca_MERS, features = features, scale.by = "radius", scale = T, cols = c("orange", "darkblue"))+coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


DotPlot(Alpaca_MERS, features = c(MERS_down_ACN4_down), scale.by = "radius", scale = T, cols = c("orange", "darkblue"))+coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

Uninf_vs_ACN4_inf.markers_up


#############################################################################################################################################

#Translate gene lists obtained by DE into the human equivalent to get GO terms with human database

#############################################################################################################################################

Uninf_vs_ACN4_inf.markers_down_list <- bitr(Uninf_vs_ACN4_inf.markers_down, fromType = "SYMBOL",
                                            toType = c( "ENTREZID"),
                                            OrgDb = org.Hs.eg.db)

Uninf_vs_ACN4_inf.markers_up_list <- bitr(Uninf_vs_ACN4_inf.markers_up, fromType = "SYMBOL",
                                          toType = c("ENTREZID"),
                                          OrgDb = org.Hs.eg.db)

Uninf_vs_ACN4_inf.markers_wilcox_list <- bitr(Uninf_vs_ACN4_inf.markers_wilcox, fromType = "SYMBOL",
                                              toType = c("ENTREZID"),
                                              OrgDb = org.Hs.eg.db)


Uninf_vs_MERS_inf.markers_down_list <- bitr(Uninf_vs_MERS_inf.markers_down, fromType = "SYMBOL",
                                            toType = c("ENTREZID"),
                                            OrgDb = org.Hs.eg.db)

Uninf_vs_MERS_inf.markers_up_list <- bitr(Uninf_vs_MERS_inf.markers_up, fromType = "SYMBOL",
                                          toType = c("ENTREZID"),
                                          OrgDb = org.Hs.eg.db)

Uninf_vs_MERS_inf.markers_wilcox_list <- bitr(Uninf_vs_MERS_inf.markers_wilcox, fromType = "SYMBOL",
                                              toType = c("ENTREZID"),
                                              OrgDb = org.Hs.eg.db)

#################################################################################################################################################

#GO Annotation

#################################################################################################################################################
#Creating geneList file of genes as rownames and FC as a value
#ACN4
names <- rownames(Uninf_vs_ACN4_inf.markers_wilcox_df)
data <- cbind(names,Uninf_vs_ACN4_inf.markers_wilcox_df)

data_fc <- data %>% dplyr::select(1,3)
colnames(data_fc)[1]  <- "SYMBOL"

#for entrezid as rownames
geneList_entrezid <- merge(Uninf_vs_ACN4_inf.markers_wilcox_list, data_fc, by = "SYMBOL", no.dups = T)
geneList_entrezid <- geneList_entrezid %>% dplyr::select(2,3)
## feature 1: numeric vector
geneList_ACN4 = geneList_entrezid[,2]
## feature 2: named vector
names(geneList_ACN4) = as.character(geneList_entrezid[,1])
## feature 3: decreasing orde
geneList_ACN4 = sort(geneList_ACN4, decreasing = TRUE)
################################################################################

?enrichGO



enrichGO <- enrichGO(gene = Uninf_vs_ACN4_inf.markers_wilcox_list$ENTREZID,
                        OrgDb = "org.Hs.eg.db",
                        keyType = "ENTREZID",
                        ont = "BP",
                        qvalueCutoff = 0.05,
                        readable = T
)

enrichGO_down <- enrichGO(gene = Uninf_vs_ACN4_inf.markers_down_list$ENTREZID,
                     OrgDb = "org.Hs.eg.db",
                     keyType = "ENTREZID",
                     ont = "BP",
                     qvalueCutoff = 0.05,
                     readable = T
)

enrichGO_up <- enrichGO(gene = Uninf_vs_ACN4_inf.markers_up_list$ENTREZID,
                     OrgDb = "org.Hs.eg.db",
                     keyType = "ENTREZID",
                     ont = "BP",
                     qvalueCutoff = 0.05,
                     readable = T
)

enrichGO_ACN4 <- enrichGO@result
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Enrich_terms")
write.csv(enrichGO_ACN4,"enrichGO_llama_ACN4_inf_terms.csv")


p1 <- barplot(
  enrichGO_down,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 5,
  font.size = 12,
  title = "",
  label_format = 30
)

p2 <- barplot(
  enrichGO_up,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 5,
  font.size = 12,
  title = "",
  label_format = 30
)

#filtering and ordering the enrichgo results to find out how many terms were significantly enriched
df <- dplyr::filter(enrichGO@result, p.adjust < 0.05)


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/barplot_overview_llama_ACN4.pdf",
    width=6.5, height=7)
barplot(
  enrichGO,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 5,
  font.size = 12,
  title = "",
  label_format = 30
)
dev.off()

?groupGO

ggo <- groupGO(gene     = Uninf_vs_ACN4_inf.markers_wilcox_list$ENTREZID,
                      OrgDb    = org.Hs.eg.db,
                      ont      = "BP",
                      level    = 3,
                      readable = TRUE)

ggo_immune <- groupGO(gene     = Uninf_vs_ACN4_inf.markers_wilcox_list$ENTREZID,
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 3,
               readable = TRUE)


ggo_viral <- groupGO(gene     = Uninf_vs_ACN4_inf.markers_wilcox_list$ENTREZID,
                      OrgDb    = org.Hs.eg.db,
                      ont      = "BP",
                      level    = 3,
                      readable = TRUE)



ggo_immune@result <- filter(ggo_immune@result, grepl('immun', ggo_immune@result$Description))
ggo_immune@result <- ggo_immune@result[order(-ggo_immune@result$Count),]


ggo_viral@result <- filter(ggo_viral@result, grepl('vir', ggo_viral@result$Description))
ggo_viral@result <- ggo_viral@result[order(-ggo_viral@result$Count),]



df1 <- ggo_immune@result
df2 <- ggo_viral@result

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Enrich_terms")
write.csv(df1,"GGO_camel_ACN4_inf_immune_terms.csv")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Enrich_terms")
write.csv(df2,"GGO_camel_ACN4_inf_viral_terms.csv")



#ggo@result <- dplyr::filter(ggo@result, Count > 50)
ggo@result <- ggo@result[order(-ggo@result$Count),]


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/barplot_overview_llama_ACN4.pdf",
    width=5, height=6)
barplot(
  ggo,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 5,
  font.size = 12,
  title = "",
  label_format = 30
)
dev.off()


ggo_immune@result <- filter(ggo_immune@result, grepl('immun', ggo_immune@result$Description))
ggo_immune@result <- ggo_immune@result[order(-ggo_immune@result$Count),]

ggo_viral@result <- filter(ggo_viral@result, grepl('vir', ggo_viral@result$Description))
ggo_viral@result <- ggo_viral@result[order(-ggo_viral@result$Count),]


heatplot(ggo,
         foldChange=geneList_ACN4) + ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red")



pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/cnet_viral_Uninf_ACN4_inf.pdf",
    width=10, height=5)
cnetplot(ggo_viral, foldChange=geneList_ACN4, colorEdge = TRUE, circular = TRUE, layout = 'gem', repel = T)
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/cnet_immune_Uninf_ACN4_inf.pdf",
    width=10, height=5)
cnetplot(ggo_immune, foldChange=geneList_ACN4, colorEdge = TRUE, circular = TRUE, layout = 'gem', repel = T)
dev.off()


#################################################################################################################################################
#Creating geneList file of genes as rownames and FC as a value
#MERS
names <- rownames(Uninf_vs_MERS_inf.markers_wilcox_df)
data <- cbind(names,Uninf_vs_MERS_inf.markers_wilcox_df)

data_fc <- data %>% dplyr::select(1,3)
colnames(data_fc)[1]  <- "SYMBOL"

#for entrezid as rownames
geneList_entrezid <- merge(Uninf_vs_MERS_inf.markers_wilcox_list, data_fc, by = "SYMBOL", no.dups = T)
geneList_entrezid <- geneList_entrezid %>% dplyr::select(2,3)
## feature 1: numeric vector
geneList_MERS = geneList_entrezid[,2]
## feature 2: named vector
names(geneList_MERS) = as.character(geneList_entrezid[,1])
## feature 3: decreasing order
geneList_MERS = sort(geneList_MERS, decreasing = TRUE)
################################################################################

#Uninfected vs MERS infected (all cells) downregulated
enrichGO <- enrichGO(gene = Uninf_vs_MERS_inf.markers_wilcox_list$ENTREZID,
                     OrgDb = "org.Hs.eg.db",
                     keyType = "ENTREZID",
                     ont = "BP",
                     qvalueCutoff = 0.05,
                     readable = T
)

enrichGO_MERS <- enrichGO@result
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Enrich_terms")
write.csv(enrichGO_MERS,"enrichGO_llama_MERS_inf_terms.csv")


#filtering and ordering the enrichgo results to find out how many terms were significantly enriched
df <- dplyr::filter(enrichGO_MERS, p.adjust < 0.05)



pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/barplot_overview_llama_MERS.pdf",
    width=6.5, height=7)
barplot(
  enrichGO,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 5,
  font.size = 12,
  title = "",
  label_format = 30
)
dev.off()







ggo <- groupGO(gene     = Uninf_vs_MERS_inf.markers_wilcox_list$ENTREZID,
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 3,
               readable = TRUE)

ggo_immune <- groupGO(gene     = Uninf_vs_MERS_inf.markers_wilcox_list$ENTREZID,
                      OrgDb    = org.Hs.eg.db,
                      ont      = "BP",
                      level    = 3,
                      readable = TRUE)


ggo_viral <- groupGO(gene     = Uninf_vs_MERS_inf.markers_wilcox_list$ENTREZID,
                     OrgDb    = org.Hs.eg.db,
                     ont      = "BP",
                     level    = 3,
                     readable = TRUE)




df <- ggo@result
df_enr <- filter(df, Count >= 1)

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Enrich_terms")
write.csv(df_enr,"GGO_llama_MERS_inf_terms.csv")


ggo_immune@result <- filter(ggo_immune@result, grepl('immun', ggo_immune@result$Description))
ggo_immune@result <- ggo_immune@result[order(-ggo_immune@result$Count),]


ggo_viral@result <- filter(ggo_viral@result, grepl('vir', ggo_viral@result$Description))
ggo_viral@result <- ggo_viral@result[order(-ggo_viral@result$Count),]




df1 <- ggo_immune@result
df2 <- ggo_viral@result

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Enrich_terms")
write.csv(df1,"GGO_llama_MERS_inf_immune_terms.csv")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Enrich_terms")
write.csv(df2,"GGO_llama_MERS_inf_viral_terms.csv")





#ggo@result <- dplyr::filter(ggo@result, Count > 50)
ggo@result <- ggo@result[order(-ggo@result$Count),]


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/barplot_overview_llama_MERS.pdf",
    width=5, height=6)
barplot(
  ggo,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 5,
  font.size = 12,
  title = "",
  label_format = 30
)
dev.off()


ggo_immune@result <- filter(ggo_immune@result, grepl('immun', ggo_immune@result$Description))
ggo_immune@result <- ggo_immune@result[order(-ggo_immune@result$Count),]

ggo_viral@result <- filter(ggo_viral@result, grepl('vir', ggo_viral@result$Description))
ggo_viral@result <- ggo_viral@result[order(-ggo_viral@result$Count),]


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/cnet_viral_Uninf_MERS_inf.pdf",
    width=10, height=5)
cnetplot(ggo_viral, foldChange=geneList_MERS, colorEdge = TRUE, circular = TRUE, layout = 'gem', repel = T)
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/cnet_immune_Uninf_MERS_inf.pdf",
    width=10, height=5)
cnetplot(ggo_immune, foldChange=geneList_MERS, colorEdge = TRUE, circular = TRUE, layout = 'gem', repel = T)
dev.off()

#################################################################################################################################################
#Uninfected vs ACN4 infected (all cells) downregulated
enrichGO_results_Uninf_vs_ACN4_inf.markers_down <- enrichGO(gene = Uninf_vs_ACN4_inf.markers_down_list$SYMBOL,
                                                            OrgDb = "org.Hs.eg.db",
                                                            keyType = "SYMBOL",
                                                            ont = "BP",
                                                            qvalueCutoff = 0.05,
                                                            minGSSize    = 100,
                                                            maxGSSize    = 500
)

# pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_Alpaca_both_virus/enrich_Uninf_MERS_inf_down.pdf",
#     width=7, height=5)
p5 <- barplot(
  enrichGO_results_Uninf_vs_ACN4_inf.markers_down,
  x = "Count",
  color = "p.adjust",
  showCategory = 8,
  minGSSize    = 100,
  maxGSSize    = 500,
  font.size = 12,
  title = "",
  label_format = 30
)

p6 <- enrichGO_results_Uninf_vs_ACN4_inf.markers_down%>%
  pairwise_termsim() %>%
  emapplot(showCategory = 8, repel = TRUE)

ggo <- groupGO(gene     = Uninf_vs_ACN4_inf.markers_down_list$ENTREZID,
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 3,
               readable = TRUE)

barplot(
  ggo,
  x = "Count",
  color = "p.adjust",
  showCategory = 10,
  minGSSize    = 100,
  maxGSSize    = 500,
  font.size = 12,
  title = "",
  label_format = 30
)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/enrich_Uninf_ACN4_inf_down.pdf",
    width=14, height=9)
p5|p6
dev.off()
#################################################################################################################################################
#Uninfected vs ACN4 infected (all cells) upregulated
enrichGO_results_Uninf_vs_ACN4_inf.markers_up <- enrichGO(gene = Uninf_vs_ACN4_inf.markers_up_list$SYMBOL,
                                                          OrgDb = "org.Hs.eg.db",
                                                          keyType = "SYMBOL",
                                                          ont = "BP",
                                                          qvalueCutoff = 0.05,
                                                          minGSSize    = 100,
                                                          maxGSSize    = 500
)


p7 <- barplot(
  enrichGO_results_Uninf_vs_ACN4_inf.markers_up,
  x = "Count",
  color = "p.adjust",
  showCategory = 8,
  minGSSize    = 100,
  maxGSSize    = 500,
  font.size = 12,
  title = "",
  label_format = 30
)


p8 <- enrichGO_results_Uninf_vs_ACN4_inf.markers_up %>%
  pairwise_termsim() %>%
  emapplot(showCategory = 8, repel = TRUE)

ggo <- groupGO(gene     = Uninf_vs_ACN4_inf.markers_up_list$ENTREZID,
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 2,
               readable = TRUE)

barplot(
  ggo,
  x = "Count",
  color = "p.adjust",
  showCategory = 10,
  minGSSize    = 100,
  maxGSSize    = 500,
  font.size = 12,
  title = "",
  label_format = 30
)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/enrich_Uninf_ACN4_inf_up.pdf",
    width=14, height=9)
p7|p8
dev.off()
#################################################################################################################################################

# Difference between Celltypes in Status (Infected, Uninfected, Bystander) for each Virus (dcCoV-ACN4, MERS-CoV)

##################################################################################################################################################
#Increase memory limit
memory.limit(24000)
# Find differentially expressed features between Infected and uninfected for one Treatment group
# MERS-CoV
Alpaca@meta.data
Idents(Alpaca) <- "celltype"

#Make a fused new identity of each treatment and status combination
Alpaca$cell.treat.status <- paste(Idents(Alpaca), Alpaca$Treat, Alpaca$Status, sep = "_")

#make it the new identity
Idents(Alpaca) <- "cell.treat.status"
#check if levels are correct
levels(Alpaca)


#Volcano Plots
#Secretory#######################################################################################################################################
# Find differentially expressed features between Mock and MERS-CoV Bystanders
x7 <- FindMarkers(Alpaca, ident.1 = "Secretory_Mock_Uninfected", ident.2 = "Secretory_MERS-CoV_Bystander", test.use = "wilcox")
p7 <- EnhancedVolcano(x7,
                      title = 'Secretory Uninfected versus Secretory MERS-CoV Bystanders',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x7),
                      col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 1.5,
                      pCutoff = 0.05,
                      colAlpha = 1)


# Find differentially expressed features between Mock and MERS-CoV Infected
x8 <- FindMarkers(Alpaca, ident.1 = "Secretory_Mock_Uninfected", ident.2 = "Secretory_MERS-CoV_Infected", test.use = "wilcox")
p8 <- EnhancedVolcano(x8,
                      title = 'Secretory Uninfected versus Secretory MERS-CoV Infected',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x8),
                      col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 1.5,
                      pCutoff = 0.05,
                      colAlpha = 1)

# Find differentially expressed features between Mock and ACN4 Bystanders
x9 <- FindMarkers(Alpaca, ident.1 = "Secretory_Mock_Uninfected", ident.2 = "Secretory_dcCoV-ACN4_Bystander", test.use = "wilcox")
p9 <- EnhancedVolcano(x9,
                      title = 'Secretory Uninfected versus Secretory dcCoV-ACN4 Bystanders',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x9),
                      col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 1.5,
                      colAlpha = 1)

# Find differentially expressed features between Mock and ACN4 Infected
x10 <- FindMarkers(Alpaca, ident.1 = "Secretory_Mock_Uninfected", ident.2 = "Secretory_dcCoV-ACN4_Infected", test.use = "wilcox")
p10 <- EnhancedVolcano(x10,
                       title = 'Secretory Uninfected versus Secretory dcCoV-ACN4 Infected',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x10),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 1.5,
                       pCutoff = 0.05,
                       colAlpha = 1)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_Alpaca_both_virus/DEG_secretory_mock_v_MERS_by_inf.pdf",
    width=15, height=7)
p7|p8
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_Alpaca_both_virus/DEG_secretory_mock_v_ACN4_by_inf.pdf",
    width=15, height=7)
p9|p10
dev.off()

#Ciliated#######################################################################################################################################
# Find differentially expressed features between Mock and MERS-CoV Bystanders
x11 <- FindMarkers(Alpaca, ident.1 = "Ciliated_Mock_Uninfected", ident.2 = "Ciliated_MERS-CoV_Bystander", test.use = "wilcox")
p11 <- EnhancedVolcano(x11,
                       title = 'Ciliated Uninfected versus Ciliated MERS-CoV Bystanders',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x11),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 0.5,
                       pCutoff = 0.05,
                       colAlpha = 1)


# Find differentially expressed features between Mock and MERS-CoV Infected
x12 <- FindMarkers(Alpaca, ident.1 = "Ciliated_Mock_Uninfected", ident.2 = "Ciliated_MERS-CoV_Infected", test.use = "wilcox")
p12 <- EnhancedVolcano(x12,
                       title = 'Ciliated Uninfected versus Ciliated MERS-CoV Infected',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x12),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 0.5,
                       pCutoff = 0.05,
                       colAlpha = 1)

# Find differentially expressed features between Mock and ACN4 Bystanders
x13 <- FindMarkers(Alpaca, ident.1 = "Ciliated_Mock_Uninfected", ident.2 = "Ciliated_dcCoV-ACN4_Bystander", test.use = "wilcox")
p13 <- EnhancedVolcano(x13,
                       title = 'Ciliated Uninfected versus Ciliated dcCoV-ACN4 Bystanders',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x13),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 0.5,
                       pCutoff = 0.05,
                       colAlpha = 1)

# Find differentially expressed features between Mock and ACN4 Infected
x14 <- FindMarkers(Alpaca, ident.1 = "Ciliated_Mock_Uninfected", ident.2 = "Ciliated_dcCoV-ACN4_Infected", test.use = "wilcox")
p14 <- EnhancedVolcano(x14,
                       title = 'Ciliated Uninfected versus Ciliated dcCoV-ACN4 Infected',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x14),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 0.5,
                       pCutoff = 0.05,
                       colAlpha = 1)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/DEG_ciliated_mock_v_MERS_by_inf.pdf",
    width=15, height=7)
p11|p12
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/DEG_ciliated_mock_v_ACN4_by_inf.pdf",
    width=15, height=7)
p13|p14
dev.off()


#Llama Cluster 1#######################################################################################################################################
# Find differentially expressed features between Mock and MERS-CoV Bystanders
x15 <- FindMarkers(Alpaca, ident.1 = "Llama Cluster 1_Mock_Uninfected", ident.2 = "Llama Cluster 1_MERS-CoV_Bystander", test.use = "wilcox")
p15 <- EnhancedVolcano(x15,
                       title = 'Llama Cluster 1 Uninfected versus Llama Cluster 1 MERS-CoV Bystanders',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x15),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 0.5,
                       pCutoff = 0.05,
                       colAlpha = 1)


# Find differentially expressed features between Mock and MERS-CoV Infected
x16 <- FindMarkers(Alpaca, ident.1 = "Llama Cluster 1_Mock_Uninfected", ident.2 = "Llama Cluster 1_MERS-CoV_Infected", test.use = "wilcox")
p16 <- EnhancedVolcano(x16,
                       title = 'Llama Cluster 1 Uninfected versus Llama Cluster 1 MERS-CoV Infected',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x16),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 0.5,
                       pCutoff = 0.05,
                       colAlpha = 1)

# Find differentially expressed features between Mock and ACN4 Bystanders
x17 <- FindMarkers(Alpaca, ident.1 = "Llama Cluster 1_Mock_Uninfected", ident.2 = "Llama Cluster 1_dcCoV-ACN4_Bystander", test.use = "wilcox")
p17 <- EnhancedVolcano(x17,
                       title = 'Llama Cluster 1 Uninfected versus Llama Cluster 1 dcCoV-ACN4 Bystanders',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x17),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 0.5,
                       pCutoff = 0.05,
                       colAlpha = 1)

# Find differentially expressed features between Mock and ACN4 Infected
x18 <- FindMarkers(Alpaca, ident.1 = "Llama Cluster 1_Mock_Uninfected", ident.2 = "Llama Cluster 1_dcCoV-ACN4_Infected", test.use = "wilcox")
p18 <- EnhancedVolcano(x18,
                       title = 'Llama Cluster 1 Uninfected versus Llama Cluster 1 dcCoV-ACN4 Infected',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x18),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 0.5,
                       pCutoff = 0.05,
                       colAlpha = 1)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_Alpaca_both_virus/DEG_Llama Cluster 1_mock_v_MERS_by_inf.pdf",
    width=15, height=7)
p15|p16
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_Alpaca_both_virus/DEG_Llama Cluster 1_mock_v_ACN4_by_inf.pdf",
    width=15, height=7)
p17|p18
dev.off()

#Llama Cluster 2#######################################################################################################################################
# Find differentially expressed features between Mock and MERS-CoV Bystanders
x19 <- FindMarkers(Alpaca, ident.1 = "Llama Cluster 2_Mock_Uninfected", ident.2 = "Llama Cluster 2_MERS-CoV_Bystander", test.use = "wilcox")
p19 <- EnhancedVolcano(x19,
                       title = 'Llama Cluster 2 Uninfected versus Llama Cluster 2 MERS-CoV Bystanders',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x19),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 1.5,
                       pCutoff = 0.05,
                       colAlpha = 1)


# Find differentially expressed features between Mock and MERS-CoV Infected
x20 <- FindMarkers(Alpaca, ident.1 = "Llama Cluster 2_Mock_Uninfected", ident.2 = "Llama Cluster 2_MERS-CoV_Infected", test.use = "wilcox")
p20 <- EnhancedVolcano(x20,
                       title = 'Llama Cluster 2 Uninfected versus Llama Cluster 2 MERS-CoV Infected',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x20),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 1.5,
                       pCutoff = 0.05,
                       colAlpha = 1)



pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_Alpaca_both_virus/DEG_Cluster 10_mock_v_MERS_by_inf.pdf",
    width=15, height=7)
p19|p20
dev.off()



#Find DEG features#############################################################################################################################################################
#Mock vs Infected
Secretory_Uninf_vs_MERS_inf.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Secretory_Mock_Uninfected", ident.1 = "Secretory_MERS-CoV_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Secretory_Uninf_vs_MERS_inf.markers_wilcox$SYMBOL <- rownames(Secretory_Uninf_vs_MERS_inf.markers_wilcox)
Secretory_Uninf_vs_ACN4_inf.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Secretory_Mock_Uninfected", ident.1 = "Secretory_dcCoV-ACN4_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Secretory_Uninf_vs_ACN4_inf.markers_wilcox$SYMBOL <- rownames(Secretory_Uninf_vs_ACN4_inf.markers_wilcox)

Ciliated_Uninf_vs_MERS_inf.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Ciliated_Mock_Uninfected", ident.1 = "Ciliated_MERS-CoV_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Ciliated_Uninf_vs_MERS_inf.markers_wilcox$SYMBOL <- rownames(Ciliated_Uninf_vs_MERS_inf.markers_wilcox)
Ciliated_Uninf_vs_ACN4_inf.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Ciliated_Mock_Uninfected", ident.1 = "Ciliated_dcCoV-ACN4_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Ciliated_Uninf_vs_ACN4_inf.markers_wilcox$SYMBOL <- rownames(Ciliated_Uninf_vs_ACN4_inf.markers_wilcox)

Club_Uninf_vs_MERS_inf.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Club_Mock_Uninfected", ident.1 = "Club_MERS-CoV_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Club_Uninf_vs_MERS_inf.markers_wilcox$SYMBOL <- rownames(Club_Uninf_vs_MERS_inf.markers_wilcox)
Club_Uninf_vs_ACN4_inf.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Club_Mock_Uninfected", ident.1 = "Club_dcCoV-ACN4_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Club_Uninf_vs_ACN4_inf.markers_wilcox$SYMBOL <- rownames(Club_Uninf_vs_ACN4_inf.markers_wilcox)

Basal_Uninf_vs_MERS_inf.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Basal_Mock_Uninfected", ident.1 = "Basal_MERS-CoV_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Basal_Uninf_vs_MERS_inf.markers_wilcox$SYMBOL <- rownames(Basal_Uninf_vs_MERS_inf.markers_wilcox)
Basal_Uninf_vs_ACN4_inf.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Basal_Mock_Uninfected", ident.1 = "Basal_dcCoV-ACN4_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Basal_Uninf_vs_ACN4_inf.markers_wilcox$SYMBOL <- rownames(Basal_Uninf_vs_ACN4_inf.markers_wilcox)

Llama_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Llama Cluster 1_Mock_Uninfected", ident.1 = "Llama Cluster 1_MERS-CoV_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Llama_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox$SYMBOL <- rownames(Llama_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox)
Llama_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Llama Cluster 1_Mock_Uninfected", ident.1 = "Llama Cluster 1_dcCoV-ACN4_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Llama_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox$SYMBOL <- rownames(Llama_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox)

Llama_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Llama Cluster 2_Mock_Uninfected", ident.1 = "Llama Cluster 2_MERS-CoV_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Llama_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox$SYMBOL <- rownames(Llama_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox)
Llama_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Llama Cluster 2_Mock_Uninfected", ident.1 = "Llama Cluster 2_dcCoV-ACN4_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Llama_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox$SYMBOL <- rownames(Llama_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox)


#Mock vs Bystander
Secretory_Uninf_vs_MERS_bys.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Secretory_Mock_Uninfected", ident.1 = "Secretory_MERS-CoV_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Secretory_Uninf_vs_MERS_bys.markers_wilcox$SYMBOL <- rownames(Secretory_Uninf_vs_MERS_bys.markers_wilcox)
Secretory_Uninf_vs_ACN4_bys.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Secretory_Mock_Uninfected", ident.1 = "Secretory_dcCoV-ACN4_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Secretory_Uninf_vs_ACN4_bys.markers_wilcox$SYMBOL <- rownames(Secretory_Uninf_vs_ACN4_bys.markers_wilcox)

Ciliated_Uninf_vs_MERS_bys.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Ciliated_Mock_Uninfected", ident.1 = "Ciliated_MERS-CoV_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Ciliated_Uninf_vs_MERS_bys.markers_wilcox$SYMBOL <- rownames(Ciliated_Uninf_vs_MERS_bys.markers_wilcox)
Ciliated_Uninf_vs_ACN4_bys.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Ciliated_Mock_Uninfected", ident.1 = "Ciliated_dcCoV-ACN4_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Ciliated_Uninf_vs_ACN4_bys.markers_wilcox$SYMBOL <- rownames(Ciliated_Uninf_vs_ACN4_bys.markers_wilcox)

Club_Uninf_vs_MERS_bys.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Club_Mock_Uninfected", ident.1 = "Club_MERS-CoV_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Club_Uninf_vs_MERS_bys.markers_wilcox$SYMBOL <- rownames(Club_Uninf_vs_MERS_bys.markers_wilcox)
Club_Uninf_vs_ACN4_bys.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Club_Mock_Uninfected", ident.1 = "Club_dcCoV-ACN4_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Club_Uninf_vs_ACN4_bys.markers_wilcox$SYMBOL <- rownames(Club_Uninf_vs_ACN4_bys.markers_wilcox)

Basal_Uninf_vs_MERS_bys.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Basal_Mock_Uninfected", ident.1 = "Basal_MERS-CoV_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Basal_Uninf_vs_MERS_bys.markers_wilcox$SYMBOL <- rownames(Basal_Uninf_vs_MERS_bys.markers_wilcox)
Basal_Uninf_vs_ACN4_bys.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Basal_Mock_Uninfected", ident.1 = "Basal_dcCoV-ACN4_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Basal_Uninf_vs_ACN4_bys.markers_wilcox$SYMBOL <- rownames(Basal_Uninf_vs_ACN4_bys.markers_wilcox)

Llama_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Llama Cluster 1_Mock_Uninfected", ident.1 = "Llama Cluster 1_MERS-CoV_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Llama_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox$SYMBOL <- rownames(Llama_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox)
Llama_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Llama Cluster 1_Mock_Uninfected", ident.1 = "Llama Cluster 1_dcCoV-ACN4_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Llama_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox$SYMBOL <- rownames(Llama_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox)

Llama_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Llama Cluster 2_Mock_Uninfected", ident.1 = "Llama Cluster 2_MERS-CoV_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Llama_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox$SYMBOL <- rownames(Llama_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox)
Llama_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox <- FindMarkers(Alpaca, ident.2 = "Llama Cluster 2_Mock_Uninfected", ident.1 = "Llama Cluster 2_dcCoV-ACN4_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Llama_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox$SYMBOL <- rownames(Llama_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox)




Secretory_Uninf_vs_MERS_inf.markers_wilcox %>% 
  dplyr::filter(Secretory_Uninf_vs_MERS_inf.markers_wilcox$avg_log2FC > 0) %>%
  count()
Secretory_Uninf_vs_MERS_inf.markers_wilcox %>% 
  dplyr::filter(Secretory_Uninf_vs_MERS_inf.markers_wilcox$avg_log2FC < 0) %>%
  count()


Secretory_Uninf_vs_MERS_bys.markers_wilcox %>% 
  dplyr::filter(Secretory_Uninf_vs_MERS_bys.markers_wilcox$avg_log2FC > 0) %>%
  count()
Secretory_Uninf_vs_MERS_bys.markers_wilcox %>% 
  dplyr::filter(Secretory_Uninf_vs_MERS_bys.markers_wilcox$avg_log2FC < 0) %>%
  count()


Secretory_Uninf_vs_ACN4_inf.markers_wilcox %>% 
  dplyr::filter(Secretory_Uninf_vs_ACN4_inf.markers_wilcox$avg_log2FC > 0) %>%
  count()
Secretory_Uninf_vs_ACN4_inf.markers_wilcox %>% 
  dplyr::filter(Secretory_Uninf_vs_ACN4_inf.markers_wilcox$avg_log2FC < 0) %>%
  count()


Secretory_Uninf_vs_ACN4_bys.markers_wilcox %>% 
  dplyr::filter(Secretory_Uninf_vs_ACN4_bys.markers_wilcox$avg_log2FC > 0) %>%
  count()
Secretory_Uninf_vs_ACN4_bys.markers_wilcox %>% 
  dplyr::filter(Secretory_Uninf_vs_ACN4_bys.markers_wilcox$avg_log2FC < 0) %>%
  count()







Ciliated_Uninf_vs_MERS_inf.markers_wilcox %>% 
  dplyr::filter(Ciliated_Uninf_vs_MERS_inf.markers_wilcox$avg_log2FC > 0) %>%
  count()
Ciliated_Uninf_vs_MERS_inf.markers_wilcox %>% 
  dplyr::filter(Ciliated_Uninf_vs_MERS_inf.markers_wilcox$avg_log2FC < 0) %>%
  count()


Ciliated_Uninf_vs_MERS_bys.markers_wilcox %>% 
  dplyr::filter(Ciliated_Uninf_vs_MERS_bys.markers_wilcox$avg_log2FC > 0) %>%
  count()
Ciliated_Uninf_vs_MERS_bys.markers_wilcox %>% 
  dplyr::filter(Ciliated_Uninf_vs_MERS_bys.markers_wilcox$avg_log2FC < 0) %>%
  count()


Ciliated_Uninf_vs_ACN4_inf.markers_wilcox %>% 
  dplyr::filter(Ciliated_Uninf_vs_ACN4_inf.markers_wilcox$avg_log2FC > 0) %>%
  count()
Ciliated_Uninf_vs_ACN4_inf.markers_wilcox %>% 
  dplyr::filter(Ciliated_Uninf_vs_ACN4_inf.markers_wilcox$avg_log2FC < 0) %>%
  count()


Ciliated_Uninf_vs_ACN4_bys.markers_wilcox %>% 
  dplyr::filter(Ciliated_Uninf_vs_ACN4_bys.markers_wilcox$avg_log2FC > 0) %>%
  count()
Ciliated_Uninf_vs_ACN4_bys.markers_wilcox %>% 
  dplyr::filter(Ciliated_Uninf_vs_ACN4_bys.markers_wilcox$avg_log2FC < 0) %>%
  count()



Basal_Uninf_vs_MERS_inf.markers_wilcox %>% 
  dplyr::filter(Basal_Uninf_vs_MERS_inf.markers_wilcox$avg_log2FC > 0) %>%
  count()
Basal_Uninf_vs_MERS_inf.markers_wilcox %>% 
  dplyr::filter(Basal_Uninf_vs_MERS_inf.markers_wilcox$avg_log2FC < 0) %>%
  count()

Basal_Uninf_vs_MERS_bys.markers_wilcox %>% 
  dplyr::filter(Basal_Uninf_vs_MERS_bys.markers_wilcox$avg_log2FC > 0) %>%
  count()
Basal_Uninf_vs_MERS_bys.markers_wilcox %>% 
  dplyr::filter(Basal_Uninf_vs_MERS_bys.markers_wilcox$avg_log2FC < 0) %>%
  count()

Basal_Uninf_vs_ACN4_inf.markers_wilcox %>% 
  dplyr::filter(Basal_Uninf_vs_ACN4_inf.markers_wilcox$avg_log2FC > 0) %>%
  count()
Basal_Uninf_vs_ACN4_inf.markers_wilcox %>% 
  dplyr::filter(Basal_Uninf_vs_ACN4_inf.markers_wilcox$avg_log2FC < 0) %>%
  count()

Basal_Uninf_vs_ACN4_bys.markers_wilcox %>% 
  dplyr::filter(Basal_Uninf_vs_ACN4_bys.markers_wilcox$avg_log2FC > 0) %>%
  count()
Basal_Uninf_vs_ACN4_bys.markers_wilcox %>% 
  dplyr::filter(Basal_Uninf_vs_ACN4_bys.markers_wilcox$avg_log2FC < 0) %>%
  count()




Club_Uninf_vs_MERS_inf.markers_wilcox %>% 
  dplyr::filter(Club_Uninf_vs_MERS_inf.markers_wilcox$avg_log2FC > 0) %>%
  count()
Club_Uninf_vs_MERS_inf.markers_wilcox %>% 
  dplyr::filter(Club_Uninf_vs_MERS_inf.markers_wilcox$avg_log2FC < 0) %>%
  count()

Club_Uninf_vs_MERS_bys.markers_wilcox %>% 
  dplyr::filter(Club_Uninf_vs_MERS_bys.markers_wilcox$avg_log2FC > 0) %>%
  count()
Club_Uninf_vs_MERS_bys.markers_wilcox %>% 
  dplyr::filter(Club_Uninf_vs_MERS_bys.markers_wilcox$avg_log2FC < 0) %>%
  count()

Club_Uninf_vs_ACN4_inf.markers_wilcox %>% 
  dplyr::filter(Club_Uninf_vs_ACN4_inf.markers_wilcox$avg_log2FC > 0) %>%
  count()
Club_Uninf_vs_ACN4_inf.markers_wilcox %>% 
  dplyr::filter(Club_Uninf_vs_ACN4_inf.markers_wilcox$avg_log2FC < 0) %>%
  count()

Club_Uninf_vs_ACN4_bys.markers_wilcox %>% 
  dplyr::filter(Club_Uninf_vs_ACN4_bys.markers_wilcox$avg_log2FC > 0) %>%
  count()
Club_Uninf_vs_ACN4_bys.markers_wilcox %>% 
  dplyr::filter(Club_Uninf_vs_ACN4_bys.markers_wilcox$avg_log2FC < 0) %>%
  count()




Llama_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox %>% 
  dplyr::filter(Llama_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox$avg_log2FC > 0) %>%
  count()
Llama_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox %>% 
  dplyr::filter(Llama_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox$avg_log2FC < 0) %>%
  count()

Llama_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox %>% 
  dplyr::filter(Llama_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox$avg_log2FC > 0) %>%
  count()
Llama_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox %>% 
  dplyr::filter(Llama_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox$avg_log2FC < 0) %>%
  count()

Llama_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox %>% 
  dplyr::filter(Llama_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox$avg_log2FC > 0) %>%
  count()
Llama_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox %>% 
  dplyr::filter(Llama_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox$avg_log2FC < 0) %>%
  count()

Llama_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox %>% 
  dplyr::filter(Llama_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox$avg_log2FC > 0) %>%
  count()
Llama_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox %>% 
  dplyr::filter(Llama_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox$avg_log2FC < 0) %>%
  count()




Llama_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox %>% 
  dplyr::filter(Llama_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox$avg_log2FC > 0) %>%
  count()
Llama_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox %>% 
  dplyr::filter(Llama_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox$avg_log2FC < 0) %>%
  count()

Llama_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox %>% 
  dplyr::filter(Llama_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox$avg_log2FC > 0) %>%
  count()
Llama_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox %>% 
  dplyr::filter(Llama_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox$avg_log2FC < 0) %>%
  count()

Llama_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox %>% 
  dplyr::filter(Llama_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox$avg_log2FC > 0) %>%
  count()
Llama_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox %>% 
  dplyr::filter(Llama_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox$avg_log2FC < 0) %>%
  count()

Llama_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox %>% 
  dplyr::filter(Llama_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox$avg_log2FC > 0) %>%
  count()
Llama_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox %>% 
  dplyr::filter(Llama_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox$avg_log2FC < 0) %>%
  count()


#save the lists
#Infected
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_both")
write.csv(Secretory_Uninf_vs_MERS_inf.markers_wilcox,"Secretory_Llama_MERS_Infected.csv")
write.csv(Secretory_Uninf_vs_ACN4_inf.markers_wilcox,"Secretory_Llama_ACN4_Infected.csv")

write.csv(Ciliated_Uninf_vs_MERS_inf.markers_wilcox,"Ciliated_Llama_MERS_Infected.csv")
write.csv(Ciliated_Uninf_vs_ACN4_inf.markers_wilcox,"Ciliated_Llama_ACN4_Infected.csv")

write.csv(Club_Uninf_vs_MERS_inf.markers_wilcox,"Club_Llama_MERS_Infected.csv")
write.csv(Club_Uninf_vs_ACN4_inf.markers_wilcox,"Club_Llama_ACN4_Infected.csv")

write.csv(Basal_Uninf_vs_MERS_inf.markers_wilcox,"Basal_Llama_MERS_Infected.csv")
write.csv(Basal_Uninf_vs_ACN4_inf.markers_wilcox,"Basal_Llama_ACN4_Infected.csv")

write.csv(Llama_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox,"Cluster1_Llama_MERS_Infected.csv")
write.csv(Llama_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox,"Cluster1_Llama_ACN4_Infected.csv")

write.csv(Llama_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox,"Cluster2_Llama_MERS_Infected.csv")
write.csv(Llama_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox,"Cluster2_Llama_ACN4_Infected.csv")

#Bystander
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_both")
write.csv(Secretory_Uninf_vs_MERS_bys.markers_wilcox,"Secretory_Llama_MERS_Bystander.csv")
write.csv(Secretory_Uninf_vs_ACN4_bys.markers_wilcox,"Secretory_Llama_ACN4_Bystander.csv")

write.csv(Ciliated_Uninf_vs_MERS_bys.markers_wilcox,"Ciliated_Llama_MERS_Bystander.csv")
write.csv(Ciliated_Uninf_vs_ACN4_bys.markers_wilcox,"Ciliated_Llama_ACN4_Bystander.csv")

write.csv(Club_Uninf_vs_MERS_bys.markers_wilcox,"Club_Llama_MERS_Bystander.csv")
write.csv(Club_Uninf_vs_ACN4_bys.markers_wilcox,"Club_Llama_ACN4_Bystander.csv")

write.csv(Basal_Uninf_vs_MERS_bys.markers_wilcox,"Basal_Llama_MERS_Bystander.csv")
write.csv(Basal_Uninf_vs_ACN4_bys.markers_wilcox,"Basal_Llama_ACN4_Bystander.csv")

write.csv(Llama_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox,"Cluster1_Llama_MERS_Bystander.csv")
write.csv(Llama_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox,"Cluster1_Llama_ACN4_Bystander.csv")

write.csv(Llama_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox,"Cluster2_Llama_MERS_Bystander.csv")
write.csv(Llama_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox,"Cluster2_Llama_ACN4_Bystander.csv")



##############################################################################################################################################################################################
#Retrieving the gene lists as vectors insteadt of the dataframe
#Infected
Secretory_Uninf_vs_MERS_inf_rownames <- rownames(Secretory_Uninf_vs_MERS_inf.markers_wilcox)
Secretory_Uninf_vs_ACN4_inf_rownames <- rownames(Secretory_Uninf_vs_ACN4_inf.markers_wilcox)

Club_Uninf_vs_MERS_inf_rownames <- rownames(Club_Uninf_vs_MERS_inf.markers_wilcox)
Club_Uninf_vs_ACN4_inf_rownames <- rownames(Club_Uninf_vs_ACN4_inf.markers_wilcox)

Basal_Uninf_vs_MERS_inf_rownames <- rownames(Basal_Uninf_vs_MERS_inf.markers_wilcox)
Basal_Uninf_vs_ACN4_inf_rownames <- rownames(Basal_Uninf_vs_ACN4_inf.markers_wilcox)

Ciliated_Uninf_vs_MERS_inf_rownames <- rownames(Ciliated_Uninf_vs_MERS_inf.markers_wilcox)
Ciliated_Uninf_vs_ACN4_inf_rownames <- rownames(Ciliated_Uninf_vs_ACN4_inf.markers_wilcox)

Llama_Cluster_1_Uninf_vs_MERS_inf_rownames <- rownames(Llama_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox)
Llama_Cluster_1_Uninf_vs_ACN4_inf_rownames <- rownames(Llama_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox)

Llama_Cluster_2_Uninf_vs_MERS_inf_rownames <- rownames(Llama_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox)
Llama_Cluster_2_Uninf_vs_ACN4_inf_rownames <- rownames(Llama_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox)

#Bystander
Secretory_Uninf_vs_MERS_bys_rownames <- rownames(Secretory_Uninf_vs_MERS_bys.markers_wilcox)
Secretory_Uninf_vs_ACN4_bys_rownames <- rownames(Secretory_Uninf_vs_ACN4_bys.markers_wilcox)

Club_Uninf_vs_MERS_bys_rownames <- rownames(Club_Uninf_vs_MERS_bys.markers_wilcox)
Club_Uninf_vs_ACN4_bys_rownames <- rownames(Club_Uninf_vs_ACN4_bys.markers_wilcox)

Basal_Uninf_vs_MERS_bys_rownames <- rownames(Basal_Uninf_vs_MERS_bys.markers_wilcox)
Basal_Uninf_vs_ACN4_bys_rownames <- rownames(Basal_Uninf_vs_ACN4_bys.markers_wilcox)

Ciliated_Uninf_vs_MERS_bys_rownames <- rownames(Ciliated_Uninf_vs_MERS_bys.markers_wilcox)
Ciliated_Uninf_vs_ACN4_bys_rownames <- rownames(Ciliated_Uninf_vs_ACN4_bys.markers_wilcox)

Llama_Cluster_1_Uninf_vs_MERS_bys_rownames <- rownames(Llama_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox)
Llama_Cluster_1_Uninf_vs_ACN4_bys_rownames <- rownames(Llama_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox)

Llama_Cluster_2_Uninf_vs_MERS_bys_rownames <- rownames(Llama_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox)
Llama_Cluster_2_Uninf_vs_ACN4_bys_rownames <- rownames(Llama_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox)

#############################################################################################################################################

#Translate gene lists obtained by DE into the human equivalent to get GO terms with human database

#############################################################################################################################################


Secretory_Llama_MERS_Infected <- bitr(Secretory_Uninf_vs_MERS_inf_rownames, fromType = "SYMBOL",
                                      toType = "ENTREZID",
                                      OrgDb = org.Hs.eg.db)
Secretory_Llama_MERS_Bystander <- bitr(Secretory_Uninf_vs_MERS_bys_rownames, fromType = "SYMBOL",
                                       toType = "ENTREZID",
                                       OrgDb = org.Hs.eg.db)
Secretory_Llama_ACN4_Infected <- bitr(Secretory_Uninf_vs_ACN4_inf_rownames, fromType = "SYMBOL",
                                      toType = "ENTREZID",
                                      OrgDb = org.Hs.eg.db)
Secretory_Llama_ACN4_Bystander <- bitr(Secretory_Uninf_vs_ACN4_bys_rownames, fromType = "SYMBOL",
                                       toType = "ENTREZID",
                                       OrgDb = org.Hs.eg.db)



Ciliated_Llama_MERS_Infected <- bitr(Ciliated_Uninf_vs_MERS_inf_rownames, fromType = "SYMBOL",
                                     toType = "ENTREZID",
                                     OrgDb = org.Hs.eg.db)
Ciliated_Llama_MERS_Bystander <- bitr(Ciliated_Uninf_vs_MERS_bys_rownames, fromType = "SYMBOL",
                                      toType = "ENTREZID",
                                      OrgDb = org.Hs.eg.db)
Ciliated_Llama_ACN4_Infected <- bitr(Ciliated_Uninf_vs_ACN4_inf_rownames, fromType = "SYMBOL",
                                     toType = "ENTREZID",
                                     OrgDb = org.Hs.eg.db)
Ciliated_Llama_ACN4_Bystander <- bitr(Ciliated_Uninf_vs_ACN4_bys_rownames, fromType = "SYMBOL",
                                      toType = "ENTREZID",
                                      OrgDb = org.Hs.eg.db)


Club_Llama_MERS_Infected <- bitr(Club_Uninf_vs_MERS_inf_rownames, fromType = "SYMBOL",
                                 toType = "ENTREZID",
                                 OrgDb = org.Hs.eg.db)
Club_Llama_MERS_Bystander <- bitr(Club_Uninf_vs_MERS_bys_rownames, fromType = "SYMBOL",
                                  toType = "ENTREZID",
                                  OrgDb = org.Hs.eg.db)
Club_Llama_ACN4_Infected <- bitr(Club_Uninf_vs_ACN4_inf_rownames, fromType = "SYMBOL",
                                 toType = "ENTREZID",
                                 OrgDb = org.Hs.eg.db)
Club_Llama_ACN4_Bystander <- bitr(Club_Uninf_vs_ACN4_bys_rownames, fromType = "SYMBOL",
                                  toType = "ENTREZID",
                                  OrgDb = org.Hs.eg.db)


Basal_Llama_MERS_Infected <- bitr(Basal_Uninf_vs_MERS_inf_rownames, fromType = "SYMBOL",
                                  toType = "ENTREZID",
                                  OrgDb = org.Hs.eg.db)
Basal_Llama_MERS_Bystander <- bitr(Basal_Uninf_vs_MERS_bys_rownames, fromType = "SYMBOL",
                                   toType = "ENTREZID",
                                   OrgDb = org.Hs.eg.db)
Basal_Llama_ACN4_Infected <- bitr(Basal_Uninf_vs_ACN4_inf_rownames, fromType = "SYMBOL",
                                  toType = "ENTREZID",
                                  OrgDb = org.Hs.eg.db)
Basal_Llama_ACN4_Bystander <- bitr(Basal_Uninf_vs_ACN4_bys_rownames, fromType = "SYMBOL",
                                   toType = "ENTREZID",
                                   OrgDb = org.Hs.eg.db)


Cluster1_Llama_MERS_Infected <- bitr(Llama_Cluster_1_Uninf_vs_MERS_inf_rownames, fromType = "SYMBOL",
                                     toType = "ENTREZID",
                                     OrgDb = org.Hs.eg.db)
Cluster1_Llama_MERS_Bystander <- bitr(Llama_Cluster_1_Uninf_vs_MERS_bys_rownames, fromType = "SYMBOL",
                                      toType = "ENTREZID",
                                      OrgDb = org.Hs.eg.db)
Cluster1_Llama_ACN4_Infected <- bitr(Llama_Cluster_1_Uninf_vs_ACN4_inf_rownames, fromType = "SYMBOL",
                                     toType = "ENTREZID",
                                     OrgDb = org.Hs.eg.db)
Cluster1_Llama_ACN4_Bystander <- bitr(Llama_Cluster_1_Uninf_vs_ACN4_bys_rownames, fromType = "SYMBOL",
                                      toType = "ENTREZID",
                                      OrgDb = org.Hs.eg.db)


Cluster2_Llama_MERS_Infected <- bitr(Llama_Cluster_2_Uninf_vs_MERS_inf_rownames, fromType = "SYMBOL",
                                     toType = "ENTREZID",
                                     OrgDb = org.Hs.eg.db)
Cluster2_Llama_MERS_Bystander <- bitr(Llama_Cluster_2_Uninf_vs_MERS_bys_rownames, fromType = "SYMBOL",
                                      toType = "ENTREZID",
                                      OrgDb = org.Hs.eg.db)
Cluster2_Llama_ACN4_Infected <- bitr(Llama_Cluster_2_Uninf_vs_ACN4_inf_rownames, fromType = "SYMBOL",
                                     toType = "ENTREZID",
                                     OrgDb = org.Hs.eg.db)
Cluster2_Llama_ACN4_Bystander <- bitr(Llama_Cluster_2_Uninf_vs_ACN4_bys_rownames, fromType = "SYMBOL",
                                      toType = "ENTREZID",
                                      OrgDb = org.Hs.eg.db)




#Merge the two datasets together
#Infected
Secretory_Llama_MERS_Infected <- merge(Secretory_Llama_MERS_Infected, Secretory_Uninf_vs_MERS_inf.markers_wilcox, by = "SYMBOL")
Secretory_Llama_ACN4_Infected <- merge(Secretory_Llama_ACN4_Infected, Secretory_Uninf_vs_ACN4_inf.markers_wilcox, by = "SYMBOL")

Ciliated_Llama_MERS_Infected <- merge(Ciliated_Llama_MERS_Infected, Ciliated_Uninf_vs_MERS_inf.markers_wilcox, by = "SYMBOL")
Ciliated_Llama_ACN4_Infected <- merge(Ciliated_Llama_ACN4_Infected, Ciliated_Uninf_vs_ACN4_inf.markers_wilcox, by = "SYMBOL")

Club_Llama_MERS_Infected <- merge(Club_Llama_MERS_Infected, Club_Uninf_vs_MERS_inf.markers_wilcox, by = "SYMBOL")
Club_Llama_ACN4_Infected <- merge(Club_Llama_ACN4_Infected, Club_Uninf_vs_ACN4_inf.markers_wilcox, by = "SYMBOL")

Basal_Llama_MERS_Infected <- merge(Basal_Llama_MERS_Infected, Basal_Uninf_vs_MERS_inf.markers_wilcox, by = "SYMBOL")
Basal_Llama_ACN4_Infected <- merge(Basal_Llama_ACN4_Infected, Basal_Uninf_vs_ACN4_inf.markers_wilcox, by = "SYMBOL")

Cluster1_Llama_MERS_Infected <- merge(Cluster1_Llama_MERS_Infected, Llama_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox, by = "SYMBOL")
Cluster1_Llama_ACN4_Infected <- merge(Cluster1_Llama_ACN4_Infected, Llama_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox, by = "SYMBOL")

Cluster2_Llama_MERS_Infected <- merge(Cluster2_Llama_MERS_Infected, Llama_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox, by = "SYMBOL")
Cluster2_Llama_ACN4_Infected <- merge(Cluster2_Llama_ACN4_Infected, Llama_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox, by = "SYMBOL")

#Bystander
Secretory_Llama_MERS_Bystander <- merge(Secretory_Llama_MERS_Bystander, Secretory_Uninf_vs_MERS_bys.markers_wilcox, by = "SYMBOL")
Secretory_Llama_ACN4_Bystander <- merge(Secretory_Llama_ACN4_Bystander, Secretory_Uninf_vs_ACN4_bys.markers_wilcox, by = "SYMBOL")

Ciliated_Llama_MERS_Bystander <- merge(Ciliated_Llama_MERS_Bystander, Ciliated_Uninf_vs_MERS_bys.markers_wilcox, by = "SYMBOL")
Ciliated_Llama_ACN4_Bystander <- merge(Ciliated_Llama_ACN4_Bystander, Ciliated_Uninf_vs_ACN4_bys.markers_wilcox, by = "SYMBOL")

Club_Llama_MERS_Bystander <- merge(Club_Llama_MERS_Bystander, Club_Uninf_vs_MERS_bys.markers_wilcox, by = "SYMBOL")
Club_Llama_ACN4_Bystander <- merge(Club_Llama_ACN4_Bystander, Club_Uninf_vs_ACN4_bys.markers_wilcox, by = "SYMBOL")

Basal_Llama_MERS_Bystander <- merge(Basal_Llama_MERS_Bystander, Basal_Uninf_vs_MERS_bys.markers_wilcox, by = "SYMBOL")
Basal_Llama_ACN4_Bystander <- merge(Basal_Llama_ACN4_Bystander, Basal_Uninf_vs_ACN4_bys.markers_wilcox, by = "SYMBOL")

Cluster1_Llama_MERS_Bystander <- merge(Cluster1_Llama_MERS_Bystander, Llama_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox, by = "SYMBOL")
Cluster1_Llama_ACN4_Bystander <- merge(Cluster1_Llama_ACN4_Bystander, Llama_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox, by = "SYMBOL")

Cluster2_Llama_MERS_Bystander <- merge(Cluster2_Llama_MERS_Bystander, Llama_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox, by = "SYMBOL")
Cluster2_Llama_ACN4_Bystander <- merge(Cluster2_Llama_ACN4_Bystander, Llama_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox, by = "SYMBOL")


#Dotplot########################################################################################################################################
#save the lists
#Infected
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_both")
write.csv(Secretory_Llama_MERS_Infected,"Secretory_Llama_MERS_Infected.csv")
write.csv(Secretory_Llama_ACN4_Infected,"Secretory_Llama_ACN4_Infected.csv")

write.csv(Ciliated_Llama_MERS_Infected,"Ciliated_Llama_MERS_Infected.csv")
write.csv(Ciliated_Llama_ACN4_Infected,"Ciliated_Llama_ACN4_Infected.csv")

write.csv(Club_Llama_MERS_Infected,"Club_Llama_MERS_Infected.csv")
write.csv(Club_Llama_ACN4_Infected,"Club_Llama_ACN4_Infected.csv")

write.csv(Basal_Llama_MERS_Infected,"Basal_Llama_MERS_Infected.csv")
write.csv(Basal_Llama_ACN4_Infected,"Basal_Llama_ACN4_Infected.csv")

write.csv(Cluster1_Llama_MERS_Infected,"Cluster1_Llama_MERS_Infected.csv")
write.csv(Cluster1_Llama_ACN4_Infected,"Cluster1_Llama_ACN4_Infected.csv")

write.csv(Cluster2_Llama_MERS_Infected,"Cluster2_Llama_MERS_Infected.csv")
write.csv(Cluster2_Llama_ACN4_Infected,"Cluster2_Llama_ACN4_Infected.csv")

#Bystander
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_both")
write.csv(Secretory_Llama_MERS_Bystander,"Secretory_Llama_MERS_Bystander.csv")
write.csv(Secretory_Llama_ACN4_Bystander,"Secretory_Llama_ACN4_Bystander.csv")

write.csv(Ciliated_Llama_MERS_Bystander,"Ciliated_Llama_MERS_Bystander.csv")
write.csv(Ciliated_Llama_ACN4_Bystander,"Ciliated_Llama_ACN4_Bystander.csv")

write.csv(Club_Llama_MERS_Bystander,"Club_Llama_MERS_Bystander.csv")
write.csv(Club_Llama_ACN4_Bystander,"Club_Llama_ACN4_Bystander.csv")

write.csv(Basal_Llama_MERS_Bystander,"Basal_Llama_MERS_Bystander.csv")
write.csv(Basal_Llama_ACN4_Bystander,"Basal_Llama_ACN4_Bystander.csv")

write.csv(Cluster1_Llama_MERS_Bystander,"Cluster1_Llama_MERS_Bystander.csv")
write.csv(Cluster1_Llama_ACN4_Bystander,"Cluster1_Llama_ACN4_Bystander.csv")

write.csv(Cluster2_Llama_MERS_Bystander,"Cluster2_Llama_MERS_Bystander.csv")
write.csv(Cluster2_Llama_ACN4_Bystander,"Cluster2_Llama_ACN4_Bystander.csv")








