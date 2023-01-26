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
library(GO.db)
library(AnnotationDbi)
library(EnhancedVolcano)
library(ggrepel)
library(DESeq2, quietly=T)
library(ComplexHeatmap, quietly = T)
library(circlize, quietly = T)
library(DEGreport, quietly=T)
library(tidyverse, quietly = T)
library(data.table, quietly = T)
library(clusterProfiler, quietly = T)
library(enrichplot, quietly = T)
library(enrichR)
library(ReactomePA, quietly = T)
library(ggVennDiagram,quietly = T)
library(PCAtools, quietly = T)
library(gprofiler2)
library(biomaRt)
library(gprofiler2)
library(ggVennDiagram)
library(ggupset)
library(MAST)
library(stringr)
library(ggvenn)

#Setting the color schemes for the script:

gene <- names(geneList)
#Treatment colors(camel-229E, MERS, Mock)
cols_treat = c("#F39C6B","#A5668B", "#96ADC8")
show_col(cols_treat)
#Status colors
cols_stat = c("#CD4631","#DEA47E")
show_col(cols_stat)
#Celltype colors for initial UMAP
cols = hcl.colors(7, "Temps")


show_col(viridis(50))
# The celltypes and their assigned color:
#Secretory = "#089392"
#Ciliated = "#40B48B"
#Club = "#9DCD84"
#Basal = "#EAE29C"
#Suprabasal = "#EAB672"
#Deuterosomal = "#E6866A"
#Ionocytes = "#CF597E"

memory.limit(24000)


##########################################################################################################################
#Loading filtered and merged Seurat Object into Script
#Infected and Mock Samples merged and cell clusters assigned

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.ferus")

Ferus <- LoadH5Seurat("Ferus.clustered.h5seurat")

# head(Ferus@meta.data)
# str(Ferus_metadata)

#DPP4 counts
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_both_dpp4_counts.pdf",
    width=5, height=5)
FeaturePlot(Ferus, features = "DPP4", cols = c("#F8A251", "darkblue"), pt.size = 1, label.size = 5, label = T, repel = T) + NoAxes()
dev.off()

#ANPEP
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/UMAP_both_anpep_counts.pdf",
    width=5, height=5)
FeaturePlot(Ferus, features = "ANPEP", cols = c("#F8A251", "#034B56"), pt.size = 1, label.size = 5, label = T, repel = T) + NoAxes()
dev.off()


##################################################################################################################
#Plot relationship between stimulated and unstimulated before to see outliers
Idents(Ferus) <- "celltype" #Parameters "Treat" needs to be marked as a Idents in the object
DefaultAssay(Ferus) <- "RNA"
##################################################################################################################################################

# Difference between Treatment (Mock, HCoV-229E, MERS-CoV) for all cells - "Bulk"

##################################################################################################################################################
# "Bulk" differential expression in the Treatment groups Mock, HcoV-229E, MERS-CoV
DefaultAssay(Ferus) <- "RNA"
#Ferus_inf.markers <- FindAllMarkers(Ferus, only.pos = TRUE, min.pct = 0.5, grouping.var = "Treat", logfc.threshold = 0.5)

# Find differentially expressed features between Mock and treated with ACN4
Mock_vs_ACN4.markers_wilcox <- FindMarkers(Ferus, slot = "data", ident.2 = "Mock", ident.1 = "dcCoV-ACN4", test.use = "wilcox")

# Find differentially expressed features between Mock and treated with MERS-CoV
Mock_vs_MERS.markers_wilcox <- FindMarkers(Ferus, slot = "data", ident.2 = "Mock", ident.1 = "MERS-CoV", test.use = "wilcox")


#save the lists
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Ferus_both")
write.csv(Mock_vs_MERS.markers_wilcox,"Camel_Mock_vs_MERS_markers_wilcox.csv")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Ferus_both")
write.csv(Mock_vs_ACN4.markers_wilcox,"Camel_Mock_vs_ACN4_markers_wilcox.csv")

#Volcano Plots####################################################################################################################################################


x1 <- FindMarkers(Ferus, ident.2 = "Mock", ident.1 = "dcCoV-ACN4", test.use = "wilcox")
p1 <- EnhancedVolcano(x1,
                  title = 'Mock versus dcCoV-ACN4 treated',
                  x = "avg_log2FC",
                  y = "p_val_adj",
                  pointSize = 3.0,
                  labSize = 6.0,
                  lab = rownames(x1),
                  col=c('black', 'black', 'black', 'red3'),
                  FCcutoff = 0.5,
                  colAlpha = 1)

x2 <- FindMarkers(Ferus, ident.2 = "Mock", ident.1 = "MERS-CoV", test.use = "wilcox")
p2 <- EnhancedVolcano(x2,
                  title = 'Mock versus MERS-CoV treated',
                  x = "avg_log2FC",
                  y = "p_val_adj",
                  pointSize = 3.0,
                  labSize = 6.0,
                  lab = rownames(x2),
                  col=c('black', 'black', 'black', 'red3'),
                  FCcutoff = 0.5,
                  colAlpha = 1)
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/DEG_Mock_v_Virus.pdf",
    width=15, height=7)
p1|p2
dev.off()

##################################################################################################################################################

# Difference between Status (Infected, Uninfected, Bystander) for each Virus (dcCoV-ACN4, MERS-CoV)

##################################################################################################################################################
DefaultAssay(Ferus) <- "RNA"
#Increase memory limit
memory.limit(24000)
# Find differentially expressed features between Infected and uninfected for one Treatment group
# MERS-CoV
str(Ferus)
Ferus@meta.data
Idents(Ferus) <- "Treat"


#Make a fused new identity of each treatment and status combination
Ferus$treat.status <- paste(Idents(Ferus), Ferus$Status, sep = "_")

#make it the new identity
Idents(Ferus) <- "treat.status"
#check if levels are correct
levels(Ferus)

#Volcano Plots####################################################################################################################################################

# Find differentially expressed features between Mock and MERS-CoV Bystanders
x3 <- FindMarkers(Ferus, ident.2 = "Mock_Uninfected", ident.1 = "MERS-CoV_Bystander", test.use = "wilcox")
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
x4 <- FindMarkers(Ferus, ident.2 = "Mock_Uninfected", ident.1 = "MERS-CoV_Infected", test.use = "wilcox")
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
x5 <- FindMarkers(Ferus, ident.2 = "Mock_Uninfected", ident.1 = "dcCoV-ACN4_Bystander", test.use = "wilcox")
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
x6 <- FindMarkers(Ferus, ident.2 = "Mock_Uninfected", ident.1 = "dcCoV-ACN4_Infected", test.use = "wilcox")
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

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/DEG_mock_v_MERS_by_inf.pdf",
    width=15, height=7)
p3|p4
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/DEG_mock_v_ACN4_by_inf.pdf",
    width=15, height=7)
p5|p6
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/DEG_MERS_inf_vs_ACN4_inf.pdf",
    width=15, height=7)
p4|p6
dev.off()

#DEG Analysis################################################################################################################################################

#Ferus <- PrepSCTFindMarkers(Ferus)

# Find differentially expressed features between Mock and MERS-CoV Bystanders
Uninf_vs_MERS_bys.markers_wilcox_df <- FindMarkers(Ferus, ident.1 = "MERS-CoV_Bystander",ident.2 = "Mock_Uninfected",test.use = "wilcox",logfc.threshold = 0.5)

# Find differentially expressed features between Mock and MERS-CoV Infected
Uninf_vs_MERS_inf.markers_wilcox_df <- FindMarkers(Ferus, ident.2 = "Mock_Uninfected", ident.1 = "MERS-CoV_Infected", test.use = "wilcox", logfc.threshold = 0.5)

# Find differentially expressed features between Mock and ACN4 Bystanders
Uninf_vs_ACN4_bys.markers_wilcox_df <- FindMarkers(Ferus, ident.2 = "Mock_Uninfected", ident.1 = "dcCoV-ACN4_Bystander", test.use = "wilcox", logfc.threshold = 0.5)

# Find differentially expressed features between Mock and ACN4 Infected
Uninf_vs_ACN4_inf.markers_wilcox_df <- FindMarkers(Ferus, ident.2 = "Mock_Uninfected", ident.1 = "dcCoV-ACN4_Infected", test.use = "wilcox", logfc.threshold = 0.5)




#save the lists
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Ferus_both")
write.csv(Uninf_vs_MERS_bys.markers_wilcox_df,"Camel_Uninf_vs_MERS_bys_markers_wilcox.csv")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Ferus_both")
write.csv(Uninf_vs_MERS_inf.markers_wilcox_df,"Camel_Uninf_vs_MERS_inf_markers_wilcox.csv")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Ferus_both")
write.csv(Uninf_vs_ACN4_bys.markers_wilcox_df,"Camel_Uninf_vs_ACN4_bys_markers_wilcox.csv")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Ferus_both")
write.csv(Uninf_vs_ACN4_inf.markers_wilcox_df,"Camel_Uninf_vs_ACN4_inf_markers_wilcox.csv")


#####################################
#up-and-down regulated mers infected
Uninf_vs_MERS_inf.markers_down <- Uninf_vs_MERS_inf.markers_wilcox_df %>%
  filter(Uninf_vs_MERS_inf.markers_wilcox_df$avg_log2FC < 0) %>%
  rownames()

Uninf_vs_MERS_inf.markers_up <- Uninf_vs_MERS_inf.markers_wilcox_df %>%
  filter(Uninf_vs_MERS_inf.markers_wilcox_df$avg_log2FC > 0) %>%
  rownames()

#count how many are LOC unannotated camel genes in up-and down regulated genes
q.data <- data.frame(string=Uninf_vs_MERS_inf.markers_down, stringsAsFactors = F)
q.data$number.of.LOC <- str_count(q.data$string, "LOC")
table(q.data['number.of.LOC'])

q.data <- data.frame(string=Uninf_vs_MERS_inf.markers_up, stringsAsFactors = F)
q.data$number.of.LOC <- str_count(q.data$string, "LOC")
table(q.data['number.of.LOC'])


#####################################
#up-and-down regulated acn4 infected
Uninf_vs_ACN4_inf.markers_down <- Uninf_vs_ACN4_inf.markers_wilcox_df %>%
  filter(Uninf_vs_ACN4_inf.markers_wilcox_df$avg_log2FC < 0) %>%
  rownames()

Uninf_vs_ACN4_inf.markers_up <- Uninf_vs_ACN4_inf.markers_wilcox_df %>%
  filter(Uninf_vs_ACN4_inf.markers_wilcox_df$avg_log2FC > 0) %>%
  rownames()

#count how many are LOC unannotated camel genes in up-and down regulated genes
q.data <- data.frame(string=Uninf_vs_ACN4_inf.markers_down, stringsAsFactors = F)
q.data$number.of.LOC <- str_count(q.data$string, "LOC")
table(q.data['number.of.LOC'])

q.data <- data.frame(string=Uninf_vs_ACN4_inf.markers_up, stringsAsFactors = F)
q.data$number.of.LOC <- str_count(q.data$string, "LOC")
table(q.data['number.of.LOC'])

#together
Uninf_vs_MERS_inf.markers_wilcox <- Uninf_vs_MERS_inf.markers_wilcox_df %>%
  rownames()
Uninf_vs_ACN4_inf.markers_wilcox <- Uninf_vs_ACN4_inf.markers_wilcox_df %>%
  rownames()

#Venn plots###############################################################################

All_genes <- list(`Camel MERS-CoV downregulated` = Uninf_vs_MERS_inf.markers_down,
                  `Camel MERS-CoV upregulated` = Uninf_vs_MERS_inf.markers_up,
                  `Camel dcCoV-ACN4 upregulated` = Uninf_vs_ACN4_inf.markers_up,
                  `Camel dcCoV-ACN4 downregulated` = Uninf_vs_ACN4_inf.markers_down)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_camel_up_down.pdf",
    width=10, height=9)
ggvenn(All_genes, 
       fill_color = c("#A3CCCC", "#00A6A6","#41521F","#A89F68", "blue"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)
dev.off()

features = c("HSPA6", "RRAD", "RSAD2", "SDCBP2", "DNAJB1", "DNAJA4")
features2 = c("CD74", "LRRQ1", "AK7", "BAIAP3", "TPPP3", "LOC102512064", "LOC102520356", "FOXJ1", "DYNLRB2", "NME9")

Ferus_MERS <- subset(x = Ferus, idents = c("Mock_Uninfected", "MERS-CoV_Infected", "dcCoV-ACN4_Infected"))

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/dotplot_shared_genes_both_virus.pdf",
    width=5, height=8)
DotPlot(Ferus_MERS, features = features, scale.by = "radius", scale = T, cols = c("orange", "darkblue"))+ coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                                                                                                                      legend.position = "top")
dev.off()


DotPlot(Ferus_MERS, features = features2, scale.by = "radius", scale = T, cols = c("orange", "darkblue"))+ coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                                                                                                                      legend.position = "top")

Uninf_vs_MERS_inf.markers_up

#Translate gene lists obtained by DE into the human equivalent to get GO terms with human database##########################################

Uninf_vs_ACN4_inf.markers_down_list <- bitr(Uninf_vs_ACN4_inf.markers_down, fromType = "SYMBOL",
                                              toType = c("ENTREZID"),
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

#GO Annotation##################################################################################################################################

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
## feature 3: decreasing orde
geneList_MERS = sort(geneList_MERS, decreasing = TRUE)


#filtering and ordering the enrichgo results to find out how many terms were significantly enriched
df <- dplyr::filter(enrichGO@result, p.adjust < 0.05)



#Making the enrichGO result for barplot
enrichGO <- enrichGO(gene = Uninf_vs_MERS_inf.markers_wilcox_list$ENTREZID,
                     OrgDb = "org.Hs.eg.db",
                     keyType = "ENTREZID",
                     ont = "BP",
                     qvalueCutoff = 0.05,
                     readable = T
)


#Making the enrichGO result for barplot
enrichGO_down <- enrichGO(gene = Uninf_vs_MERS_inf.markers_down_list$ENTREZID,
                     OrgDb = "org.Hs.eg.db",
                     keyType = "ENTREZID",
                     ont = "BP",
                     qvalueCutoff = 0.05,
                     readable = T
)


#Making the enrichGO result for barplot
enrichGO_up <- enrichGO(gene = Uninf_vs_MERS_inf.markers_up_list$ENTREZID,
                     OrgDb = "org.Hs.eg.db",
                     keyType = "ENTREZID",
                     ont = "BP",
                     qvalueCutoff = 0.05,
                     readable = T
)

barplot(
  enrichGO_up,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 5,
  font.size = 12,
  title = "",
  label_format = 30
)

barplot(
  enrichGO_down,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 5,
  font.size = 12,
  title = "",
  label_format = 30
)


enrichGO_MERS <- enrichGO@result
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Enrich_terms")
write.csv(enrichGO_MERS,"enrichGO_camel_MERS_inf_terms.csv")



#making the barplot
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/barplot_overview_camel_MERS.pdf",
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

?group

ggo <- groupGO(gene     = Uninf_vs_MERS_inf.markers_wilcox_list$ENTREZID,
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 3,
               readable = TRUE,
               qvalueCutoff = 0.05)


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
write.csv(df_enr,"GGO_camel_MERS_inf_terms.csv")


ggo_immune@result <- filter(ggo_immune@result, grepl('immun', ggo_immune@result$Description))
ggo_immune@result <- ggo_immune@result[order(-ggo_immune@result$Count),]


ggo_viral@result <- filter(ggo_viral@result, grepl('vir', ggo_viral@result$Description))
ggo_viral@result <- ggo_viral@result[order(-ggo_viral@result$Count),]




df1 <- ggo_immune@result
df2 <- ggo_viral@result

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Enrich_terms")
write.csv(df1,"GGO_camel_MERS_inf_immune_terms.csv")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Enrich_terms")
write.csv(df2,"GGO_camel_MERS_inf_viral_terms.csv")


ggo@result <- dplyr::filter(ggo@result, Count > 50)
ggo@result <- ggo@result[order(-ggo@result$Count),]



heatplot(ggo,
         foldChange=geneList_MERS) + ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red")

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/heatplot_Uninf_MERS_inf.pdf",
    width=11, height=15)
heatplot(ggo_5,
         foldChange=geneList_MERS) + ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red")
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/cnet_viral_Uninf_MERS_inf.pdf",
    width=10, height=5)
cnetplot(df, foldChange=geneList_MERS, colorEdge = TRUE, circular = TRUE, layout = 'gem', repel = T)
dev.off()

df <- enrichGO %>%
  pairwise_termsim()

filter(df@result, grepl('viral|virus', df@result$Description)) %>%
  cnetplot(foldChange=geneList_MERS, colorEdge = TRUE, circular = TRUE, layout = 'gem', repel = T)


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/cnet_immune_Uninf_MERS_inf.pdf",
    width=10, height=5)
cnetplot(ggo_immune, foldChange=geneList_MERS, colorEdge = TRUE, circular = TRUE, layout = 'gem', repel = T)
dev.off()

#################################################################################################################################################
#ACN4
#Creating geneList file of genes as rownames and FC as a value
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





#Making the enrichGO result for barplot
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


barplot(
  enrichGO_down,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 5,
  font.size = 12,
  title = "",
  label_format = 30
)


barplot(
  enrichGO_up,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 5,
  font.size = 12,
  title = "",
  label_format = 30
)

enrichGO_ACN4 <- enrichGO@result
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Enrich_terms")
write.csv(enrichGO_ACN4,"enrichGO_camel_ACN4_inf_terms.csv")


#filtering and ordering the enrichgo results to find out how many terms were significantly enriched
df <- dplyr::filter(enrichGO@result, p.adjust < 0.05)



#making the barplot
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/barplot_overview_camel_ACN4.pdf",
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


#filtering and ordering the enrichgo results to find out how many terms were significantly enriched
df <- dplyr::filter(enrichGO@result, p.adjust > 0.05)






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


df <- ggo_ACN4@result
df_enr <- filter(df, Count >= 1)

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Enrich_terms")
write.csv(df_enr,"GGO_llama_ACN4_inf_terms.csv")



ggo@result <- dplyr::filter(ggo@result, Count > 50)
ggo@result <- ggo@result[order(-ggo@result$Count),]



ggo_immune@result <- filter(ggo_immune@result, grepl('immun', ggo_immune@result$Description))
ggo_immune@result <- ggo_immune@result[order(-ggo_immune@result$Count),]


ggo_viral@result <- filter(ggo_viral@result, grepl('vir', ggo_viral@result$Description))
ggo_viral@result <- ggo_viral@result[order(-ggo_viral@result$Count),]

df1 <- ggo_immune@result
df2 <- ggo_viral@result

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Enrich_terms")
write.csv(df1,"GGO_llama_ACN4_inf_immune_terms.csv")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Enrich_terms")
write.csv(df2,"GGO_llama_ACN4_inf_viral_terms.csv")



#Overview, ordered and filtered for > 50 gene counts GO terms
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/barplot_overview_camel_ACN4.pdf",
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

heatplot(ggo,
         foldChange=geneList_ACN4) + coord_flip()

ggo2 <- pairwise_termsim(ggo)


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/cnet_viral_Uninf_ACN4_inf.pdf",
    width=10, height=5)
cnetplot(ggo_viral, foldChange=geneList_ACN4, colorEdge = TRUE, circular = TRUE, layout = 'gem', repel = T)
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/cnet_immune_Uninf_ACN4_inf.pdf",
    width=10, height=5)
cnetplot(ggo_immune, foldChange=geneList_ACN4, colorEdge = TRUE, circular = TRUE, layout = 'gem', repel = T)
dev.off()


##################################################################################################################################################

# Difference between Celltypes in Status (Infected, Uninfected, Bystander) for each Virus (dcCoV-ACN4, MERS-CoV)

##################################################################################################################################################
DefaultAssay(Ferus) <- "RNA"
#Increase memory limit
memory.limit(24000)
# Find differentially expressed features between Infected and uninfected for one Treatment group
# MERS-CoV
Ferus@meta.data
Idents(Ferus) <- "celltype"

#Make a fused new identity of each treatment and status combination
Ferus$cell.treat.status <- paste(Idents(Ferus), Ferus$Treat, Ferus$Status, sep = "_")

#make it the new identity
Idents(Ferus) <- "cell.treat.status"
#check if levels are correct
levels(Ferus)


#Volcano Plots
#Secretory#######################################################################################################################################
# Find differentially expressed features between Mock and MERS-CoV Bystanders
x7 <- FindMarkers(Ferus, ident.2 = "Secretory_Mock_Uninfected", ident.1 = "Secretory_MERS-CoV_Bystander", test.use = "wilcox")
p7 <- EnhancedVolcano(x7,
                      title = 'Secretory Uninfected versus Secretory MERS-CoV Bystanders',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x7),
                      col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 0.5,
                      pCutoff = 0.05,
                      colAlpha = 1)


# Find differentially expressed features between Mock and MERS-CoV Infected
x8 <- FindMarkers(Ferus, ident.2 = "Secretory_Mock_Uninfected", ident.1 = "Secretory_MERS-CoV_Infected", test.use = "wilcox")
p8 <- EnhancedVolcano(x8,
                      title = 'Secretory Uninfected versus Secretory MERS-CoV Infected',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x8),
                      col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 0.5,
                      pCutoff = 0.05,
                      colAlpha = 1)

# Find differentially expressed features between Mock and ACN4 Bystanders
x9 <- FindMarkers(Ferus, ident.2 = "Secretory_Mock_Uninfected", ident.1 = "Secretory_dcCoV-ACN4_Bystander", test.use = "wilcox")
p9 <- EnhancedVolcano(x9,
                      title = 'Secretory Uninfected versus Secretory dcCoV-ACN4 Bystanders',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x9),
                      col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 0.5,
                      pCutoff = 0.05,
                      colAlpha = 1)

# Find differentially expressed features between Mock and ACN4 Infected
x10 <- FindMarkers(Ferus, ident.12 = "Secretory_Mock_Uninfected", ident.1 = "Secretory_dcCoV-ACN4_Infected", test.use = "wilcox")
p10 <- EnhancedVolcano(x10,
                      title = 'Secretory Uninfected versus Secretory dcCoV-ACN4 Infected',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x10),
                      col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 0.5,
                      pCutoff = 0.05,
                      colAlpha = 1)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/DEG_secretory_mock_v_MERS_by_inf.pdf",
    width=15, height=7)
p7|p8
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/DEG_secretory_mock_v_ACN4_by_inf.pdf",
    width=15, height=7)
p9|p10
dev.off()

#Ciliated#######################################################################################################################################
# Find differentially expressed features between Mock and MERS-CoV Bystanders
x11 <- FindMarkers(Ferus, ident.1 = "Ciliated_Mock_Uninfected", ident.2 = "Ciliated_MERS-CoV_Bystander", test.use = "wilcox")
p11 <- EnhancedVolcano(x11,
                      title = 'Ciliated Uninfected versus Ciliated MERS-CoV Bystanders',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x11),
                      col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 1.5,
                      pCutoff = 0.05,
                      colAlpha = 1)


# Find differentially expressed features between Mock and MERS-CoV Infected
x12 <- FindMarkers(Ferus, ident.1 = "Ciliated_Mock_Uninfected", ident.2 = "Ciliated_MERS-CoV_Infected", test.use = "wilcox")
p12 <- EnhancedVolcano(x12,
                      title = 'Ciliated Uninfected versus Ciliated MERS-CoV Infected',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x12),
                      col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 1.5,
                      pCutoff = 0.05,
                      colAlpha = 1)

# Find differentially expressed features between Mock and ACN4 Bystanders
x13 <- FindMarkers(Ferus, ident.1 = "Ciliated_Mock_Uninfected", ident.2 = "Ciliated_dcCoV-ACN4_Bystander", test.use = "wilcox")
p13 <- EnhancedVolcano(x13,
                      title = 'Ciliated Uninfected versus Ciliated dcCoV-ACN4 Bystanders',
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      pointSize = 3.0,
                      labSize = 6.0,
                      lab = rownames(x13),
                      col=c('black', 'black', 'black', 'red3'),
                      FCcutoff = 1.5,
                      pCutoff = 0.05,
                      colAlpha = 1)

# Find differentially expressed features between Mock and ACN4 Infected
x14 <- FindMarkers(Ferus, ident.1 = "Ciliated_Mock_Uninfected", ident.2 = "Ciliated_dcCoV-ACN4_Infected", test.use = "wilcox")
p14 <- EnhancedVolcano(x14,
                       title = 'Ciliated Uninfected versus Ciliated dcCoV-ACN4 Infected',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x14),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 1.5,
                       pCutoff = 0.05,
                       colAlpha = 1)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/DEG_ciliated_mock_v_MERS_by_inf.pdf",
    width=15, height=7)
p11|p12
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/DEG_ciliated_mock_v_ACN4_by_inf.pdf",
    width=15, height=7)
p13|p14
dev.off()


#Cluster 9#######################################################################################################################################
# Find differentially expressed features between Mock and MERS-CoV Bystanders
x15 <- FindMarkers(Ferus, ident.1 = "Cluster 9_Mock_Uninfected", ident.2 = "Cluster 9_MERS-CoV_Bystander", test.use = "wilcox")
p15 <- EnhancedVolcano(x15,
                       title = 'Cluster 9 Uninfected versus Cluster 9 MERS-CoV Bystanders',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x15),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 1.5,
                       pCutoff = 0.05,
                       colAlpha = 1)


# Find differentially expressed features between Mock and MERS-CoV Infected
x16 <- FindMarkers(Ferus, ident.1 = "Cluster 9_Mock_Uninfected", ident.2 = "Cluster 9_MERS-CoV_Infected", test.use = "wilcox")
p16 <- EnhancedVolcano(x16,
                       title = 'Cluster 9 Uninfected versus Cluster 9 MERS-CoV Infected',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x16),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 1.5,
                       pCutoff = 0.05,
                       colAlpha = 1)

# Find differentially expressed features between Mock and ACN4 Bystanders
x17 <- FindMarkers(Ferus, ident.1 = "Cluster 9_Mock_Uninfected", ident.2 = "Cluster 9_dcCoV-ACN4_Bystander", test.use = "wilcox")
p17 <- EnhancedVolcano(x17,
                       title = 'Cluster 9 Uninfected versus Cluster 9 dcCoV-ACN4 Bystanders',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x17),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 1.5,
                       pCutoff = 0.05,
                       colAlpha = 1)

# Find differentially expressed features between Mock and ACN4 Infected
x18 <- FindMarkers(Ferus, ident.1 = "Cluster 9_Mock_Uninfected", ident.2 = "Cluster 9_dcCoV-ACN4_Infected", test.use = "wilcox")
p18 <- EnhancedVolcano(x18,
                       title = 'Cluster 9 Uninfected versus Cluster 9 dcCoV-ACN4 Infected',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x18),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 1.5,
                       pCutoff = 0.05,
                       colAlpha = 1)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/DEG_Cluster 9_mock_v_MERS_by_inf.pdf",
    width=15, height=7)
p15|p16
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/DEG_Cluster 9_mock_v_ACN4_by_inf.pdf",
    width=15, height=7)
p17
dev.off()

#Cluster 10#######################################################################################################################################
# Find differentially expressed features between Mock and MERS-CoV Bystanders
x19 <- FindMarkers(Ferus, ident.1 = "Cluster 10_Mock_Uninfected", ident.2 = "Cluster 10_MERS-CoV_Bystander", test.use = "wilcox")
p19 <- EnhancedVolcano(x19,
                       title = 'Cluster 10 Uninfected versus Cluster 10 MERS-CoV Bystanders',
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
x20 <- FindMarkers(Ferus, ident.1 = "Cluster 10_Mock_Uninfected", ident.2 = "Cluster 10_MERS-CoV_Infected", test.use = "wilcox")
p20 <- EnhancedVolcano(x20,
                       title = 'Cluster 10 Uninfected versus Cluster 10 MERS-CoV Infected',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x20),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 1.5,
                       pCutoff = 0.05,
                       colAlpha = 1)

# Find differentially expressed features between Mock and ACN4 Bystanders
x21 <- FindMarkers(Ferus, ident.1 = "Cluster 10_Mock_Uninfected", ident.2 = "Cluster 10_dcCoV-ACN4_Bystander", test.use = "wilcox")
p21 <- EnhancedVolcano(x21,
                       title = 'Cluster 10 Uninfected versus Cluster 10 dcCoV-ACN4 Bystanders',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x21),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 1.5,
                       pCutoff = 0.05,
                       colAlpha = 1)

# Find differentially expressed features between Mock and ACN4 Infected
x22 <- FindMarkers(Ferus, ident.1 = "Cluster 10_Mock_Uninfected", ident.2 = "Cluster 10_dcCoV-ACN4_Infected", test.use = "wilcox")
p22 <- EnhancedVolcano(x22,
                       title = 'Cluster 10 Uninfected versus Cluster 10 dcCoV-ACN4 Infected',
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pointSize = 3.0,
                       labSize = 6.0,
                       lab = rownames(x22),
                       col=c('black', 'black', 'black', 'red3'),
                       FCcutoff = 1.5,
                       pCutoff = 0.05,
                       colAlpha = 1)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/DEG_Cluster 10_mock_v_MERS_by_inf.pdf",
    width=15, height=7)
p20
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/DEG_Cluster 10_mock_v_ACN4_by_inf.pdf",
    width=15, height=7)
p21|p22
dev.off()





#Find DEG features#############################################################################################################################################################
#Mock vs Infected
Secretory_Uninf_vs_MERS_inf.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Secretory_Mock_Uninfected", ident.1 = "Secretory_MERS-CoV_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Secretory_Uninf_vs_MERS_inf.markers_wilcox$SYMBOL <- rownames(Secretory_Uninf_vs_MERS_inf.markers_wilcox)
Secretory_Uninf_vs_ACN4_inf.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Secretory_Mock_Uninfected", ident.1 = "Secretory_dcCoV-ACN4_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Secretory_Uninf_vs_ACN4_inf.markers_wilcox$SYMBOL <- rownames(Secretory_Uninf_vs_ACN4_inf.markers_wilcox)

Ciliated_Uninf_vs_MERS_inf.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Ciliated_Mock_Uninfected", ident.1 = "Ciliated_MERS-CoV_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Ciliated_Uninf_vs_MERS_inf.markers_wilcox$SYMBOL <- rownames(Ciliated_Uninf_vs_MERS_inf.markers_wilcox)
Ciliated_Uninf_vs_ACN4_inf.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Ciliated_Mock_Uninfected", ident.1 = "Ciliated_dcCoV-ACN4_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Ciliated_Uninf_vs_ACN4_inf.markers_wilcox$SYMBOL <- rownames(Ciliated_Uninf_vs_ACN4_inf.markers_wilcox)

Club_Uninf_vs_MERS_inf.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Club_Mock_Uninfected", ident.1 = "Club_MERS-CoV_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Club_Uninf_vs_MERS_inf.markers_wilcox$SYMBOL <- rownames(Club_Uninf_vs_MERS_inf.markers_wilcox)
Club_Uninf_vs_ACN4_inf.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Club_Mock_Uninfected", ident.1 = "Club_dcCoV-ACN4_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Club_Uninf_vs_ACN4_inf.markers_wilcox$SYMBOL <- rownames(Club_Uninf_vs_ACN4_inf.markers_wilcox)

Basal_Uninf_vs_MERS_inf.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Basal_Mock_Uninfected", ident.1 = "Basal_MERS-CoV_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Basal_Uninf_vs_MERS_inf.markers_wilcox$SYMBOL <- rownames(Basal_Uninf_vs_MERS_inf.markers_wilcox)
Basal_Uninf_vs_ACN4_inf.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Basal_Mock_Uninfected", ident.1 = "Basal_dcCoV-ACN4_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Basal_Uninf_vs_ACN4_inf.markers_wilcox$SYMBOL <- rownames(Basal_Uninf_vs_ACN4_inf.markers_wilcox)

Camel_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Camel Cluster 1_Mock_Uninfected", ident.1 = "Camel Cluster 1_MERS-CoV_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Camel_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox$SYMBOL <- rownames(Camel_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox)
Camel_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Camel Cluster 1_Mock_Uninfected", ident.1 = "Camel Cluster 1_dcCoV-ACN4_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Camel_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox$SYMBOL <- rownames(Camel_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox)

Camel_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Camel Cluster 2_Mock_Uninfected", ident.1 = "Camel Cluster 2_MERS-CoV_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Camel_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox$SYMBOL <- rownames(Camel_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox)
Camel_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Camel Cluster 2_Mock_Uninfected", ident.1 = "Camel Cluster 2_dcCoV-ACN4_Infected", test.use = "wilcox", logfc.threshold = 0.5)
Camel_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox$SYMBOL <- rownames(Camel_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox)


#Mock vs Bystander
Secretory_Uninf_vs_MERS_bys.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Secretory_Mock_Uninfected", ident.1 = "Secretory_MERS-CoV_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Secretory_Uninf_vs_MERS_bys.markers_wilcox$SYMBOL <- rownames(Secretory_Uninf_vs_MERS_bys.markers_wilcox)
Secretory_Uninf_vs_ACN4_bys.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Secretory_Mock_Uninfected", ident.1 = "Secretory_dcCoV-ACN4_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Secretory_Uninf_vs_ACN4_bys.markers_wilcox$SYMBOL <- rownames(Secretory_Uninf_vs_ACN4_bys.markers_wilcox)

Ciliated_Uninf_vs_MERS_bys.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Ciliated_Mock_Uninfected", ident.1 = "Ciliated_MERS-CoV_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Ciliated_Uninf_vs_MERS_bys.markers_wilcox$SYMBOL <- rownames(Ciliated_Uninf_vs_MERS_bys.markers_wilcox)
Ciliated_Uninf_vs_ACN4_bys.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Ciliated_Mock_Uninfected", ident.1 = "Ciliated_dcCoV-ACN4_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Ciliated_Uninf_vs_ACN4_bys.markers_wilcox$SYMBOL <- rownames(Ciliated_Uninf_vs_ACN4_bys.markers_wilcox)

Club_Uninf_vs_MERS_bys.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Club_Mock_Uninfected", ident.1 = "Club_MERS-CoV_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Club_Uninf_vs_MERS_bys.markers_wilcox$SYMBOL <- rownames(Club_Uninf_vs_MERS_bys.markers_wilcox)
Club_Uninf_vs_ACN4_bys.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Club_Mock_Uninfected", ident.1 = "Club_dcCoV-ACN4_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Club_Uninf_vs_ACN4_bys.markers_wilcox$SYMBOL <- rownames(Club_Uninf_vs_ACN4_bys.markers_wilcox)

Basal_Uninf_vs_MERS_bys.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Basal_Mock_Uninfected", ident.1 = "Basal_MERS-CoV_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Basal_Uninf_vs_MERS_bys.markers_wilcox$SYMBOL <- rownames(Basal_Uninf_vs_MERS_bys.markers_wilcox)
Basal_Uninf_vs_ACN4_bys.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Basal_Mock_Uninfected", ident.1 = "Basal_dcCoV-ACN4_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Basal_Uninf_vs_ACN4_bys.markers_wilcox$SYMBOL <- rownames(Basal_Uninf_vs_ACN4_bys.markers_wilcox)

Camel_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Camel Cluster 1_Mock_Uninfected", ident.1 = "Camel Cluster 1_MERS-CoV_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Camel_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox$SYMBOL <- rownames(Camel_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox)
Camel_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Camel Cluster 1_Mock_Uninfected", ident.1 = "Camel Cluster 1_dcCoV-ACN4_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Camel_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox$SYMBOL <- rownames(Camel_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox)

Camel_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Camel Cluster 2_Mock_Uninfected", ident.1 = "Camel Cluster 2_MERS-CoV_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Camel_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox$SYMBOL <- rownames(Camel_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox)
Camel_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox <- FindMarkers(Ferus, ident.2 = "Camel Cluster 2_Mock_Uninfected", ident.1 = "Camel Cluster 2_dcCoV-ACN4_Bystander", test.use = "wilcox", logfc.threshold = 0.5)
Camel_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox$SYMBOL <- rownames(Camel_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox)




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




Camel_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox %>% 
  dplyr::filter(Camel_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox$avg_log2FC > 0) %>%
  count()
Camel_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox %>% 
  dplyr::filter(Camel_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox$avg_log2FC < 0) %>%
  count()

Camel_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox %>% 
  dplyr::filter(Camel_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox$avg_log2FC > 0) %>%
  count()
Camel_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox %>% 
  dplyr::filter(Camel_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox$avg_log2FC < 0) %>%
  count()

Camel_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox %>% 
  dplyr::filter(Camel_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox$avg_log2FC > 0) %>%
  count()
Camel_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox %>% 
  dplyr::filter(Camel_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox$avg_log2FC < 0) %>%
  count()

Camel_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox %>% 
  dplyr::filter(Camel_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox$avg_log2FC > 0) %>%
  count()
Camel_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox %>% 
  dplyr::filter(Camel_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox$avg_log2FC < 0) %>%
  count()




Camel_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox %>% 
  dplyr::filter(Camel_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox$avg_log2FC > 0) %>%
  count()
Camel_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox %>% 
  dplyr::filter(Camel_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox$avg_log2FC < 0) %>%
  count()

Camel_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox %>% 
  dplyr::filter(Camel_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox$avg_log2FC > 0) %>%
  count()
Camel_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox %>% 
  dplyr::filter(Camel_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox$avg_log2FC < 0) %>%
  count()

Camel_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox %>% 
  dplyr::filter(Camel_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox$avg_log2FC > 0) %>%
  count()
Camel_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox %>% 
  dplyr::filter(Camel_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox$avg_log2FC < 0) %>%
  count()

Camel_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox %>% 
  dplyr::filter(Camel_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox$avg_log2FC > 0) %>%
  count()
Camel_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox %>% 
  dplyr::filter(Camel_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox$avg_log2FC < 0) %>%
  count()


#save the lists
#Infected
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Ferus_both")
write.csv(Secretory_Uninf_vs_MERS_inf.markers_wilcox,"Secretory_Camel_MERS_Infected.csv")
write.csv(Secretory_Uninf_vs_ACN4_inf.markers_wilcox,"Secretory_Camel_ACN4_Infected.csv")

write.csv(Ciliated_Uninf_vs_MERS_inf.markers_wilcox,"Ciliated_Camel_MERS_Infected.csv")
write.csv(Ciliated_Uninf_vs_ACN4_inf.markers_wilcox,"Ciliated_Camel_ACN4_Infected.csv")

write.csv(Club_Uninf_vs_MERS_inf.markers_wilcox,"Club_Camel_MERS_Infected.csv")
write.csv(Club_Uninf_vs_ACN4_inf.markers_wilcox,"Club_Camel_ACN4_Infected.csv")

write.csv(Basal_Uninf_vs_MERS_inf.markers_wilcox,"Basal_Camel_MERS_Infected.csv")
write.csv(Basal_Uninf_vs_ACN4_inf.markers_wilcox,"Basal_Camel_ACN4_Infected.csv")

write.csv(Camel_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox,"Cluster1_Camel_MERS_Infected.csv")
write.csv(Camel_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox,"Cluster1_Camel_ACN4_Infected.csv")

write.csv(Camel_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox,"Cluster2_Camel_MERS_Infected.csv")
write.csv(Camel_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox,"Cluster2_Camel_ACN4_Infected.csv")

#Bystander
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Ferus_both")
write.csv(Secretory_Uninf_vs_MERS_bys.markers_wilcox,"Secretory_Camel_MERS_Bystander.csv")
write.csv(Secretory_Uninf_vs_ACN4_bys.markers_wilcox,"Secretory_Camel_ACN4_Bystander.csv")

write.csv(Ciliated_Uninf_vs_MERS_bys.markers_wilcox,"Ciliated_Camel_MERS_Bystander.csv")
write.csv(Ciliated_Uninf_vs_ACN4_bys.markers_wilcox,"Ciliated_Camel_ACN4_Bystander.csv")

write.csv(Club_Uninf_vs_MERS_bys.markers_wilcox,"Club_Camel_MERS_Bystander.csv")
write.csv(Club_Uninf_vs_ACN4_bys.markers_wilcox,"Club_Camel_ACN4_Bystander.csv")

write.csv(Basal_Uninf_vs_MERS_bys.markers_wilcox,"Basal_Camel_MERS_Bystander.csv")
write.csv(Basal_Uninf_vs_ACN4_bys.markers_wilcox,"Basal_Camel_ACN4_Bystander.csv")

write.csv(Camel_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox,"Cluster1_Camel_MERS_Bystander.csv")
write.csv(Camel_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox,"Cluster1_Camel_ACN4_Bystander.csv")

write.csv(Camel_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox,"Cluster2_Camel_MERS_Bystander.csv")
write.csv(Camel_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox,"Cluster2_Camel_ACN4_Bystander.csv")


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

Camel_Cluster_1_Uninf_vs_MERS_inf_rownames <- rownames(Camel_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox)
Camel_Cluster_1_Uninf_vs_ACN4_inf_rownames <- rownames(Camel_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox)

Camel_Cluster_2_Uninf_vs_MERS_inf_rownames <- rownames(Camel_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox)
Camel_Cluster_2_Uninf_vs_ACN4_inf_rownames <- rownames(Camel_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox)

#Bystander
Secretory_Uninf_vs_MERS_bys_rownames <- rownames(Secretory_Uninf_vs_MERS_bys.markers_wilcox)
Secretory_Uninf_vs_ACN4_bys_rownames <- rownames(Secretory_Uninf_vs_ACN4_bys.markers_wilcox)

Club_Uninf_vs_MERS_bys_rownames <- rownames(Club_Uninf_vs_MERS_bys.markers_wilcox)
Club_Uninf_vs_ACN4_bys_rownames <- rownames(Club_Uninf_vs_ACN4_bys.markers_wilcox)

Basal_Uninf_vs_MERS_bys_rownames <- rownames(Basal_Uninf_vs_MERS_bys.markers_wilcox)
Basal_Uninf_vs_ACN4_bys_rownames <- rownames(Basal_Uninf_vs_ACN4_bys.markers_wilcox)

Ciliated_Uninf_vs_MERS_bys_rownames <- rownames(Ciliated_Uninf_vs_MERS_bys.markers_wilcox)
Ciliated_Uninf_vs_ACN4_bys_rownames <- rownames(Ciliated_Uninf_vs_ACN4_bys.markers_wilcox)

Camel_Cluster_1_Uninf_vs_MERS_bys_rownames <- rownames(Camel_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox)
Camel_Cluster_1_Uninf_vs_ACN4_bys_rownames <- rownames(Camel_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox)

Camel_Cluster_2_Uninf_vs_MERS_bys_rownames <- rownames(Camel_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox)
Camel_Cluster_2_Uninf_vs_ACN4_bys_rownames <- rownames(Camel_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox)

#############################################################################################################################################

#Translate gene lists obtained by DE into the human equivalent to get GO terms with human database

#############################################################################################################################################


Secretory_Camel_MERS_Infected <- bitr(Secretory_Uninf_vs_MERS_inf_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)
Secretory_Camel_MERS_Bystander <- bitr(Secretory_Uninf_vs_MERS_bys_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)
Secretory_Camel_ACN4_Infected <- bitr(Secretory_Uninf_vs_ACN4_inf_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)
Secretory_Camel_ACN4_Bystander <- bitr(Secretory_Uninf_vs_ACN4_bys_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)



Ciliated_Camel_MERS_Infected <- bitr(Ciliated_Uninf_vs_MERS_inf_rownames, fromType = "SYMBOL",
                                          toType = "ENTREZID",
                                          OrgDb = org.Hs.eg.db)
Ciliated_Camel_MERS_Bystander <- bitr(Ciliated_Uninf_vs_MERS_bys_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)
Ciliated_Camel_ACN4_Infected <- bitr(Ciliated_Uninf_vs_ACN4_inf_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)
Ciliated_Camel_ACN4_Bystander <- bitr(Ciliated_Uninf_vs_ACN4_bys_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)


Club_Camel_MERS_Infected <- bitr(Club_Uninf_vs_MERS_inf_rownames, fromType = "SYMBOL",
                                          toType = "ENTREZID",
                                          OrgDb = org.Hs.eg.db)
Club_Camel_MERS_Bystander <- bitr(Club_Uninf_vs_MERS_bys_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)
Club_Camel_ACN4_Infected <- bitr(Club_Uninf_vs_ACN4_inf_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)
Club_Camel_ACN4_Bystander <- bitr(Club_Uninf_vs_ACN4_bys_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)


Basal_Camel_MERS_Infected <- bitr(Basal_Uninf_vs_MERS_inf_rownames, fromType = "SYMBOL",
                                          toType = "ENTREZID",
                                          OrgDb = org.Hs.eg.db)
Basal_Camel_MERS_Bystander <- bitr(Basal_Uninf_vs_MERS_bys_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)
Basal_Camel_ACN4_Infected <- bitr(Basal_Uninf_vs_ACN4_inf_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)
Basal_Camel_ACN4_Bystander <- bitr(Basal_Uninf_vs_ACN4_bys_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)


Cluster1_Camel_MERS_Infected <- bitr(Camel_Cluster_1_Uninf_vs_MERS_inf_rownames, fromType = "SYMBOL",
                                          toType = "ENTREZID",
                                          OrgDb = org.Hs.eg.db)
Cluster1_Camel_MERS_Bystander <- bitr(Camel_Cluster_1_Uninf_vs_MERS_bys_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)
Cluster1_Camel_ACN4_Infected <- bitr(Camel_Cluster_1_Uninf_vs_ACN4_inf_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)
Cluster1_Camel_ACN4_Bystander <- bitr(Camel_Cluster_1_Uninf_vs_ACN4_bys_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)


Cluster2_Camel_MERS_Infected <- bitr(Camel_Cluster_2_Uninf_vs_MERS_inf_rownames, fromType = "SYMBOL",
                                          toType = "ENTREZID",
                                          OrgDb = org.Hs.eg.db)
Cluster2_Camel_MERS_Bystander <- bitr(Camel_Cluster_2_Uninf_vs_MERS_bys_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)
Cluster2_Camel_ACN4_Infected <- bitr(Camel_Cluster_2_Uninf_vs_ACN4_inf_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)
Cluster2_Camel_ACN4_Bystander <- bitr(Camel_Cluster_2_Uninf_vs_ACN4_bys_rownames, fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db)



#Merge the two datasets together
#Infected
Secretory_Camel_MERS_Infected <- merge(Secretory_Camel_MERS_Infected, Secretory_Uninf_vs_MERS_inf.markers_wilcox, by = "SYMBOL")
Secretory_Camel_ACN4_Infected <- merge(Secretory_Camel_ACN4_Infected, Secretory_Uninf_vs_ACN4_inf.markers_wilcox, by = "SYMBOL")

Ciliated_Camel_MERS_Infected <- merge(Ciliated_Camel_MERS_Infected, Ciliated_Uninf_vs_MERS_inf.markers_wilcox, by = "SYMBOL")
Ciliated_Camel_ACN4_Infected <- merge(Ciliated_Camel_ACN4_Infected, Ciliated_Uninf_vs_ACN4_inf.markers_wilcox, by = "SYMBOL")

Club_Camel_MERS_Infected <- merge(Club_Camel_MERS_Infected, Club_Uninf_vs_MERS_inf.markers_wilcox, by = "SYMBOL")
Club_Camel_ACN4_Infected <- merge(Club_Camel_ACN4_Infected, Club_Uninf_vs_ACN4_inf.markers_wilcox, by = "SYMBOL")

Basal_Camel_MERS_Infected <- merge(Basal_Camel_MERS_Infected, Basal_Uninf_vs_MERS_inf.markers_wilcox, by = "SYMBOL")
Basal_Camel_ACN4_Infected <- merge(Basal_Camel_ACN4_Infected, Basal_Uninf_vs_ACN4_inf.markers_wilcox, by = "SYMBOL")

Cluster1_Camel_MERS_Infected <- merge(Cluster1_Camel_MERS_Infected, Camel_Cluster_1_Uninf_vs_MERS_inf.markers_wilcox, by = "SYMBOL")
Cluster1_Camel_ACN4_Infected <- merge(Cluster1_Camel_ACN4_Infected, Camel_Cluster_1_Uninf_vs_ACN4_inf.markers_wilcox, by = "SYMBOL")

Cluster2_Camel_MERS_Infected <- merge(Cluster2_Camel_MERS_Infected, Camel_Cluster_2_Uninf_vs_MERS_inf.markers_wilcox, by = "SYMBOL")
Cluster2_Camel_ACN4_Infected <- merge(Cluster2_Camel_ACN4_Infected, Camel_Cluster_2_Uninf_vs_ACN4_inf.markers_wilcox, by = "SYMBOL")

#Bystander
Secretory_Camel_MERS_Bystander <- merge(Secretory_Camel_MERS_Bystander, Secretory_Uninf_vs_MERS_bys.markers_wilcox, by = "SYMBOL")
Secretory_Camel_ACN4_Bystander <- merge(Secretory_Camel_ACN4_Bystander, Secretory_Uninf_vs_ACN4_bys.markers_wilcox, by = "SYMBOL")

Ciliated_Camel_MERS_Bystander <- merge(Ciliated_Camel_MERS_Bystander, Ciliated_Uninf_vs_MERS_bys.markers_wilcox, by = "SYMBOL")
Ciliated_Camel_ACN4_Bystander <- merge(Ciliated_Camel_ACN4_Bystander, Ciliated_Uninf_vs_ACN4_bys.markers_wilcox, by = "SYMBOL")

Club_Camel_MERS_Bystander <- merge(Club_Camel_MERS_Bystander, Club_Uninf_vs_MERS_bys.markers_wilcox, by = "SYMBOL")
Club_Camel_ACN4_Bystander <- merge(Club_Camel_ACN4_Bystander, Club_Uninf_vs_ACN4_bys.markers_wilcox, by = "SYMBOL")

Basal_Camel_MERS_Bystander <- merge(Basal_Camel_MERS_Bystander, Basal_Uninf_vs_MERS_bys.markers_wilcox, by = "SYMBOL")
Basal_Camel_ACN4_Bystander <- merge(Basal_Camel_ACN4_Bystander, Basal_Uninf_vs_ACN4_bys.markers_wilcox, by = "SYMBOL")

Cluster1_Camel_MERS_Bystander <- merge(Cluster1_Camel_MERS_Bystander, Camel_Cluster_1_Uninf_vs_MERS_bys.markers_wilcox, by = "SYMBOL")
Cluster1_Camel_ACN4_Bystander <- merge(Cluster1_Camel_ACN4_Bystander, Camel_Cluster_1_Uninf_vs_ACN4_bys.markers_wilcox, by = "SYMBOL")

Cluster2_Camel_MERS_Bystander <- merge(Cluster2_Camel_MERS_Bystander, Camel_Cluster_2_Uninf_vs_MERS_bys.markers_wilcox, by = "SYMBOL")
Cluster2_Camel_ACN4_Bystander <- merge(Cluster2_Camel_ACN4_Bystander, Camel_Cluster_2_Uninf_vs_ACN4_bys.markers_wilcox, by = "SYMBOL")


#Dotplot########################################################################################################################################
#save the lists
#Infected
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Ferus_both")
write.csv(Secretory_Camel_MERS_Infected,"Secretory_Camel_MERS_Infected.csv")
write.csv(Secretory_Camel_ACN4_Infected,"Secretory_Camel_ACN4_Infected.csv")

write.csv(Ciliated_Camel_MERS_Infected,"Ciliated_Camel_MERS_Infected.csv")
write.csv(Ciliated_Camel_ACN4_Infected,"Ciliated_Camel_ACN4_Infected.csv")

write.csv(Club_Camel_MERS_Infected,"Club_Camel_MERS_Infected.csv")
write.csv(Club_Camel_ACN4_Infected,"Club_Camel_ACN4_Infected.csv")

write.csv(Basal_Camel_MERS_Infected,"Basal_Camel_MERS_Infected.csv")
write.csv(Basal_Camel_ACN4_Infected,"Basal_Camel_ACN4_Infected.csv")

write.csv(Cluster1_Camel_MERS_Infected,"Cluster1_Camel_MERS_Infected.csv")
write.csv(Cluster1_Camel_ACN4_Infected,"Cluster1_Camel_ACN4_Infected.csv")

write.csv(Cluster2_Camel_MERS_Infected,"Cluster2_Camel_MERS_Infected.csv")
write.csv(Cluster2_Camel_ACN4_Infected,"Cluster2_Camel_ACN4_Infected.csv")

#Bystander
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Ferus_both")
write.csv(Secretory_Camel_MERS_Bystander,"Secretory_Camel_MERS_Bystander.csv")
write.csv(Secretory_Camel_ACN4_Bystander,"Secretory_Camel_ACN4_Bystander.csv")

write.csv(Ciliated_Camel_MERS_Bystander,"Ciliated_Camel_MERS_Bystander.csv")
write.csv(Ciliated_Camel_ACN4_Bystander,"Ciliated_Camel_ACN4_Bystander.csv")

write.csv(Club_Camel_MERS_Bystander,"Club_Camel_MERS_Bystander.csv")
write.csv(Club_Camel_ACN4_Bystander,"Club_Camel_ACN4_Bystander.csv")

write.csv(Basal_Camel_MERS_Bystander,"Basal_Camel_MERS_Bystander.csv")
write.csv(Basal_Camel_ACN4_Bystander,"Basal_Camel_ACN4_Bystander.csv")

write.csv(Cluster1_Camel_MERS_Bystander,"Cluster1_Camel_MERS_Bystander.csv")
write.csv(Cluster1_Camel_ACN4_Bystander,"Cluster1_Camel_ACN4_Bystander.csv")

write.csv(Cluster2_Camel_MERS_Bystander,"Cluster2_Camel_MERS_Bystander.csv")
write.csv(Cluster2_Camel_ACN4_Bystander,"Cluster2_Camel_ACN4_Bystander.csv")


#################################################################################################################################################
# #GO Annotation
# #MERS
# GO_results_Secretory_Uninf_vs_MERS_inf_wilcox <- enrichGO(gene = Ciliated_Uninf_vs_MERS_inf.markers_wilcox_human_gene_names,
#                                                 OrgDb = "org.Hs.eg.db",
#                                                 keyType = "SYMBOL",
#                                                 ont = "BP")
# 
# GO_results_Ciliated_Uninf_vs_MERS_inf_wilcox <- enrichGO(gene = Cluster_9_Uninf_vs_MERS_inf.markers_wilcox_human_gene_names,
#                                                 OrgDb = "org.Hs.eg.db",
#                                                 keyType = "SYMBOL",
#                                                 ont = "BP")
# 
# GO_results_Secretory_Uninf_vs_MERS_inf_wilcox <- enrichGO(gene = Cluster_10_Uninf_vs_MERS_inf.markers_wilcox_human_gene_names,
#                                                 OrgDb = "org.Hs.eg.db",
#                                                 keyType = "SYMBOL",
#                                                 ont = "BP")
# 
# #ACN4
# GO_results_Secretory_Uninf_vs_ACN4_inf_wilcox <- enrichGO(gene = Ciliated_Uninf_vs_ACN4_inf.markers_wilcox_human_gene_names,
#                                                           OrgDb = "org.Hs.eg.db",
#                                                           keyType = "SYMBOL",
#                                                           ont = "BP")
# 
# GO_results_Ciliated_Uninf_vs_ACN4_inf_wilcox <- enrichGO(gene = Cluster_9_Uninf_vs_ACN4_inf.markers_wilcox_human_gene_names,
#                                                          OrgDb = "org.Hs.eg.db",
#                                                          keyType = "SYMBOL",
#                                                          ont = "BP")
# 
# GO_results_Secretory_Uninf_vs_ACN4_inf_wilcox <- enrichGO(gene = Cluster_10_Uninf_vs_ACN4_inf.markers_wilcox_human_gene_names,
#                                                           OrgDb = "org.Hs.eg.db",
#                                                           keyType = "SYMBOL",
#                                                           ont = "BP")




