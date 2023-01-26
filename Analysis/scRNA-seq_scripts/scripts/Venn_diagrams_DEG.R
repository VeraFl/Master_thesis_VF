
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
library(ComplexHeatmap, quietly = T)
library(circlize, quietly = T)
library(DEGreport, quietly=T)
library(tidyverse, quietly = T)
library(data.table, quietly = T)
library(clusterProfiler, quietly = T)
library(enrichplot, quietly = T)
library(ReactomePA, quietly = T)
library(ggVennDiagram,quietly = T)
library(PCAtools, quietly = T)
library(gprofiler2)
library(biomaRt)
library(gprofiler2)
library(ggVennDiagram)
library(ggupset)
library(MAST)
library(ggvenn)


#Treatment colors(camel-229E, MERS, Mock) camel
cols_treat1 = c("#F39C6B","#A5668B", "#96ADC8")
show_col(cols_treat1)

#Treatment colors (HCoV-229E, MERS-CoV, Mock) human
cols_treat2 = c("#84A59D","#A5668B", "#96ADC8")
show_col(cols_treat2)

#Status colors
cols_stat = c("#CD4631","#DEA47E")
show_col(cols_stat)
#Celltype colors
cols = hcl.colors(7, "Temps")
memory.limit(24000)

?FindConservedMarkers


#Global Analysis #########################################################################################################################################
#Import lists of marker genes that are differentially expressed
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Ferus_both")
Camel_Mock_vs_ACN4_markers_wilcox <- read_csv("Camel_Mock_vs_ACN4_markers_wilcox.csv")
Camel_Mock_vs_ACN4_markers <- as.vector(Camel_Mock_vs_ACN4_markers_wilcox$SYMBOL)

Camel_Mock_vs_MERS_markers_wilcox <- read_csv("Camel_Mock_vs_MERS_markers_wilcox.csv")
Camel_Mock_vs_MERS_markers <- as.vector(Camel_Mock_vs_MERS_markers_wilcox$SYMBOL)

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_both")
Llama_Mock_vs_ACN4_markers_wilcox <- read_csv("Llama_Mock_vs_ACN4_markers_wilcox.csv")
Llama_Mock_vs_ACN4_markers <- as.vector(Llama_Mock_vs_ACN4_markers_wilcox$SYMBOL)

Llama_Mock_vs_MERS_markers_wilcox <- read_csv("Llama_Mock_vs_MERS_markers_wilcox.csv")
Llama_Mock_vs_MERS_markers <- as.vector(Llama_Mock_vs_MERS_markers_wilcox$SYMBOL)

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Human_both")
Human_Mock_vs_HCoV_markers_wilcox <- read_csv("Human_Mock_vs_HCoV_markers_wilcox.csv")
Human_Mock_vs_HCoV_markers <- as.vector(Human_Mock_vs_HCoV_markers_wilcox$SYMBOL)

Human_Mock_vs_MERS_markers_wilcox <- read_csv("Human_Mock_vs_MERS_markers_wilcox.csv")
Human_Mock_vs_MERS_markers <- as.vector(Human_Mock_vs_MERS_markers_wilcox$SYMBOL)

Heatmap()


#venn diagrams
MERS_genes <- list(`Camel MERS-CoV` = Camel_Mock_vs_MERS_markers,
              `Llama MERS-CoV` = Llama_Mock_vs_MERS_markers,
              `Human MERS-CoV` = Human_Mock_vs_MERS_markers)

c_h_229E_genes <- list(`Camel dcCoV-ACN4` = Camel_Mock_vs_ACN4_markers,
                   `Llama dcCoV-ACN4` = Llama_Mock_vs_ACN4_markers,
                   `Human HCoV-229E` = Human_Mock_vs_HCoV_markers)

All <- list(`Camel dcCoV-ACN4` = Camel_Mock_vs_ACN4_markers,
            `Camel MERS-CoV` = Camel_Mock_vs_MERS_markers,
            `Llama MERS-CoV` = Llama_Mock_vs_MERS_markers,
            `Llama dcCoV-ACN4` = Llama_Mock_vs_ACN4_markers)

Camel_genes <- list(`Camel dcCoV-ACN4` = Camel_Mock_vs_ACN4_markers,
                    `Camel MERS-CoV` = Camel_Mock_vs_MERS_markers)

Human_genes <- list(`Human HCoV-229E` = Human_Mock_vs_HCoV_markers,
                    `Human MERS-CoV` = Human_Mock_vs_MERS_markers)

Llama_genes <- list(`Llama dcCoV-ACN4` = Llama_Mock_vs_ACN4_markers,
                    `Llama MERS-CoV` = Llama_Mock_vs_MERS_markers)


camel_genes_intersect <- intersect(Camel_Mock_vs_MERS_markers,Camel_Mock_vs_ACN4_markers)
camel_genes_intersect
human_genes_intersect <- intersect(Human_Mock_vs_HCoV_markers,Human_Mock_vs_MERS_markers)
human_genes_intersect

mers_genes_h_c_intersect <- intersect(Camel_Mock_vs_MERS_markers,Human_Mock_vs_MERS_markers)
mers_genes_l_c_intersect <- intersect(Llama_Mock_vs_MERS_markers,Camel_Mock_vs_MERS_markers)
mers_genes_h_l_intersect <- intersect(Human_Mock_vs_MERS_markers,Llama_Mock_vs_MERS_markers)
mers_genes_h_c_intersect
mers_genes_l_c_intersect
mers_genes_h_l_intersect

#Camel-camel229E, camel-MERS-CoV, Llama-MERS-COV, Llama-camel229E, Human-HCOV-229E, Human-MERS-COV
#"#999999", "#E69F00", "#56B4E9", "#009E73", "#B8336A", "#7768AE"

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_all_Virus_camelids.pdf",
    width=10, height=9)
ggvenn(All, 
       fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_229E_all_species.pdf",
    width=7, height=6)
ggvenn(c_h_229E_genes, 
       fill_color = c("#999999", "#009E73", "#B8336A"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 6)
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_mers_all_species.pdf",
    width=7, height=6)
ggvenn(MERS_genes,
       fill_color = c("#E69F00", "#56B4E9", "#7768AE"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 6)
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_all_Virus_camel.pdf",
    width=9, height=6)
ggvenn(Camel_genes, 
       fill_color = c("#999999","#E69F00"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 6)
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_all_Virus_human.pdf",
    width=9, height=6)
ggvenn(Human_genes, 
       fill_color = c("#B8336A", "#7768AE"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 6)
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_all_Virus_llama.pdf",
    width=9, height=6)
ggvenn(Llama_genes, 
       fill_color = c("#009E73", "#56B4E9"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 6)
dev.off()


#Infected vs Non-infected #########################################################################################################################################

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Ferus_both")
Camel_Uninf_vs_MERS_bys_markers_wilcox <- read_csv("Camel_Uninf_vs_MERS_bys_markers_wilcox.csv")
Camel_Uninf_vs_MERS_bys_markers <- as.vector(Camel_Uninf_vs_MERS_bys_markers_wilcox$SYMBOL)

Camel_Uninf_vs_MERS_inf_markers_wilcox <- read_csv("Camel_Uninf_vs_MERS_inf_markers_wilcox.csv")
Camel_Uninf_vs_MERS_inf_markers <- as.vector(Camel_Uninf_vs_MERS_inf_markers_wilcox$SYMBOL)

Camel_Uninf_vs_ACN4_bys_markers_wilcox <- read_csv("Camel_Uninf_vs_ACN4_bys_markers_wilcox.csv")
Camel_Uninf_vs_ACN4_bys_markers <- as.vector(Camel_Uninf_vs_ACN4_bys_markers_wilcox$SYMBOL)

Camel_Uninf_vs_ACN4_inf_markers_wilcox <- read_csv("Camel_Uninf_vs_ACN4_inf_markers_wilcox.csv")
Camel_Uninf_vs_ACN4_inf_markers <- as.vector(Camel_Uninf_vs_ACN4_inf_markers_wilcox$SYMBOL)


setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_both")
Llama_Uninf_vs_MERS_bys_markers_wilcox <- read_csv("Llama_Uninf_vs_MERS_bys_markers_wilcox.csv")
Llama_Uninf_vs_MERS_bys_markers <- as.vector(Llama_Uninf_vs_MERS_bys_markers_wilcox$SYMBOL)

Llama_Uninf_vs_MERS_inf_markers_wilcox <- read_csv("Llama_Uninf_vs_MERS_inf_markers_wilcox.csv")
Llama_Uninf_vs_MERS_inf_markers <- as.vector(Llama_Uninf_vs_MERS_inf_markers_wilcox$SYMBOL)

Llama_Uninf_vs_ACN4_bys_markers_wilcox <- read_csv("Llama_Uninf_vs_ACN4_bys_markers_wilcox.csv")
Llama_Uninf_vs_ACN4_bys_markers <- as.vector(Llama_Uninf_vs_ACN4_bys_markers_wilcox$SYMBOL)

Llama_Uninf_vs_ACN4_inf_markers_wilcox <- read_csv("Llama_Uninf_vs_ACN4_inf_markers_wilcox.csv")
Llama_Uninf_vs_ACN4_inf_markers <- as.vector(Llama_Uninf_vs_ACN4_inf_markers_wilcox$SYMBOL)

#venn diagrams
All_camel_genes <- list(`Camel MERS-CoV infected` = Camel_Uninf_vs_MERS_inf_markers,
                   `Camel MERS-CoV bystander` = Camel_Uninf_vs_MERS_bys_markers,
                   `Camel dcCoV-ACN4 bystander` = Camel_Uninf_vs_ACN4_bys_markers,
                   `Camel dcCoV-ACN4 infected` = Camel_Uninf_vs_ACN4_inf_markers)

All_llama_genes <- list(`Llama MERS-CoV infected` = Llama_Uninf_vs_MERS_inf_markers,
                  `Llama MERS-CoV bystander` = Llama_Uninf_vs_MERS_bys_markers,
                  `Llama dcCoV-ACN4 bystander` = Llama_Uninf_vs_ACN4_bys_markers,
                  `Llama dcCoV-ACN4 infected` = Llama_Uninf_vs_ACN4_inf_markers)


Bystander_camel_genes_both_viruses <- list(`Camel MERS-CoV bystander` = Camel_Uninf_vs_MERS_bys_markers,
                                     `Camel dcCoV-ACN4 bystander` = Camel_Uninf_vs_ACN4_bys_markers)

Bystander_llama_genes_both_viruses <- list(`Llama MERS-CoV bystander` = Llama_Uninf_vs_MERS_bys_markers,
                                     `Llama dcCoV-ACN4 bystander` = Llama_Uninf_vs_ACN4_bys_markers)

Infected_camel_genes_both_viruses <- list(`Camel MERS-CoV infected` = Camel_Uninf_vs_MERS_inf_markers,
                                    `Camel dcCoV-ACN4 infected` = Camel_Uninf_vs_ACN4_inf_markers)

Bystander_llama_genes_both_viruses <- list(`Llama MERS-CoV bystander` = Llama_Uninf_vs_MERS_bys_markers,
                                     `Llama dcCoV-ACN4 bystander` = Llama_Uninf_vs_ACN4_bys_markers)

#camel vs Llama

Camel_llama_MERS_CoV_bystander <- list(`Camel MERS-CoV bystander` = Camel_Uninf_vs_MERS_bys_markers,
                             `Llama MERS-CoV bystander` = Llama_Uninf_vs_MERS_bys_markers)

Camel_llama_MERS_CoV_infected <- list(`Camel MERS-CoV infected` = Camel_Uninf_vs_MERS_inf_markers,
                             `Llama MERS-CoV infected` = Llama_Uninf_vs_MERS_inf_markers)


Camel_llama_dcCoV_ACN4_bystander <- list(`Camel dcCoV-ACN4 bystander` = Camel_Uninf_vs_ACN4_bys_markers,
                                         `Llama dcCoV-ACN4 bystander` = Llama_Uninf_vs_ACN4_bys_markers)

Camel_llama_dcCoV_ACN4_infected <- list(`Camel dcCoV-ACN4 infected` = Camel_Uninf_vs_ACN4_inf_markers,
                                        `Llama dcCoV-ACN4 infected` = Llama_Uninf_vs_ACN4_inf_markers)


Camel_llama_both_virus_infected <- list(`Camel dcCoV-ACN4 infected` = Camel_Uninf_vs_ACN4_inf_markers,
                                        `Camel MERS-CoV infected` = Camel_Uninf_vs_MERS_inf_markers,
                                        `Llama MERS-CoV infected` = Llama_Uninf_vs_MERS_inf_markers,
                                        `Llama dcCoV-ACN4 infected` = Llama_Uninf_vs_ACN4_inf_markers)

Camel_llama_both_virus_bystander <- list(`Camel dcCoV-ACN4 bystander` = Camel_Uninf_vs_ACN4_bys_markers,
                                        `Camel MERS-CoV bystander` = Camel_Uninf_vs_MERS_bys_markers,
                                        `Llama MERS-CoV bystander` = Llama_Uninf_vs_MERS_bys_markers,
                                        `Llama dcCoV-ACN4 bystander` = Llama_Uninf_vs_ACN4_bys_markers)


camel_inf_genes_both_virus <- intersect(Camel_Uninf_vs_MERS_inf_markers,Camel_Uninf_vs_ACN4_inf_markers)
llama_inf_genes_both_virus <- intersect(Llama_Uninf_vs_MERS_inf_markers,Llama_Uninf_vs_ACN4_inf_markers)

#Camel##################################################################################################################

All_camel_genes <- list(`Camel MERS-CoV infected` = Camel_Uninf_vs_MERS_inf_markers,
                        `Camel MERS-CoV bystander` = Camel_Uninf_vs_MERS_bys_markers,
                        `Camel dcCoV-ACN4 bystander` = Camel_Uninf_vs_ACN4_bys_markers,
                        `Camel dcCoV-ACN4 infected` = Camel_Uninf_vs_ACN4_inf_markers)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_camel_status_.pdf",
    width=10, height=9)
ggvenn(All_camel_genes, 
       fill_color = c("#F3A712", "#A8C686","#E4572E","#669BBC"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)
dev.off()

#Llama##################################################################################################################
All_llama_genes <- list(`Llama MERS-CoV infected` = Llama_Uninf_vs_MERS_inf_markers,
                        `Llama MERS-CoV bystander` = Llama_Uninf_vs_MERS_bys_markers,
                        `Llama dcCoV-ACN4 bystander` = Llama_Uninf_vs_ACN4_bys_markers,
                        `Llama dcCoV-ACN4 infected` = Llama_Uninf_vs_ACN4_inf_markers)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_llama_status_.pdf",
    width=10, height=9)
ggvenn(All_llama_genes, 
       fill_color = c("#224429", "#4F86C6","#744FC6","#C83E4D"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)
dev.off()


#Llama%Camel##################################################################################################################
Camel_llama_both_virus_infected <- list(`Camel dcCoV-ACN4 infected` = Camel_Uninf_vs_ACN4_inf_markers,
                                        `Camel MERS-CoV infected` = Camel_Uninf_vs_MERS_inf_markers,
                                        `Llama MERS-CoV infected` = Llama_Uninf_vs_MERS_inf_markers,
                                        `Llama dcCoV-ACN4 infected` = Llama_Uninf_vs_ACN4_inf_markers)


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_camel_vs_llama_infected.pdf",
    width=10, height=9)
ggvenn(Camel_llama_both_virus_infected, 
       fill_color = c("#669BBC", "#F3A712","#224429","#C83E4D"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4,
       show_elements = F)
dev.off()

?ggvenn

#Camel - ACN4 vs MERS
intersect(Camel_Uninf_vs_ACN4_inf_markers, Camel_Uninf_vs_MERS_inf_markers)

#Llama - ACN4 vs MERS
intersect(Llama_Uninf_vs_ACN4_inf_markers, Llama_Uninf_vs_MERS_inf_markers)

#Camel - ACN4 vs MERS vs Llama - MERS
intersect(intersect(Camel_Uninf_vs_ACN4_inf_markers, Camel_Uninf_vs_MERS_inf_markers), Llama_Uninf_vs_MERS_inf_markers)

#Camel - 229E vs MERS vs Llama - 229E
intersect(intersect(Camel_Uninf_vs_ACN4_inf_markers, Camel_Uninf_vs_MERS_inf_markers), Llama_Uninf_vs_ACN4_inf_markers)

#Llama - ACN4 vs MERS vs Camel - MERS
intersect(intersect(Llama_Uninf_vs_ACN4_inf_markers, Llama_Uninf_vs_MERS_inf_markers), Camel_Uninf_vs_MERS_inf_markers)

#Llama - 229E vs MERS vs Camel - 229E
intersect(intersect(Llama_Uninf_vs_ACN4_inf_markers, Llama_Uninf_vs_MERS_inf_markers), Camel_Uninf_vs_ACN4_inf_markers)

#Camel - ACN4 vs Llama - MERS
intersect(Camel_Uninf_vs_ACN4_inf_markers,Llama_Uninf_vs_MERS_inf_markers)

#Camel - MERS vs Llama - ACN4
intersect(Camel_Uninf_vs_MERS_inf_markers,Llama_Uninf_vs_ACN4_inf_markers)



#################################################

Camel_llama_both_virus_bystander <- list(`Camel dcCoV-ACN4 bystander` = Camel_Uninf_vs_ACN4_bys_markers,
                                         `Camel MERS-CoV bystander` = Camel_Uninf_vs_MERS_bys_markers,
                                         `Llama MERS-CoV bystander` = Llama_Uninf_vs_MERS_bys_markers,
                                         `Llama dcCoV-ACN4 bystander` = Llama_Uninf_vs_ACN4_bys_markers)
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_camel_vs_llama_bystander.pdf",
    width=10, height=9)
ggvenn(Camel_llama_both_virus_bystander, 
       fill_color = c("#E4572E", "#A8C686","#4F86C6","#744FC6"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)
dev.off()

Llama_Uninf_vs_ACN4_bys_markers
Camel_Uninf_vs_MERS_bys_markers


camel_bys_genes_both_virus <- intersect(Camel_Uninf_vs_MERS_bys_markers,Camel_Uninf_vs_ACN4_bys_markers)
llama_bys_genes_both_virus <- intersect(Llama_Uninf_vs_MERS_bys_markers,Llama_Uninf_vs_ACN4_bys_markers)
camel_bys_genes_both_virus
llama_bys_genes_both_virus

camel_v_llama_MERS_bys_genes_both_virus <- intersect(Camel_Uninf_vs_MERS_bys_markers,Llama_Uninf_vs_MERS_bys_markers)
camel_v_llama_ACN4_bys_genes_both_virus <- intersect(Camel_Uninf_vs_ACN4_bys_markers,Llama_Uninf_vs_ACN4_bys_markers)
camel_v_llama_MERS_bys_genes_both_virus
camel_v_llama_ACN4_bys_genes_both_virus

intersect(Camel_Uninf_vs_MERS_bys_markers,Llama_Uninf_vs_ACN4_bys_markers)

?ggvenn
#Across Celltypes#########################################################################################################################################

#save the lists
#Infected
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/DEG/Ferus_both/Celltype")
Secretory_Camel_MERS_Infected <- read_csv("Secretory_Camel_MERS_Infected.csv")
Secretory_Camel_ACN4_Infected <- read_csv("Secretory_Camel_ACN4_Infected.csv")
Secretory_Camel_MERS_Infected <- as.vector(Secretory_Camel_MERS_Infected$SYMBOL)
Secretory_Camel_ACN4_Infected <- as.vector(Secretory_Camel_ACN4_Infected$SYMBOL)

Ciliated_Camel_MERS_Infected <- read_csv("Ciliated_Camel_MERS_Infected.csv")
Ciliated_Camel_ACN4_Infected <- read_csv("Ciliated_Camel_ACN4_Infected.csv")
Ciliated_Camel_MERS_Infected <- as.vector(Ciliated_Camel_MERS_Infected$SYMBOL)
Ciliated_Camel_ACN4_Infected <- as.vector(Ciliated_Camel_ACN4_Infected$SYMBOL)

Club_Camel_MERS_Infected <- read_csv("Club_Camel_MERS_Infected.csv")
Club_Camel_ACN4_Infected <- read_csv("Club_Camel_ACN4_Infected.csv")
Club_Camel_MERS_Infected <- as.vector(Club_Camel_MERS_Infected$SYMBOL)
Club_Camel_ACN4_Infected <- as.vector(Club_Camel_ACN4_Infected$SYMBOL)

Basal_Camel_MERS_Infected <- read_csv("Basal_Camel_MERS_Infected.csv")
Basal_Camel_ACN4_Infected <- read_csv("Basal_Camel_ACN4_Infected.csv")
Basal_Camel_MERS_Infected <- as.vector(Basal_Camel_MERS_Infected$SYMBOL)
Basal_Camel_ACN4_Infected <- as.vector(Basal_Camel_ACN4_Infected$SYMBOL)

Cluster1_Camel_MERS_Infected <- read_csv("Cluster1_Camel_MERS_Infected.csv")
Cluster1_Camel_ACN4_Infected <- read_csv("Cluster1_Camel_ACN4_Infected.csv")
Cluster1_Camel_MERS_Infected <- as.vector(Cluster1_Camel_MERS_Infected$SYMBOL)
Cluster1_Camel_ACN4_Infected <- as.vector(Cluster1_Camel_ACN4_Infected$SYMBOL)

Cluster2_Camel_MERS_Infected <- read_csv("Cluster2_Camel_MERS_Infected.csv")
Cluster2_Camel_ACN4_Infected <- read_csv("Cluster2_Camel_ACN4_Infected.csv")
Cluster2_Camel_MERS_Infected <- as.vector(Cluster2_Camel_MERS_Infected$SYMBOL)
Cluster2_Camel_ACN4_Infected <- as.vector(Cluster2_Camel_ACN4_Infected$SYMBOL)

#Bystander
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/DEG/Ferus_both/Celltype")
Secretory_Camel_MERS_Bystander <- read_csv("Secretory_Camel_MERS_Bystander.csv")
Secretory_Camel_ACN4_Bystander <- read_csv("Secretory_Camel_ACN4_Bystander.csv")
Secretory_Camel_MERS_Bystander <- as.vector(Secretory_Camel_MERS_Bystander$SYMBOL)
Secretory_Camel_ACN4_Bystander <- as.vector(Secretory_Camel_ACN4_Bystander$SYMBOL)

Ciliated_Camel_MERS_Bystander <- read_csv("Ciliated_Camel_MERS_Bystander.csv")
Ciliated_Camel_ACN4_Bystander <- read_csv("Ciliated_Camel_ACN4_Bystander.csv")
Ciliated_Camel_MERS_Bystander <- as.vector(Ciliated_Camel_MERS_Bystander$SYMBOL)
Ciliated_Camel_ACN4_Bystander <- as.vector(Ciliated_Camel_ACN4_Bystander$SYMBOL)

Club_Camel_MERS_Bystander <- read_csv("Club_Camel_MERS_Bystander.csv")
Club_Camel_ACN4_Bystander <- read_csv("Club_Camel_ACN4_Bystander.csv")
Club_Camel_MERS_Bystander <- as.vector(Club_Camel_MERS_Bystander$SYMBOL)
Club_Camel_ACN4_Bystander <- as.vector(Club_Camel_ACN4_Bystander$SYMBOL)

Basal_Camel_MERS_Bystander <- read_csv("Basal_Camel_MERS_Bystander.csv")
Basal_Camel_ACN4_Bystander <- read_csv("Basal_Camel_ACN4_Bystander.csv")
Basal_Camel_MERS_Bystander <- as.vector(Basal_Camel_MERS_Bystander$SYMBOL)
Basal_Camel_ACN4_Bystander <- as.vector(Basal_Camel_ACN4_Bystander$SYMBOL)

Cluster1_Camel_MERS_Bystander <- read_csv("Cluster1_Camel_MERS_Bystander.csv")
Cluster1_Camel_ACN4_Bystander <- read_csv("Cluster1_Camel_ACN4_Bystander.csv")
Cluster1_Camel_MERS_Bystander <- as.vector(Cluster1_Camel_MERS_Bystander$SYMBOL)
Cluster1_Camel_ACN4_Bystander <- as.vector(Cluster1_Camel_ACN4_Bystander$SYMBOL)

Cluster2_Camel_MERS_Bystander <- read_csv("Cluster2_Camel_MERS_Bystander.csv")
Cluster2_Camel_ACN4_Bystander <- read_csv("Cluster2_Camel_ACN4_Bystander.csv")
Cluster2_Camel_MERS_Bystander <- as.vector(Cluster2_Camel_MERS_Bystander$SYMBOL)
Cluster2_Camel_ACN4_Bystander <- as.vector(Cluster2_Camel_ACN4_Infected$SYMBOL)


setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/DEG/Llama_both/Celltype")
#Infected
Secretory_Llama_MERS_Infected <- read_csv("Secretory_Llama_MERS_Infected.csv")
Secretory_Llama_ACN4_Infected <- read_csv("Secretory_Llama_ACN4_Infected.csv")
Secretory_Llama_MERS_Infected <- as.vector(Secretory_Llama_MERS_Infected$SYMBOL)
Secretory_Llama_ACN4_Infected <- as.vector(Secretory_Llama_ACN4_Infected$SYMBOL)

Ciliated_Llama_MERS_Infected <- read_csv("Ciliated_Llama_MERS_Infected.csv")
Ciliated_Llama_ACN4_Infected <- read_csv("Ciliated_Llama_ACN4_Infected.csv")
Ciliated_Llama_MERS_Infected <- as.vector(Ciliated_Llama_MERS_Infected$SYMBOL)
Ciliated_Llama_ACN4_Infected <- as.vector(Ciliated_Llama_ACN4_Infected$SYMBOL)

Club_Llama_MERS_Infected <- read_csv("Club_Llama_MERS_Infected.csv")
Club_Llama_ACN4_Infected <- read_csv("Club_Llama_ACN4_Infected.csv")
Club_Llama_MERS_Infected <- as.vector(Club_Llama_MERS_Infected$SYMBOL)
Club_Llama_ACN4_Infected <- as.vector(Club_Llama_ACN4_Infected$SYMBOL)

Basal_Llama_MERS_Infected <- read_csv("Basal_Llama_MERS_Infected.csv")
Basal_Llama_ACN4_Infected <- read_csv("Basal_Llama_ACN4_Infected.csv")
Basal_Llama_MERS_Infected <- as.vector(Basal_Llama_MERS_Infected$SYMBOL)
Basal_Llama_ACN4_Infected <- as.vector(Basal_Llama_ACN4_Infected$SYMBOL)

Cluster1_Llama_MERS_Infected <- read_csv("Cluster1_Llama_MERS_Infected.csv")
Cluster1_Llama_ACN4_Infected <- read_csv("Cluster1_Llama_ACN4_Infected.csv")
Cluster1_Llama_MERS_Infected <- as.vector(Cluster1_Llama_MERS_Infected$SYMBOL)
Cluster1_Llama_ACN4_Infected <- as.vector(Cluster1_Llama_ACN4_Infected$SYMBOL)

Cluster2_Llama_MERS_Infected <- read_csv("Cluster2_Llama_MERS_Infected.csv")
Cluster2_Llama_ACN4_Infected <- read_csv("Cluster2_Llama_ACN4_Infected.csv")
Cluster2_Llama_MERS_Infected <- as.vector(Cluster2_Llama_MERS_Infected$SYMBOL)
Cluster2_Llama_ACN4_Infected <- as.vector(Cluster2_Llama_ACN4_Infected$SYMBOL)

#Bystander
Secretory_Llama_MERS_Bystander <- read_csv("Secretory_Llama_MERS_Bystander.csv")
Secretory_Llama_ACN4_Bystander <- read_csv("Secretory_Llama_ACN4_Bystander.csv")
Secretory_Llama_MERS_Bystander <- as.vector(Secretory_Llama_MERS_Bystander$SYMBOL)
Secretory_Llama_ACN4_Bystander <- as.vector(Secretory_Llama_ACN4_Bystander$SYMBOL)

Ciliated_Llama_MERS_Bystander <- read_csv("Ciliated_Llama_MERS_Bystander.csv")
Ciliated_Llama_ACN4_Bystander <- read_csv("Ciliated_Llama_ACN4_Bystander.csv")
Ciliated_Llama_MERS_Bystander <- as.vector(Ciliated_Llama_MERS_Bystander$SYMBOL)
Ciliated_Llama_ACN4_Bystander <- as.vector(Ciliated_Llama_ACN4_Bystander$SYMBOL)

Club_Llama_MERS_Bystander <- read_csv("Club_Llama_MERS_Bystander.csv")
Club_Llama_ACN4_Bystander <- read_csv("Club_Llama_ACN4_Bystander.csv")
Club_Llama_MERS_Bystander <- as.vector(Club_Llama_MERS_Bystander$SYMBOL)
Club_Llama_ACN4_Bystander <- as.vector(Club_Llama_ACN4_Bystander$SYMBOL)

Basal_Llama_MERS_Bystander <- read_csv("Basal_Llama_MERS_Bystander.csv")
Basal_Llama_ACN4_Bystander <- read_csv("Basal_Llama_ACN4_Bystander.csv")
Basal_Llama_MERS_Bystander <- as.vector(Basal_Llama_MERS_Bystander$SYMBOL)
Basal_Llama_ACN4_Bystander <- as.vector(Basal_Llama_ACN4_Bystander$SYMBOL)

Cluster1_Llama_MERS_Bystander <- read_csv("Cluster1_Llama_MERS_Bystander.csv")
Cluster1_Llama_ACN4_Bystander <- read_csv("Cluster1_Llama_ACN4_Bystander.csv")
Cluster1_Llama_MERS_Bystander <- as.vector(Cluster1_Llama_MERS_Bystander$SYMBOL)
Cluster1_Llama_ACN4_Bystander <- as.vector(Cluster1_Llama_ACN4_Bystander$SYMBOL)

Cluster2_Llama_MERS_Bystander <- read_csv("Cluster2_Llama_MERS_Bystander.csv")
Cluster2_Llama_ACN4_Bystander <- read_csv("Cluster2_Llama_ACN4_Bystander.csv")
Cluster2_Llama_MERS_Bystander <- as.vector(Cluster2_Llama_MERS_Bystander$SYMBOL)
Cluster2_Llama_ACN4_Bystander <- as.vector(Cluster2_Llama_ACN4_Bystander$SYMBOL)


#Camel##########################################################################
camel_celltypePalette <- structure(c("#807EBF", "#7EA8BE", "#089392", "#40B48B", "#9DCD84", "#EAE29C"), 
                          names = c('Cluster1','Cluster2','Secretory', 'Ciliated', 'Club', 'Basal'))


#MERS infected#################################################################
Camel_MERS_inf <- list(`Ciliated` = Ciliated_Camel_MERS_Infected,
                       `Secretory` = Secretory_Camel_MERS_Infected,
                       `Cluster 1` = Cluster1_Camel_MERS_Infected,
                       `Club` = Club_Camel_MERS_Infected)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_camel_MERS_inf_celltypes.pdf",
    width=10, height=10)
ggvenn(Camel_MERS_inf, 
       fill_color = c("#40B48B","#089392","#807EBF","#9DCD84"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)
dev.off()

intersect(intersect(Ciliated_Camel_MERS_Infected, Secretory_Camel_MERS_Infected), intersect(Cluster1_Camel_MERS_Infected, Club_Camel_MERS_Infected))
intersect(intersect(Ciliated_Camel_MERS_Infected, Secretory_Camel_MERS_Infected), Club_Camel_MERS_Infected)
intersect(intersect(Ciliated_Camel_MERS_Infected, Secretory_Camel_MERS_Infected), Cluster1_Camel_MERS_Infected)
intersect(Secretory_Camel_MERS_Infected, intersect(Cluster1_Camel_MERS_Infected, Club_Camel_MERS_Infected))
intersect(Ciliated_Camel_MERS_Infected, intersect(Cluster1_Camel_MERS_Infected, Club_Camel_MERS_Infected))
intersect(Ciliated_Camel_MERS_Infected, Club_Camel_MERS_Infected)
intersect(Secretory_Camel_MERS_Infected, Club_Camel_MERS_Infected)
intersect(Club_Camel_MERS_Infected, Cluster1_Camel_MERS_Infected)
intersect(Ciliated_Camel_MERS_Infected, Cluster1_Camel_MERS_Infected)
intersect(Secretory_Camel_MERS_Infected, Cluster1_Camel_MERS_Infected)

#MERS bystander################################################################
Camel_MERS_bys <- list(`Ciliated` = Ciliated_Camel_MERS_Bystander,
                       `Secretory` = Secretory_Camel_MERS_Bystander,
                       `Cluster 1` = Cluster1_Camel_MERS_Bystander,
                       `Club` = Club_Camel_MERS_Bystander)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_camel_MERS_bys_celltypes.pdf",
    width=10, height=10)
ggvenn(Camel_MERS_bys, 
       fill_color = c("#40B48B","#089392","#807EBF","#9DCD84"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)
dev.off()

intersect(intersect(Ciliated_Camel_MERS_Bystander, Secretory_Camel_MERS_Bystander), intersect(Cluster1_Camel_MERS_Bystander, Club_Camel_MERS_Bystander))
intersect(intersect(Ciliated_Camel_MERS_Bystander, Secretory_Camel_MERS_Bystander), Club_Camel_MERS_Bystander)
intersect(intersect(Ciliated_Camel_MERS_Bystander, Secretory_Camel_MERS_Bystander), Cluster1_Camel_MERS_Bystander)
intersect(Secretory_Camel_MERS_Bystander, intersect(Cluster1_Camel_MERS_Bystander, Club_Camel_MERS_Bystander))
intersect(Ciliated_Camel_MERS_Bystander, intersect(Cluster1_Camel_MERS_Bystander, Club_Camel_MERS_Bystander))
intersect(Ciliated_Camel_MERS_Bystander, Club_Camel_MERS_Bystander)
intersect(Secretory_Camel_MERS_Bystander, Club_Camel_MERS_Bystander)
intersect(Club_Camel_MERS_Bystander, Cluster1_Camel_MERS_Bystander)
intersect(Ciliated_Camel_MERS_Bystander, Cluster1_Camel_MERS_Bystander)
intersect(Secretory_Camel_MERS_Bystander, Cluster1_Camel_MERS_Bystander)

#ACN4 infected#################################################################
Camel_ACN4_inf <- list(`Ciliated` = Ciliated_Camel_ACN4_Infected,
                       `Secretory` = Secretory_Camel_ACN4_Infected,
                       `Cluster 2` = Cluster2_Camel_ACN4_Infected)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_camel_ACN4_inf_celltypes.pdf",
    width=10, height=10)
ggvenn(Camel_ACN4_inf, 
       fill_color = c("#40B48B","#089392","#7EA8BE"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)
dev.off()

#ACN4 bystander################################################################
Camel_ACN4_bys <- list(`Basal` = Basal_Camel_ACN4_Bystander,
                       `Secretory` = Secretory_Camel_ACN4_Bystander,
                       `Cluster 1` = Cluster1_Camel_ACN4_Bystander,
                       `Club` = Club_Camel_ACN4_Bystander)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_camel_ACN4_bys_celltypes.pdf",
    width=10, height=10)
ggvenn(Camel_ACN4_bys, 
       fill_color = c("#EAE29C","#089392","#807EBF","#9DCD84"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)
dev.off()

#Cluster1 vs Cluster2
Camel_Cluster1_v_cluster2 <- list(`Camel_Cluster 1` = Cluster1_Camel_MERS_Infected,
                                  `Camel_Cluster 2` = Cluster2_Camel_MERS_Infected)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_camel_MERS_inf_cluster1_v_cluster2.pdf",
    width=10, height=10)
ggvenn(Camel_Cluster1_v_cluster2, 
       fill_color = c("#807EBF", "#7EA8BE"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)
dev.off()

intersect(Cluster1_Camel_MERS_Infected, Cluster2_Camel_MERS_Infected)


#Camel_vs_Llama##################################################################

Secretory_llama_camel <- list(`Llama Secretory` = Secretory_Llama_MERS_Infected,
                                  `Camel Secretory` = Secretory_Camel_MERS_Infected)

v1 <- ggvenn(Secretory_llama_camel, 
       fill_color = c("black", "white"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)

Club_llama_camel <- list(`Llama Club` = Club_Llama_MERS_Infected,
                              `Camel Club` = Club_Camel_MERS_Infected)

v2 <- ggvenn(Club_llama_camel, 
       fill_color = c("black", "white"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)

Ciliated_llama_camel <- list(`Llama Ciliated` = Ciliated_Llama_MERS_Infected,
                              `Camel Ciliated` = Ciliated_Camel_MERS_Infected)

v3 <- ggvenn(Ciliated_llama_camel, 
       fill_color = c("black", "white"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)

Basal_llama_camel <- list(`Llama Basal` = Basal_Llama_MERS_Infected,
                             `Camel Basal` = Basal_Camel_MERS_Infected)

v4 <- ggvenn(Basal_llama_camel, 
       fill_color = c("black", "white"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)

(v1/v2)|(v3/v4)


intersect(Secretory_Llama_MERS_Infected, Secretory_Camel_MERS_Infected)
intersect(Ciliated_Llama_MERS_Infected, Ciliated_Camel_MERS_Infected)
intersect(Club_Llama_MERS_Infected, Club_Camel_MERS_Infected)
intersect(Basal_Llama_MERS_Infected, Basal_Camel_MERS_Infected)

#Llama#########################################################################
llama_celltypePalette <- structure(c("#8A6E7E", "#40679B", "#089392", "#40B48B", "#9DCD84", "#EAE29C"), 
                          names = c("Cluster1", "Cluster2", "Secretory", "Ciliated", "Club", "Basal"))


#MERS infected#################################################################
Llama_MERS_inf <- list(`Ciliated` = Ciliated_Llama_MERS_Infected,
                       `Secretory` = Secretory_Llama_MERS_Infected,
                       `Cluster 2` = Cluster2_Llama_MERS_Infected,
                       `Club` = Club_Llama_MERS_Infected)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_llama_MERS_inf_celltypes.pdf",
    width=10, height=10)
ggvenn(Llama_MERS_inf, 
       fill_color = c("#40B48B","#089392","#40679B","#9DCD84"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)
dev.off()

intersect(intersect(Ciliated_Llama_MERS_Infected, Secretory_Llama_MERS_Infected), intersect(Cluster2_Llama_MERS_Infected, Club_Llama_MERS_Infected))
intersect(intersect(Ciliated_Llama_MERS_Infected, Secretory_Llama_MERS_Infected), Club_Llama_MERS_Infected)
intersect(intersect(Ciliated_Llama_MERS_Infected, Secretory_Llama_MERS_Infected), Cluster1_Llama_MERS_Infected)
intersect(Secretory_Llama_MERS_Infected, intersect(Cluster1_Llama_MERS_Infected, Club_Llama_MERS_Infected))
intersect(Ciliated_Llama_MERS_Infected, intersect(Cluster1_Llama_MERS_Infected, Club_Llama_MERS_Infected))
intersect(Ciliated_Llama_MERS_Infected, Club_Llama_MERS_Infected)
intersect(Secretory_Llama_MERS_Infected, Club_Llama_MERS_Infected)
intersect(Club_Llama_MERS_Infected, Cluster1_Llama_MERS_Infected)
intersect(Ciliated_Llama_MERS_Infected, Secretory_Llama_MERS_Infected)
intersect(Ciliated_Llama_MERS_Infected, Cluster1_Llama_MERS_Infected)
intersect(Secretory_Llama_MERS_Infected, Cluster1_Llama_MERS_Infected)

#MERS bystander################################################################
Llama_MERS_bys <- list(`Ciliated` = Ciliated_Llama_MERS_Bystander,
                       `Secretory` = Secretory_Llama_MERS_Bystander,
                       `Cluster 2` = Cluster2_Llama_MERS_Bystander,
                       `Club` = Club_Llama_MERS_Bystander)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_llama_MERS_bys_celltypes.pdf",
    width=10, height=10)
ggvenn(Llama_MERS_bys, 
       fill_color = c("#40B48B","#089392","#40679B","#9DCD84"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)
dev.off()

intersect(Ciliated_Llama_MERS_Bystander, Cluster2_Llama_MERS_Bystander)

#ACN4 infected#################################################################
Llama_ACN4_inf <- list(`Ciliated` = Ciliated_Llama_ACN4_Infected,
                       `Secretory` = Secretory_Llama_ACN4_Infected,
                       `Cluster 1` = Cluster1_Llama_ACN4_Infected,
                       `Basal` = Basal_Llama_ACN4_Infected)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_llama_ACN4_inf_celltypes.pdf",
    width=10, height=10)
ggvenn(Llama_ACN4_inf, 
       fill_color = c("#40B48B","#089392","#8A6E7E", "#EAE29C"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)
dev.off()

#ACN4 bystander################################################################
Llama_ACN4_bys <- list(`Ciliated` = Ciliated_Llama_ACN4_Bystander,
                       `Secretory` = Secretory_Llama_ACN4_Bystander,
                       `Cluster 1` = Cluster1_Llama_ACN4_Bystander,
                       `Basal` = Basal_Llama_ACN4_Bystander)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_llama_ACN4_bys_celltypes.pdf",
    width=10, height=10)
ggvenn(Llama_ACN4_bys, 
       fill_color = c("#40B48B","#089392","#8A6E7E", "#EAE29C"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)
dev.off()

######################
#Cluster1 vs Cluster2
Llama_Cluster1_v_cluster2 <- list(`Llama_Cluster 1` = Cluster1_Llama_MERS_Infected,
                                  `Llama_Cluster 2` = Cluster2_Llama_MERS_Infected)

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Venn_llama_MERS_inf_cluster1_v_cluster2.pdf",
    width=10, height=10)
ggvenn(Llama_Cluster1_v_cluster2, 
       fill_color = c("#8A6E7E","#40679B"),
       show_percentage = FALSE, text_size = 7,
       stroke_size = 0.5, set_name_size = 4)
dev.off()

intersect(Cluster1_Llama_MERS_Infected,Cluster2_Llama_MERS_Infected)
#Dotplot####################################################################################################################################

#Camel#################################################################################
# Import DEG files
# Define viral genes
virus_genes <- c("mers-3UTR","mers-orf8b","mers-N","mers-M","mers-E","mers-orf5","mers-orf4b","mers-orf4a","mers-orf3","mers-S","mers-orf1ab","mers-5UTR",
                 "camel229E-3UTR","camel229E-ORF8","camel229E-N","camel229E-M","camel229E-E","camel229E-ORF4","camel229E-S","camel229E-ORF1a", "camel229E-ORF1ab","camel229E-5UTR")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Ferus_both")
# import files
filelist = list.files(pattern="*.csv", full.names = TRUE)
datalist <- lapply(filelist, FUN= read.csv, na.strings="")

# Modify datalist
names(datalist) <- gsub(".csv", "", basename(filelist))
datalist <- lapply(datalist, function(x) x[!x$SYMBOL %in% virus_genes, ]) # Remove viral genes from files
datalist_entrezid <- lapply(datalist, function(x) {x$ENTREZID <- as.character(x$ENTREZID); return(x)})

#deleting the empty lists
datalist_entrezid[c(1, 7, 13, 16)] = NULL
filelist <- filelist[-c(1, 7, 13, 16)]

# Add variable information
datalist_entrezid <- lapply(names(datalist_entrezid), function(x) {datalist_entrezid[[x]]$name = x; datalist_entrezid[[x]]}) # Add column "name" with dataframe name
datalist_entrezid <- lapply(datalist_entrezid, tidyr::separate, name, c("CellType", "Species", "Condition", "Status"), remove = FALSE, sep = "_")
datalist_entrezid <- lapply(datalist_entrezid, function(x) {x$expr <- ifelse(x$avg_log2FC > 0, "upregulated", "downregulated"); return(x)}) # Annotate if gene is up or downregulated

names(datalist_entrezid) <- gsub(".csv", "", basename(filelist))

# #get group GO results
# groupGO_result <- compareCluster(ENTREZID~expr+CellType+Condition+Status, data=rbindlist(datalist_entrezid),
#                                  keyType = "ENTREZID",
#                                  fun = "groupGO",
#                                  ont      = "BP",
#                                  level    = 3,
#                                  OrgDb='org.Hs.eg.db',
#                                  readable = TRUE)
# 
# #dotplot of groupgo results
# p.groupGO_result <- dotplot(groupGO_result, x="CellType", showCategory = 1, font.size = 16) +
#   facet_grid(.~Condition + Status ~expr) + theme_classic() + 
#   theme(axis.text = element_text(size = 14, family = "sans", ),
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         legend.text = element_text(size = 14, family = "sans"),
#         legend.title = element_text(size = 14, family = "sans"),
#         axis.title = element_text(size = 14, family = "sans", face = "bold"),
#         axis.line = element_line(linetype = "solid", colour = "black", size=1),
#         strip.text = element_text(size = 16)) +
#   scale_colour_gradientn(colours = c("red", "lightblue"))
# p.groupGO_result_result


# data_entrezid <- lapply(datalist_entrezid, function(x) {x$ENTREZID})
# 
# ck <- compareCluster(geneCluster = data_entrezid, fun = enrichGO, OrgDb='org.Hs.eg.db')
# ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
# head(ck) 
# cnetplot(ck)


#get enrichgo results
enrichGO_result <- compareCluster(ENTREZID~expr+CellType+Condition+Status, data=rbindlist(datalist_entrezid), fun="enrichGO", OrgDb='org.Hs.eg.db', pvalueCutoff=0.01)
p.enrichGO_result <- dotplot(enrichGO_result, x="CellType", showCategory = 1, font.size = 16) +
  facet_grid(.~expr ~Condition + Status ) + theme_classic() + 
  theme(axis.text = element_text(size = 14, family = "sans", ),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = 14, family = "sans"),
        legend.title = element_text(size = 14, family = "sans"),
        axis.title = element_text(size = 14, family = "sans", face = "bold"),
        axis.line = element_line(linetype = "solid", colour = "black", size=1),
        strip.text = element_text(size = 16)) +
  scale_colour_gradientn(colours = c("red", "lightblue"))



pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/dotplot_across_all_conditions.pdf",
    width=10, height=10)
p.enrichGO_result +theme(legend.position = "top")
dev.off()

p1 <- cnetplot(ck, repel = T,
               cex_label_category = 1.5,
               cex_label_gene = 1.5) + NoLegend()
?cnetplot
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/cnetplot_across_all_conditions.pdf",
    width=10, height=10)
p1
dev.off()


#Llama#################################################################################
# Import DEG files
# Define viral genes
virus_genes <- c("mers-3UTR","mers-orf8b","mers-N","mers-M","mers-E","mers-orf5","mers-orf4b","mers-orf4a","mers-orf3","mers-S","mers-orf1ab","mers-5UTR",
                 "camel229E-3UTR","camel229E-ORF8","camel229E-N","camel229E-M","camel229E-E","camel229E-ORF4","camel229E-S","camel229E-ORF1a", "camel229E-ORF1ab","camel229E-5UTR")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Llama_both")
# import files
filelist = list.files(pattern="*.csv", full.names = TRUE)
datalist <- lapply(filelist, FUN= read.csv, na.strings="")

# Modify datalist
names(datalist) <- gsub(".csv", "", basename(filelist))
datalist <- lapply(datalist, function(x) x[!x$SYMBOL %in% virus_genes, ]) # Remove viral genes from files
datalist_entrezid <- lapply(datalist, function(x) {x$ENTREZID <- as.character(x$ENTREZID); return(x)})


# Add variable information
datalist_entrezid <- lapply(names(datalist_entrezid), function(x) {datalist_entrezid[[x]]$name = x; datalist_entrezid[[x]]}) # Add column "name" with dataframe name
datalist_entrezid <- lapply(datalist_entrezid, tidyr::separate, name, c("CellType", "Species", "Condition", "Status"), remove = FALSE, sep = "_")
datalist_entrezid <- lapply(datalist_entrezid, function(x) {x$expr <- ifelse(x$avg_log2FC > 0, "upregulated", "downregulated"); return(x)}) # Annotate if gene is up or downregulated

names(datalist_entrezid) <- gsub(".csv", "", basename(filelist))

# #get group GO results
# groupGO_result <- compareCluster(ENTREZID~expr+CellType+Condition+Status, data=rbindlist(datalist_entrezid),
#                                  keyType = "ENTREZID",
#                                  fun = "groupGO",
#                                  ont      = "BP",
#                                  level    = 3,
#                                  OrgDb='org.Hs.eg.db',
#                                  readable = TRUE)
# 
# #dotplot of groupgo results
# p.groupGO_result <- dotplot(groupGO_result, x="CellType", showCategory = 1, font.size = 16) +
#   facet_grid(.~Condition + Status ~expr) + theme_classic() + 
#   theme(axis.text = element_text(size = 14, family = "sans", ),
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         legend.text = element_text(size = 14, family = "sans"),
#         legend.title = element_text(size = 14, family = "sans"),
#         axis.title = element_text(size = 14, family = "sans", face = "bold"),
#         axis.line = element_line(linetype = "solid", colour = "black", size=1),
#         strip.text = element_text(size = 16)) +
#   scale_colour_gradientn(colours = c("red", "lightblue"))
# p.groupGO_result_result

data_entrezid <- lapply(datalist_entrezid, function(x) {x$ENTREZID})

ck <- compareCluster(geneCluster = data_entrezid, fun = enrichGO, OrgDb='org.Hs.eg.db')
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck) 
cnetplot(ck)

#get enrichgo results
enrichGO_result <- compareCluster(ENTREZID~expr+CellType+Condition+Status, data=rbindlist(datalist_entrezid), fun="enrichGO", OrgDb='org.Hs.eg.db', pvalueCutoff=0.01, readable = T)
p.enrichGO_result <- dotplot(enrichGO_result, x="CellType", showCategory = 10, font.size = 16) +
  facet_grid(.~expr ~Condition + Status) + theme_classic() + 
  theme(axis.text = element_text(size = 14, family = "sans", ),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = 14, family = "sans"),
        legend.title = element_text(size = 14, family = "sans"),
        axis.title = element_text(size = 14, family = "sans", face = "bold"),
        axis.line = element_line(linetype = "solid", colour = "black", size=1),
        strip.text = element_text(size = 16)) +
  scale_colour_gradientn(colours = c("red", "lightblue"))



pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/dotplot_across_all_conditions.pdf",
    width=15, height=15)
p.enrichGO_result + theme(legend.position = "top")
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/cnet_across_all_conditions.pdf",
    width=10, height=10)
cnetplot(enrichGO_result)
dev.off()

#Llama and Camel#################################################################################
# Import DEG files
# Define viral genes
virus_genes <- c("mers-3UTR","mers-orf8b","mers-N","mers-M","mers-E","mers-orf5","mers-orf4b","mers-orf4a","mers-orf3","mers-S","mers-orf1ab","mers-5UTR",
                 "camel229E-3UTR","camel229E-ORF8","camel229E-N","camel229E-M","camel229E-E","camel229E-ORF4","camel229E-S","camel229E-ORF1a", "camel229E-ORF1ab","camel229E-5UTR")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/DEG/Ferus_both/Celltype")
# import files
filelist = list.files(pattern="*.csv", full.names = TRUE)
datalist <- lapply(filelist, FUN= read.csv, na.strings="")

# Modify datalist
names(datalist) <- gsub(".csv", "", basename(filelist))
datalist <- lapply(datalist, function(x) x[!x$SYMBOL %in% virus_genes, ]) # Remove viral genes from files
datalist_entrezid <- lapply(datalist, function(x) {x$ENTREZID <- as.character(x$ENTREZID); return(x)})

#deleting the empty lists
datalist_entrezid[c(1, 15, 29, 34)] = NULL
filelist <- filelist[-c(1, 15, 29, 34)]


# Add variable information
datalist_entrezid <- lapply(names(datalist_entrezid), function(x) {datalist_entrezid[[x]]$name = x; datalist_entrezid[[x]]}) # Add column "name" with dataframe name
datalist_entrezid <- lapply(datalist_entrezid, tidyr::separate, name, c("CellType", "Species", "Condition", "Status"), remove = FALSE, sep = "_")
datalist_entrezid <- lapply(datalist_entrezid, function(x) {x$expr <- ifelse(x$avg_log2FC > 0, "upregulated", "downregulated"); return(x)}) # Annotate if gene is up or downregulated

names(datalist_entrezid) <- gsub(".csv", "", basename(filelist))


# data_entrezid <- lapply(datalist_entrezid, function(x) {x$ENTREZID})
# 
# ck <- compareCluster(geneCluster = data_entrezid, fun = enrichGO, OrgDb='org.Hs.eg.db')
# ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
# head(ck) 
# cnetplot(ck)

#get enrichgo results
enrichGO_result <- compareCluster(ENTREZID~expr+CellType+Condition+Status+Species, data=rbindlist(datalist_entrezid), fun="enrichGO", ont = "BP", OrgDb='org.Hs.eg.db', readable = T)

enrichGO_result_1 <- enrichGO_result@compareClusterResult

number_terms <- enrichGO_result_1 %>%
  group_by(expr, CellType, Condition, Status, Species) %>%
  dplyr::count(ID) %>%
  summarize(Sum = sum(n),
            n = n,
            expr = expr,
            CellType = CellType,
            Condition = Condition,
            Status = Status,
            Species = Species)

number_terms_uniq <- number_terms[!duplicated(number_terms), ]

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Enrich_terms")
write.csv(number_terms_uniq,"EnrichGo_across_all_conditions_count_terms.csv")


enrichGO_result_1 <- enrichGO_result_1
enrichGO_result_2 <- enrichGO_result
enrichGO_result_3 <- enrichGO_result
enrichGO_result_4 <- enrichGO_result


setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/")
write.csv(enrichGO_result@compareClusterResult,"EnrichGo_across_all_conditions_dotplot_terms.csv")


enrichGO_result@compareClusterResult <- enrichGO_result@compareClusterResult[order(-enrichGO_result@compareClusterResult$p.adjust),]

p1 <- as.data.frame(unique(enrichGO_result@compareClusterResult$Description))
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Enrich_terms")
write.csv(p1,"EnrichGo_across_all_conditions_dotplot_unique_terms.csv")



enrichGO_result_1@compareClusterResult <- filter(enrichGO_result_1@compareClusterResult, grepl('immun', enrichGO_result_1@compareClusterResult$Description))
enrichGO_result_2@compareClusterResult <- filter(enrichGO_result_2@compareClusterResult, grepl('virus|viral', enrichGO_result_2@compareClusterResult$Description))
enrichGO_result_3@compareClusterResult <- filter(enrichGO_result_3@compareClusterResult, grepl('apop', enrichGO_result_3@compareClusterResult$Description))
enrichGO_result_4@compareClusterResult <- filter(enrichGO_result_4@compareClusterResult, grepl('immun|virus|viral|apop', enrichGO_result_4@compareClusterResult$Description))



p2 <- as.data.frame(unique(enrichGO_result_4@compareClusterResult$Description))
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/results/Enrich_terms")
write.csv(p2,"EnrichGo_across_all_conditions_dotplot_immune_viral_apop_terms.csv")

enrichGO_result@compareClusterResult <- enrichGO_result@compareClusterResult[order(-enrichGO_result@compareClusterResult$p.adjust),]

#enrich go all terms
p.enrichGO_result <- dotplot(enrichGO_result, x="CellType", showCategory = 1, font.size = 16) +
  facet_grid(.~expr ~Species + Condition + Status) + theme_classic() + 
  theme(axis.text = element_text(size = 14, family = "sans", ),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = 14, family = "sans"),
        legend.title = element_text(size = 14, family = "sans"),
        axis.title = element_text(size = 14, family = "sans", face = "bold"),
        axis.line = element_line(linetype = "solid", colour = "black", size=1),
        strip.text = element_text(size = 16)) +
  scale_colour_gradientn(colours = c("red", "lightblue"))


p.enrichGO_result_1 <- dotplot(enrichGO_result_1, x="CellType", showCategory = 1, font.size = 16) +
  facet_grid(.~expr ~Species + Condition + Status) + theme_classic() + 
  theme(axis.text = element_text(size = 14, family = "sans", ),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = 14, family = "sans"),
        legend.title = element_text(size = 14, family = "sans"),
        axis.title = element_text(size = 14, family = "sans", face = "bold"),
        axis.line = element_line(linetype = "solid", colour = "black", size=1),
        strip.text = element_text(size = 16)) +
  scale_colour_gradientn(colours = c("red", "lightblue"))

p.enrichGO_result_2 <- dotplot(enrichGO_result_2, x="CellType", showCategory = 1, font.size = 16) +
  facet_grid(.~expr ~Species + Condition + Status) + theme_classic() + 
  theme(axis.text = element_text(size = 14, family = "sans", ),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = 14, family = "sans"),
        legend.title = element_text(size = 14, family = "sans"),
        axis.title = element_text(size = 14, family = "sans", face = "bold"),
        axis.line = element_line(linetype = "solid", colour = "black", size=1),
        strip.text = element_text(size = 16)) +
  scale_colour_gradientn(colours = c("red", "lightblue"))

p.enrichGO_result_3 <- dotplot(enrichGO_result_3, x="CellType", showCategory = 1, font.size = 16) +
  facet_grid(.~expr ~Species + Condition + Status) + theme_classic() + 
  theme(axis.text = element_text(size = 14, family = "sans", ),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = 14, family = "sans"),
        legend.title = element_text(size = 14, family = "sans"),
        axis.title = element_text(size = 14, family = "sans", face = "bold"),
        axis.line = element_line(linetype = "solid", colour = "black", size=1),
        strip.text = element_text(size = 16)) +
  scale_colour_gradientn(colours = c("red", "lightblue"))

p.enrichGO_result_4 <- dotplot(enrichGO_result_4, x="CellType", showCategory = 1, font.size = 16) +
  facet_grid(.~expr ~Species + Condition + Status) + theme_classic() + 
  theme(axis.text = element_text(size = 14, family = "sans", ),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = 14, family = "sans"),
        legend.title = element_text(size = 14, family = "sans"),
        axis.title = element_text(size = 14, family = "sans", face = "bold"),
        axis.line = element_line(linetype = "solid", colour = "black", size=1),
        strip.text = element_text(size = 16)) +
  scale_colour_gradientn(colours = c("red", "lightblue"))

p.enrichGO_result_1
?compareCluster

?dotplot

p.enrichGO_result_2 + theme(legend.position = "top")
p.enrichGO_result_3 + theme(legend.position = "top")


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/dotplot_across_all_conditions.pdf",
    width=12, height=35)
p.enrichGO_result + theme(legend.position = "top")
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/dotplot_across_all_conditions_immune.pdf",
    width=10, height=35)
p.enrichGO_result_1 + theme(legend.position = "top")
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/dotplot_across_all_conditions_vir.pdf",
    width=10, height=35)
p.enrichGO_result_2 + theme(legend.position = "top")
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/dotplot_across_all_conditions_apop.pdf",
    width=10, height=35)
p.enrichGO_result_3 + theme(legend.position = "top")
dev.off()

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/dotplot_across_all_conditions_vir_immune_apop.pdf",
    width=12, height=35)
p.enrichGO_result_4 + theme(legend.position = "top")
dev.off()


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/cnet_across_all_conditions.pdf",
    width=10, height=10)
cnetplot(enrichGO_result)
dev.off()

