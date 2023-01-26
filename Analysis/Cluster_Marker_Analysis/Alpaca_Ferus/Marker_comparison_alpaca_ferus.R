#Marker gene comparison

library(ggvenn)
library(dplyr)
library(tidyr)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(patchwork)
library(readr)
library(gdata)
library(Matrix)
library(writexl)
library(grid)
library(viridis)
library(gridExtra)
library(superheat)
library(heatmap3)

cols = viridis(5)

#########################################################################################################################
#Marker comparison Alpaca-Human

#Import the Alpaca_Human shared genes per cluster from python script
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Alpaca_Ferus_Comparison")
temp = list.files(pattern="*.txt")
names(temp) <- temp
df <- lapply(temp, FUN=read.table, header=FALSE)

name <- strsplit(names(df), '[_.]')
species_A <- lapply(name, `[`, 1) %>% unlist()
species_B <- lapply(name, `[`, 3) %>% unlist()
cluster_A <- lapply(name, `[`, 2) %>% unlist()
cluster_B <- lapply(name, `[`, 4) %>% unlist()

filelist <- mapply(cbind, df, "Species_A"=species_A, SIMPLIFY=F)
filelist <- mapply(cbind, filelist, "Species_B"=species_B, SIMPLIFY=F)
filelist <- mapply(cbind, filelist, "Cluster_A"=cluster_A, SIMPLIFY=F)
filelist <- mapply(cbind, filelist, "Cluster_B"=cluster_B, SIMPLIFY=F)

# Bind all txt files together
Alpaca_Ferus <- bind_rows(filelist, .id = "column_label")

DF <-Alpaca_Ferus %>% mutate(Alpaca = case_when(Cluster_A == 0 ~ "Cluster 0",
                                               Cluster_A == 1 ~ "Cluster 1",
                                               Cluster_A == 2 ~ "Cluster 2",
                                               Cluster_A == 3 ~ "Cluster 3",
                                               Cluster_A == 4 ~ "Cluster 4",
                                               Cluster_A == 5 ~ "Cluster 5"),
                             Ferus = case_when(Cluster_B == 0 ~ "Cluster 0",
                                               Cluster_B == 1 ~ "Cluster 1",
                                               Cluster_B == 2 ~ "Cluster 2",
                                               Cluster_B == 3 ~ "Cluster 3",
                                               Cluster_B == 4 ~ "Cluster 4",
                                               Cluster_B == 5 ~ "Cluster 5"))


#Having the data as a dataframe
Heatmap <- DF %>%
  group_by(Alpaca, Ferus) %>%
  count(column_label) %>%
  select(Alpaca, n)

Heatmap[is.na(Heatmap)] = 0

#Having the data as a matrix
Plot <- spread(Heatmap, Ferus, n)
Plot[is.na(Plot)] = 0
Plot <- data.matrix(Plot, rownames.force = TRUE)
rownames(Plot) <- c("Cluster 0", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")
#colnames(Plot) <- c("trash", "Cluster 0", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6")
Plot<-Plot[,-1]



#Base-R Heatmap with dendrogram
heatmap(Plot, scale = "column", margins = c(9,9), xlab = "Ferus", ylab = "Alpaca", col= viridis(5))
legend(x="bottomright", legend=c("low", "high"), 
       fill=viridis(5))


#GGplot Heatmap
Heatmap %>% ggplot(aes(x = Alpaca, y = Ferus, fill = n))+
  geom_tile()+
  scale_fill_viridis()

#Heatmap3 Heatmap
heatmap3(Plot, scale = "column")



#Superheat Heatmap with dendrogram
pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/cluster_comparison_ferus.pdf",
    width=11, height=6)
superheat(X = Plot,
          X.text = Plot,
          X.text.size = 6,
          heat.na.col = "grey",
          row.dendrogram = TRUE, 
          col.dendrogram = TRUE, 
          grid.hline = FALSE, 
          grid.vline = FALSE,
          column.title = "Ferus",
          column.title.size = 7,
          row.title = "Alpaca",
          row.title.size = 7
)
dev.off()


#########################################################################################################################
#Import the Alpaca genes
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Single_cluster_data")
temp = list.files(pattern="*.txt")
names(temp) <- temp
df <- lapply(temp, FUN=read.table, header=TRUE)

name <- strsplit(names(df), '[_.]')
species <- lapply(name, `[`, 1)
virus <- lapply(name, `[`, 2)
cluster <- lapply(name, `[`, 4) %>% unlist()
filelist <- mapply(cbind, df, "Cluster"=cluster, SIMPLIFY=F)
filelist <- mapply(cbind, filelist, "Virus"=virus, SIMPLIFY=F)
filelist <- mapply(cbind, filelist, "Species"=species, SIMPLIFY=F)

# Bind all txt files together
Alpaca_Ferus <- bind_rows(filelist, .id = "column_label")
#########################################################################################################################

cluster0_alpaca <- Alpaca_Ferus %>%
  subset(Species == "Alpaca") %>%
  subset(Cluster == 0) %>%
  select(2,)
cluster1_alpaca <- Alpaca_Ferus %>%
  subset(Species == "Alpaca") %>%
  subset(Cluster == 1) %>%
  select(2,)
cluster2_alpaca <- Alpaca_Ferus %>%
  subset(Species == "Alpaca") %>%
  subset(Cluster == 2) %>%
  select(2,)
cluster3_alpaca <- Alpaca_Ferus %>%
  subset(Species == "Alpaca") %>%
  subset(Cluster == 3) %>%
  select(2,)
cluster4_alpaca <- Alpaca_Ferus %>%
  subset(Species == "Alpaca") %>%
  subset(Cluster == 4) %>%
  select(2,)
# cluster5_alpaca <- Alpaca_Ferus %>%
#   subset(Species == "Alpaca") %>%
#   subset(Cluster == 5) %>%
#   select(2,)
  
alpaca_cluster0 <- cluster0_alpaca$gene
alpaca_cluster1 <- cluster1_alpaca$gene
alpaca_cluster2 <- cluster2_alpaca$gene
alpaca_cluster3 <- cluster3_alpaca$gene
alpaca_cluster4 <- cluster4_alpaca$gene
# alpaca_cluster5 <- cluster5_alpaca$gene

cluster0_Ferus <- Alpaca_Ferus %>%
  subset(Species == "Ferus") %>%
  subset(Cluster == 0) %>%
  select(2,)
cluster1_Ferus <- Alpaca_Ferus %>%
  subset(Species == "Ferus") %>%
  subset(Cluster == 1) %>%
  select(2,)
cluster2_Ferus <- Alpaca_Ferus %>%
  subset(Species == "Ferus") %>%
  subset(Cluster == 2) %>%
  select(2,)
cluster3_Ferus <- Alpaca_Ferus %>%
  subset(Species == "Ferus") %>%
  subset(Cluster == 3) %>%
  select(2,)
cluster4_Ferus <- Alpaca_Ferus %>%
  subset(Species == "Ferus") %>%
  subset(Cluster == 4) %>%
  select(2,)
cluster5_Ferus <- Alpaca_Ferus %>%
  subset(Species == "Ferus") %>%
  subset(Cluster == 5) %>%
  select(2,)


Ferus_cluster0 <- cluster0_Ferus$gene
Ferus_cluster1 <- cluster1_Ferus$gene
Ferus_cluster2 <- cluster2_Ferus$gene
Ferus_cluster3 <- cluster3_Ferus$gene
Ferus_cluster4 <- cluster4_Ferus$gene
Ferus_cluster5 <- cluster5_Ferus$gene


?ggvenn


#venn diagrams
genes <- list(alpaca_cluster0 = alpaca_cluster0,
              alpaca_cluster1 = alpaca_cluster1,
              alpaca_cluster2 = alpaca_cluster2,
              alpaca_cluster3 = alpaca_cluster3,
              alpaca_cluster4 = alpaca_cluster4,
              # alpaca_cluster5 = alpaca_cluster5,
              human_cluster0 = human_cluster0,
              human_cluster1 = human_cluster1,
              human_cluster2 = human_cluster2,
              human_cluster3 = human_cluster3,
              human_cluster4 = human_cluster4,
              human_cluster5 = human_cluster5)
                        


#alpaca Cluster 0

p1 <- ggvenn(genes, 
             c("alpaca_cluster0", "human_cluster0"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             show_percentage = FALSE, text_size = 5,
             stroke_size = 0.5, set_name_size = 2.5)

p2 <- ggvenn(genes, 
             c("alpaca_cluster0", "human_cluster1"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             show_percentage = FALSE, text_size = 5,
             stroke_size = 0.5, set_name_size = 2.5)

p3 <- ggvenn(genes, 
             c("alpaca_cluster0", "human_cluster2"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             show_percentage = FALSE, text_size = 5,
             stroke_size = 0.5, set_name_size = 2.5)

p4 <- ggvenn(genes, 
             c("alpaca_cluster0", "human_cluster3"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             show_percentage = FALSE, text_size = 5,
             stroke_size = 0.5, set_name_size = 2.5)

p5 <- ggvenn(genes, 
             c("alpaca_cluster0", "human_cluster4"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             show_percentage = FALSE, text_size = 5,
             stroke_size = 0.5, set_name_size = 2.5)

p6 <- ggvenn(genes, 
             c("alpaca_cluster0", "human_cluster5"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             show_percentage = FALSE, text_size = 5,
             stroke_size = 0.5, set_name_size = 2.5)

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)

########################################################################################################
#alpaca cluster 1
p1 <- ggvenn(genes, 
             c("alpaca_cluster1", "human_cluster0"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p2 <- ggvenn(genes, 
             c("alpaca_cluster1", "human_cluster1"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p3 <- ggvenn(genes, 
             c("alpaca_cluster1", "human_cluster2"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p4 <- ggvenn(genes, 
             c("alpaca_cluster1", "human_cluster3"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p5 <- ggvenn(genes, 
             c("alpaca_cluster1", "human_cluster4"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p6 <- ggvenn(genes, 
             c("alpaca_cluster1", "human_cluster5"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)

######################################################################################################
#alpaca cluster 2

p1 <- ggvenn(genes, 
             c("alpaca_cluster2", "human_cluster0"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p2 <- ggvenn(genes, 
             c("alpaca_cluster2", "human_cluster1"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p3 <- ggvenn(genes, 
             c("alpaca_cluster2", "human_cluster2"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p4 <- ggvenn(genes, 
             c("alpaca_cluster2", "human_cluster3"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p5 <- ggvenn(genes, 
             c("alpaca_cluster2", "human_cluster4"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p6 <- ggvenn(genes, 
             c("alpaca_cluster2", "human_cluster5"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)

######################################################################################################
#alpaca cluster 3

p1 <- ggvenn(genes, 
             c("alpaca_cluster3", "human_cluster0"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p2 <- ggvenn(genes, 
             c("alpaca_cluster3", "human_cluster1"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p3 <- ggvenn(genes, 
             c("alpaca_cluster3", "human_cluster2"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p4 <- ggvenn(genes, 
             c("alpaca_cluster3", "human_cluster3"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p5 <- ggvenn(genes, 
             c("alpaca_cluster3", "human_cluster4"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p6 <- ggvenn(genes, 
             c("alpaca_cluster3", "human_cluster5"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)

######################################################################################################
#alpaca cluster 4

p1 <- ggvenn(genes, 
             c("alpaca_cluster4", "human_cluster0"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p2 <- ggvenn(genes, 
             c("alpaca_cluster4", "human_cluster1"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p3 <- ggvenn(genes, 
             c("alpaca_cluster4", "human_cluster2"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p4 <- ggvenn(genes, 
             c("alpaca_cluster4", "human_cluster3"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p5 <- ggvenn(genes, 
             c("alpaca_cluster4", "human_cluster4"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p6 <- ggvenn(genes, 
             c("alpaca_cluster4", "human_cluster5"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)


##########################################################################################
#combinatorial comparison of all the human clusters with each other

p1 <- ggvenn(genes, 
             c("human_cluster0", "human_cluster1"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p2 <- ggvenn(genes, 
             c("human_cluster0", "human_cluster2"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p3 <- ggvenn(genes, 
             c("human_cluster0", "human_cluster3"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p4 <- ggvenn(genes, 
             c("human_cluster0", "human_cluster4"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p5 <- ggvenn(genes, 
             c("human_cluster0", "human_cluster5"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p6 <- ggvenn(genes, 
             c("human_cluster1", "human_cluster2"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p7 <- ggvenn(genes, 
             c("human_cluster1", "human_cluster3"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p8 <- ggvenn(genes, 
             c("human_cluster1", "human_cluster4"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p9 <- ggvenn(genes, 
             c("human_cluster1", "human_cluster5"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p10 <- ggvenn(genes, 
             c("human_cluster2", "human_cluster3"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p11 <- ggvenn(genes, 
             c("human_cluster2", "human_cluster4"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p12 <- ggvenn(genes, 
             c("human_cluster2", "human_cluster5"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p13 <- ggvenn(genes, 
              c("human_cluster3", "human_cluster4"),
              fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
              stroke_size = 0.5, set_name_size = 4)

p14 <- ggvenn(genes, 
              c("human_cluster3", "human_cluster5"),
              fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
              stroke_size = 0.5, set_name_size = 4)

p15 <- ggvenn(genes, 
              c("human_cluster4", "human_cluster5"),
              fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
              stroke_size = 0.5, set_name_size = 4)

grid.arrange(p1, p2, p3, p4, p5, nrow = 3)
grid.arrange(p6, p7, p8, p9, p10, nrow = 3)
grid.arrange(p11, p12, p13, p14, p15, nrow = 3)

##########################################################################################
#combinatorial comparison of all the alpaca clusters with each other

p1 <- ggvenn(genes, 
             c("alpaca_cluster0", "alpaca_cluster1"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p2 <- ggvenn(genes, 
             c("alpaca_cluster0", "alpaca_cluster2"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p3 <- ggvenn(genes, 
             c("alpaca_cluster0", "alpaca_cluster3"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p4 <- ggvenn(genes, 
             c("alpaca_cluster0", "alpaca_cluster4"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p5 <- ggvenn(genes, 
             c("alpaca_cluster1", "alpaca_cluster2"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p6 <- ggvenn(genes, 
             c("alpaca_cluster1", "alpaca_cluster3"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p7 <- ggvenn(genes, 
             c("alpaca_cluster1", "alpaca_cluster4"),
             fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             stroke_size = 0.5, set_name_size = 4)

p8 <- ggvenn(genes, 
              c("alpaca_cluster2", "alpaca_cluster3"),
              fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
              stroke_size = 0.5, set_name_size = 4)

p9 <- ggvenn(genes, 
              c("alpaca_cluster2", "alpaca_cluster4"),
              fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
              stroke_size = 0.5, set_name_size = 4)

p10 <- ggvenn(genes, 
              c("alpaca_cluster3", "alpaca_cluster4"),
              fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
              stroke_size = 0.5, set_name_size = 4)

grid.arrange(p1, p2, p3, p4, p5, nrow = 3)
grid.arrange(p6, p7, p8, p9, p10, nrow = 3)

