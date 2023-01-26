#Plots

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
library(ComplexHeatmap)
library(ggplot2)
library(ggcorrplot)
library(scales)
library(grDevices)
library(colorspace)
library(forcats)

#Treatment colors
cols_treat = c("#E09992","#3DBCC2")
show_col(cols_treat)

#Status colors
#Infected, Uninfected
cols_stat = c('#9DCD84', '#089392')
show_col(cols_stat)
#Infected, Uninfected, Bystander
cols_stat_2 = c("#519872", "#CD4631","#DEA47E")
show_col(cols_stat_2)

#Celltype colors
cols_cell = hcl.colors(7, "Temps")
show_col(cols_cell)

################################################################################################################
#Import data from all three species
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens")
Human_MERS_combined <- LoadH5Seurat("MERS.combined.h5seurat")
Human_HCoV_combined <- LoadH5Seurat("HCoV.combined.h5seurat")
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/V.pacos")
Llama_MERS_combined <- LoadH5Seurat("MERS.combined.h5seurat")
Llama_camel229E_combined <- LoadH5Seurat("camel229E.combined.h5seurat")
setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.ferus")
Ferus_MERS_combined <- LoadH5Seurat("MERS.combined.h5seurat")
Ferus_camel229E_combined <- LoadH5Seurat("camel229E.combined.h5seurat")
#################################################################################################################
head(Human_MERS_combined@meta.data)
head(Llama_MERS_combined@meta.data)
head(Ferus_MERS_combined@meta.data)


Human_MERS_combined <- Human_MERS_combined@meta.data %>% select(c("Status", "Treat", "Phase", "Total.Viral.Counts", "Percent.Viral", "nCount_RNA")) 
Human_HCoV_combined <- Human_HCoV_combined@meta.data %>% select(c("Status", "Treat", "Phase", "Total.Viral.Counts", "Percent.Viral", "nCount_RNA"))
Llama_MERS_combined <- Llama_MERS_combined@meta.data %>% select(c("Status", "Treat", "Phase", "Total.Viral.Counts", "Percent.Viral", "nCount_RNA"))
Llama_camel229E_combined <- Llama_camel229E_combined@meta.data %>% select(c("Status", "Treat", "Phase", "Total.Viral.Counts", "Percent.Viral", "nCount_RNA"))
Ferus_MERS_combined <- Ferus_MERS_combined@meta.data %>% select(c("Status", "Treat", "Phase", "Total.Viral.Counts", "Percent.Viral", "nCount_RNA"))
Ferus_camel229E_combined <- Ferus_camel229E_combined@meta.data %>% select(c("Status", "Treat", "Phase", "Total.Viral.Counts", "Percent.Viral", "nCount_RNA"))


################################################################################
#Plot1 - #How many cells to we have in each treatment group?
Human_MERS_combined_groups <- Human_MERS_combined %>%
  group_by(Status, Treat) %>%
  count(Status, Treat) %>% 
  mutate(Species="Human", Virus="MERS") 

Human_HCoV_combined_groups <- Human_HCoV_combined %>%
  group_by(Status, Treat) %>%
  count(Status, Treat) %>% 
  mutate(Species="Human", Virus="HCoV") 

Llama_MERS_combined_groups <- Llama_MERS_combined %>%
  group_by(Status, Treat) %>%
  count(Status, Treat) %>% 
  mutate(Species="Llama", Virus="MERS")

Llama_camel229E_combined_groups <- Llama_camel229E_combined %>%
  group_by(Status, Treat) %>%
  count(Status, Treat) %>% 
  mutate(Species="Llama", Virus="camel229E")

Ferus_MERS_combined_groups <- Ferus_MERS_combined %>%
  group_by(Status, Treat) %>%
  count(Status, Treat) %>% 
  mutate(Species="Ferus", Virus="MERS")

Ferus_camel229E_combined_groups <- Ferus_camel229E_combined %>%
  group_by(Status, Treat) %>%
  count(Status, Treat) %>% 
  mutate(Species="Ferus", Virus="camel229E")

plot1 <- bind_rows(Human_MERS_combined_groups,Human_HCoV_combined_groups, 
                         Ferus_MERS_combined_groups, Ferus_camel229E_combined_groups,
                         Llama_MERS_combined_groups, Llama_camel229E_combined_groups)

human_plot1 <- bind_rows(Human_MERS_combined_groups,Human_HCoV_combined_groups)
human_plot1 <- human_plot1 %>%
  group_by(Virus, Treat) %>%
  summarize(Sum = sum(n),
            n = n,
            Treat = Treat,
            Virus = Virus,
            Status = Status) %>%
  mutate("Percent" = n/Sum*100) %>%
  mutate(across(where(is.numeric), ~ round(., 0)))



ferus_plot1 <- bind_rows(Ferus_MERS_combined_groups, Ferus_camel229E_combined_groups)
ferus_plot1 <- ferus_plot1 %>%
  group_by(Virus, Treat) %>%
  summarize(Sum = sum(n),
            n = n,
            Treat = Treat,
            Virus = Virus,
            Status = Status) %>%
  mutate("Percent" = n/Sum*100) %>%
  mutate(across(where(is.numeric), ~ round(., 0)))


llama_plot1 <- bind_rows(Llama_MERS_combined_groups, Llama_camel229E_combined_groups)
llama_plot1 <- llama_plot1 %>%
  group_by(Virus, Treat) %>%
  summarize(Sum = sum(n),
            n = n,
            Treat = Treat,
            Virus = Virus,
            Status = Status) %>%
  mutate("Percent" = n/Sum*100) %>%
  mutate(across(where(is.numeric), ~ round(., 0)))


p1 <- human_plot1 %>%
  ggplot(aes(x= Treat, y = n, fill = Status))+
  geom_col()+
  geom_text(aes(label = Percent), size = 4, vjust = -1, colour = "brown", position = "stack")+
  geom_text(aes(label = n), size = 4, vjust = -2.5, position = "stack")+
  ggtitle("Homo sapiens") +
  ylim(0, 5500)+
  ylab("Number of cells")+
  xlab("Treatment")+
  scale_fill_manual(values = cols_stat)+
  theme(#axis.text.x=element_blank(), #remove x axis labels
    axis.ticks.x=element_blank(),
    plot.title = element_text(size = rel(2)),
    legend.position="right") +
  facet_grid(cols = vars(Virus)) +
  NoLegend()



p2 <- ferus_plot1 %>% ggplot(aes(x= Treat, y = n, fill = Status))+
  geom_col()+
  geom_text(aes(label = Percent), size = 4, vjust = -1, colour = "brown", position = "stack")+
  geom_text(aes(label = n), size = 4, vjust = -2.5, position = "stack")+
  ggtitle("Camelus bactrianus") +
  ylim(0, 5500)+
  ylab("Number of cells")+
  xlab("Treatment")+
  scale_fill_manual(values = cols_stat)+
  theme(#axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = rel(2)),
        legend.position="right") +
  facet_grid(cols = vars(Virus))+
  NoLegend()

p3 <- llama_plot1 %>% ggplot(aes(x= Treat, y = n, fill = Status))+
  geom_col()+
  geom_text(aes(label = Percent), size = 4, vjust = 1, colour = "brown", position = "stack")+
  geom_text(aes(label = n), size = 4, vjust = -2.5, position = "stack")+
  ggtitle("Llama glama") +
  ylim(0, 5500)+
  ylab("Number of cells")+
  xlab("Treatment")+
  scale_fill_manual(values = cols_stat)+
  theme(#axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = rel(2)),
        legend.position="right") +
  NoLegend()


grid.arrange(p1, p2, p3, nrow = 1)

############################################################################################################
#Human_both_virus
#What celltype clusters to we have in each sample?

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens")
Human_clustered <- LoadH5Seurat("Human.clustered.h5seurat")

head(Human_clustered@meta.data)

Idents(Human_clustered) <- "celltype"
Human_clustered$celltype <- Idents(Human_clustered)

Human_clustered <- Human_clustered@meta.data %>%  select(c("Percent.MERS", "Percent.hcov", "Status", "Treat", "celltype"))
Human_clustered <- rename(Human_clustered, Celltype = celltype)

Human_clustered <- Human_clustered %>%
  mutate(Cilia_status = case_when(Celltype == "Secretory" ~ "Non-ciliated",
                                  Celltype == "Basal" ~ "Non-ciliated",
                                  Celltype == "Club" ~ "Non-ciliated",
                                  Celltype == "Deuterosomal" ~ "Non-ciliated",
                                  Celltype == "Ionocytes" ~ "Non-ciliated",
                                  Celltype == "Suprabasal" ~ "Non-ciliated",
                                  Celltype == "Ciliated" ~ "Ciliated"))


#Grouping by celltype and treatment allos wo count the cells in each celltype per treat
Human_celltype <- Human_clustered %>%
  group_by(Cilia_status, Treat) %>%
  dplyr::count(Status) %>%
  summarize(Sum = sum(n),
            n = n,
            Cilia_status = Cilia_status,
            Treat = Treat,
            Status = Status) %>%
  mutate("Percentage" = n/Sum*100)


# Reordering Levels
Human_celltype$Treat <- factor(Human_celltype$Treat, levels=c('Mock','HCoV-229E','MERS-CoV'))

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Human_both_virus/fraction_of_cells_per_celltype.pdf",
    width=6, height=5)
Human_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = Cilia_status, y = n, fill = Status))+
  geom_bar(position = "stack", stat = "identity", width = 0.8) +
  #scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  #geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_stat_2)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="right") +
  facet_grid(cols = vars(Treat))
dev.off()


#Viral Reads
Human_viral_celltype <- Human_clustered %>%
  group_by(Celltype, Treat, Status) %>%
  count(Status)%>%
  filter(Treat == "MERS" | Treat == "HCOV") %>%
  group_by(Celltype, Treat) %>%  
  summarize(Sum = sum(n),
            n = n,
            Celltype = Celltype,
            Treat = Treat,
            Status = Status) %>%
  mutate("Percentage" = n/Sum*100)

Human_infected_celltype_table <- Human_viral_celltype %>% 
  filter(Status == "Infected") %>%
  select(c("Celltype", "Treat", "n"))

#Infected cells Human
Human_infected_celltype_table %>%
  ggplot()+
  geom_col(mapping = aes(x = Celltype, y = n, fill = Treat), width = 0.8)+
  scale_fill_manual(values = cols3)+
  ylab("Infected cells")+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 10),
        legend.position = "bottom")

Human_viral_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = Treat, y = Percentage, fill = Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Percentage")+
  xlab("Condition")+
  facet_grid(cols = vars(Celltype))+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 4, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols7)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 10),
        legend.position = "bottom")



############################################################################################################
#Alpaca_both_virus
#What celltype clusters to we have in each sample?

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/V.pacos")
Alpaca.clustered <- LoadH5Seurat("Alpaca.clustered.h5seurat")

head(Alpaca.clustered@meta.data)

#Alpaca.clustered$celltype <- Idents(Alpaca.clustered)

Alpaca.clustered <- Alpaca.clustered@meta.data %>%  dplyr::select(c("Percent.MERS", "Percent.camel229E", "Status", "Treat", "celltype"))

cols = c("#8A6E7E", "#40679B", "#089392", "#40B48B", "#9DCD84", "#EAE29C")
show_col(cols)
Alpaca.clustered <- rename(Alpaca.clustered, Celltype = celltype)

Alpaca.clustered <- Alpaca.clustered %>%
  mutate(Cilia_status = case_when(Celltype == "Secretory" ~ "Non-ciliated",
                                  Celltype == "Basal" ~ "Non-ciliated",
                                  Celltype == "Club" ~ "Non-ciliated",
                                  Celltype == "Llama Cluster 1" ~ "Non-ciliated",
                                  Celltype == "Llama Cluster 2" ~ "Non-ciliated",
                                  Celltype == "Ciliated" ~ "Ciliated"))

#Grouping by celltype and treatment allos wo count the cells in each celltype per treat
Alpaca_celltype <- Alpaca.clustered %>%
  group_by(Cilia_status, Treat) %>%
  dplyr::count(Status) %>%
  summarize(Sum = sum(n),
            n = n,
            Cilia_status = Cilia_status,
            Treat = Treat,
            Status = Status) %>%
  mutate("Percentage" = n/Sum*100)


# Reordering Levels
Alpaca_celltype$Treat <- factor(Alpaca_celltype$Treat, levels=c('Mock','dcCoV-ACN4','MERS-CoV'))
#Alpaca_celltype$celltype <- factor(Alpaca_celltype$celltype, levels=c('Llama Cluster 1','Llama Cluster 2','Secretory','Ciliated','Club', 'Basal'))

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Llama_alpaca_both_virus/fraction_of_cells_per_celltype.pdf",
    width=6, height=5)
Alpaca_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = Cilia_status, y = n, fill = Status))+
  geom_bar(position = "stack", stat = "identity", width = 0.8) +
  #scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  #geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_stat_2)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="right") +
  facet_grid(cols = vars(Treat))
dev.off()






Alpaca_viral_celltype <- Alpaca.clustered %>%
  group_by(Celltype, Treat, Status) %>%
  count(Status)%>%
  filter(Treat == "MERS-CoV" | Treat == "dcCoV-ACN4") %>%
  group_by(Celltype, Treat) %>%
  summarize(Sum = sum(n),
            n = n,
            Celltype = Celltype,
            Treat = Treat,
            Status = Status) %>%
  mutate("Percentage" = n/Sum*100)



Alpaca_viral_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = Treat, y = Percentage, fill = Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Percentage")+
  xlab("Condition")+
  facet_grid(cols = vars(Celltype))+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 4, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_stat)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 10),
        legend.position = "bottom")

############################################################################################################
#Ferus_both_virus
#What celltype clusters to we have in each sample?

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.ferus")
Ferus.clustered <- LoadH5Seurat("Ferus.clustered.h5seurat")

head(Ferus.clustered@meta.data)


Ferus_clustered <- Ferus.clustered@meta.data %>%  dplyr::select(c("Percent.MERS", "Percent.camel229E", "Status", "Treat", "celltype"))
levels(Ferus.clustered)

cols_cells = c("#807EBF", "#7EA8BE", "#089392", "#40B48B", "#9DCD84", "#EAE29C")
show_col(cols)
Ferus_clustered <- rename(Ferus_clustered, Celltype = celltype)


#Grouping by celltype and treatment allos wo count the cells in each celltype per treat

Ferus_clustered <- Ferus_clustered %>%
  mutate(Cilia_status = case_when(Celltype == "Secretory" ~ "Non-ciliated",
                                  Celltype == "Basal" ~ "Non-ciliated",
                                  Celltype == "Club" ~ "Non-ciliated",
                                  Celltype == "Camel Cluster 1" ~ "Non-ciliated",
                                  Celltype == "Camel Cluster 2" ~ "Non-ciliated",
                                  Celltype == "Ciliated" ~ "Ciliated"))

#Grouping by celltype and treatment allos wo count the cells in each celltype per treat
Ferus_clustered <- Ferus_clustered %>%
  group_by(Cilia_status, Treat) %>%
  dplyr::count(Status) %>%
  summarize(Sum = sum(n),
            n = n,
            Cilia_status = Cilia_status,
            Treat = Treat,
            Status = Status) %>%
  mutate("Percentage" = n/Sum*100)


# Reordering Levels
Ferus_clustered$Treat <- factor(Ferus_clustered$Treat, levels=c('Mock','dcCoV-ACN4','MERS-CoV'))
#Alpaca_celltype$celltype <- factor(Alpaca_celltype$celltype, levels=c('Llama Cluster 1','Llama Cluster 2','Secretory','Ciliated','Club', 'Basal'))

pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/fraction_of_inf_per_celia_status.pdf",
    width=6, height=5)
Ferus_clustered %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = Cilia_status, y = n, fill = Status))+
  geom_bar(position = "stack", stat = "identity", width = 0.8) +
  #scale_y_continuous(labels = scales::percent_format())+
  ylab("Fraction of cells")+
  xlab("Condition")+
  #geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 5, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = cols_stat_2)+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 12, angle = 45),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, vjust = -0.5, size = 20, margin = margin(r = 15)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position="right") +
  facet_grid(cols = vars(Treat))
dev.off()


#Viral reads in the celltypes?


Ferus_viral_celltype <- Ferus_clustered %>%
  group_by(Celltype, Treat, Status) %>%
  count(Status)%>%
  filter(Treat == "MERS-CoV" | Treat == "dcCoV-ACN4") %>%
  group_by(Celltype, Treat) %>%  
  summarize(Sum = sum(n),
            n = n,
            Celltype = Celltype,
            Treat = Treat,
            Status = Status) %>%
  mutate("Percentage" = n/Sum*100)


pdf(file="C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/figures/Camel_ferus_both_virus/Proportion_Inf_celltypes.pdf",
    width=5, height=5)
Ferus_viral_celltype %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  ggplot(aes(x = Celltype, y = Percentage, fill = Status))+
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  ylab("Percentage")+
  xlab("Condition")+
  facet_grid(cols = vars(Treat))+
  geom_text(aes(label = Percentage), position = position_fill(vjust = 0.5), col = "black", size = 4, hjust = 0.5, check_overlap=TRUE)+
  scale_fill_manual(values = c("#519872", "#CD4631"))+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 10),
        legend.position = "bottom")
dev.off()

############################################################################################################


setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/H.sapiens")
Human_MERS_clustered <- LoadH5Seurat("MERS.clustered.h5seurat")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/V.pacos")
Llama_MERS_clustered <- LoadH5Seurat("MERS.clustered.h5seurat")

setwd("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/scRNA_seq/seurat_objects/C.ferus")
Ferus_MERS_clustered <- LoadH5Seurat("MERS.clustered.h5seurat")

head(Human_MERS_clustered@meta.data)
head(Llama_MERS_clustered@meta.data)
head(Ferus_MERS_clustered@meta.data)

Human_MERS_clustered$celltype <- Idents(Human_MERS_clustered)
Llama_MERS_clustered$celltype <- Idents(Llama_MERS_clustered)
Ferus_MERS_clustered$celltype <- Idents(Ferus_MERS_clustered)


Human_MERS_clustered<- Human_MERS_clustered@meta.data %>% select(c("Percent.Viral", "Status", "Treat", "celltype")) 
Llama_MERS_clustered <- Llama_MERS_clustered@meta.data %>% select(c("Percent.Viral","Status", "Treat", "celltype"))
Ferus_MERS_clustered <- Ferus_MERS_clustered@meta.data %>% select(c("Percent.Viral","Status", "Treat", "celltype"))


Human_MERS_combined_celltype <- Human_MERS_clustered %>%
  group_by(celltype, Treat) %>%
  count(celltype) %>% 
  mutate(Species="Human")


Llama_MERS_combined_celltype <- Llama_MERS_clustered %>%
  group_by(celltype, Treat) %>%
  count(celltype) %>% 
  mutate(Species="Llama")


Ferus_MERS_combined_celltype <- Ferus_MERS_clustered %>%
  group_by(celltype, Treat) %>%
  count(celltype) %>% 
  mutate(Species="Ferus")

plot3 <- bind_rows(Human_MERS_combined_celltype, Llama_MERS_combined_celltype, Ferus_MERS_combined_celltype)

plotX <- plot3 %>%
  group_by(Species, Treat) %>%
  summarize(Sum = sum(n),
            n = n,
            celltype = celltype,
            Treat = Treat,
            Species = Species) %>%
  mutate("Percent" = n/Sum*100)

plotX <- plotX %>% mutate("pos" = cumsum(Percent)-0.5*Percent)

plotX <- plotX %>% mutate(across(where(is.numeric), ~ round(., 0)))
plotX %>% 
  ggplot(aes(x = Treat, y = Percent, fill = celltype))+
  geom_col() +
  geom_text(aes(label = Percent, y = pos), size = 4)+
  scale_fill_manual(values = cols)+
  facet_grid(cols = vars(Species))+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = rel(2)),
        legend.position="bottom")

?position_stack()

?geom_col()
#How are the Viral reads distributed in the different cell types per species?

Human_MERS_combined_groups <- Human_MERS_clustered %>%
  group_by(celltype) %>%
  summarise_at(vars(Percent.Viral), list(average_viral_percentage = mean)) %>%
  mutate(Species="Human")


Llama_MERS_combined_groups <- Llama_MERS_clustered %>%
  group_by(celltype) %>%
  summarise_at(vars(Percent.Viral), list(average_viral_percentage = mean)) %>%
  mutate(Species="Llama")


Ferus_MERS_combined_groups <- Ferus_MERS_clustered %>%
  group_by(celltype) %>%
  summarise_at(vars(Percent.Viral), list(average_viral_percentage = mean)) %>%
  mutate(Species="Ferus")


plot4 <- bind_rows(Human_MERS_combined_groups, Llama_MERS_combined_groups, Ferus_MERS_combined_groups)


plot4 %>% 
  ggplot(aes(x = celltype, y = average_viral_percentage, fill = celltype))+
  geom_col() +
  scale_fill_manual(values = cols)+
  facet_grid(cols = vars(Species))+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = rel(2)),
        legend.position="bottom")



