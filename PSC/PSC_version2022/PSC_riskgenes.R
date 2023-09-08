###################################################
# Title: PSC riskgenes in PSC-I vs NI
# Date: 25-07-2023
###################################################


######################
# Library
######################

library(ggplot2) 
library(dplyr) 
library(tidyr) 
library(stringr)
library(patchwork)
library(grid)
library(gridExtra) 
library(ggpubr)
library(RColorBrewer) 
library(Seurat)
library(readr)
library(MAST)
library(readxl)
library(openxlsx)
library(enrichR)
require(tidyverse)

######################
# Main codes
######################

# read in dataframes
epi <- readRDS("Nieuw/epi_azimuth_duox2.rds")
imm <- readRDS("Nieuw/imm_azimuth_with_plasma.rds")
str <- readRDS("Nieuw/stro_azimuth.rds")

DefaultAssay(epi) = "RNA"
DefaultAssay(imm) = "RNA"
DefaultAssay(str) = "RNA"
Idents(epi) <- "celltype.final"
Idents(imm) <- "predicted.cell_type.pred"
Idents(str) <- "predicted.cell_type.pred"

# create DE lists
data = str # Repeat for all three datasets (epi, stromal, immune)
celltypes <- as.character(unique(Idents(data)))
for (x in celltypes) {
  # safe pathnames
  path <- paste("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/DE_PSCINI_all/",x,".csv",sep="")
  # find markers
  x <- FindMarkers(data, subset.ident = x, group.by = "state", test.use = "MAST", ident.1 = "PSC-I", ident.2 = "PSC-NI")
  if(nrow(x) == 0)
  {
    #skip iteration because of zero significant genes
    next
  }
  #rest of iteration for case of no error
  # create a 'state' column
  x$state = NA
  x = as.data.frame(x)
  x$state = ifelse(x$avg_log2FC > 0,"PSC_I","PSC_NI")
  # create a 'gene' column
  x$gene = NA
  x$gene = rownames(x)
  # save dataframes
  write.csv(x, path)
}

## Repeat for all three datasets (epi, stromal, immune)


# create dataframe with all DE PSC riskgenes from PSC-I vs PSC-NI
PSC_riskgenes <- c("ATXN2","BACH2", "BCL2", "CCDC88B", "CCL20", 
                   "CD226","CD28","CLEC16A","CTLA4","CXCR2","FOXP1",
                   "FUT2","GPR35","HDAC7","HHEX","HLA-DRB1","IL2",
                   "IL21","IL2RA","MMEL1","MST1","NFKB1","NKX2-3",
                   "PRDX5","PRKCB","PRKD2","PSMG1","PTPN2","PUS10",
                   "PVT1","RFX4","RIC8B","SH2B3","SIK2","SOCS1",
                   "STRN4","TCF4", "TNFRSF14","UBASH3A")

setwd("Results/DE_2023/DE_PSCINI_all/")
file_list <- list.files(path="/Users/amberbangma/Documents/R/PSC/Results/DE_2023/DE_PSCINI_all/")
dataset <- data.frame()
for (i in 1:length(file_list)){
  temp_data <- read_csv(file_list[i]) #each file will be read in, specify which columns you need read in to avoid any errors
  temp_data$celltype <- sapply(strsplit(gsub(".csv", "", file_list[i]), "/"), function(x){x[1]}) #clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
  dataset <- rbind(dataset, temp_data) #for each iteration, bind the new data to the building dataset
}
dataset$...1 <- NULL

dataset <- dataset[dataset$gene %in% PSC_riskgenes, ]
PSCriskgenes_05 <- filter(dataset, dataset$p_val < 0.05)

# make celltype a factor variable
all <- PSCriskgenes_05
all <- all[,c(1,2,7,8)]
all$celltype <- as.factor(all$celltype)
# check the current levels
levels(all$celltype)
# reorder the levels
#all$celltype <- factor(all$celltype, levels = levels(all$celltype)[c(0-46)])
all$celltype <- factor(all$celltype, levels = levels(all$celltype)[c(23,20,21,10,2,3,4,5,40,7,8,9,11,32,24,28,30,13,14,6,42,43,44,
                                                                     45,35,27,31,33,34,19,16,26,22,41,29,37,12,36,38,39,25,17,18,1,15)])
levels(all$celltype)

# make figure
all$group <- cut(all$avg_log2FC, breaks = c(-Inf, -1, -0.5, 0, 0.5, 1, Inf),
                 labels = c("<-1", "(-1,-0.5]", "(-0.5,0]", "(0,0.5]", "(0.5,1]", ">1"))
#all <- subset(all, !(group %in% c("(-0.5,0]", "(0,0.5]")))

#plot the heatmap with p<0.05
all %>% 
  ggplot(aes(x=gene, y=celltype, fill = factor(group))) +
  coord_fixed(ratio = 1) +
  geom_tile() +
  scale_fill_manual(
    breaks = levels(all$group),
    values = c(
      "#225EA8",   # <-1
      "#8CC4F4",   # (-1,-0.5]
      "#C6DBEF",   # (-0.5,0]
      "#FF9999",   # (0,0.5]
      "#e82525",   # (0.5,1]
      "#9E1717"   # >1
    )) +
  ggtitle('')  +
  geom_tile( colour = 'white' ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y=NULL, x=NULL) +
  guides(fill=guide_legend(title="avg_log2FC")) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"),
        legend.spacing.y = unit(0.5, 'cm'))
ggsave("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/Riskgenes/PSCHeatmap_uncorr.pdf",  width=6, height=6)

#plot the heatmap with p<0.00128
all <- filter(all, all$p_val < 0.00128)
all %>% 
  ggplot(aes(x=gene, y=celltype, fill = factor(group))) +
  coord_fixed(ratio = 1) +
  geom_tile() +
  scale_fill_manual(
    breaks = levels(all$group),
    values = c(
    "#225EA8",   # <-1
    "#8CC4F4",   # (-1,-0.5]
    "#C6DBEF",   # (-0.5,0]
    "#FF9999",   # (0,0.5]
    "#e82525",   # (0.5,1]
    "#9E1717"   # >1
    )) +
  ggtitle('')  +
  geom_tile( colour = 'white' ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y=NULL, x=NULL) +
  guides(fill=guide_legend(title="avg_log2FC")) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"),
        legend.spacing.y = unit(0.5, 'cm'))
ggsave("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/Riskgenes/PSCHeatmap_corr.pdf",  width=6, height=6)

############ same for UC #################
# create DE lists
data = str # Repeat for all three datasets (epi, stromal, immune)
celltypes <- as.character(unique(Idents(data)))
for (x in celltypes) {
  # safe pathnames
  path <- paste("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/DE_UCINI_all/",x,".csv",sep="")
  # find markers
  x <- FindMarkers(data, subset.ident = x, group.by = "state", test.use = "MAST", ident.1 = "UC-I", ident.2 = "UC-NI")
  if(nrow(x) == 0)
  {
    #skip iteration because of zero significant genes
    next
  }
  #rest of iteration for case of no error
  # create a 'state' column
  x$state = NA
  x = as.data.frame(x)
  x$state = ifelse(x$avg_log2FC > 0,"UC_I","UC_NI")
  # create a 'gene' column
  x$gene = NA
  x$gene = rownames(x)
  # save dataframes
  write.csv(x, path)
}

setwd("Results/DE_2023/DE_UCINI_all/")
file_list <- list.files(path="/Users/amberbangma/Documents/R/PSC/Results/DE_2023/DE_UCINI_all/")
dataset <- data.frame()
for (i in 1:length(file_list)){
  temp_data <- read_csv(file_list[i]) #each file will be read in, specify which columns you need read in to avoid any errors
  temp_data$celltype <- sapply(strsplit(gsub(".csv", "", file_list[i]), "/"), function(x){x[1]}) #clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
  dataset <- rbind(dataset, temp_data) #for each iteration, bind the new data to the building dataset
}
dataset$...1 <- NULL

dataset <- dataset[dataset$gene %in% PSC_riskgenes, ]
PSCriskgenes_05 <- filter(dataset, dataset$p_val < 0.05)

# make celltype a factor variable
all <- PSCriskgenes_05
all <- all[,c(1,2,7,8)]
all$celltype <- as.factor(all$celltype)
# check the current levels
levels(all$celltype)
# reorder the levels
#all$celltype <- factor(all$celltype, levels = levels(all$celltype)[c(0-46)])
all$celltype <- factor(all$celltype, levels = levels(all$celltype)[c(23, 24, 25, 19, 20, 9, 2, 3, 4, 40, 6, 7, 8, 11, 26, 10, 30, 32, 13, 5, 41, 42, 43, 44, 45, 36, 29, 33, 34, 35, 18, 15, 21, 28, 22, 31, 12, 37, 38, 39, 27, 16, 17, 1, 14)])
levels(all$celltype)

# make figure
all$group <- cut(all$avg_log2FC, breaks = c(-Inf, -1, -0.5, 0, 0.5, 1, Inf),
                 labels = c("<-1", "(-1,-0.5]", "(-0.5,0]", "(0,0.5]", "(0.5,1]", ">1"))
#all <- subset(all, !(group %in% c("(-0.5,0]", "(0,0.5]")))

#plot the heatmap with p<0.05
all %>% 
  ggplot(aes(x=gene, y=celltype, fill = factor(group))) +
  coord_fixed(ratio = 1) +
  geom_tile() +
  scale_fill_manual(
    breaks = levels(all$group),
    values = c(
      "#225EA8",   # <-1
      "#8CC4F4",   # (-1,-0.5]
      "#C6DBEF",   # (-0.5,0]
      "#FF9999",   # (0,0.5]
      "#e82525",   # (0.5,1]
      "#9E1717"   # >1
    )) +
  ggtitle('')  +
  geom_tile( colour = 'white' ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y=NULL, x=NULL) +
  guides(fill=guide_legend(title="avg_log2FC")) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"),
        legend.spacing.y = unit(0.5, 'cm'))
ggsave("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/Riskgenes/UCHeatmap_uncorr.pdf",  width=6, height=6)

#plot the heatmap with p<0.00128
all <- filter(all, all$p_val < 0.00128)
all %>% 
  ggplot(aes(x=gene, y=celltype, fill = factor(group))) +
  coord_fixed(ratio = 1) +
  geom_tile() +
  scale_fill_manual(
    breaks = levels(all$group),
    values = c(
      "#225EA8",   # <-1
      "#8CC4F4",   # (-1,-0.5]
      "#C6DBEF",   # (-0.5,0]
      "#FF9999",   # (0,0.5]
      "#e82525",   # (0.5,1]
      "#9E1717"   # >1
    )) +
  ggtitle('')  +
  geom_tile( colour = 'white' ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y=NULL, x=NULL) +
  guides(fill=guide_legend(title="avg_log2FC")) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"),
        legend.spacing.y = unit(0.5, 'cm'))
ggsave("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/Riskgenes/UCHeatmap_corr.pdf",  width=6, height=6)




