#######################################
# Generate DUOX2 markers     #
#######################################

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)
library(tidyr)
library(MAST)
library(readxl)
library(openxlsx)
library(enrichR)

#load in datasets
epi <- readRDS("Nieuw/epi_azimuth_duox2.rds")
DefaultAssay(epi) = "RNA"
Idents(epi) <- "celltype.final"

#Create DUOX2 marker pathways
DUOX2markers <- FindMarkers(epi, group.by = "celltype.final", test.use = "MAST", ident.1 = "DUOX2 enterocytes")
DUOX2markers$Gene <- rownames(DUOX2markers)
DUOX2markers <- filter(DUOX2markers, DUOX2markers$p_val_adj < 0.05)
DUOX2pathways <- enrichr(DUOX2markers$Gene, databases = "GO_Biological_Process_2018")$GO_Biological_Process
DUOX2pathways <- filter(DUOX2pathways, DUOX2pathways$Adjusted.P.value < 0.05)

#Create enterocyte marker pathways
Enteromarkers <- FindMarkers(epi, group.by = "celltype.final", test.use = "MAST", ident.1 = "Enterocytes")
Enteromarkers$Gene <- rownames(Enteromarkers)
Enteromarkers <- filter(Enteromarkers, Enteromarkers$p_val_adj < 0.05)
Enteropathways <- enrichr(Enteromarkers$Gene, databases = "GO_Biological_Process_2018")$GO_Biological_Process
Enteropathways <- filter(Enteropathways, Enteropathways$Adjusted.P.value < 0.05)

#Create specific marker pathways list and write df's
DUOX2specific <- subset(DUOX2pathways,!(Term%in%Enteropathways$Term))
DUOX2nonspecific <- subset(DUOX2pathways,(Term%in%Enteropathways$Term))
write.csv(DUOX2specific, "Results/DE_2023/DUOX2specific.csv")
write.csv(DUOX2nonspecific, "Results/DE_2023/DUOX2nonspecific.csv")
write.csv(DUOX2pathways, "Results/DE_2023/DUOX2pathways.csv")
write.csv(Enteropathways, "Results/DE_2023/Enteropathways.csv")


