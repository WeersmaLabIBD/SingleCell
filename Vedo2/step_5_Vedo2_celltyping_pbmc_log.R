library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)
library(tidyr)
library(MAST)
library(readxl)
library(openxlsx)

# Load in database
data <- readRDS("Data/batch1_integrated/vedo2_batch1_PBMC_integrated_noribo.rds")
data = UpdateSeuratObject(object = data)
DefaultAssay(data) = "RNA"

#not normalized dataset
DimPlot(data)
FeaturePlot(data, features = c("CD3E", "CD14", "CD4", "CD8A", "MS4A1", "IL7R"))

#normalization
data <- NormalizeData(data, verbose = FALSE)
DimPlot(data, label = T) +NoLegend()

#add metadata
data@meta.data$patient = sapply(strsplit(data$Final_HTO,"-"), `[`, 1)
data@meta.data$timepoint = sapply(strsplit(data$Final_HTO,"-"), `[`, 2)

table(data$patient, data$Final_HTO)

#Setting cluster resolution (used on "integrated" assay)
DefaultAssay(data) <- "integrated"
data <- FindClusters(data, resolution = 0.7)
DimPlot(data, label = T)

#Marker genes pbmc
DefaultAssay(data) <- "RNA"
markers_pbmc <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25)
write.csv(markers_pbmc, "Results/markergenes/markers_pbmc_log.csv")
selected_markers_pbmc <- markers_pbmc %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(selected_markers_pbmc, "Results/markergenes/selected_markers_pbmc_log.csv")
