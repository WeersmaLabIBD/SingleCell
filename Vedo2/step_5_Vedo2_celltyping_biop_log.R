# general stuff for sc data

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
data <- readRDS("Data/batch1_integrated/vedo2_batch1_biop_integrated_noribo_nonorm_withVEDO01.rds")
data = UpdateSeuratObject(object = data)
DefaultAssay(data) = "RNA"

#not normalized dataset
DimPlot(data)

#normalization
data <- NormalizeData(data, verbose = FALSE)
DimPlot(data, label = T) +NoLegend()

#add metadata
data$Final_HTO <- sub("V001_T0_ascendens_I", "V001-T0-ascendens-I-3191", data$Final_HTO)
data$Final_HTO <- sub("V001_T0_transversum_NI", "V001-T0-transversum-NI-3191", data$Final_HTO)
data$Final_HTO <- sub("V002-T0-colonasc-NI-3160", "V002-T0-ascendens-NI-3160", data$Final_HTO)
data$Final_HTO <- sub("V005-T0-transvers-NI-3228", "V005-T0-transversum-NI-3228", data$Final_HTO)
data$Final_HTO <- sub("V005-T0-sigmoid-NI-3228", "V005-T0-sigmoid-I-3228", data$Final_HTO)
data@meta.data$patient = sapply(strsplit(data$Final_HTO,"-"), `[`, 1)
data@meta.data$timepoint = sapply(strsplit(data$Final_HTO,"-"), `[`, 2)
data@meta.data$location = sapply(strsplit(data$Final_HTO,"-"), `[`, 3)
    #data@meta.data$state <- paste(data$disease, data$inflammation, sep='-')
data@meta.data$inflammation = sapply(strsplit(data$Final_HTO,"-"), `[`, 4)
data@meta.data$sample = sapply(strsplit(data$Final_HTO,"-"), `[`, 5)

table(data$inflammation, data$Final_HTO)

#Setting cluster resolution (used on "integrated" assay)
DefaultAssay(data) <- "integrated"
data <- FindClusters(data, resolution = 0.7)
DimPlot(data, label = T)
data <- RunTSNE(data, reduction = "pca", dims = 1:5, seed.use = 1, tsne.method = "Rtsne", dim.embed = 2, reduction.name = "tsne", reduction.key = "tSNE_")
DimPlot(data, reduction = "tsne", label = T)

# recluster subset
# epi = average expression EPCAM > 4
AverageExpression(data, features = "EPCAM", assays = "RNA")
epi<-subset(data, idents = c("3", "5", "7", "12", "14", "16", "19", "27", "28"))
DefaultAssay(epi)<-"integrated"
epi<- RunPCA(epi)
epi <- RunUMAP(epi, dims = 1:30)
epi<-FindNeighbors(epi, dims = 1:30)
epi<-FindClusters(epi, resolution = 0.6)
DimPlot(epi, label=T)
DefaultAssay(epi)<-"RNA"
FeaturePlot(epi, features = c("BEST4", "DUOX2", "MUC2"))

#epi markers
markers_epi <- FindAllMarkers(epi, only.pos = TRUE, min.pct = 0.25)
write.csv(markers_epi, "Results/markergenes/markers_epi_log.csv")
#to export only top (10/20/etc) gene markers
selected_markers_epi <- markers_epi %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(selected_markers_epi, "Results/markergenes/selected_markers_epi_log.csv")

# leuko = average expression PTPRC (CD45) > 3
AverageExpression(data, features = "PTPRC", assays = "RNA")
leuko<-subset(data, idents = c("6", "8", "11", "13","15","17","22","23","24","31"))
DefaultAssay(leuko)<-"integrated"
leuko<- RunPCA(leuko)
leuko <- RunUMAP(leuko, dims = 1:30)
leuko<-FindNeighbors(leuko, dims = 1:30)
leuko<-FindClusters(leuko, resolution = 0.8)
DimPlot(leuko, label=T)
DefaultAssay(leuko)<-"RNA"
FeaturePlot(leuko, features = c("CTSG", "IL17A"))

#leuko markers
markers_leuko <- FindAllMarkers(leuko, only.pos = TRUE, min.pct = 0.25)
write.csv(markers_leuko, "Results/markergenes/markers_leuko_log.csv")
selected_markers_leuko <- markers_leuko %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(selected_markers_leuko, "Results/markergenes/selected_markers_leuko_log.csv")


# stromal =  average expression THY1 > 1 (fibroblast), SOX10 > 2 (glia), MADCAM1 > 1 (endotheel)
AverageExpression(data, features = "THY1", assays = "RNA")
AverageExpression(data, features = "SOX10", assays = "RNA")
AverageExpression(data, features = "MADCAM1", assays = "RNA")
stromal<-subset(data, idents = c("1", "9", "20", "25", "26", "29"))
DefaultAssay(stromal)<-"integrated"
stromal<- RunPCA(stromal)
stromal <- RunUMAP(stromal, dims = 1:30)
stromal <-FindNeighbors(stromal, dims = 1:30)
stromal <-FindClusters(stromal, resolution = 0.4)
DimPlot(stromal, label=T)
DefaultAssay(stromal)<-"RNA"
FeaturePlot(stromal, features = c("WNT2B", "WNT5B", "MADCAM1", "RSPO3", "SOX6"))

#subcluster again (only clusters 1,5 and 6)to find IAFs
fibro <- subset(stromal, idents = c("1", "5", "6"))
DefaultAssay(fibro)<-"integrated"
fibro <- RunPCA(fibro)
fibro <- RunUMAP(fibro, dims = 1:30)
fibro <-FindNeighbors(fibro, dims = 1:30)
fibro <-FindClusters(fibro, resolution = 0.4)
DimPlot(fibro, label=T)
DefaultAssay(fibro)<-"RNA"
FeaturePlot(fibro, features = c("WNT2B", "WNT5B", "WNT5A", "SOX6", "RSPO3", "TNFRSF11B", "IL11", "IL13RA2", "CHI3L1", "IL6","CXCL3", "IL33", "MMP3", "CCL19"), cols = c("grey", "red"), min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(fibro, features = c("IL11", "IL13RA2", "CHI3L1"), cols = c("grey", "red"), split.by = "inflammation",min.cutoff = "q10", max.cutoff = "q90")
DotPlot(fibro, features = c("WNT2B", "WNT5B", "WNT5A", "SOX6", "RSPO3", "TNFRSF11B", "IL11", "IL13RA2", "CHI3L1", "IL6","CXCL3", "IL33", "MMP3", "CCL19"), cols = c("grey", "red"), dot.scale = 6) + RotatedAxis() + coord_flip()
DotPlot(fibro, features = c("WNT2B", "WNT5B", "WNT5A", "SOX6", "RSPO3", "TNFRSF11B", "IL11", "IL13RA2", "CHI3L1", "IL6","CXCL3", "IL33", "MMP3", "CCL19"), cols = c("green", "blue"),split.by = "timepoint", dot.scale = 6) + RotatedAxis() + coord_flip()

#stromal markers
markers_stromal <- FindAllMarkers(stromal, only.pos = TRUE, min.pct = 0.25)
write.csv(markers_stromal, "Results/markergenes/markers_stromal_log.csv")
selected_markers_stromal <- markers_stromal %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(selected_markers_stromal, "Results/markergenes/selected_markers_stromal_log.csv")


# plasma = remaining cells
plasma<-subset(data, idents = c("0", "2", "4", "10", "18", "21", "30"))
DefaultAssay(plasma)<-"integrated"
plasma<- RunPCA(plasma)
plasma <- RunUMAP(plasma, dims = 1:30)
plasma <-FindNeighbors(plasma, dims = 1:30)
plasma <-FindClusters(plasma, resolution = 0.1)
DimPlot(plasma, label=T)
DefaultAssay(plasma)<-"RNA"
FeaturePlot(plasma, features = c("IGHG1", "IGHA1", "IGHM"))

#plasma markers
markers_plasma <- FindAllMarkers(plasma, only.pos = TRUE, min.pct = 0.25)
write.csv(markers_plasma, "Results/markergenes//markers_plasma_log.csv")
selected_markers_plasma <- markers_plasma %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(selected_markers_plasma, "Results/markergenes/selected_markers_plasma_log.csv")

