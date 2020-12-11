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
data@meta.data$patient = sapply(strsplit(data$Final_HTO,"-"), `[`, 1)
data@meta.data$timepoint = sapply(strsplit(data$Final_HTO,"-"), `[`, 2)
data@meta.data$location = sapply(strsplit(data$Final_HTO,"-"), `[`, 3)
    #data@meta.data$state <- paste(data$disease, data$inflammation, sep='-')
data@meta.data$inflammation = sapply(strsplit(data$Final_HTO,"-"), `[`, 4)
data@meta.data$sample = sapply(strsplit(data$Final_HTO,"-"), `[`, 5)

table(data$inflammation, data$Final_HTO)

# recluster subset
# epi = average expression EPCAM > 5
AverageExpression(data, features = "EPCAM", assays = "RNA")
epi<-subset(data, idents = c("0","2","8","15","17","18","21","25","26"))
DefaultAssay(epi)<-"integrated"
epi<- RunPCA(epi)
epi <- RunUMAP(epi, dims = 1:30)
epi<-FindNeighbors(epi, dims = 1:30)
epi<-FindClusters(epi, resolution = 0.5)
DimPlot(epi, label=T)
DefaultAssay(epi)<-"RNA"
FeaturePlot(epi, features = c("BEST4", "DUOX2", "MUC2"))

#epi markers
markers_epi <- FindAllMarkers(epi, only.pos = TRUE, min.pct = 0.25)
write.csv(markers_epi, "Results/markergenes/markers_epi_log.csv")
selected_markers_epi <- markers_epi %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(selected_markers_epi, "Results/markergenes/selected_markers_epi_log.csv")

#after celltyping!
new.ident=c("Immature_enterocyte", "Cycling_TA", "Absorptive_enterocyte", "Ribo_TA","Stem", "MT-Hi-enterocyte", "BEST4_enterocyte", "Immature_goblet", "DUOX2_enterocyte", "PLCG2_TA", "Absorptive_TA", "Tuft", "REG_TA","Enteroendocrine", "Goblet", "MAST")
names(new.ident) <- levels(epi)
epi=RenameIdents(epi,new.ident)
DimPlot(epi, label=T)
epi[["celltypes"]] <- Idents(object = epi)
meta_epi<-epi@meta.data

# leuko = average expression PTPRC (CD45) > 3
AverageExpression(data, features = "PTPRC", assays = "RNA")
leuko<-subset(data, idents = c("5","6","11","12","13","14","22","23", "28"))
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

#after celltyping leuko
new.idents<-c("Follicular_B", "T_activated_FOS_low", "Treg", "CD4_T_memory", "MT_Hi_T", "CD8T", "APC", "CD8T", "T_activated_FOS_high", "Follicular_B", "APC", "MT_Hi_B", "Doublets", "GC", "MT_Hi_T", "Cycling_B", "APC", "Cycling_T", "Doublets")
names(new.idents) <- levels(leuko)
leuko=RenameIdents(leuko,new.idents)
DimPlot(leuko, label=T)
leuko[["celltypes"]] <- Idents(object = leuko)
meta_leuko<-leuko@meta.data

# stromal =  average expression THY1 > 1 (fibroblast), SOX10 > 2 (glia), MADCAM1 > 1 (endotheel)
AverageExpression(data, features = "THY1", assays = "RNA")
AverageExpression(data, features = "SOX10", assays = "RNA")
AverageExpression(data, features = "MADCAM1", assays = "RNA")
stromal<-subset(data, idents = c("7", "4", "20", "24", "27"))
DefaultAssay(stromal)<-"integrated"
stromal<- RunPCA(stromal)
stromal <- RunUMAP(stromal, dims = 1:30)
stromal <-FindNeighbors(stromal, dims = 1:30)
stromal <-FindClusters(stromal, resolution = 0.4)
DimPlot(stromal, label=T)
DefaultAssay(stromal)<-"RNA"
FeaturePlot(stromal, features = c("WNT2B", "WNT5B", "MADCAM1", "RSPO3", "SOX6"))

#stromal markers
markers_stromal <- FindAllMarkers(stromal, only.pos = TRUE, min.pct = 0.25)
write.csv(markers_stromal, "Results/markergenes/markers_stromal_log.csv")
selected_markers_stromal <- markers_stromal %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(selected_markers_stromal, "Results/markergenes/selected_markers_stromal_log.csv")

#after celltyping stromal
new.idents<-c("WNT5B+", "Inflammatory_fibroblast", "RSPO3+", "WNT2B+", "endothelial", "endothelial", "Mt-Hi_stromal", "pericytes", "glia", "myofibroblasts", "endothelial", "Doublets")
names(new.idents) <- levels(stromal)
stromal=RenameIdents(stromal,new.idents)
DimPlot(stromal, label=T)
stromal[["celltypes"]] <- Idents(object = stromal)
DimPlot(stromal)
meta_stromal<-stromal@meta.data

# plasma = remaining cells
plasma<-subset(data, idents = c("1", "3", "9", "10", "16", "19"))
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

#after celltyping plasma
new.idents<-c("IgA_plasma", "IgA_plasma", "IgG_plasma", "IgA_plasma", "IgG_plasma", "MT-Hi_plasma", "IgM_plasma", "IgA_plasma", "IgM_plasma")
names(new.idents) <- levels(plasma)
plasma=RenameIdents(plasma,new.idents)
DimPlot(plasma, label=T)
plasma[["celltypes"]] <- Idents(object = plasma)
DimPlot(plasma)
meta_plasma<-plasma@meta.data

# create meta_all
colnames(meta_stromal)
colnames(meta_leuko)
colnames(meta_plasma)
colnames(meta_epi)
meta_plasma<-meta_plasma[c(1,29)]
meta_leuko<-meta_leuko[c(1,29)]
meta_stromal<-meta_stromal[c(1,29)]
meta_epi <- meta_epi[c(1,29)]
x<-rbind(meta_stromal, meta_leuko)
x<-rbind(x, meta_plasma)
x<-rbind(x, meta_epi)
write.csv(x, "Data/meta_all.csv")

# add celltypes
CellsMeta<-data.frame(data@meta.data)
x$NAME<-rownames(x)
CellsMeta$NAME<-row.names(CellsMeta)
x<-x[,c(2,3)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
CellsMeta<-CellsMeta[,c(1,29)]
rownames(CellsMeta)<-CellsMeta$NAME
data<-AddMetaData(data, CellsMeta)
Idents(data)=data$celltypes
DimPlot(data, label = T, repel = T)

data = subset(data, subset = celltypes != "MT_Hi_T")
data = subset(data, subset = celltypes != "Doublets")
data = subset(data, subset = celltypes != "MT_Hi_B")
data = subset(data, subset = celltypes != "MT-Hi_plasma")
data = subset(data, subset = celltypes != "Mt-Hi_stromal")
data = subset(data, subset = celltypes != "MT-Hi-enterocyte")
data = subset(data, subset = sample != "3296")

Idents(data)=data$celltypes
DimPlot(data, label = T, repel = T)

DefaultAssay(data)<-"integrated"
data <- RunUMAP(data, dims = 1:30)
data <-FindNeighbors(data, dims = 1:30)
data<-FindClusters(data, resolution = 0.4)

Idents(data)=data$celltypes
DimPlot(data, label = T, repel = T, order = c("T_activated_FOS_low", "Absorptive_enterocyte", "IgA_plasma", "WNT5B+", "Aborptive_TA"))

saveRDS(data, "Data/PSC_processed_oct.rds")

# create figures
# DimPlot(stromal, label=T, pt.size = 0.3, label.size = 6, repel = T) + NoLegend()
# ggsave("Data/clusters_stromal_oct.pdf", height = 9, width = 16)
# DimPlot(data1, label = T, repel = T, order = c("T_activated_FOS_low", "Absorptive_enterocyte", "IgA_plasma", "WNT5B+", "Aborptive_TA", "REG_TA"), pt.size = 0.25, label.size = 5) + NoLegend()