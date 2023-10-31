library(Seurat)
library(patchwork)
library(ggplot2)
library(scran)
library(scater)
library(dplyr)

#Read dataframe 08 lane1
dataset_10x<-Read10X("Data/Souporcell/08_lane1/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "PSC")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.9999)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/Souporcell/clusters_08_lane1.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="UC-NI-3249"
All$assignment[All$assignment == "1"]="PSC-I-3086"
All$assignment[All$assignment == "2"]="PSC-NI-3047"
All$assignment[All$assignment == "3"]="PSC-I-3316"
All$status[All$status == "doublet"]="Doublet"
All$status[All$status == "singlet"]="Singlet"
All$status[All$status == "unassigned"]="Negative"
for (i in 1:nrow(All)){ 
  if (All[i,2] == "Negative"){
    All[i,2]<-All[i,5] 
  } else {
    All[i,2]<-All[i,2]
  }
}
for (i in 1:nrow(All)){ 
  if (All[i,3] == "Negative"){
    All[i,3]<-All[i,4] 
  } else {
    All[i,3]<-All[i,3]
  }
}

Final_HTO <- All$HTO
names(Final_HTO) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO,
  col.name = 'Final_HTO')
Final_HTO_status <- All$HTO_status
names(Final_HTO_status) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO_status,
  col.name = 'Final_HTO_status')

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_08lane1 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_08lane1 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="lane08_1")
dataset@meta.data$lane<-"08_1"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "IGHA2")
#saveRDS(dataset, file="~/Documents/R/SingleCell/Data/200108_lane1_sct.rds")

#Read dataframe 08 lane2
dataset_10x<-Read10X("Data/Souporcell/08_lane2/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "PSC")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.965)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/Souporcell/clusters_08_lane2.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="UC-NI-3006"
All$assignment[All$assignment == "2"]="PSC-NI-3267"
All$assignment[All$assignment == "3"]="UC-I-3249"
All$assignment[All$assignment == "4"]="PSC-NI-3191"
All$status[All$status == "doublet"]="Doublet"
All$status[All$status == "singlet"]="Singlet"
All$status[All$status == "unassigned"]="Negative"
for (i in 1:nrow(All)){ 
  if (All[i,2] == "Negative" & All[i,5] != "1"){
    All[i,2]<-All[i,5] 
  } else {
    All[i,2]<-All[i,2]
  }
}
for (i in 1:nrow(All)){ 
  if (All[i,3] == "Negative" & All[i,5] != "1"){
    All[i,3]<-All[i,4] 
  } else {
    All[i,3]<-All[i,3]
  }
}

Final_HTO <- All$HTO
names(Final_HTO) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO,
  col.name = 'Final_HTO')
Final_HTO_status <- All$HTO_status
names(Final_HTO_status) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO_status,
  col.name = 'Final_HTO_status')
table(dataset@meta.data[["Final_HTO_status"]], dataset@meta.data[["hash.ID"]])
table(dataset@meta.data[["Final_HTO"]], dataset@meta.data[["hash.ID"]])

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_08lane2 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_08lane2 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="lane08_2")
dataset@meta.data$lane<-"08_2"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "IGHA2")
FeaturePlot(dataset, "BEST4")
FeaturePlot(dataset, "CD3E")
saveRDS(dataset, file="~/Documents/R/SingleCell/Data/200108_lane2_sct.rds")

#Read dataframe 09 lane1
dataset_10x<-Read10X("Data/Souporcell/09_lane1/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "PSC")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.91)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/Souporcell/clusters_09_lane1.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="UC-NI-3107"
All$assignment[All$assignment == "1"]="HC-NI-3083"
All$assignment[All$assignment == "2"]="UC-NI-3195"
All$assignment[All$assignment == "3"]="PSC-I-3191"
All$status[All$status == "doublet"]="Doublet"
All$status[All$status == "singlet"]="Singlet"
All$status[All$status == "unassigned"]="Negative"
for (i in 1:nrow(All)){ 
  if (All[i,2] == "Negative"){
    All[i,2]<-All[i,5] 
  } else {
    All[i,2]<-All[i,2]
  }
}
for (i in 1:nrow(All)){ 
  if (All[i,3] == "Negative"){
    All[i,3]<-All[i,4] 
  } else {
    All[i,3]<-All[i,3]
  }
}

Final_HTO <- All$HTO
names(Final_HTO) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO,
  col.name = 'Final_HTO')
Final_HTO_status <- All$HTO_status
names(Final_HTO_status) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO_status,
  col.name = 'Final_HTO_status')
table(dataset@meta.data[["Final_HTO_status"]], dataset@meta.data[["hash.ID"]])
table(dataset@meta.data[["Final_HTO"]], dataset@meta.data[["hash.ID"]])

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_09lane1 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_09lane1 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
table(dataset@meta.data[["Final_HTO"]])
dataset<-RenameCells(dataset, add.cell.id="lane09_1")
dataset@meta.data$lane<-"09_1"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "IGHG4")
FeaturePlot(dataset, "BEST4")
saveRDS(dataset, file="~/Documents/R/SingleCell/Data/200109_lane1_sct.rds")

#Read dataframe 09 lane1
dataset_10x<-Read10X("Data/Souporcell/09_lane2/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "PSC")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.985)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/Souporcell/clusters_09_lane2.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="PSC-NI-3086"
All$assignment[All$assignment == "1"]="UC-NI-3044"
All$assignment[All$assignment == "2"]="PSC-NI-3041"
All$assignment[All$assignment == "3"]="UC-I-3107"
All$status[All$status == "doublet"]="Doublet"
All$status[All$status == "singlet"]="Singlet"
All$status[All$status == "unassigned"]="Negative"
for (i in 1:nrow(All)){ 
  if (All[i,2] == "Negative"){
    All[i,2]<-All[i,5] 
  } else {
    All[i,2]<-All[i,2]
  }
}
for (i in 1:nrow(All)){ 
  if (All[i,3] == "Negative"){
    All[i,3]<-All[i,4] 
  } else {
    All[i,3]<-All[i,3]
  }
}

Final_HTO <- All$HTO
names(Final_HTO) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO,
  col.name = 'Final_HTO')
Final_HTO_status <- All$HTO_status
names(Final_HTO_status) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO_status,
  col.name = 'Final_HTO_status')
table(dataset@meta.data[["Final_HTO_status"]], dataset@meta.data[["hash.ID"]])
table(dataset@meta.data[["Final_HTO"]], dataset@meta.data[["hash.ID"]])

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_09lane2 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_09lane2 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
table(dataset@meta.data[["Final_HTO"]])
dataset<-RenameCells(dataset, add.cell.id="lane09_2")
dataset@meta.data$lane<-"09_2"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "IGHA2")
FeaturePlot(dataset, "BEST4")
saveRDS(dataset, file="~/Documents/R/SingleCell/Data/200109_lane2_sct.rds")

#Read dataframe 13 lane1
dataset_10x<-Read10X("Data/Souporcell/13_lane1/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "PSC")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.98)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/Souporcell/clusters_13_lane1.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="PSC-NI-3015"
All$assignment[All$assignment == "1"]="HC-NI-3129"
All$assignment[All$assignment == "2"]="UC-I-3023"
All$assignment[All$assignment == "3"]="UC-NI-3085"
All$status[All$status == "doublet"]="Doublet"
All$status[All$status == "singlet"]="Singlet"
All$status[All$status == "unassigned"]="Negative"
for (i in 1:nrow(All)){ 
  if (All[i,2] == "Negative"){
    All[i,2]<-All[i,5] 
  } else {
    All[i,2]<-All[i,2]
  }
}
for (i in 1:nrow(All)){ 
  if (All[i,3] == "Negative"){
    All[i,3]<-All[i,4] 
  } else {
    All[i,3]<-All[i,3]
  }
}

Final_HTO <- All$HTO
names(Final_HTO) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO,
  col.name = 'Final_HTO')
Final_HTO_status <- All$HTO_status
names(Final_HTO_status) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO_status,
  col.name = 'Final_HTO_status')
table(dataset@meta.data[["Final_HTO_status"]], dataset@meta.data[["hash.ID"]])
table(dataset@meta.data[["Final_HTO"]], dataset@meta.data[["hash.ID"]])

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_13lane1 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_13lane1 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="lane13_1")
dataset@meta.data$lane<-"13_1"
table(dataset@meta.data[["Final_HTO"]])
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "IGHG4")
FeaturePlot(dataset, "BEST4")
saveRDS(dataset, file="~/Documents/R/SingleCell/Data/200113_lane1_sct.rds")

#Read dataframe 13 lane2
dataset_10x<-Read10X("Data/Souporcell/13_lane2/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "PSC")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.96)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/Souporcell/clusters_13_lane2.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="UC-I-3006"
All$assignment[All$assignment == "2"]="PSC-NI-3005"
All$status[All$status == "doublet"]="Doublet"
All$status[All$status == "singlet"]="Singlet"
All$status[All$status == "unassigned"]="Negative"
for (i in 1:nrow(All)){ 
  if (All[i,2] == "Negative" & (All[i,5] == "UC-I-3006" | All[i,5] == "PSC-NI-3005")){
    All[i,2]<-All[i,5] 
  } else {
    All[i,2]<-All[i,2]
  }
}
for (i in 1:nrow(All)){ 
  if (All[i,3] == "Negative" & (All[i,5] == "UC-I-3006" | All[i,5] == "PSC-NI-3005")){
    All[i,3]<-All[i,4] 
  } else {
    All[i,3]<-All[i,3]
  }
}

Final_HTO <- All$HTO
names(Final_HTO) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO,
  col.name = 'Final_HTO')
Final_HTO_status <- All$HTO_status
names(Final_HTO_status) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO_status,
  col.name = 'Final_HTO_status')
table(dataset@meta.data[["Final_HTO_status"]], dataset@meta.data[["hash.ID"]])
table(dataset@meta.data[["Final_HTO"]], dataset@meta.data[["hash.ID"]])

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_13lane2 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_13lane2 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="lane13_2")
dataset@meta.data$lane<-"13_2"
table(dataset@meta.data[["Final_HTO"]])
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "IGHA2")
FeaturePlot(dataset, "BEST4")
saveRDS(dataset, file="~/Documents/R/SingleCell/Data/200113_lane2_sct.rds")

#Read dataframe 14 lane1
dataset_10x<-Read10X("Data/Souporcell/14_lane1/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "PSC")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.995)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/Souporcell/clusters_14_lane1.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="HC-NI-3002"
All$assignment[All$assignment == "1"]="PSC-NI-3176"
All$assignment[All$assignment == "2"]="UC-NI-3287"
All$assignment[All$assignment == "3"]="UC-I-3125"
All$status[All$status == "doublet"]="Doublet"
All$status[All$status == "singlet"]="Singlet"
All$status[All$status == "unassigned"]="Negative"
for (i in 1:nrow(All)){ 
  if (All[i,2] == "Negative"){
    All[i,2]<-All[i,5] 
  } else {
    All[i,2]<-All[i,2]
  }
}
for (i in 1:nrow(All)){ 
  if (All[i,3] == "Negative"){
    All[i,3]<-All[i,4] 
  } else {
    All[i,3]<-All[i,3]
  }
}

Final_HTO <- All$HTO
names(Final_HTO) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO,
  col.name = 'Final_HTO')
Final_HTO_status <- All$HTO_status
names(Final_HTO_status) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO_status,
  col.name = 'Final_HTO_status')
table(dataset@meta.data[["Final_HTO_status"]], dataset@meta.data[["hash.ID"]])
table(dataset@meta.data[["Final_HTO"]], dataset@meta.data[["hash.ID"]])

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_14lane1 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_14lane1 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="lane14_1")
dataset@meta.data$lane<-"14_1"
table(dataset@meta.data[["Final_HTO"]])
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "IGHA2")
FeaturePlot(dataset, "BEST4")
saveRDS(dataset, file="~/Documents/R/SingleCell/Data/200114_lane1_sct.rds")

#Read dataframe 14 lane2
dataset_10x<-Read10X("Data/Souporcell/14_lane2/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "PSC")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.995)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/Souporcell/clusters_14_lane2.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="PSC-NI-3162"
All$assignment[All$assignment == "1"]="HC-NI-3049"
All$assignment[All$assignment == "2"]="PSC-I-3019"
All$assignment[All$assignment == "3"]="PSC-NI-3263"
All$status[All$status == "doublet"]="Doublet"
All$status[All$status == "singlet"]="Singlet"
All$status[All$status == "unassigned"]="Negative"
for (i in 1:nrow(All)){ 
  if (All[i,2] == "Negative"){
    All[i,2]<-All[i,5] 
  } else {
    All[i,2]<-All[i,2]
  }
}
for (i in 1:nrow(All)){ 
  if (All[i,3] == "Negative"){
    All[i,3]<-All[i,4] 
  } else {
    All[i,3]<-All[i,3]
  }
}

Final_HTO <- All$HTO
names(Final_HTO) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO,
  col.name = 'Final_HTO')
Final_HTO_status <- All$HTO_status
names(Final_HTO_status) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO_status,
  col.name = 'Final_HTO_status')
table(dataset@meta.data[["Final_HTO_status"]], dataset@meta.data[["hash.ID"]])
table(dataset@meta.data[["Final_HTO"]], dataset@meta.data[["hash.ID"]])

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_14lane2 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_14lane2 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="lane14_1")
dataset@meta.data$lane<-"14_2"
table(dataset@meta.data[["Final_HTO"]])
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "IGHG4")
FeaturePlot(dataset, "BEST4")
saveRDS(dataset, file="~/Documents/R/SingleCell/Data/200114_lane2_sct.rds")

#Read dataframe 15 lane1
dataset_10x<-Read10X("Data/Souporcell/15_lane1/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "PSC")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.99)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/Souporcell/clusters_15_lane1.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="UC-NI-3076"
All$assignment[All$assignment == "1"]="UC-I-3059"
All$assignment[All$assignment == "2"]="PSC-NI-3325"
All$assignment[All$assignment == "3"]="PSC-NI-3069"
All$status[All$status == "doublet"]="Doublet"
All$status[All$status == "singlet"]="Singlet"
All$status[All$status == "unassigned"]="Negative"
for (i in 1:nrow(All)){ 
  if (All[i,2] == "Negative"){
    All[i,2]<-All[i,5] 
  } else {
    All[i,2]<-All[i,2]
  }
}
for (i in 1:nrow(All)){ 
  if (All[i,3] == "Negative"){
    All[i,3]<-All[i,4] 
  } else {
    All[i,3]<-All[i,3]
  }
}

Final_HTO <- All$HTO
names(Final_HTO) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO,
  col.name = 'Final_HTO')
Final_HTO_status <- All$HTO_status
names(Final_HTO_status) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO_status,
  col.name = 'Final_HTO_status')
table(dataset@meta.data[["Final_HTO_status"]], dataset@meta.data[["hash.ID"]])
table(dataset@meta.data[["Final_HTO"]], dataset@meta.data[["hash.ID"]])

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_15lane1 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_15lane1 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="lane15_1")
dataset@meta.data$lane<-"15_1"
table(dataset@meta.data[["Final_HTO"]])
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "IGHA2")
FeaturePlot(dataset, "BEST4")
saveRDS(dataset, file="~/Documents/R/SingleCell/Data/200115_lane1_sct.rds")

#Read dataframe 15 lane2
dataset_10x<-Read10X("Data/Souporcell/15_lane2/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "PSC")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.96)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/Souporcell/clusters_15_lane2.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="HC-NI-3030"
All$assignment[All$assignment == "1"]="UC-NI-3125"
All$assignment[All$assignment == "2"]="PSC-I-3325"
All$assignment[All$assignment == "3"]="PSC-I-3069"
All$status[All$status == "doublet"]="Doublet"
All$status[All$status == "singlet"]="Singlet"
All$status[All$status == "unassigned"]="Negative"
for (i in 1:nrow(All)){ 
  if (All[i,2] == "Negative"){
    All[i,2]<-All[i,5] 
  } else {
    All[i,2]<-All[i,2]
  }
}
for (i in 1:nrow(All)){ 
  if (All[i,3] == "Negative"){
    All[i,3]<-All[i,4] 
  } else {
    All[i,3]<-All[i,3]
  }
}

Final_HTO <- All$HTO
names(Final_HTO) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO,
  col.name = 'Final_HTO')
Final_HTO_status <- All$HTO_status
names(Final_HTO_status) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO_status,
  col.name = 'Final_HTO_status')
table(dataset@meta.data[["Final_HTO_status"]], dataset@meta.data[["hash.ID"]])
table(dataset@meta.data[["Final_HTO"]], dataset@meta.data[["hash.ID"]])

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_15lane2 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_15lane2 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="lane15_2")
dataset@meta.data$lane<-"15_2"
table(dataset@meta.data[["Final_HTO"]])
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "IGHG4")
FeaturePlot(dataset, "BEST4")
saveRDS(dataset, file="~/Documents/R/SingleCell/Data/200115_lane2_sct.rds")

#Read dataframe 16 lane1
dataset_10x<-Read10X("Data/Souporcell/16_lane1/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "PSC")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Custom`)
dataset[['ADT']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- NormalizeData(dataset, assay = "ADT", normalization.method = "CLR")
dataset <- ScaleData(dataset, assay = "ADT")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.975)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/Souporcell/clusters_16_lane1.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="UC-I-3076"
All$assignment[All$assignment == "1"]="UC-I-3085"
All$assignment[All$assignment == "2"]="PSC-I-3317"
All$assignment[All$assignment == "3"]="PSC-NI-3147"
All$status[All$status == "doublet"]="Doublet"
All$status[All$status == "singlet"]="Singlet"
All$status[All$status == "unassigned"]="Negative"
for (i in 1:nrow(All)){ 
  if (All[i,2] == "Negative"){
    All[i,2]<-All[i,5] 
  } else {
    All[i,2]<-All[i,2]
  }
}
for (i in 1:nrow(All)){ 
  if (All[i,3] == "Negative"){
    All[i,3]<-All[i,4] 
  } else {
    All[i,3]<-All[i,3]
  }
}

Final_HTO <- All$HTO
names(Final_HTO) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO,
  col.name = 'Final_HTO')
Final_HTO_status <- All$HTO_status
names(Final_HTO_status) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO_status,
  col.name = 'Final_HTO_status')
table(dataset@meta.data[["Final_HTO_status"]], dataset@meta.data[["hash.ID"]])
table(dataset@meta.data[["Final_HTO"]], dataset@meta.data[["hash.ID"]])

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_16lane1 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_16lane1 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="lane16_1")
dataset@meta.data$lane<-"16_1"
table(dataset@meta.data[["Final_HTO"]])
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "IGHG4")
FeaturePlot(dataset, "BEST4")
FeaturePlot(dataset, "CD103-TotalA")
saveRDS(dataset, file="~/Documents/R/SingleCell/Data/200116_lane1_sct.rds")

#Read dataframe 16 lane2
dataset_10x<-Read10X("Data/Souporcell/16_lane2/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "PSC")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Custom`)
dataset[['ADT']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- NormalizeData(dataset, assay = "ADT", normalization.method = "CLR")
dataset <- ScaleData(dataset, assay = "ADT")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.98)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/Souporcell/clusters_16_lane2.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="UC-NI-3059"
All$assignment[All$assignment == "1"]="HC-NI-3296"
All$assignment[All$assignment == "2"]="PSC-NI-3312"
All$assignment[All$assignment == "3"]="PSC-I-3147"
All$status[All$status == "doublet"]="Doublet"
All$status[All$status == "singlet"]="Singlet"
All$status[All$status == "unassigned"]="Negative"
for (i in 1:nrow(All)){ 
  if (All[i,2] == "Negative"){
    All[i,2]<-All[i,5] 
  } else {
    All[i,2]<-All[i,2]
  }
}
for (i in 1:nrow(All)){ 
  if (All[i,3] == "Negative"){
    All[i,3]<-All[i,4] 
  } else {
    All[i,3]<-All[i,3]
  }
}

Final_HTO <- All$HTO
names(Final_HTO) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO,
  col.name = 'Final_HTO')
Final_HTO_status <- All$HTO_status
names(Final_HTO_status) <- rownames(All)
dataset <- AddMetaData(
  object = dataset,
  metadata = Final_HTO_status,
  col.name = 'Final_HTO_status')
table(dataset@meta.data[["Final_HTO_status"]], dataset@meta.data[["hash.ID"]])
table(dataset@meta.data[["Final_HTO"]], dataset@meta.data[["hash.ID"]])

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_16lane2 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_16lane2 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="lane16_2")
dataset@meta.data$lane<-"16_2"
table(dataset@meta.data[["Final_HTO"]])
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "IGHA2")
FeaturePlot(dataset, "BEST4")
saveRDS(dataset, file="~/Documents/R/SingleCell/Data/200116_lane2_sct.rds")