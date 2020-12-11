library(Seurat)
library(patchwork)
library(ggplot2)
library(scran)
library(scater)
library(dplyr)

#Read dataframe 11 lane1
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200611_lane1/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
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
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200611_lane1.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="V004-T4-rectum-NI-3300"
All$assignment[All$assignment == "1"]="V050-T0-sigmoid-I-3131"
All$assignment[All$assignment == "2"]="V008-T4-rectum-I-3356"
All$assignment[All$assignment == "3"]="V007-T0-sigmoid-NI-3298"
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

table(All$HTO)
table(All$HTO_status)
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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_11lane1 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_11lane1 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="11_lane1")
dataset@meta.data$lane<-"11_1"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "EPCAM")
FeaturePlot(dataset, "BEST4")
saveRDS(dataset, file="Data/preprocessed_lanes/200611_lane1_sct.rds")


#Read dataframe 200611_lane2
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200611_lane2/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.894)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200611_lane2.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="V005-T0-ileum-NI-3228"
All$assignment[All$assignment == "1"]="V006-T4-ileum-I-3351"
All$assignment[All$assignment == "2"]="V002-T4-sigmoid-I-3222"
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

table(All$HTO)
table(All$HTO_status)
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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_11lane2 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_11lane2 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="11_lane2")
dataset@meta.data$lane<-"11_2"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "EPCAM")
FeaturePlot(dataset, "AOC3")
saveRDS(dataset, file="Data/preprocessed_lanes/200611_lane2_sct.rds")


#Read dataframe 200611_lane3
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200611_lane3/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.898)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200611_lane3.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="V050-T0-sigmoid-I-3131"
All$assignment[All$assignment == "1"]="V005-T0-ileum-NI-3228"
All$assignment[All$assignment == "2"]="V004-T4-rectum-NI-3300"
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

table(All$HTO)
table(All$HTO_status)
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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_11lane3 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_11lane3 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="11_lane3")
dataset@meta.data$lane<-"11_3"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "EPCAM")
FeaturePlot(dataset, "IL11")
saveRDS(dataset, file="Data/preprocessed_lanes/200611_lane3_sct.rds")


#Read dataframe 200611_lane4
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200611_lane4/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
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
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200611_lane4.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="V007-T0-sigmoid-NI-3298"
All$assignment[All$assignment == "1"]="V008-T4-rectum-I-3356"
All$assignment[All$assignment == "2"]="V006-T4-ileum-I-3351"
All$assignment[All$assignment == "3"]="V002-T4-sigmoid-I-3222"
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

table(All$HTO)
table(All$HTO_status)

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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_11lane4 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_11lane4 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="11_lane4")
dataset@meta.data$lane<-"11_4"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "EPCAM")
FeaturePlot(dataset, "IL13RA2")
saveRDS(dataset, file="Data/preprocessed_lanes/200611_lane4_sct.rds")


#Read dataframe 200612_lane1
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200612_lane1/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.947)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200612_lane1.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="V007-T4-sigmoid-NI-3339"
All$assignment[All$assignment == "1"]="V008-T0-ileum-NI-3245"
All$assignment[All$assignment == "2"]="V001-T4-sigmoid-NI-3271"
All$assignment[All$assignment == "3"]="V004-T0-sigmoid-NI-3188"
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

table(All$HTO)
table(All$HTO_status)

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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_12lane1 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_12lane1 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="12_lane1")
dataset@meta.data$lane<-"12_1"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "EPCAM")
saveRDS(dataset, file="Data/preprocessed_lanes/200612_lane1_sct.rds")


#Read dataframe 200612_lane2
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200612_lane2/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.852)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200612_lane2.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="V002-T0-rectum-I-3160"
All$assignment[All$assignment == "1"]="V006-T0-ileum-NI-3214"
All$assignment[All$assignment == "2"]="V005-T4-sigmoid-NI-3392"
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

table(All$HTO)
table(All$HTO_status)

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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_12lane2 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_12lane2 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="12_lane2")
dataset@meta.data$lane<-"12_2"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "EPCAM")
FeaturePlot(dataset, "IGHA2")
saveRDS(dataset, file="Data/preprocessed_lanes/200612_lane2_sct.rds")


#Read dataframe 200612_lane3
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200612_lane3/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.95)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200612_lane3.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="V007-T4-sigmoid-NI-3339"
All$assignment[All$assignment == "1"]="V004-T0-sigmoid-NI-3188"
All$assignment[All$assignment == "2"]="V005-T4-sigmoid-NI-3392"
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

table(All$HTO)
table(All$HTO_status)

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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_12lane3 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_12lane3 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="12_lane3")
dataset@meta.data$lane<-"12_3"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "EPCAM")
saveRDS(dataset, file="Data/preprocessed_lanes/200612_lane3_sct.rds")


#Read dataframe 200612_lane4
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200612_lane4/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.915)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200612_lane4.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="V006-T0-ileum-NI-3214"
All$assignment[All$assignment == "1"]="V008-T0-ileum-NI-3245"
All$assignment[All$assignment == "2"]="V001-T4-sigmoid-NI-3271"
All$assignment[All$assignment == "3"]="V002-T0-rectum-I-3160"
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

table(All$HTO)
table(All$HTO_status)

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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_12lane4 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_12lane4 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="12_lane4")
dataset@meta.data$lane<-"12_4"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "EPCAM")
saveRDS(dataset, file="Data/preprocessed_lanes/200612_lane4_sct.rds")


#Read dataframe 200625_lane1
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200625_lane1/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.951)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200625_lane1.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="V004-T4-sigmoid-NI-3300"
All$assignment[All$assignment == "1"]="V008-T4-sigmoid-NI-3356"
All$assignment[All$assignment == "2"]="V007-T0-ileum-NI-3298"
All$assignment[All$assignment == "3"]="V050-T4-sigmoid-I-3215"
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

table(All$HTO)
table(All$HTO_status)

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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_25lane1 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_25lane1 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="25_lane1")
dataset@meta.data$lane<-"25_1"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "EPCAM")
saveRDS(dataset, file="Data/preprocessed_lanes/200625_lane1_sct.rds")


#Read dataframe 200625_lane2
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200625_lane2/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.865)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200625_lane2.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="V006-T4-sigmoid-I-3351"
All$assignment[All$assignment == "1"]="V002-T0-colonasc-NI-3160"
All$assignment[All$assignment == "2"]="V005-T0-transvers-NI-3228"
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

table(All$HTO)
table(All$HTO_status)

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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_25lane2 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_25lane2 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="25_lane2")
dataset@meta.data$lane<-"25_2"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "EPCAM")
saveRDS(dataset, file="Data/preprocessed_lanes/200625_lane2_sct.rds")


#Read dataframe 200625_lane3
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200625_lane3/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.88)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200625_lane3.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="V007-T0-ileum-NI-3298"
All$assignment[All$assignment == "1"]="V005-T0-transvers-NI-3228"
All$assignment[All$assignment == "2"]="V004-T4-sigmoid-NI-3300"
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

table(All$HTO)
table(All$HTO_status)

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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_25lane3 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_25lane3 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="25_lane3")
dataset@meta.data$lane<-"25_3"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "EPCAM")
saveRDS(dataset, file="Data/preprocessed_lanes/200625_lane3_sct.rds")


#Read dataframe 200625_lane4
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200625_lane4/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.955)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200625_lane4.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="V002-T0-colonasc-NI-3160"
All$assignment[All$assignment == "1"]="V008-T4-sigmoid-NI-3356"
All$assignment[All$assignment == "2"]="V050-T4-sigmoid-I-3215"
All$assignment[All$assignment == "3"]="V006-T4-sigmoid-I-3351"
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

table(All$HTO)
table(All$HTO_status)

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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_25lane4 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_25lane4 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="25_lane4")
dataset@meta.data$lane<-"25_4"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "EPCAM")
saveRDS(dataset, file="Data/preprocessed_lanes/200625_lane4_sct.rds")


#Read dataframe 200626_lane1
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200626_lane1/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.86)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200626_lane1.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="V004-T0-rectum-I-3188"
All$assignment[All$assignment == "1"]="V005-T4-rectum-NI-3392"
All$assignment[All$assignment == "2"]="V006-T0-sigmoid-NI-3214"
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

table(All$HTO)
table(All$HTO_status)

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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_26lane1 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_26lane1 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="26_lane1")
dataset@meta.data$lane<-"26_1"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "EPCAM")
saveRDS(dataset, file="Data/preprocessed_lanes/200626_lane1_sct.rds")


#Read dataframe 200626_lane2
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200626_lane2/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.86)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200626_lane2.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="V005-T0-sigmoid-NI-3228"
All$assignment[All$assignment == "1"]="V007-T4-ileum-NI-3339"
All$assignment[All$assignment == "2"]="V008-T0-sigmoid-I-3245"
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

table(All$HTO)
table(All$HTO_status)

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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_26lane2 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_26lane2 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="26_lane2")
dataset@meta.data$lane<-"26_2"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "BEST4")
FeaturePlot(dataset, "EPCAM")
saveRDS(dataset, file="Data/preprocessed_lanes/200626_lane2_sct.rds")


#Read dataframe 200626_lane3
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200626_lane3/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.869)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200626_lane3.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="V005-T4-rectum-NI-3392"
All$assignment[All$assignment == "1"]="V008-T0-sigmoid-I-3245"
All$assignment[All$assignment == "2"]="V006-T0-sigmoid-NI-3214"
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

table(All$HTO)
table(All$HTO_status)

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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_26lane3 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_26lane3 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="26_lane3")
dataset@meta.data$lane<-"26_3"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "BEST4")
FeaturePlot(dataset, "EPCAM")
saveRDS(dataset, file="Data/preprocessed_lanes/200626_lane3_sct.rds")


#Read dataframe 200626_lane4
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200626_lane4/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.869)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200626_lane4.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="V007-T4-ileum-NI-3339"
All$assignment[All$assignment == "1"]="V005-T0-sigmoid-NI-3228"
All$assignment[All$assignment == "2"]="V004-T0-rectum-I-3188"
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

table(All$HTO)
table(All$HTO_status)

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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 60)
HTO_number_26lane4 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_26lane4 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="26_lane4")
dataset@meta.data$lane<-"26_4"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "BEST4")
FeaturePlot(dataset, "EPCAM")
saveRDS(dataset, file="Data/preprocessed_lanes/200626_lane4_sct.rds")


#Read dataframe 200618_lane1 (PBMC)
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200618_lane1/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
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
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200618_lane1.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
options(max.print=9999)
table(All$HTO, All$assignment)
All$assignment[All$assignment == "0"]="V050-T1-180619"
All$assignment[All$assignment == "1"]="V008-T0-151019"
All$assignment[All$assignment == "2"]="V006-T4-100220"
All$assignment[All$assignment == "3"]="V002-T5-291019"
All$assignment[All$assignment == "4"]="V007-T2-241219"
All$assignment[All$assignment == "5"]="V005-T4-190220"
All$assignment[All$assignment == "6"]="V004-T3-251019"
All$assignment[All$assignment == "7"]="V001-T0-140819"
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

table(All$HTO)
table(All$HTO_status)

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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 25)
HTO_number_18lane1 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_18lane1 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="18_lane1")
dataset@meta.data$lane<-"18_1"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "CD3E")
saveRDS(dataset, file="Data/preprocessed_lanes/200618_lane1_sct.rds")


#Read dataframe 200618_lane2 (PBMC)
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200618_lane2/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.9905)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200618_lane2.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
options(max.print=9999)
table(All$HTO, All$assignment)
All$assignment[All$assignment == "1"]="V002-T2-050819"
All$assignment[All$assignment == "2"]="V005-T3-241219"
All$assignment[All$assignment == "3"]="V050-T3-300719"
All$status[All$status == "doublet"]="Doublet"
All$status[All$status == "singlet"]="Singlet"
All$status[All$status == "unassigned"]="Negative"

for (i in 1:nrow(All)){ 
  if ((All[i,2] == "Negative") & (All[i,5] == "0")) {
    All[i,2]<-All[i,2] 
  } else if ((All[i,2] == "Negative") & (All[i,5] == "4")) {
    All[i,2]<-All[i,2]
  } else if (All[i,2] == "Negative") {
    All[i,2]<-All[i,5]
  } else {
    All[i,2]<-All[i,2]
  }
}

for (i in 1:nrow(All)){ 
  if (All[i,2] == "Negative") {
    All[i,4]<-All[i,3] 
  } else if (All[i,3] == "Negative") {
    All[i,3]<-All[i,4]
  } else {
    All[i,3]<-All[i,3]
  }
}

table(All$HTO)
table(All$HTO_status)

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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 25)
HTO_number_18lane2 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_18lane2 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="18_lane2")
dataset@meta.data$lane<-"18_2"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "CD3E")
saveRDS(dataset, file="Data/preprocessed_lanes/200618_lane2_sct.rds")


#Read dataframe 200618_lane3 (PBMC)
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200618_lane3/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.9922)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200618_lane3.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
options(max.print=9999)
table(All$HTO, All$assignment)
All$assignment[All$assignment == "3"]="V050-T1-180619"
All$status[All$status == "doublet"]="Doublet"
All$status[All$status == "singlet"]="Singlet"
All$status[All$status == "unassigned"]="Negative"

for (i in 1:nrow(All)){ 
  if ((All[i,2] == "Negative") & (All[i,5] == "V050-T1-180619")) {
    All[i,2]<-All[i,5] 
  } else {
    All[i,2]<-All[i,2]
  }
}

for (i in 1:nrow(All)){ 
  if (All[i,2] == "Negative") {
    All[i,4]<-All[i,3] 
  } else if (All[i,3] == "Negative") {
    All[i,3]<-All[i,4]
  } else {
    All[i,3]<-All[i,3]
  }
}

table(All$HTO)
table(All$HTO_status)

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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 25)
HTO_number_18lane3 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_18lane3 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="18_lane3")
dataset@meta.data$lane<-"18_3"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "CD3E")
saveRDS(dataset, file="Data/preprocessed_lanes/200618_lane3_sct.rds")


#Read dataframe 200618_lane4 (PBMC)
dataset_10x<-Read10X("Data/filtered_feature_bc_matrix/200618_lane4/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Vedo2")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.9815)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = 'Data/souporcell_clusters/clusters_200618_lane4.tsv', sep = '\t', header = TRUE)
souporcell <- select(souporcell, c("barcode", "status", "assignment"))
All <- merge(seurat_call, souporcell, by = "barcode", all = F)
rownames(All) <- All$barcode
options(max.print=9999)
table(All$HTO, All$assignment)
All$assignment[All$assignment == "1"]="V001-T4-071119"
All$assignment[All$assignment == "2"]="V004-T1-130919"
All$assignment[All$assignment == "3"]="V050-T3-300719"
All$assignment[All$assignment == "4"]="V007-T2-241219"
All$assignment[All$assignment == "5"]="V008-T0-151019"
All$assignment[All$assignment == "6"]="V006-T4-100220"
All$status[All$status == "doublet"]="Doublet"
All$status[All$status == "singlet"]="Singlet"
All$status[All$status == "unassigned"]="Negative"

for (i in 1:nrow(All)){ 
  if ((All[i,2] == "Negative") & (All[i,5] == "0")) {
    All[i,2]<-All[i,2]
  } else if (All[i,2] == "Negative") {
    All[i,2]<-All[i,5]
  } else {
    All[i,2]<-All[i,2]
  }
}

for (i in 1:nrow(All)){ 
  if (All[i,2] == "Negative") {
    All[i,4]<-All[i,3] 
  } else if (All[i,3] == "Negative") {
    All[i,3]<-All[i,4]
  } else {
    All[i,3]<-All[i,3]
  }
}

table(All$HTO)
table(All$HTO_status)

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
table(dataset$Final_HTO, dataset$HTO_classification)

#qc and preprocess
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
VlnPlot(dataset, "percent.mt")
dataset<-subset(dataset, subset=percent.mt < 25)
HTO_number_18lane4 <- sum(dataset$HTO_classification.global == "Singlet")
souporcell_number_18lane4 <- sum(dataset$Final_HTO_status == "Singlet")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
dataset<-RenameCells(dataset, add.cell.id="18_lane4")
dataset@meta.data$lane<-"18_4"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "CD3E")
saveRDS(dataset, file="Data/preprocessed_lanes/200618_lane4_sct.rds")

