#Read dataframe 08 lane2
dataset_10x<-Read10X("../VEDO2/Data/08_lane2/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "PSC")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- HTODemux(dataset, assay = "HTO", positive.quantile = 0.993)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)

#demultiplex
HTO <- dataset@meta.data[["HTO_classification"]]
HTO_status <- dataset@meta.data[["HTO_classification.global"]]
barcode <- dataset@assays[["RNA"]]@data@Dimnames[[2]]
seurat_call <- data.frame(HTO, HTO_status, barcode)
souporcell <- read.table(file = '../VEDO2/Data/clusters_08_lane2.tsv', sep = '\t', header = TRUE)
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
Idents(dataset) <- "Final_HTO"
dataset <- subset(dataset, idents = "PSC-NI-3191")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
table(dataset@meta.data[["Final_HTO"]])
dataset<-RenameCells(dataset, add.cell.id="VEDO01_T0_NI")
dataset@meta.data$lane<-"VEDO01_T0_NI"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "IGHA2")
FeaturePlot(dataset, "EPCAM")
FeaturePlot(dataset, "CD3E")
dataset@meta.data[["Final_HTO"]] <- "V001_T0_transversum_NI"
table(dataset@meta.data[["Final_HTO"]])
saveRDS(dataset, file="../VEDO2/Data/VEDO01_T0_NI_sct.rds")

#Read dataframe 09 lane1
dataset_10x<-Read10X("../VEDO2/Data/09_lane1/")
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
souporcell <- read.table(file = '../VEDO2/Data/clusters_09_lane1.tsv', sep = '\t', header = TRUE)
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
Idents(dataset) <- "Final_HTO"
dataset <- subset(dataset, idents = "PSC-I-3191")
Idents(dataset) <- "Final_HTO_status"
dataset <- subset(dataset, idents = "Singlet")
table(dataset@meta.data[["Final_HTO"]])
dataset<-RenameCells(dataset, add.cell.id="VEDO01_T0_I")
dataset@meta.data$lane<-"VEDO01_T0_I"
dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "IGHG4")
FeaturePlot(dataset, "EPCAM")
dataset@meta.data[["Final_HTO"]] <- "V001_T0_ascendens_I"
table(dataset@meta.data[["Final_HTO"]])
saveRDS(dataset, file="../VEDO2/Data/VEDO01_T0_I_sct.rds")

