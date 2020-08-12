# for all of the datasets, pre-process data

### Proteasedata

dataset<-Read10X("~/filtered_feature_bc_matrix/")
dataset<-suppressWarnings(CreateSeuratObject(counts = dataset, project = "Protease"))

# process ADT/HTO now if present, see below script

DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
dataset<-suppressWarnings(SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F))
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.1, verbose = FALSE)

saveRDS(dataset, file="~/dataset_sct.rds")

dataset<-Read10X("~/filtered_feature_bc_matrix/")
dataset<-suppressWarnings(CreateSeuratObject(counts = dataset, project = "Protease"))

# process ADT/HTO now if present, see below script

DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
dataset<-suppressWarnings(SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F))
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.1, verbose = FALSE)

saveRDS(dataset, file="~/dataset_sct.rds")


dataset<-Read10X("~/Desktop/methods_paper/Sangerdata_202004/OTARscRNA8356111/filtered_feature_bc_matrix/")
dataset<-suppressWarnings(CreateSeuratObject(counts = dataset, project = "Protease"))

# process ADT/HTO now if present, see below script

DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
dataset<-suppressWarnings(SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F))
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.1, verbose = FALSE)

saveRDS(dataset, file="~/dataset_sct.rds")

dataset<-Read10X("~/filtered_feature_bc_matrix/")
dataset<-suppressWarnings(CreateSeuratObject(counts = dataset, project = "Protease"))

# process ADT/HTO now if present, see below script

DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
dataset<-suppressWarnings(SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F))
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.1, verbose = FALSE)

saveRDS(dataset, file="~/dataset_sct.rds")

dataset<-Read10X("~/filtered_feature_bc_matrix/")
dataset<-suppressWarnings(CreateSeuratObject(counts = dataset, project = "Protease"))

# process ADT/HTO now if present, see below script

DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
dataset<-suppressWarnings(SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F))
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.1, verbose = FALSE)

saveRDS(dataset, file="~/dataset_sct.rds")

dataset<-Read10X("~/filtered_feature_bc_matrix/")
dataset<-suppressWarnings(CreateSeuratObject(counts = dataset, project = "Protease"))

# process ADT/HTO now if present, see below script

DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
dataset<-suppressWarnings(SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F))
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.1, verbose = FALSE)

saveRDS(dataset, file="~/dataset_sct.rds")

dataset<-Read10X("~/filtered_feature_bc_matrix/")
dataset<-suppressWarnings(CreateSeuratObject(counts = dataset, project = "Protease"))

# process ADT/HTO now if present, see below script

DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
dataset<-suppressWarnings(SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F))
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.1, verbose = FALSE)

saveRDS(dataset, file="~/dataset_sct.rds")

dataset<-Read10X("~/filtered_feature_bc_matrix/")
dataset<-suppressWarnings(CreateSeuratObject(counts = dataset, project = "Protease"))

# process ADT/HTO now if present, see below script

DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
dataset<-suppressWarnings(SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F))
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.1, verbose = FALSE)

saveRDS(dataset, file="~/dataset_sct.rds")

dataset<-Read10X("~/filtered_feature_bc_matrix/")
dataset<-suppressWarnings(CreateSeuratObject(counts = dataset, project = "Protease"))

# process ADT/HTO now if present, see below script

DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)
dataset<-suppressWarnings(SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F))
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.1, verbose = FALSE)

saveRDS(dataset, file="~/dataset_sct.rds")


DimPlot(dataset)

### smillie data - previously done

epi_fresh_split_smillie<-Read10X("/epi_smillie_2019/", unique.features = TRUE,  gene.column=1)
imm_fresh_split_smillie<-Read10X("/imm_smillie_2019/", unique.features = TRUE,  gene.column=1)
fib_fresh_split_smillie<-Read10X("/fibro_smillie_2019/", unique.features = TRUE,  gene.column=1)
epi_fresh_split <- CreateSeuratObject(counts = epi_fresh_split_smillie, project = "Smillie", min.cells = 3, min.features = 200)
imm_fresh_split <- CreateSeuratObject(counts = imm_fresh_split_smillie, project = "Smillie", min.cells = 3, min.features = 200)
fib_fresh_split <- CreateSeuratObject(counts = fib_fresh_split_smillie, project = "Smillie", min.cells = 3, min.features = 200)

epi_fresh_split[["percent.mt"]] <- PercentageFeatureSet(epi_fresh_split, pattern = "^MT-")
epi_fresh_split<-subset(epi_fresh_split, subset=percent.mt < 60)

epi_fresh_split<-SCTransform(epi_fresh_split,vars.to.regress = "percent.mt", verbose = F)
epi_fresh_split <- RunPCA(epi_fresh_split, verbose = FALSE)
epi_fresh_split <- RunUMAP(epi_fresh_split, dims = 1:30, verbose = FALSE)
epi_fresh_split <- FindNeighbors(epi_fresh_split, dims = 1:30, verbose = FALSE)
epi_fresh_split <- FindClusters(epi_fresh_split, resolution = 0.2, verbose = FALSE)
save(epi_fresh_split, file="/epi_fresh_split_sct.Rd")

imm_fresh_split[["percent.mt"]] <- PercentageFeatureSet(imm_fresh_split, pattern = "^MT-")
imm_fresh_split<-subset(imm_fresh_split, subset=percent.mt < 60)

imm_fresh_split<-SCTransform(imm_fresh_split,vars.to.regress = "percent.mt", verbose = F)
imm_fresh_split <- RunPCA(imm_fresh_split, verbose = FALSE)
imm_fresh_split <- RunUMAP(imm_fresh_split, dims = 1:30, verbose = FALSE)
imm_fresh_split <- FindNeighbors(imm_fresh_split, dims = 1:30, verbose = FALSE)
imm_fresh_split <- FindClusters(imm_fresh_split, resolution = 0.2, verbose = FALSE)
save(imm_fresh_split, file="/imm_fresh_split_sct.Rd")

fib_fresh_split[["percent.mt"]] <- PercentageFeatureSet(fib_fresh_split, pattern = "^MT-")
fib_fresh_split<-subset(fib_fresh_split, subset=percent.mt < 60)

fib_fresh_split<-SCTransform(fib_fresh_split,vars.to.regress = "percent.mt", verbose = F)
fib_fresh_split <- RunPCA(fib_fresh_split, verbose = FALSE)
fib_fresh_split <- RunUMAP(fib_fresh_split, dims = 1:30, verbose = FALSE)
fib_fresh_split <- FindNeighbors(fib_fresh_split, dims = 1:30, verbose = FALSE)
fib_fresh_split <- FindClusters(fib_fresh_split, resolution = 0.2, verbose = FALSE)
save(fib_fresh_split, file="/fib_fresh_split_sct.Rd")

# cluster on gearshift and extract healthy cells only

## whole collagenase samples
# 3083		1
# 3002		1
# 3049		B
# 3030		A
# 3296		C
library(Seurat)
dataset_10x<-Read10X("~/filtered_feature_bc_matrix/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Whole_collagenase")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Custom`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset<-HTODemux(dataset, assay = "HTO", positive.quantile = 0.95)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)
Idents(dataset) <- "HTO_classification.global"
dataset <- subset(dataset, idents = "Singlet")
Idents(dataset) <- "HTO_maxID"
dataset <- subset(dataset, idents = "HC-NI-3002")
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)

dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
saveRDS(dataset, file="~/HC-NI-3002_sct.rds")

dataset_10x<-Read10X("~/filtered_feature_bc_matrix/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Whole_collagenase")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Custom`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset<-HTODemux(dataset, assay = "HTO", positive.quantile = 0.97)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)
Idents(dataset) <- "HTO_classification.global"
dataset <- subset(dataset, idents = "Singlet")
Idents(dataset) <- "HTO_maxID"
dataset <- subset(dataset, idents = "HC-NI-3129")
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)

dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
saveRDS(dataset, file="~/HC-NI-3129_sct.rds")
DimPlot(dataset)

####
dataset_10x<-Read10X("~/filtered_feature_bc_matrix/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Whole_collagenase")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Custom`)
dataset[['ADT']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- NormalizeData(dataset, assay = "ADT", normalization.method = "CLR")
dataset<-HTODemux(dataset, assay = "HTO", positive.quantile = 0.95)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)
Idents(dataset) <- "HTO_classification.global"
dataset <- subset(dataset, idents = "Singlet")
Idents(dataset) <- "HTO_maxID"
dataset <- subset(dataset, idents = "HC-NI-3296")
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)

dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "CD103-TotalA")
saveRDS(dataset, file="~/HC-NI-3296_sct.rds")

###
dataset_10x<-Read10X("~/filtered_feature_bc_matrix/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Whole_collagenase")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Custom`)
#dataset[['ADT']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
#dataset <- NormalizeData(dataset, assay = "ADT", normalization.method = "CLR")
dataset<-HTODemux(dataset, assay = "HTO", positive.quantile = 0.93)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
Idents(dataset) <- "HTO_classification.global"
dataset <- subset(dataset, idents = "Singlet")
Idents(dataset) <- "HTO_maxID"
dataset <- subset(dataset, idents = "HC-NI-3030")
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)

dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
FeaturePlot(dataset, "CD103-TotalA")
saveRDS(dataset, file="~/HC-NI-3030_sct.rds")


#######
dataset_10x<-Read10X("~/filtered_feature_bc_matrix/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Whole_collagenase")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
#dataset[['ADT']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
#dataset <- NormalizeData(dataset, assay = "ADT", normalization.method = "CLR")
#dataset <- ScaleData(dataset, assay = "ADT")
dataset<-HTODemux(dataset, assay = "HTO", positive.quantile = 0.87)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)
Idents(dataset) <- "HTO_classification.global"
dataset <- subset(dataset, idents = "Singlet")
Idents(dataset) <- "HTO_maxID"
dataset <- subset(dataset, idents = "HC-NI-3083")
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)

dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
#FeaturePlot(dataset, "CD103-TotalA")
saveRDS(dataset, file="~/HC-NI-3083_sct.rds")


#######
dataset_10x<-Read10X("~/filtered_feature_bc_matrix/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Whole_collagenase")
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
#dataset[['ADT']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
#dataset <- NormalizeData(dataset, assay = "ADT", normalization.method = "CLR")
#dataset <- ScaleData(dataset, assay = "ADT")
dataset<-HTODemux(dataset, assay = "HTO", positive.quantile = 0.93)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)
Idents(dataset) <- "HTO_classification.global"
dataset <- subset(dataset, idents = "Singlet")
Idents(dataset) <- "HTO_maxID"
dataset <- subset(dataset, idents = "HC-NI-3049")
DefaultAssay(dataset)<-"RNA"
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset<-subset(dataset, subset=percent.mt < 60)

dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
#FeaturePlot(dataset, "CD103-TotalA")
saveRDS(dataset, file="~/HC-NI-3049_sct.rds")

#######
dataset_10x<-Read10X("~/filtered_feature_bc_matrix/")
dataset<-CreateSeuratObject(counts = dataset_10x$`Gene Expression`, project = "Whole_collagenase")
dataset_10x$`Custom`<-dataset_10x$`Antibody Capture`[c(1:15),]
dataset_10x$`Antibody Capture`<-dataset_10x$`Antibody Capture`[c(16,17,18,19),]
dataset[['HTO']]<-CreateAssayObject(counts = dataset_10x$`Antibody Capture`)
dataset[['ADT']]<-CreateAssayObject(counts = dataset_10x$`Custom`)
dataset<-subset(dataset, nCount_RNA > 200)
dataset <- NormalizeData(dataset, assay = "HTO", normalization.method = "CLR")
dataset <- NormalizeData(dataset, assay = "ADT", normalization.method = "CLR")
dataset <- ScaleData(dataset, assay = "ADT")
dataset<-HTODemux(dataset, assay = "HTO", positive.quantile = 0.999)
table(dataset$HTO_classification.global)
table(dataset$HTO_maxID)
HTOHeatmap(dataset, assay = "HTO", ncells = 1000)
Idents(dataset) <- "HTO_classification.global"
dataset <- subset(dataset, idents = "Singlet")
Idents(dataset) <- "HTO_maxID"
dataset_1 <- subset(dataset, idents = "C-3034-HC")
dataset_2 <- subset(dataset, idents = "D-3037-HC")

DefaultAssay(dataset_1)<-"RNA"
dataset_1[["percent.mt"]] <- PercentageFeatureSet(dataset_1, pattern = "^MT-")
dataset<-subset(dataset_1, subset=percent.mt < 60)

dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
#FeaturePlot(dataset, "CD103-TotalA")
saveRDS(dataset, file="~/Desktop/CITEseq/C-3034-HC_sct.rds")
dataset_1<-dataset

DefaultAssay(dataset_2)<-"RNA"
dataset_2[["percent.mt"]] <- PercentageFeatureSet(dataset_2, pattern = "^MT-")
dataset<-subset(dataset_2, subset=percent.mt < 60)

dataset<-SCTransform(dataset,vars.to.regress = "percent.mt", verbose = F)
dataset <- RunPCA(dataset, verbose = FALSE)
dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:30, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = 0.2, verbose = FALSE)
DimPlot(dataset)
saveRDS(dataset, file="~/D-3037-HC_sct.rds")


