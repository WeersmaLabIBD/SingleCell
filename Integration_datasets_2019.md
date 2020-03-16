---
# "Integration of multiple datasets on himem cluster"
author: W. Uniken Venema
---
## Prepare files
```
library(Seurat)
install.packages("devtools")
library(devtools)
install.packages("hdf5r")
library(hdf5r) <!--install libraries-->
```
open files
```
TI_frozen_coldprotease<-Read10X_h5("~/x/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
TI_frozen_collagenase<-Read10X_h5("~/x/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
TI_fresh_collagenase<-Read10X_h5("~/x/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
Rectum_fresh_coldprotease<-Read10X_h5("~/x/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
TI_fresh_coldprotease<-Read10X("~/x/filtered_feature_bc_matrix/")
Colon_frozen_collagenase<-Read10X("~/x/filtered_feature_bc_matrix/", unique.features = TRUE)
```
create seurat objects
```
TI_fresh_prot <- CreateSeuratObject(counts = TI_fresh_coldprotease, project = "protocols", min.cells = 3, min.features = 200)
TI_frozen_prot <- CreateSeuratObject(counts = TI_frozen_coldprotease, project = "protocols", min.cells = 3, min.features = 200)
TI_fresh_coll <- CreateSeuratObject(counts = TI_fresh_collagenase, project = "protocols", min.cells = 3, min.features = 200)
TI_frozen_coll <- CreateSeuratObject(counts = TI_frozen_collagenase, project = "protocols", min.cells = 3, min.features = 200)
Re_fresh_prot<-CreateSeuratObject(counts = Rectum_fresh_coldprotease, project = "protocols", min.cells = 3, min.features = 200)
Col_frozen_coll<-CreateSeuratObject(counts = Colon_frozen_collagenase$`Gene Expression`, project = "protocols")
```
add hashtag data
```
Colon_frozen_collagenase$`HTO Capture` <- Colon_frozen_collagenase$`Antibody Capture`
Colon_frozen_collagenase$`HTO Capture`<-Colon_frozen_collagenase$`HTO Capture`[c(16:19),]
Colon_frozen_collagenase$`Antibody Capture`<-Colon_frozen_collagenase$`Antibody Capture`[c(1:15),]
Col_frozen_coll[['HTO']]<-CreateAssayObject(counts = Colon_frozen_collagenase$`HTO Capture`)
Col_frozen_coll[['ADT']]<-CreateAssayObject(counts = Colon_frozen_collagenase$`Antibody Capture`)
Col_frozen_coll<-subset(Col_frozen_coll, nCount_RNA > 200)
```
set percentage mitochondrial reads
```
TI_fresh_prot[["percent.mt"]] <- PercentageFeatureSet(TI_fresh_prot, pattern = "^MT-")
TI_fresh_coll[["percent.mt"]] <- PercentageFeatureSet(TI_fresh_coll, pattern = "^MT-")
TI_frozen_prot[["percent.mt"]] <- PercentageFeatureSet(TI_frozen_prot, pattern = "^MT-")
TI_frozen_coll[["percent.mt"]] <- PercentageFeatureSet(TI_frozen_coll, pattern = "^MT-")
Re_fresh_prot[["percent.mt"]] <- PercentageFeatureSet(Re_fresh_prot, pattern = "^MT-")
Col_frozen_coll[["percent.mt"]] <- PercentageFeatureSet(Col_frozen_coll, pattern = "^MT-")
```
appoint hashtags
```
Col_frozen_coll <- NormalizeData(Col_frozen_coll, assay = "HTO", normalization.method = "CLR")
Col_frozen_coll<-MULTIseqDemux(Col_frozen_coll, assay = "HTO", quantile = 0.99)
Col_frozen_coll<-HTODemux(Col_frozen_coll, assay = "HTO", positive.quantile = 0.99)
```
select only singlets
```
table(Col_frozen_coll$HTO_classification.global)
table(Col_frozen_coll$HTO_maxID)
Idents(Col_frozen_coll) <- "HTO_classification.global"
Col_frozen_coll_singlets <- subset(Col_frozen_coll, idents = "Singlet")
DefaultAssay(Col_frozen_coll_singlets)<-"RNA"
```
remove cells with >60% mito counts as these will influence clustering too much
```
TI_fresh_prot <- subset(TI_fresh_prot, subset = percent.mt < 60)
dim(TI_fresh_prot@assays$RNA)
TI_fresh_coll <- subset(TI_fresh_coll, subset = percent.mt < 60)
dim(TI_fresh_coll@assays$RNA)
TI_frozen_prot <- subset(TI_frozen_prot, subset = percent.mt < 60)
dim(TI_frozen_prot@assays$RNA)
TI_frozen_coll <- subset(TI_frozen_coll, subset = percent.mt < 60)
dim(TI_frozen_coll@assays$RNA)
Re_fresh_prot <- subset(Re_fresh_prot, subset = percent.mt < 60)
dim(Re_fresh_prot@assays$RNA)
Col_frozen_coll_singlets <- subset(Col_frozen_coll_singlets, subset = percent.mt < 60)
dim(Col_frozen_coll_singlets@assays$RNA)
```
SCtransform and regress mitochondrial %, run PCA, UMAP, find neighbors and clusters
```
Col_frozen_coll_singlets<-SCTransform(Col_frozen_coll_singlets,vars.to.regress = "percent.mt", verbose = F)
Col_frozen_coll_singlets <- RunPCA(Col_frozen_coll_singlets, verbose = FALSE)
Col_frozen_coll_singlets <- RunUMAP(Col_frozen_coll_singlets, dims = 1:30, verbose = FALSE)
Col_frozen_coll_singlets <- FindNeighbors(Col_frozen_coll_singlets, dims = 1:30, verbose = FALSE)
VlnPlot(Col_frozen_coll, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
Col_frozen_coll_singlets <- FindClusters(Col_frozen_coll_singlets, resolution = 0.2, verbose = FALSE)
DimPlot(Col_frozen_coll_singlets)

TI_fresh_prot<-SCTransform(TI_fresh_prot,vars.to.regress = "percent.mt", verbose = F)
TI_fresh_prot <- RunPCA(TI_fresh_prot, verbose = FALSE)
TI_fresh_prot <- RunUMAP(TI_fresh_prot, dims = 1:30, verbose = FALSE)
TI_fresh_prot <- FindNeighbors(TI_fresh_prot, dims = 1:30, verbose = FALSE)
TI_fresh_prot <- FindClusters(TI_fresh_prot, resolution = 0.2, verbose = FALSE)
DimPlot(TI_fresh_prot)

TI_fresh_coll<-SCTransform(TI_fresh_coll,vars.to.regress = "percent.mt", verbose = F)
TI_fresh_coll <- RunPCA(TI_fresh_coll, verbose = FALSE)
TI_fresh_coll <- RunUMAP(TI_fresh_coll, dims = 1:30, verbose = FALSE)
TI_fresh_coll <- FindNeighbors(TI_fresh_coll, dims = 1:30, verbose = FALSE)
TI_fresh_coll <- FindClusters(TI_fresh_coll, resolution = 0.2, verbose = FALSE)
DimPlot(TI_fresh_coll)

TI_frozen_coll<-SCTransform(TI_frozen_coll,vars.to.regress = "percent.mt", verbose = F)
TI_frozen_coll <- RunPCA(TI_frozen_coll, verbose = FALSE)
TI_frozen_coll <- RunUMAP(TI_frozen_coll, dims = 1:30, verbose = FALSE)
TI_frozen_coll <- FindNeighbors(TI_frozen_coll, dims = 1:30, verbose = FALSE)
TI_frozen_coll <- FindClusters(TI_frozen_coll, resolution = 0.2, verbose = FALSE)
DimPlot(TI_frozen_coll)

TI_frozen_prot<-SCTransform(TI_frozen_prot,vars.to.regress = "percent.mt", verbose = F)
TI_frozen_prot <- RunPCA(TI_frozen_prot, verbose = FALSE)
TI_frozen_prot <- RunUMAP(TI_frozen_prot, dims = 1:30, verbose = FALSE)
TI_frozen_prot <- FindNeighbors(TI_frozen_prot, dims = 1:30, verbose = FALSE)
TI_frozen_prot <- FindClusters(TI_frozen_prot, resolution = 0.2, verbose = FALSE)
DimPlot(TI_frozen_prot)

Re_fresh_prot<-SCTransform(Re_fresh_prot,vars.to.regress = "percent.mt", verbose = F)
Re_fresh_prot <- RunPCA(Re_fresh_prot, verbose = FALSE)
Re_fresh_prot <- RunUMAP(Re_fresh_prot, dims = 1:30, verbose = FALSE)
Re_fresh_prot <- FindNeighbors(Re_fresh_prot, dims = 1:30, verbose = FALSE)
Re_fresh_prot <- FindClusters(Re_fresh_prot, resolution = 0.2, verbose = FALSE)
DimPlot(Re_fresh_prot)
```
view number of cells left after QC
```
dim(TI_fresh_prot@assays$RNA)
dim(TI_fresh_coll@assays$RNA)
dim(TI_frozen_prot@assays$RNA)
dim(TI_frozen_coll@assays$RNA)
dim(Re_fresh_prot@assays$RNA)
dim(Col_frozen_coll_singlets@assays$RNA)
```



