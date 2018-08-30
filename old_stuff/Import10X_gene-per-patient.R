## 10X import for Takeda RIPK2
## Cytotoxic cells
## 10X import for Takeda RIPK2
## Cytotoxic cells

library(Seurat)

## Cytotoxic T cells

Cytotox <- Read10X(data.dir = "/Users/festeneam/Downloads/filtered_matrices_mex_cytotox/hg19")
Cytotox_T <- CreateSeuratObject(raw.data = Cytotox, min.cells = 3,
                                min.genes = 200,
                                project = "10X Project",
                                names.delim = "-",
                                names.field = 2)

Cytotox_T <- SetIdent(Cytotox_T, ident.use = c("Cyto"))

mito.genes <- grep(pattern = "^MT-", x = rownames(x = Cytotox_T@data), value = TRUE) 

# Get TSS normalized mitochodrial counts
col.total.counts <- Matrix::colSums(expm1(Cytotox_T@data))
percent.mito <- Matrix::colSums(Cytotox_T@raw.data[mito.genes,])/Matrix::colSums(Cytotox_T@raw.data)

Cytotox_T <- AddMetaData(object = Cytotox_T, metadata = percent.mito, col.name = "percent.mito")

Cytotox_T <- NormalizeData(object = Cytotox_T, normalization.method = "LogNormalize",
                           scale.factor = 1e4)

Cytotox_T <- ScaleData(Cytotox_T)
Cytotox_T <- RunPCA(object = Cytotox_T, pc.genes=rownames(Cytotox_T@data),
                    do.print = TRUE,
                    pcs.print = 1:5,
                    genes.print = 5)

PrintPCA(object = Cytotox_T, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

Cytotox_T <- RunTSNE(object = Cytotox_T, dims.use = 1:10, do.fast = TRUE) 
TSNEPlot(object = Cytotox_T)

VlnPlot(Cytotox_T, c("ZBTB7B","CD8B", "GZMA","RIPK2"), nCol = 4)

FeaturePlot(Cytotox_T, c("ZBTB7B","CD8B", "GZMA","RIPK2"), cols.use = c("lightgrey","deeppink"), 
            pt.size = 1)

Cytotox_T <- FindClusters(object = Cytotox_T, reduction.type = "pca",
                          dims.use = 1:10,
                          resolution = 0.6,
                          print.output = 0,
                          save.SNN = TRUE)
## Monocytes

CD14 <- Read10X(data.dir = "/Users/festeneam/Downloads/filtered_matrices_mex_CD14/hg19")
CD14_mon <- CreateSeuratObject(raw.data = CD14, min.cells = 3,
                               min.genes = 200,
                               project = "10X Project",
                               names.delim = "-",
                               names.field = 2)

CD14_mon <- SetIdent(CD14_mon, ident.use = c("CD14"))

mito.genes <- grep(pattern = "^MT-", x = rownames(x = CD14_mon@data), value = TRUE) 

# Get TSS normalized mitochodrial counts
col.total.counts <- Matrix::colSums(expm1(CD14_mon@data))
percent.mito <- Matrix::colSums(CD14_mon@raw.data[mito.genes,])/Matrix::colSums(CD14_mon@raw.data)

CD14_mon <- AddMetaData(object = CD14_mon, metadata = percent.mito, col.name = "percent.mito")

CD14_mon <- NormalizeData(object = CD14_mon, normalization.method = "LogNormalize",
                          scale.factor = 1e4)

CD14_mon <- ScaleData(CD14_mon)
CD14_mon <- RunPCA(object = CD14_mon, pc.genes=rownames(CD14_mon@data),
                   do.print = TRUE,
                   pcs.print = 1:5,
                   genes.print = 5)

PrintPCA(object = CD14_mon, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

CD14_mon <- RunTSNE(object = CD14_mon, dims.use = 1:10, do.fast = TRUE) 
TSNEPlot(object = CD14_mon)

VlnPlot(CD14_mon, c("CD14","SELL", "HLA-DRA","RIPK2"), nCol = 4)

FeaturePlot(CD14_mon, c("CD14","SELL", "HLA-DRA","RIPK2"), cols.use = c("lightgrey","deeppink"), 
            pt.size = 1)

CD14_mon <- FindClusters(object = CD14_mon, reduction.type = "pca",
                         dims.use = 1:10,
                         resolution = 0.6,
                         print.output = 0,
                         save.SNN = TRUE)

## NK cells

CD56 <- Read10X(data.dir = "/Users/festeneam/Downloads/filtered_matrices_mex_CD56/hg19")
CD56_mon <- CreateSeuratObject(raw.data = CD56, min.cells = 3,
                               min.genes = 200,
                               project = "10X Project",
                               names.delim = "-",
                               names.field = 2)

CD56_mon <- SetIdent(CD56_mon, ident.use = c("NK"))

mito.genes <- grep(pattern = "^MT-", x = rownames(x = CD56_mon@data), value = TRUE) 

# Get TSS normalized mitochodrial counts
col.total.counts <- Matrix::colSums(expm1(CD56_mon@data))
percent.mito <- Matrix::colSums(CD56_mon@raw.data[mito.genes,])/Matrix::colSums(CD56_mon@raw.data)

CD56_mon <- AddMetaData(object = CD56_mon, metadata = percent.mito, col.name = "percent.mito")

CD56_mon <- NormalizeData(object = CD56_mon, normalization.method = "LogNormalize",
                          scale.factor = 1e4)

CD56_mon <- ScaleData(CD56_mon)
CD56_mon <- RunPCA(object = CD56_mon, pc.genes=rownames(CD56_mon@data),
                   do.print = TRUE,
                   pcs.print = 1:5,
                   genes.print = 5)

PrintPCA(object = CD56_mon, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

CD56_mon <- RunTSNE(object = CD56_mon, dims.use = 1:10, do.fast = TRUE) 
TSNEPlot(object = CD56_mon)

VlnPlot(CD56_mon, c("NCAM1","FCGR3A", "ITGAM","RIPK2"), nCol = 4)

FeaturePlot(CD56_mon, c("NCAM1","FCGR3A", "ITGAM","RIPK2"), cols.use = c("lightgrey","deeppink"), 
            pt.size = 1)

CD56_mon <- FindClusters(object = CD56_mon, reduction.type = "pca",
                         dims.use = 1:10,
                         resolution = 0.6,
                         print.output = 0,
                         save.SNN = TRUE)

total <- MergeSeurat(NaiveCD4_T, CD14_mon, add.cell.id1 = "CD4_Naive", add.cell.id2 = "CD14")
total2 <- MergeSeurat(total, CD56_mon, add.cell.id2 = "NK")
TenX_cells <- MergeSeurat(total2, Cytotox_T, add.cell.id2 = "Cytotox")

TenX_cells <- ScaleData(TenX_cells)
TenX_cells <- RunPCA(object = TenX_cells, pc.genes=rownames(TenX_cells@data),
                   do.print = TRUE,
                   pcs.print = 1:5,
                   genes.print = 5)

PrintPCA(object = TenX_cells, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

TenX_cells <- RunTSNE(object = TenX_cells, dims.use = 1:10, do.fast = TRUE) 
TSNEPlot(object = TenX_cells)

VlnPlot(TenX_cells, c("OSM","NOD2","OPRL1","RIPK2"), nCol = 4)

FeaturePlot(TenX_cells, c("CCR7","CD160","GZMA","RIPK2"), cols.use = c("lightgrey","deeppink"), 
            pt.size = 1)

## RIPK2 positivity
# Table positive cells for a gene per patient

RIPK2.pos <- TenX_cells@data["RIPK2",]>0
RIPK2.pos <- 1*RIPK2.pos

TenX_cells <- AddMetaData(TenX_cells, metadata = RIPK2.pos, col.name = "RIPK2.pos")
TenX_cells <- AddMetaData(TenX_cells, metadata = TenX_cells@ident, col.name = "ident")


