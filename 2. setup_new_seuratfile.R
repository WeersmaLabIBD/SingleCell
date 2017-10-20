## create new seurat (v2.0) object of raw data file
datafile <- CreateSeuratObject(raw.data = x, min.cells = 3, min.genes = 200, project = "datafile_BROAD")
# discover mitochondrial genes (high percentage of total expressed genes can mark cell breakdown)
mito.genes <- grep(pattern = "^MT\\.", x = rownames(x = datafile@data), value = TRUE)
percent.mito <- colSums(datafile@raw.data[mito.genes, ])/colSums(datafile@raw.data)
datafile <- AddMetaData(object = datafile, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = datafile, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = datafile, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = datafile, gene1 = "nUMI", gene2 = "nGene")
# for ngene: 200-2500, percent mito: <5%
datafile <- FilterCells(object = datafile, subset.names = "nGene", low.thresholds = 200, high.thresholds = 2500)
datafile <- FilterCells(object = datafile, subset.names = "percent.mito", low.thresholds = -Inf, high.thresholds =  0.05)

# log normalize
datafile <- NormalizeData(object = datafile, normalization.method = "LogNormalize", scale.factor = 10000)

# add metadata
metadata<-read.csv("..")
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
## add location
#extract file that is a copy from @data.info
CellsMeta = datafile@meta.data
head(CellsMeta)
CellsMeta$cellen=row.names(CellsMeta)
row.names(CellsMeta)=NULL
dim(CellsMeta)
CellsMeta<-keeping.order(CellsMeta, merge, y=metadata, by = "cellen", all=FALSE)
row.names(CellsMeta)<-CellsMeta$cellen
dim(CellsMeta)
CellsMeta$cellen<-NULL
head(CellsMeta)
datafile <- AddMetaData(datafile, CellsMeta)

# scale, perform PCAs and JackStraw analysis to define significance of PCAs
datafile <- ScaleData(object = datafile, vars.to.regress = c("nUMI", "percent.mito"))
datafile<-FindVariableGenes(datafile, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 2, y.cutoff = 2)
datafile <- RunPCA(object = datafile, pc.genes = row.names(datafile@data), do.print = TRUE, pcs.print = 1:5, genes.print = 5)
#PrintPCA(object = datafile, pcs.print = 1:7, genes.print = 5, use.full = FALSE)
#VizPCA(object = datafile, pcs.use = 1:2)
#PCAPlot(object = datafile, dim.1 = 1, dim.2 = 2)
datafile = JackStraw(datafile, num.replicate = 5, do.print = T)
datafile <- ProjectPCA(object = datafile, do.print = T)
#datafile <- RunPCA(object = datafile, pc.genes = row.names(datafile@data), do.print = TRUE, pcs.print = 1:5, genes.print = 5)
#PCHeatmap(object = datafile, pc.use = 6, cells.use = 500, do.balanced = TRUE,label.columns = FALSE)
#PCHeatmap(object = datafile, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
datafile_siggenes<-PCASigGenes(datafile, 1:20, pval.cut=1e-5, max.per.pc = 250)
datafile_sig<-RunPCA(datafile, pc.genes = datafile_siggenes, do.print=T)
datafile_sig <- JackStraw(object = datafile_sig, num.replicate = 100, do.print = FALSE)
# visualize significance of PCAs
JackStrawPlot(object = datafile, PCs = 1:15)
PCElbowPlot(object = datafile)

# Find clusters using relevant PCAs, adapt resolution to amount of clusters that is wishful/makes sense
datafile_sig <- FindClusters(object = datafile_sig, reduction.type = "pca", dims.use = 1:18, resolution = 0.6, print.output = 0, save.SNN = TRUE)
#PrintFindClustersParams(object = datafile)

# run TSNE and visualize
datafile_sig <- RunTSNE(object = datafile_sig, dims.use = 1:18, do.fast = TRUE)
TSNEPlot(object = datafile_sig)

