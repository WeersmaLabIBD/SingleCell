# Author: WTC
# Date: July 2018
# Script for sampling 251 cells from 10x CTL dataset, 
#   combine those (using Seurat's CCA function) with 251 CTLs from CD dataset, randomly assign source IDs 
#   and perform differential expression analyses on these randomly made groups of cells, each existing of 251 cells. 

# download Cytotoxic T cell gene/cell matrix (raw) from https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/cytotoxic_t
load("~/Desktop/Single_cell/CTL_10x_CD/CTL_10x.Rda")
load("~/Desktop/Single_cell/CTL_10x_CD/CTL_CD.Rda") # import CD CTL dataset
risk_use<-read.csv("~/Desktop/Single_cell/final_data_paper/DE/risk_genes_in_cd_and_10x.csv")
genes_use<-read.csv("~/Desktop/Single_cell/final_data_paper/DE/genes_in_cd_and_10x.csv")

a<-"number_DEgenes"
b<-"number_DEriskgenes"
number_DE_genes<-data.frame(a)
number_DE_riskgenes<-data.frame(b)
for(i in (1:1000)){
cat("iteration", i, "\n")
# select subset of 251 cells
samp_col_idxs_1  <- sample(ncol(CTL_10x@data), 251)
samp_col_names_1 <- colnames(CTL_10x@data) [samp_col_idxs_1]
CTL_10x_subset<-CTL_10x@data[,samp_col_names_1]
dim(CTL_10x_subset)
CTL_10x_subset<-CreateSeuratObject(raw.data = CTL_10x_subset, min.cells = 3, project = "subset")

mito.genes <- grep(pattern = "^MT-", x = rownames(x = CTL_10x_subset@data), value = TRUE)
percent.mito <- Matrix::colSums(CTL_10x_subset@raw.data[mito.genes, ])/Matrix::colSums(CTL_10x_subset@raw.data)
CTL_10x_subset <- AddMetaData(object = CTL_10x_subset, metadata = percent.mito, col.name = "percent.mito")

# for ngene: 200-2500, percent mito: <5%
CTL_10x_subset <- FilterCells(object = CTL_10x_subset, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

# log normalize
CTL_10x_subset <- NormalizeData(object = CTL_10x_subset, normalization.method = "LogNormalize", scale.factor = 10000)
CTL_10x_subset <- ScaleData(object = CTL_10x_subset, vars.to.regress = c("nUMI", "percent.mito"))
dim(CTL_10x_subset@data)

CTL_10x_subset@meta.data$dataset <- "subset10x"

save(CTL_10x_subset, file=paste0("CTL_10x_subset", i, ".Rda"))


## same for CD
CTL_CD@meta.data$dataset <- "cd"

## Combine the cdulated and ctrlulated cells into a single object
CTL_subset.combined <- RunCCA(CTL_10x_subset, CTL_CD, genes.use = intersect(rownames(CTL_10x_subset@data), rownames(CTL_CD@data)), num.cc = 30, scale.data=T)

###############################################################
# Process the data after combining the two datasets
###############################################################

#CTL_subset.combined <- RunPCA(object = CTL_subset.combined, pc.genes = CTL_subset.combined@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
#CTL_subset.combined = RunTSNE(CTL_subset.combined, dims.use = 1:15, do.fast = T)

## CCA equivalent of the PCElbowPlot
#MetageneBicorPlot(CTL_subset.combined, grouping.var = "dataset", dims.eval = 1:30, display.progress = FALSE)

## Use the CCA reduction to better overlap the two original datasetes
#CTL_subset.combined <- AlignSubspace(CTL_subset.combined, reduction.type = "cca", grouping.var = "dataset", dims.align = 1:15)

#CTL_subset.combined <- RunTSNE(CTL_subset.combined, reduction.use = "cca.aligned", dims.use = 1:15, do.fast = T)
#CTL_subset.combined <- FindClusters(CTL_subset.combined, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:15)
#save(CTL_subset.combined, file=paste0("CTL_subset_", i, "combined_aligned.Rda"))

### add metadata to randomize cell ids


## add location
#extract file that is a copy from @data.info
CellsMeta = CTL_subset.combined@meta.data
CellsMeta$random_asssignment <- "A" 
CellsMeta$random_asssignment[sample(nrow(CellsMeta), abs(nrow(CellsMeta)/2), replace=FALSE)] <- "B"
CTL_subset.combined <- AddMetaData(CTL_subset.combined, CellsMeta)

CTL_subset.combined<-SetAllIdent(CTL_subset.combined, "random_asssignment")
subset_A<-WhichCells(CTL_subset.combined, "A")
subset_B<-WhichCells(CTL_subset.combined, "B")

# calculate DE

# genes present in the two datasets:
# DE with genes present in the two datasets, in >1%, MAST
CTL_subset.combined_DE_markers = FindAllMarkers(CTL_subset.combined, only.pos = T, test.use = "MAST", genes.use = genes_use$gene)
write.table(CTL_subset.combined_DE_markers, file = paste0("CTL_subset", i, "_DEmarkers_MAST.txt"))

# and the significant ones (p<0.05)
CTL_subset.combined_DE_markers<-CTL_subset.combined_DE_markers[CTL_subset.combined_DE_markers$p_val_adj <0.05,]
number_DE_genes[i,a]<-nrow(CTL_subset.combined_DE_markers)

# merge with risk genes
# merge all DE results with risk genes
risk_all_filtered_subset<-merge(risk_use, CTL_subset.combined_DE_markers, by="gene", all=F)

number_DE_riskgenes[i,b]<-nrow(risk_all_filtered_subset)}

write.csv(number_DE_genes, "number_DE_genes.csv")
write.csv(number_DE_riskgenes, "number_DE_riskgenes.csv")



