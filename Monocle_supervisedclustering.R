
library(monocle)

## The first step in working with Monocle is to load up your data into Monocle's main class, 
## CellDataSet:
        
#samplesheet = 
#read.csv("/Users/festeneam/Artikelen/ScRNAseq/run02/allcellsallptsextra/Monocle/allcells_allptsextra_info_masterlegend_4monocle.csv", 
#         header =T)
#samplesheet[is.na(samplesheet)] <- 0
#expr_matrix = read.csv("/Users/festeneam/Artikelen/ScRNAseq/run02/allcellsallptsextra/Monocle/allcellallpts_basematrix.csv")
#expr_matrix[is.na(expr_matrix)] <- 0
#exp_matrix = as.matrix(expr_matrix)

#pd <- new("AnnotatedDataFrame", data = samplesheet)
## Rownames become nrs, convert samplenames to rownames to avoid error
#rownames(pd) <- pd$sampleNames
## Create featuredata to avoid errors
#fd <- rownames(expr_matrix)
#fd <- as.data.frame(fd)
#fd$gene_short_name <- fd$fd
#rownames(fd) <- fd$fd
#fd <- new("AnnotatedDataFrame", data = fd)

##fd <- new("AnnotatedDataFrame", data = gene_annotation)
#cds <- newCellDataSet(as(as.matrix(exp_matrix), "sparseMatrix"), 
#                      phenoData = pd, 
#                      featureData = fd,
#                      lowerDetectionLimit = 0.1,
#                      expressionFamily=negbinomial())

## Load Seuratfiles and convert them into Monocle
load("/Users/festeneam/Artikelen/ScRNAseq/FinalAnalysisSet/seuratfile_bloodcells.Rda")
load("/Users/festeneam/Artikelen/ScRNAseq/FinalAnalysisSet/seuratfile_mucosacells.Rda")
load("/Users/festeneam/Artikelen/ScRNAseq/FinalAnalysisSet/seuratfile_allcells.Rda")
allpts_blood_CDS <- importCDS(seuratfile_bloodcells, import_all = T)
allpts_mucosa_CDS <- importCDS(seuratfile_mucosacells, import_all = T)
allcells_allpts_CDS <- importCDS(seuratfile_allcells, import_all = T)

## Size factors help to normalize differences in mRNA recovered across cells
## Dispersion values help to perform differential expression analysis
allpts_blood_CDS <- estimateSizeFactors(allpts_blood_CDS)
allpts_blood_CDS <- estimateDispersions(allpts_blood_CDS)
allpts_mucosa_CDS <- estimateSizeFactors(allpts_mucosa_CDS)
allpts_mucosa_CDS <- estimateDispersions(allpts_mucosa_CDS)
allcells_allpts_CDS <- estimateSizeFactors(allcells_allpts_CDS)
allcells_allpts_CDS <- estimateDispersions(allcells_allpts_CDS)

# Quick QC filtering out genes expressed in less than 5 cells or lower than 0.1
allpts_blood_CDS <- detectGenes(allpts_blood_CDS, min_expr = 0.1)
print(head(fData(allpts_blood_CDS)))
expressed_genes <- row.names(subset(fData(allpts_blood_CDS),
                                    num_cells_expressed >= 5))

allpts_mucosa_CDS <- detectGenes(allpts_mucosa_CDS, min_expr = 0.1)
print(head(fData(allpts_mucosa_CDS)))
expressed_genes <- row.names(subset(fData(allpts_mucosa_CDS),
                                    num_cells_expressed >= 5))

allcells_allpts_CDS <- detectGenes(allcells_allpts_CDS, min_expr = 0.1)
print(head(fData(allcells_allpts_CDS)))
expressed_genes <- row.names(subset(fData(allcells_allpts_CDS),
                                    num_cells_expressed >= 5))

## Classify cells in blood set
## Make CellTypeHierarchy object
Th17_id1 <- row.names(subset(fData(allpts_blood_CDS), gene_short_name == c("IL17A")))
Th17_id2 <- row.names(subset(fData(allpts_blood_CDS),gene_short_name == c("IL17F")))
Th17_id3 <- row.names(subset(fData(allpts_blood_CDS), gene_short_name == c("IL22")))
Treg_id1 <- row.names(subset(fData(allpts_blood_CDS), gene_short_name == c("TNFRSF18")))
Treg_id2 <- row.names(subset(fData(allpts_blood_CDS), gene_short_name == c("TNFRSF4")))
Treg_id3 <- row.names(subset(fData(allpts_blood_CDS), gene_short_name == c("FOXP3")))
Quies_id1 <- row.names(subset(fData(allpts_blood_CDS), gene_short_name == c("NOG")))
Quies_id2 <- row.names(subset(fData(allpts_blood_CDS), gene_short_name == c("CCR7")))
Th1_id <- row.names(subset(fData(allpts_blood_CDS), gene_short_name == c("IFNG")))
Cytotox_id1 <- row.names(subset(fData(allpts_blood_CDS), gene_short_name == c("NKG7")))
Cytotox_id2 <- row.names(subset(fData(allpts_blood_CDS), gene_short_name == c("GNLY")))

## Define the celltype calling tree for blood
cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Th17", classify_func = function(x) { x[Th17_id1,] >= 1 | x[Th17_id2,] >= 1})
cth <- addCellType(cth, "Treg", classify_func = function(x)
{ x[Th17_id1,] < 1 & x[Th17_id2,] < 1 & x[Treg_id1,] > 1 | x[Treg_id2,] > 1 | x[Treg_id3,] > 1 })
cth <- addCellType(cth, "Th1", classify_func = function(x)
{ x[Th17_id1,] < 1 & x[Th17_id2,] < 1 & x[Treg_id1,] < 1 & x[Treg_id2,] < 1 & 
    x[Treg_id3,] < 1 & x[Th1_id,] > 1 })
cth <- addCellType(cth, "Quies", classify_func = function(x)
{ x[Th17_id1,] < 1 & x[Th17_id2,] < 1 & x[Treg_id1,] < 1 & x[Treg_id2,] < 1 & 
    x[Treg_id3,] < 1 & x[Th1_id,] < 1 & x[Quies_id1,] > 1 | x[Quies_id2,] > 1 })
cth <- addCellType(cth, "Cytotox", classify_func = function(x)
{ x[Th17_id1,] < 1 & x[Th17_id2,] < 1 & x[Treg_id1,] < 1 & x[Treg_id2,] < 1 & 
    x[Treg_id3,] < 1 & x[Th1_id,] < 1 & x[Quies_id1,] < 1 & x[Quies_id2,] < 1 &
    x[Cytotox_id1,] > 1 | x[Cytotox_id2,] > 1 })

## Use CellTypeHierarchy object, to classify all the cells in the blood set:
allpts_blood_CDS <- classifyCells(allpts_blood_CDS, cth, 0.1)



### Clustering cells using marker genes
## Pick genes that co-vary with marker genes. In a sense, building a large list of genes to use 
## as markers, so that even if a cell doesn't have the main marker, it might be recognizable
marker_diff <- markerDiffTable(allpts_blood_CDS[expressed_genes,],
                               cth,
                               residualModelFormulaStr = "~patient + num_genes_expressed",
                               cores = 1)
## Identifies genes that are differentially expressed between the cell types (with provisional 
## residual model of effects to exclude from this test). Returns a data frame of test results, 
## to use to pick the genes for clustering
## Pick the top 10 or 20 genes that are most specific for each cell type. This ensures that the 
## clustering genes aren't dominated by markers for one cell type:
## Ranking genes by how restricted expression is for cell type - score between 0-1, 1 most specific
candidate_clustering_genes <-
  row.names(subset(marker_diff, qval < 0.01))
marker_spec <-
  calculateMarkerSpecificity(allpts_blood_CDS[candidate_clustering_genes,], cth)
head(selectTopMarkers(marker_spec, 3))
## To cluster the cells, choose the top 200 markers for each of the cell types
semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
allpts_blood_CDS <- setOrderingFilter(allpts_blood_CDS, semisup_clustering_genes)
plot_ordering_genes(allpts_blood_CDS)

plot_pc_variance_explained(allpts_blood_CDS, return_all = F)

allpts_blood_CDS <- reduceDimension(allpts_blood_CDS, max_components = 2, num_dim = 40,
                        norm_method = 'log',
                        reduction_method = 'tSNE',
                        residualModelFormulaStr = "~patient + num_genes_expressed",
                        verbose = T)
allpts_blood_CDS <- clusterCells(allpts_blood_CDS, num_clusters = 2)
plot_cell_clusters(allpts_blood_CDS, 1, 2, color = "CellType")

## Imputing cell type
allpts_blood_CDS <- clusterCells(allpts_blood_CDS,
                     num_clusters = 8,
                     frequency_thresh = 0.1,
                     cell_type_hierarchy = cth)
plot_cell_clusters(allpts_blood_CDS, 1, 2, color = "CellType",
                   markers = c("CCR7", "NKG7"))

### MUCOSA
## Classify cells in mucosa set (finish blood first!)
## Make CellTypeHierarchy object
Th17_id1 <- row.names(subset(fData(allpts_mucosa_CDS), gene_short_name == c("IL17A")))
Th17_id2 <- row.names(subset(fData(allpts_mucosa_CDS),gene_short_name == c("IL17F")))
Th17_id3 <- row.names(subset(fData(allpts_mucosa_CDS), gene_short_name == c("IL22")))
Treg_id1 <- row.names(subset(fData(allpts_mucosa_CDS), gene_short_name == c("TNFRSF18")))
Treg_id2 <- row.names(subset(fData(allpts_mucosa_CDS), gene_short_name == c("TNFRSF4")))
Treg_id3 <- row.names(subset(fData(allpts_mucosa_CDS), gene_short_name == c("FOXP3")))
Quies_id1 <- row.names(subset(fData(allpts_mucosa_CDS), gene_short_name == c("NOG")))
Quies_id2 <- row.names(subset(fData(allpts_mucosa_CDS), gene_short_name == c("CCR7")))
Th1_id <- row.names(subset(fData(allpts_mucosa_CDS), gene_short_name == c("IFNG")))
Cytotox_id1 <- row.names(subset(fData(allpts_mucosa_CDS), gene_short_name == c("NKG7")))
Cytotox_id2 <- row.names(subset(fData(allpts_mucosa_CDS), gene_short_name == c("GNLY")))

## Define the celltype calling tree for mucosa
cthm <- newCellTypeHierarchy()
cthm <- addCellType(cthm, "Th17", classify_func = function(x) { x[Th17_id1,] >= 1 | x[Th17_id2,] >= 1})
cthm <- addCellType(cthm, "Treg", classify_func = function(x)
{ x[Th17_id1,] < 1 & x[Th17_id2,] < 1 & x[Treg_id1,] > 1 | x[Treg_id2,] > 1 | x[Treg_id3,] > 1 })
cthm <- addCellType(cthm, "Th1", classify_func = function(x)
{ x[Th17_id1,] < 1 & x[Th17_id2,] < 1 & x[Treg_id1,] < 1 & x[Treg_id2,] < 1 & 
    x[Treg_id3,] < 1 & x[Th1_id,] > 1 })
cthm <- addCellType(cthm, "Quies", classify_func = function(x)
{ x[Th17_id1,] < 1 & x[Th17_id2,] < 1 & x[Treg_id1,] < 1 & x[Treg_id2,] < 1 & 
    x[Treg_id3,] < 1 & x[Th1_id,] < 1 & x[Quies_id1,] > 1 | x[Quies_id2,] > 1 })
cthm <- addCellType(cthm, "Cytotox", classify_func = function(x)
{ x[Th17_id1,] < 1 & x[Th17_id2,] < 1 & x[Treg_id1,] < 1 & x[Treg_id2,] < 1 & 
    x[Treg_id3,] < 1 & x[Th1_id,] < 1 & x[Quies_id1,] < 1 & x[Quies_id2,] < 1 &
    x[Cytotox_id1,] > 1 | x[Cytotox_id2,] > 1 })

## Use CellTypeHierarchy object, to classify all the cells in the mucosa set:
allpts_mucosa_CDS <- classifyCells(allpts_mucosa_CDS, cthm, 0.1)

### Clustering mucosa cells using marker genes
## Pick genes that co-vary with marker genes. In a sense, building a large list of genes to use 
## as markers, so that even if a cell doesn't have the main marker, it might be recognizable
marker_diff <- markerDiffTable(allpts_mucosa_CDS[expressed_genes,],
                               cthm,
                               residualModelFormulaStr = "~patient + num_genes_expressed",
                               cores = 1)
## Identifies genes that are differentially expressed between the cell types (with provisional 
## residual model of effects to exclude from this test). Returns a data frame of test results, 
## to use to pick the genes for clustering
## Pick the top 10 or 20 genes that are most specific for each cell type. This ensures that the 
## clustering genes aren't dominated by markers for one cell type:
## Ranking genes by how restricted expression is for cell type - score between 0-1, 1 most specific
candidate_clustering_genes <-
  row.names(subset(marker_diff, qval < 0.01))
marker_spec_m <-
  calculateMarkerSpecificity(allpts_mucosa_CDS[candidate_clustering_genes,], cthm)
head(selectTopMarkers(marker_spec_m, 3))
## To cluster the cells, choose the top 200-500 markers for each of the cell types
semisup_clustering_genes <- unique(selectTopMarkers(marker_spec_m, 500)$gene_id)
allpts_mucosa_CDS <- setOrderingFilter(allpts_mucosa_CDS, semisup_clustering_genes)
plot_ordering_genes(allpts_mucosa_CDS)

plot_pc_variance_explained(allpts_mucosa_CDS, return_all = F)

allpts_mucosa_CDS <- reduceDimension(allpts_mucosa_CDS, max_components = 2, num_dim = 20,
                                    norm_method = 'log',
                                    reduction_method = 'tSNE',
                                    #residualModelFormulaStr = "~patient + num_genes_expressed",
                                    verbose = T)
allpts_mucosa_CDS <- clusterCells(allpts_mucosa_CDS, num_clusters = 5)
plot_cell_clusters(allpts_mucosa_CDS, 1, 2, color = "CellType")

## Imputing cell type
allpts_mucosa_CDS <- clusterCells(allpts_mucosa_CDS,
                                 num_clusters = 8,
                                 frequency_thresh = 0.1,
                                 cell_type_hierarchy = cthm)
plot_cell_clusters(allpts_mucosa_CDS, 1, 2, color = "CellType",
                   markers = c("CCR7", "NKG7", "IL17A", "TNFRSF4"))

## Classify cells in allcells set
## Make CellTypeHierarchy object
Th17_id1 <- row.names(subset(fData(allcells_allpts_CDS), gene_short_name == c("IL17A")))
Th17_id2 <- row.names(subset(fData(allcells_allpts_CDS),gene_short_name == c("IL17F")))
Th17_id3 <- row.names(subset(fData(allcells_allpts_CDS), gene_short_name == c("IL22")))
Treg_id1 <- row.names(subset(fData(allcells_allpts_CDS), gene_short_name == c("TNFRSF18")))
Treg_id2 <- row.names(subset(fData(allcells_allpts_CDS), gene_short_name == c("TNFRSF4")))
Treg_id3 <- row.names(subset(fData(allcells_allpts_CDS), gene_short_name == c("FOXP3")))
Quies_id1 <- row.names(subset(fData(allcells_allpts_CDS), gene_short_name == c("NOG")))
Quies_id2 <- row.names(subset(fData(allcells_allpts_CDS), gene_short_name == c("CCR7")))
Th1_id <- row.names(subset(fData(allcells_allpts_CDS), gene_short_name == c("IFNG")))
Cytotox_id1 <- row.names(subset(fData(allcells_allpts_CDS), gene_short_name == c("NKG7")))
Cytotox_id2 <- row.names(subset(fData(allcells_allpts_CDS), gene_short_name == c("GNLY")))

## Define the celltype calling tree for blood
ctha <- newCellTypeHierarchy()
ctha <- addCellType(ctha, "Th17", classify_func = function(x) { x[Th17_id1,] >= 1 | x[Th17_id2,] >= 1})
ctha <- addCellType(ctha, "Treg", classify_func = function(x)
{ x[Th17_id1,] < 1 & x[Th17_id2,] < 1 & x[Treg_id1,] > 1 | x[Treg_id2,] > 1 | x[Treg_id3,] > 1 })
ctha <- addCellType(ctha, "Th1", classify_func = function(x)
{ x[Th17_id1,] < 1 & x[Th17_id2,] < 1 & x[Treg_id1,] < 1 & x[Treg_id2,] < 1 & 
    x[Treg_id3,] < 1 & x[Th1_id,] > 1 })
ctha <- addCellType(ctha, "Quies", classify_func = function(x)
{ x[Th17_id1,] < 1 & x[Th17_id2,] < 1 & x[Treg_id1,] < 1 & x[Treg_id2,] < 1 & 
    x[Treg_id3,] < 1 & x[Th1_id,] < 1 & x[Quies_id1,] > 1 | x[Quies_id2,] > 1 })
ctha <- addCellType(ctha, "Cytotox", classify_func = function(x)
{ x[Th17_id1,] < 1 & x[Th17_id2,] < 1 & x[Treg_id1,] < 1 & x[Treg_id2,] < 1 & 
    x[Treg_id3,] < 1 & x[Th1_id,] < 1 & x[Quies_id1,] < 1 & x[Quies_id2,] < 1 &
    x[Cytotox_id1,] > 1 | x[Cytotox_id2,] > 1 })

## Use CellTypeHierarchy object, to classify all the cells in the blood set:
allcells_allpts_CDS <- classifyCells(allcells_allpts_CDS, ctha, 0.1)



### Clustering cells using marker genes
## Pick genes that co-vary with marker genes. In a sense, building a large list of genes to use 
## as markers, so that even if a cell doesn't have the main marker, it might be recognizable
marker_diff <- markerDiffTable(allcells_allpts_CDS[expressed_genes,],
                               ctha,
                               residualModelFormulaStr = "~patient + num_genes_expressed",
                               cores = 1)
## Identifies genes that are differentially expressed between the cell types (with provisional 
## residual model of effects to exclude from this test). Returns a data frame of test results, 
## to use to pick the genes for clustering
## Pick the top 10 or 20 genes that are most specific for each cell type. This ensures that the 
## clustering genes aren't dominated by markers for one cell type:
## Ranking genes by how restricted expression is for cell type - score between 0-1, 1 most specific
candidate_clustering_genes <-
  row.names(subset(marker_diff, qval < 0.01))
marker_spec <-
  calculateMarkerSpecificity(allcells_allpts_CDS[candidate_clustering_genes,], ctha)
head(selectTopMarkers(marker_spec, 3))
## To cluster the cells, choose the top 200 markers for each of the cell types
semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
allcells_allpts_CDS <- setOrderingFilter(allcells_allpts_CDS, semisup_clustering_genes)
plot_ordering_genes(allcells_allpts_CDS)

plot_pc_variance_explained(allcells_allpts_CDS, return_all = F)

allcells_allpts_CDS <- reduceDimension(allcells_allpts_CDS, max_components = 2, num_dim = 40,
                                    norm_method = 'log',
                                    reduction_method = 'tSNE',
                                    #residualModelFormulaStr = "~patient + num_genes_expressed",
                                    verbose = T)
allcells_allpts_CDS <- clusterCells(allcells_allpts_CDS, num_clusters = 2)
plot_cell_clusters(allcells_allpts_CDS, 1, 2, color = "CellType")

## Imputing cell type
allcells_allpts_CDS <- clusterCells(allcells_allpts_CDS,
                                 num_clusters = 8,
                                 frequency_thresh = 0.1,
                                 cell_type_hierarchy = cth)
plot_cell_clusters(allcells_allpts_CDS, 1, 2, color = "CellType",
                   markers = c("CCR7", "NKG7"))

