## Integrating datasets using SCTransform integration, Seurat V3 integration and harmony
```
# first: load all QC-ed, SCTransformed and filtered datasets
```
### SCT integration
testing SCT integration, takes a long time and a lot of computational power to compute
```
DefaultAssay(Day6_lane2)<-"SCT" # set all default assays to 'SCT'
alldata_202002_integration_list <- list(Day1_lane1, Day1_lane2, Day2_lane1, Day2_lane2, Day3_lane1, Day3_lane2,Day4_lane1, Day4_lane2,Day5_lane1, Day5_lane2,Day6_lane1, Day6_lane2)

print("loaded")
features <- SelectIntegrationFeatures(object.list = alldata_202002_integration_list, nfeatures = 3000)
print("features_selected")
alldata_202002_integration_list <- PrepSCTIntegration(object.list = alldata_202002_integration_list, anchor.features = features)
print("integration_prepped")
anchors <- FindIntegrationAnchors(object.list = alldata_202002_integration_list, anchor.features = features, normalization.method = "SCT")
print("achors_found")
rm(alldata_202002_integration_list)
alldata.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
print("integrated")
rm(anchors)
alldata.integrated <- RunPCA(alldata.integrated, verbose = FALSE)
print("pca2_ran")
alldata.integrated <- RunUMAP(alldata.integrated, dims = 1:30)
print("umap_ran")
alldata.integrated <- FindNeighbors(alldata.integrated, dims = 1:30)
alldata.integrated <- FindClusters(alldata.integrated, resolution = 0.5)

DefaultAssay(alldata.integrated) <- "RNA"
alldata.integrated <- NormalizeData(alldata.integrated, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
saveRDS(alldata.integrated, file="/data_complete_sct_integrated.rds")

```
### Seurat integration
testing Seurat v3 integration
```
options(future.globals.maxSize = 300000 * 1024^2) # set to max necessary
alldata_202002_integration_anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
alldata.integrated <- IntegrateData(anchorset = alldata_202002_integration_anchors, dims = 1:30)
DefaultAssay(alldata.integrated) <- "integrated"
alldata.integrated <- ScaleData(alldata.integrated, verbose = FALSE)
alldata.integrated <- RunPCA(alldata.integrated, npcs = 30, verbose = FALSE)
alldata.integrated <- RunUMAP(alldata.integrated, reduction = "pca", dims = 1:30)
saveRDS(alldata.integrated, "/alldata_normal_integration.rds")
```
### harmony integration
super fast integration
```
# start with integrated dataset from above, or merge different datasets into one and start with that
T_cells_harmony <- RunHarmony(T_cells, "dataset", assay.use="SCT")
T_cells_harmony  <- RunUMAP(T_cells_harmony, reduction = "harmony", dims = 1:30)
DimPlot(T_cells_harmony, group.by="Sample")
DefaultAssay(T_cells_harmony)<-"SCT"
T_cells_harmony<-FindNeighbors(T_cells_harmony, reduction = "harmony", dims = 1:30)
T_cells_harmony<-FindClusters(T_cells_harmony,resolution = 0.5)
saveRDS(T_cells_harmony, "/T_cells_harmony_integration.rds")
```
