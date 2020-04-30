## Integrating datasets using SCTransform integration, Seurat V3 integration and harmony

### SCT integration
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
