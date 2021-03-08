library(Seurat)
library(ggplot2)

options(future.globals.maxSize = 150000 * 1024^2)

lane18_1<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_new/200618_lane1_sct.rds")
DefaultAssay(lane18_1)<-"SCT"
lane18_2<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_new/200618_lane2_sct.rds")  
DefaultAssay(lane18_2)<-"SCT"
lane18_3<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_new/200618_lane3_sct.rds")
DefaultAssay(lane18_3)<-"SCT"
lane18_4<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_new/200618_lane4_sct.rds")
DefaultAssay(lane18_4)<-"SCT"

lane18_1[["percent.ribo"]] <- PercentageFeatureSet(lane18_1, pattern = "^RPL|^RPS")
lane18_2[["percent.ribo"]] <- PercentageFeatureSet(lane18_2, pattern = "^RPL|^RPS")
lane18_3[["percent.ribo"]] <- PercentageFeatureSet(lane18_3, pattern = "^RPL|^RPS")
lane18_4[["percent.ribo"]] <- PercentageFeatureSet(lane18_4, pattern = "^RPL|^RPS")

alldata_vedo2_batch1_PBMC_integration_list <- list(lane18_1, lane18_2, lane18_3, lane18_4)

print("loaded")
features <- SelectIntegrationFeatures(object.list = alldata_vedo2_batch1_PBMC_integration_list, nfeatures = 3000)
print("features_selected")
alldata_vedo2_batch1_PBMC_integration_list <- PrepSCTIntegration(object.list = alldata_vedo2_batch1_PBMC_integration_list, anchor.features = features)
print("integration_prepped")
anchors <- FindIntegrationAnchors(object.list = alldata_vedo2_batch1_PBMC_integration_list, anchor.features = features, normalization.method = "SCT")
print("achors_found")
rm(alldata_vedo2_batch1_PBMC_integration_list)
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

saveRDS(alldata.integrated, file="/groups/umcg-weersma/tmp01/Emilia/batch1_integration/vedo2_batch1_PBMC_integrated_noribo_withmito.rds")
