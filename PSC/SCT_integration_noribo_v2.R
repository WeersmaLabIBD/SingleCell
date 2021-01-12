library(Seurat)
library(ggplot2)

options(future.globals.maxSize = 150000 * 1024^2)

lane08_1 <-readRDS("/groups/umcg-weersma/tmp01/Amber/PSC_preprocessed_new/200108_lane1_sct.rds")
DefaultAssay(lane08_1)<-"SCT"
lane08_2<-readRDS("/groups/umcg-weersma/tmp01/Amber/PSC_preprocessed_new/200108_lane2_sct.rds")  
DefaultAssay(lane08_2)<-"SCT"
lane09_1<-readRDS("/groups/umcg-weersma/tmp01/Amber/PSC_preprocessed_new/200109_lane1_sct.rds")
DefaultAssay(lane09_1)<-"SCT"
lane09_2<-readRDS("/groups/umcg-weersma/tmp01/Amber/PSC_preprocessed_new/200109_lane2_sct.rds")
DefaultAssay(lane09_2)<-"SCT"
lane13_1<-readRDS("/groups/umcg-weersma/tmp01/Amber/PSC_preprocessed_new/200113_lane1_sct.rds")
DefaultAssay(lane13_1)<-"SCT"
lane13_2<-readRDS("/groups/umcg-weersma/tmp01/Amber/PSC_preprocessed_new/200113_lane2_sct.rds")
DefaultAssay(lane13_2)<-"SCT"
lane14_1<-readRDS("/groups/umcg-weersma/tmp01/Amber/PSC_preprocessed_new/200114_lane1_sct.rds")
DefaultAssay(lane14_1)<-"SCT"
lane14_2<-readRDS("/groups/umcg-weersma/tmp01/Amber/PSC_preprocessed_new/200114_lane2_sct.rds")
DefaultAssay(lane14_2)<-"SCT"
lane15_1<-readRDS("/groups/umcg-weersma/tmp01/Amber/PSC_preprocessed_new/200115_lane1_sct.rds")
DefaultAssay(lane15_1)<-"SCT"
lane15_2<-readRDS("/groups/umcg-weersma/tmp01/Amber/PSC_preprocessed_new/200115_lane2_sct.rds")
DefaultAssay(lane15_2)<-"SCT"
lane16_1<-readRDS("/groups/umcg-weersma/tmp01/Amber/PSC_preprocessed_new/200116_lane1_sct.rds")
DefaultAssay(lane16_1)<-"SCT"
lane16_2<-readRDS("/groups/umcg-weersma/tmp01/Amber/PSC_preprocessed_new/200116_lane2_sct.rds")
DefaultAssay(lane16_2)<-"SCT"

lane08_1[["percent.ribo"]] <- PercentageFeatureSet(lane08_1, pattern = "^RPL|^RPS")
lane08_2[["percent.ribo"]] <- PercentageFeatureSet(lane08_2, pattern = "^RPL|^RPS")
lane09_1[["percent.ribo"]] <- PercentageFeatureSet(lane09_1, pattern = "^RPL|^RPS")
lane09_2[["percent.ribo"]] <- PercentageFeatureSet(lane09_2, pattern = "^RPL|^RPS")
lane13_1[["percent.ribo"]] <- PercentageFeatureSet(lane13_1, pattern = "^RPL|^RPS")
lane13_2[["percent.ribo"]] <- PercentageFeatureSet(lane13_2, pattern = "^RPL|^RPS")
lane14_1[["percent.ribo"]] <- PercentageFeatureSet(lane14_1, pattern = "^RPL|^RPS")
lane14_2[["percent.ribo"]] <- PercentageFeatureSet(lane14_2, pattern = "^RPL|^RPS")
lane15_1[["percent.ribo"]] <- PercentageFeatureSet(lane15_1, pattern = "^RPL|^RPS")
lane15_2[["percent.ribo"]] <- PercentageFeatureSet(lane15_2, pattern = "^RPL|^RPS")
lane16_1[["percent.ribo"]] <- PercentageFeatureSet(lane16_1, pattern = "^RPL|^RPS")
lane16_2[["percent.ribo"]] <- PercentageFeatureSet(lane16_2, pattern = "^RPL|^RPS")

alldata_PSC_202002_integration_list <- list(lane08_1, lane08_2, lane09_1, lane09_2, lane13_1, lane13_2, lane14_1, lane14_2, lane15_1, lane15_2, lane16_1, lane16_2)

print("loaded")
features <- SelectIntegrationFeatures(object.list = alldata_PSC_202002_integration_list, nfeatures = 3000)
print("features_selected")
alldata_PSC_202002_integration_list <- PrepSCTIntegration(object.list = alldata_PSC_202002_integration_list, anchor.features = features)
print("integration_prepped")
anchors <- FindIntegrationAnchors(object.list = alldata_PSC_202002_integration_list, anchor.features = features, normalization.method = "SCT")
print("achors_found")
rm(alldata_PSC_202002_integration_list)
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

saveRDS(alldata.integrated, file="/groups/umcg-weersma/tmp01/Amber/PSC_202002_integrated_v2_noribo.rds")
