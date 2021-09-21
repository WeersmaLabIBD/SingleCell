library(Seurat)
library(ggplot2)

options(future.globals.maxSize = 150000 * 1024^2)

lane11_1<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_nomito/200611_lane1_nomito_sct.rds")
DefaultAssay(lane11_1)<-"SCT"
lane11_2<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_nomito/200611_lane2_nomito_sct.rds")  
DefaultAssay(lane11_2)<-"SCT"
lane11_3<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_nomito/200611_lane3_nomito_sct.rds")
DefaultAssay(lane11_3)<-"SCT"
lane11_4<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_nomito/200611_lane4_nomito_sct.rds")
DefaultAssay(lane11_4)<-"SCT"
lane12_1<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_nomito/200612_lane1_nomito_sct.rds")
DefaultAssay(lane12_1)<-"SCT"
lane12_2<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_nomito/200612_lane2_nomito_sct.rds")
DefaultAssay(lane12_2)<-"SCT"
lane12_3<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_nomito/200612_lane3_nomito_sct.rds")
DefaultAssay(lane12_3)<-"SCT"
lane12_4<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_nomito/200612_lane4_nomito_sct.rds")
DefaultAssay(lane12_4)<-"SCT"
lane25_1<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_nomito/200625_lane1_nomito_sct.rds")
DefaultAssay(lane25_1)<-"SCT"
lane25_2<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_nomito/200625_lane2_nomito_sct.rds")
DefaultAssay(lane25_2)<-"SCT"
lane25_3<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_nomito/200625_lane3_nomito_sct.rds")
DefaultAssay(lane25_3)<-"SCT"
lane25_4<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_nomito/200625_lane4_nomito_sct.rds")
DefaultAssay(lane25_4)<-"SCT"
lane26_1<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_nomito/200626_lane1_nomito_sct.rds")
DefaultAssay(lane26_1)<-"SCT"
lane26_2<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_nomito/200626_lane2_nomito_sct.rds")
DefaultAssay(lane26_2)<-"SCT"
lane26_3<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_nomito/200626_lane3_nomito_sct.rds")
DefaultAssay(lane26_3)<-"SCT"
lane26_4<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/preprocessed_lanes_nomito/200626_lane4_nomito_sct.rds")
DefaultAssay(lane26_4)<-"SCT"
Vedo01_T0_NI<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/VEDO01_T0_NI_sct.rds")
DefaultAssay(Vedo01_T0_NI)<-"SCT"
Vedo01_T0_I<-readRDS("/groups/umcg-weersma/tmp01/Emilia/batch1_integration/VEDO01_T0_I_sct.rds")
DefaultAssay(Vedo01_T0_I)<-"SCT"

lane11_1[["percent.ribo"]] <- PercentageFeatureSet(lane11_1, pattern = "^RPL|^RPS")
lane11_2[["percent.ribo"]] <- PercentageFeatureSet(lane11_2, pattern = "^RPL|^RPS")
lane11_3[["percent.ribo"]] <- PercentageFeatureSet(lane11_3, pattern = "^RPL|^RPS")
lane11_4[["percent.ribo"]] <- PercentageFeatureSet(lane11_4, pattern = "^RPL|^RPS")
lane12_1[["percent.ribo"]] <- PercentageFeatureSet(lane12_1, pattern = "^RPL|^RPS")
lane12_2[["percent.ribo"]] <- PercentageFeatureSet(lane12_2, pattern = "^RPL|^RPS")
lane12_3[["percent.ribo"]] <- PercentageFeatureSet(lane12_3, pattern = "^RPL|^RPS")
lane12_4[["percent.ribo"]] <- PercentageFeatureSet(lane12_4, pattern = "^RPL|^RPS")
lane25_1[["percent.ribo"]] <- PercentageFeatureSet(lane25_1, pattern = "^RPL|^RPS")
lane25_2[["percent.ribo"]] <- PercentageFeatureSet(lane25_2, pattern = "^RPL|^RPS")
lane25_3[["percent.ribo"]] <- PercentageFeatureSet(lane25_3, pattern = "^RPL|^RPS")
lane25_4[["percent.ribo"]] <- PercentageFeatureSet(lane25_4, pattern = "^RPL|^RPS")
lane26_1[["percent.ribo"]] <- PercentageFeatureSet(lane26_1, pattern = "^RPL|^RPS")
lane26_2[["percent.ribo"]] <- PercentageFeatureSet(lane26_2, pattern = "^RPL|^RPS")
lane26_3[["percent.ribo"]] <- PercentageFeatureSet(lane26_3, pattern = "^RPL|^RPS")
lane26_4[["percent.ribo"]] <- PercentageFeatureSet(lane26_4, pattern = "^RPL|^RPS")
Vedo01_T0_NI[["percent.ribo"]] <- PercentageFeatureSet(Vedo01_T0_NI, pattern = "^RPL|^RPS")
Vedo01_T0_I[["percent.ribo"]] <- PercentageFeatureSet(Vedo01_T0_I, pattern = "^RPL|^RPS")

alldata_vedo2_batch1_biop_integration_list <- list(lane11_1, lane11_2, lane11_3, lane11_4, lane12_1, lane12_2, lane12_3, lane12_4, lane25_1, lane25_2, lane25_3, lane25_4, lane26_1, lane26_2, lane26_3, lane26_4, Vedo01_T0_NI, Vedo01_T0_I)

print("loaded")
features <- SelectIntegrationFeatures(object.list = alldata_vedo2_batch1_biop_integration_list, nfeatures = 3000)
print("features_selected")
alldata_vedo2_batch1_biop_integration_list <- PrepSCTIntegration(object.list = alldata_vedo2_batch1_biop_integration_list, anchor.features = features)
print("integration_prepped")
anchors <- FindIntegrationAnchors(object.list = alldata_vedo2_batch1_biop_integration_list, anchor.features = features, normalization.method = "SCT")

#when using reference lanes
#anchors <- FindIntegrationAnchors(object.list = alldata_vedo2_batch1_biop_integration_list, anchor.features = features, normalization.method = "SCT", reference = c(1,4,7,8,16))

print("achors_found")
rm(alldata_vedo2_batch1_biop_integration_list)
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

saveRDS(alldata.integrated, file="/groups/umcg-weersma/tmp01/Emilia/batch1_integration/vedo2_batch1_biop_integrated_noribo_nonorm_withmito.rds")
