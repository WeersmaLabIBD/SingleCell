
library(Seurat)
library(ggplot2)

options(future.globals.maxSize = 150000 * 1024^2)



HC_3296<-readRDS("/methodspaper_data/HC-NI-3296_sct.rds")
DefaultAssay(HC_3296)<-"SCT"
HC_3002<-readRDS("/methodspaper_data/HC-NI-3002_sct.rds")  
DefaultAssay(HC_3002)<-"SCT"
HC_3030<-readRDS("/methodspaper_data/HC-NI-3030_sct.rds")
DefaultAssay(HC_3030)<-"SCT"
HC_3049<-readRDS("/methodspaper_data/HC-NI-3049_sct.rds")
DefaultAssay(HC_3049)<-"SCT"
HC_3083<-readRDS("/methodspaper_data/HC-NI-3083_sct.rds")
DefaultAssay(HC_3083)<-"SCT"
HC_3037<-readRDS("/methodspaper_data/D-3037-HC_sct.rds")
DefaultAssay(HC_3037)<-"SCT"
HC_3034<-readRDS("/methodspaper_data/C-3034-HC_sct.rds")
DefaultAssay(HC_3034)<-"SCT"

HC_1<-readRDS("/methodspaper_data/HC_Sanger1_dataset_sct.rds")
HC_2<-readRDS("/methodspaper_data/HC_Sanger2_dataset_sct.rds")
HC_3<-readRDS("/methodspaper_data/HC_Sanger3_dataset_sct.rds")
HC_4<-readRDS("/methodspaper_data/HC_Sanger4_dataset_sct.rds")
HC_5<-readRDS("/methodspaper_data/HC_Sanger5_dataset_sct.rds")
HC_6<-readRDS("/methodspaper_data/HC_Sanger6_dataset_sct.rds")
HC_7<-readRDS("/methodspaper_data/HC_Sanger7_dataset_sct.rds")
HC_8<-readRDS("/methodspaper_data/HC_Sanger8_dataset_sct.rds")
HC_9<-readRDS("/methodspaper_data/HC_Sanger9_dataset_sct.rds")
DefaultAssay(HC_1)<-"SCT"
DefaultAssay(HC_2)<-"SCT"
DefaultAssay(HC_3)<-"SCT"
DefaultAssay(HC_4)<-"SCT"
DefaultAssay(HC_5)<-"SCT"
DefaultAssay(HC_6)<-"SCT"
DefaultAssay(HC_7)<-"SCT"
DefaultAssay(HC_8)<-"SCT"
DefaultAssay(HC_9)<-"SCT"
imm_healthy<-readRDS("/imm_healthy.rds")
epi_healthy<-readRDS("/epi_healthy.rds")
fib_healthy<-readRDS("/fib_healthy.rds")

HC_3296<-RenameCells(HC_3296, add.cell.id="HC3296")
HC_3002<-RenameCells(HC_3002, add.cell.id="HC3002")

HC_3030<-RenameCells(HC_3030, add.cell.id="HC3030")
HC_3049<-RenameCells(HC_3049, add.cell.id="HC3049")
HC_3083<-RenameCells(HC_3083, add.cell.id="HC3083")

HC_3034<-RenameCells(HC_3034, add.cell.id="HC3034")
HC_3037<-RenameCells(HC_3037, add.cell.id="HC3037")

HC_1<-RenameCells(HC_1, add.cell.id="HC1")
HC_2<-RenameCells(HC_2, add.cell.id="HC2")
HC_3<-RenameCells(HC_3, add.cell.id="HC3")
HC_4<-RenameCells(HC_4, add.cell.id="HC4")
HC_5<-RenameCells(HC_5, add.cell.id="HC5")
HC_6<-RenameCells(HC_6, add.cell.id="HC6")
HC_7<-RenameCells(HC_7, add.cell.id="HC7")
HC_8<-RenameCells(HC_8, add.cell.id="HC8")
HC_9<-RenameCells(HC_9, add.cell.id="HC9")

imm_healthy@meta.data$dataset<-"imm_healthy"
epi_healthy@meta.data$dataset<-"epi_healthy"
fib_healthy@meta.data$dataset<-"fib_healthy"

HC_3296@meta.data$dataset<-"HC_3296"
HC_3002@meta.data$dataset<-"HC_3002"

HC_3030@meta.data$dataset<-"HC_3030"
HC_3049@meta.data$dataset<-"HC_3049"
HC_3083@meta.data$dataset<-"HC_3083"
HC_3037@meta.data$dataset<-"HC_3037"
HC_3034@meta.data$dataset<-"HC_3034"

HC_1@meta.data$dataset<-"HC_1"
HC_2@meta.data$dataset<-"HC_2"
HC_3@meta.data$dataset<-"HC_3"
HC_4@meta.data$dataset<-"HC_4"
HC_5@meta.data$dataset<-"HC_5"
HC_6@meta.data$dataset<-"HC_6"
HC_7@meta.data$dataset<-"HC_7"
HC_8@meta.data$dataset<-"HC_8"
HC_9@meta.data$dataset<-"HC_9"

imm_healthy@meta.data$method<-"splitcollagenase"
epi_healthy@meta.data$method<-"splitcollagenase"
fib_healthy@meta.data$method<-"splitcollagenase"

HC_3296@meta.data$method<-"wholecollagenase"
HC_3002@meta.data$method<-"wholecollagenase"

HC_3030@meta.data$method<-"wholecollagenase"
HC_3049@meta.data$method<-"wholecollagenase"
HC_3083@meta.data$method<-"wholecollagenase"
HC_3037@meta.data$method<-"wholecollagenase"
HC_3034@meta.data$method<-"wholecollagenase"

HC_1@meta.data$method<-"splitprotease"
HC_2@meta.data$method<-"splitprotease"
HC_3@meta.data$method<-"splitprotease"
HC_4@meta.data$method<-"splitprotease"
HC_5@meta.data$method<-"splitprotease"
HC_6@meta.data$method<-"splitprotease"
HC_7@meta.data$method<-"splitprotease"
HC_8@meta.data$method<-"splitprotease"
HC_9@meta.data$method<-"splitprotease"


imm_healthy_sub[["percent.ribo"]] <- PercentageFeatureSet(imm_healthy_sub, pattern = "^RPL|^RPS")
epi_healthy_sub[["percent.ribo"]] <- PercentageFeatureSet(epi_healthy_sub, pattern = "^RPL|^RPS")
fib_healthy_sub[["percent.ribo"]] <- PercentageFeatureSet(fib_healthy_sub, pattern = "^RPL|^RPS")

imm_healthy[["percent.ribo"]] <- PercentageFeatureSet(imm_healthy, pattern = "^RPL|^RPS")
epi_healthy[["percent.ribo"]] <- PercentageFeatureSet(epi_healthy, pattern = "^RPL|^RPS")
fib_healthy[["percent.ribo"]] <- PercentageFeatureSet(fib_healthy, pattern = "^RPL|^RPS")

HC_2[["percent.ribo"]] <- PercentageFeatureSet(HC_2, pattern = "^RPL|^RPS")
HC_6[["percent.ribo"]] <- PercentageFeatureSet(HC_6, pattern = "^RPL|^RPS")
HC_9[["percent.ribo"]] <- PercentageFeatureSet(HC_9, pattern = "^RPL|^RPS")

HC_3296[["percent.ribo"]] <- PercentageFeatureSet(HC_3296, pattern = "^RPL|^RPS")
HC_3030[["percent.ribo"]] <- PercentageFeatureSet(HC_3030, pattern = "^RPL|^RPS")
HC_3083[["percent.ribo"]] <- PercentageFeatureSet(HC_3083, pattern = "^RPL|^RPS")
HC_3037[["percent.ribo"]] <- PercentageFeatureSet(HC_3037, pattern = "^RPL|^RPS")
HC_3034[["percent.ribo"]] <- PercentageFeatureSet(HC_3034, pattern = "^RPL|^RPS")


alldata_202002_integration_list <- list(HC_3296, HC_3030, HC_3083, HC_3034, HC_3037,HC_2, HC_6,HC_9, epi_healthy_sub,imm_healthy_sub,fib_healthy_sub)

for (i in 1:length(alldata_202002_integration_list)){
alldata_202002_integration_list[[i]] <- SCTransform(alldata_202002_integration_list[[i]],vars.to.regress = c("percent.mt", "percent.ribo"), verbose = FALSE)
}

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


for (i in 1:nrow(methodspaper@meta.data)){ 
  if (methodspaper@meta.data[i,21] == "unassigned"){
    methodspaper@meta.data[i,23]<-methodspaper@meta.data[i,22] 
  } else {
    methodspaper@meta.data[i,23]<-methodspaper@meta.data[i,21]
  }
}

