methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Best4+ Enterocytes"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Enterocytes"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Enterocyte Progenitors"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Immature Enterocytes 1"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Immature Enterocytes 2"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "TA 1"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "TA 2"]<-"epithelium"

# epi general
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "M cells"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Cycling TA"]<-"epithelium"

# epi secretory
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Goblet"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Enteroendocrine"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Immature Goblet"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Secretory TA"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Stem"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Tuft"]<-"epithelium"


#Tcell
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD4+ Activated Fos-lo"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD4+ Activated Fos-hi"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD4+ Memory"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD4+ PD1+"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD8+ IELs"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD8+ IL17+"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD8+ LP"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "MT-hi"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Tregs"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Cycling T"]<-"Tcell"


# ILS?NK
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "ILCs"]<-"Lymphoid"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "NKs"]<-"Lymphoid"




#B
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Cycling B"]<-"Bcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "GC"]<-"Bcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Plasma"]<-"Bcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Follicular"]<-"Bcell"

#Myeloid
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD69+ Mast"]<-"myeloid"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD69- Mast"]<-"myeloid"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Inflammatory Monocytes"]<-"myeloid"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Macrophages"]<-"myeloid"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Cycling Monocytes"]<-"myeloid"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "DC1"]<-"myeloid"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "DC2"]<-"myeloid"

#mesenchymal


methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Glia"]<-"Glia"



# fibroblasts
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Inflammatory Fibroblasts"]<-"fibroblast"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "RSPO3+"]<-"fibroblast"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT2B+ Fos-hi"]<-"fibroblast"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT2B+ Fos-lo 1"]<-"fibroblast"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT2B+ Fos-lo 2"]<-"fibroblast"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT5B+ 1"]<-"fibroblast"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT5B+ 2"]<-"fibroblast"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT2B+ Fos-lo 1"]<-"fibroblast"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Myofibroblasts"]<-"fibroblast"

#endothelial
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Post-capillary Venules"]<-"endothelial"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Pericytes"]<-"endothelial"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Endothelial"]<-"endothelial"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Microvascular"]<-"endothelial"



write.csv(methodspaper@meta.data, "/metadata_100perc_smillie_20200628.csv")

Idents(methodspaper)<-"method"

Tcells<-subset(methodspaper, subset = Predicted_all == "Tcell")
Bcells<-subset(methodspaper, subset = Predicted_all == "Bcell")
Fibroblasts<-subset(methodspaper, subset = Predicted_all == "fibroblast")
Myeloid<-subset(methodspaper, subset = Predicted_all == "myeloid")
Epithelium<-subset(methodspaper, subset = Predicted_all == "epithelium")
Endothelial<-subset(methodspaper, subset = Predicted_all == "endothelial")
Lymphoid<-subset(methodspaper, subset = Predicted_all == "Lymphoid")
Glia<-subset(methodspaper, subset = Predicted_all == "Glia")

DefaultAssay(methodspaper) # must be RNA

Tcell_whole_splitcollagenase<-FindMarkers(Tcells, test.use="MAST", ident.1="wholecollagenase", ident.2="splitcollagenase")
write.csv(Tcell_whole_splitcollagenase, "/Tcell_whole_splitcollagenase.csv")

Tcell_whole_splitprotease<-FindMarkers(Tcells, test.use="MAST", ident.1="wholecollagenase", ident.2="splitprotease")
write.csv(Tcell_whole_splitprotease, "/Tcell_whole_splitprotease.csv")

Tcell_whole_splitprotease_epi<-FindMarkers(Tcells, test.use="MAST", ident.1="wholecollagenase", ident.2="splitprotease_epi")
write.csv(Tcell_whole_splitprotease_epi, "/Tcell_whole_splitprotease_epi.csv")

Bcell_whole_splitcollagenase<-FindMarkers(Bcells, test.use="MAST", ident.1="wholecollagenase", ident.2="splitcollagenase")
write.csv(Bcell_whole_splitcollagenase, "/Bcell_whole_splitcollagenase.csv")

Bcell_whole_splitprotease<-FindMarkers(Bcells, test.use="MAST", ident.1="wholecollagenase", ident.2="splitprotease")
write.csv(Bcell_whole_splitprotease, "/Bcell_whole_splitprotease.csv")

Bcell_whole_splitprotease_epi<-FindMarkers(Bcells, test.use="MAST", ident.1="wholecollagenase", ident.2="splitprotease_epi")
write.csv(Bcell_whole_splitprotease_epi, "/Bcell_whole_splitprotease_epi.csv")


Fibroblast_whole_splitcollagenase<-FindMarkers(Fibroblasts, test.use="MAST", ident.1="wholecollagenase", ident.2="splitcollagenase")
write.csv(Fibroblast_whole_splitcollagenase, "/Fibroblast_whole_splitcollagenase.csv")

Fibroblast_whole_splitprotease<-FindMarkers(Fibroblasts, test.use="MAST", ident.1="wholecollagenase", ident.2="splitprotease")
write.csv(Fibroblast_whole_splitprotease, "/Fibroblast_whole_splitprotease.csv")

Fibroblast_whole_splitprotease_epi<-FindMarkers(Fibroblasts, test.use="MAST", ident.1="wholecollagenase", ident.2="splitprotease_epi")
write.csv(Fibroblast_whole_splitprotease_epi, "/Fibroblast_whole_splitprotease_epi.csv")


Myeloid_whole_splitcollagenase<-FindMarkers(Myeloid, test.use="MAST", ident.1="wholecollagenase", ident.2="splitcollagenase")
write.csv(Myeloid_whole_splitcollagenase, "/Myeloid_whole_splitcollagenase.csv")

Myeloid_whole_splitprotease<-FindMarkers(Myeloid, test.use="MAST", ident.1="wholecollagenase", ident.2="splitprotease")
write.csv(Myeloid_whole_splitprotease, "/Myeloid_whole_splitprotease.csv")

Myeloid_whole_splitprotease_epi<-FindMarkers(Myeloid, test.use="MAST", ident.1="wholecollagenase", ident.2="splitprotease_epi")
write.csv(Myeloid_whole_splitprotease_epi, "/Myeloid_whole_splitprotease_epi.csv")


Epithelium_whole_splitcollagenase<-FindMarkers(Epithelium, test.use="MAST", ident.1="wholecollagenase", ident.2="splitcollagenase")
write.csv(Epithelium_whole_splitcollagenase, "/Epithelium_whole_splitcollagenase.csv")

Epithelium_whole_splitprotease<-FindMarkers(Epithelium, test.use="MAST", ident.1="wholecollagenase", ident.2="splitprotease")
write.csv(Epithelium_whole_splitprotease, "/Epithelium_whole_splitprotease.csv")

Epithelium_whole_splitprotease_epi<-FindMarkers(Epithelium, test.use="MAST", ident.1="wholecollagenase", ident.2="splitprotease_epi")
write.csv(Epithelium_whole_splitprotease_epi, "/Epithelium_whole_splitprotease_epi.csv")


Lymphoid_whole_splitcollagenase<-FindMarkers(Lymphoid, test.use="MAST", ident.1="wholecollagenase", ident.2="splitcollagenase")
write.csv(Lymphoid_whole_splitcollagenase, "/Lymphoid_whole_splitcollagenase.csv")

Lymphoid_whole_splitprotease<-FindMarkers(Lymphoid, test.use="MAST", ident.1="wholecollagenase", ident.2="splitprotease")
write.csv(Lymphoid_whole_splitprotease, "/Lymphoid_whole_splitprotease.csv")

Lymphoid_whole_splitprotease_epi<-FindMarkers(Lymphoid, test.use="MAST", ident.1="wholecollagenase", ident.2="splitprotease_epi")
write.csv(Lymphoid_whole_splitprotease_epi, "/Lymphoid_whole_splitprotease_epi.csv")

Endothelial_whole_splitcollagenase<-FindMarkers(Endothelial, test.use="MAST", ident.1="wholecollagenase", ident.2="splitcollagenase")
write.csv(Endothelial_whole_splitcollagenase, "/Endothelial_whole_splitcollagenase.csv")

Endothelial_whole_splitprotease<-FindMarkers(Endothelial, test.use="MAST", ident.1="wholecollagenase", ident.2="splitprotease")
write.csv(Endothelial_whole_splitprotease, "/Endothelial_whole_splitprotease.csv")

Endothelial_whole_splitprotease_epi<-FindMarkers(Endothelial, test.use="MAST", ident.1="wholecollagenase", ident.2="splitprotease_epi")
write.csv(Endothelial_whole_splitprotease_epi, "/Endothelial_whole_splitprotease_epi.csv")

Glia_whole_splitcollagenase<-FindMarkers(Glia, test.use="MAST", ident.1="wholecollagenase", ident.2="splitcollagenase")
write.csv(Glia_whole_splitcollagenase, "/Glia_whole_splitcollagenase.csv")

Glia_whole_splitprotease<-FindMarkers(Glia, test.use="MAST", ident.1="wholecollagenase", ident.2="splitprotease")
write.csv(Glia_whole_splitprotease, "/Glia_whole_splitprotease.csv")

Glia_whole_splitprotease_epi<-FindMarkers(Glia, test.use="MAST", ident.1="wholecollagenase", ident.2="splitprotease_epi")
write.csv(Glia_whole_splitprotease_epi, "/Glia_whole_splitprotease_epi.csv")



saveRDS(methodspaper, "/SCT_integration_100perc_smillie.rds")
