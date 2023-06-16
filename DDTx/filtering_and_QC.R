library(Seurat)

#open Seuratfile
ddtx<-readRDS("/source/ddtx_merged_demultiplexed_clustered_compartment_azi_elmentaiteadultileum_demuxlet.rds")
dim(ddtx)

# filter 60% mito
data<-subset(ddtx, subset=percent.mt < 60)
dim(data)
# reduction of 90553 to 88167 cells

# exclude T8. patient 3 (wrongly sampled, only recipient cells)
data@meta.data$pt_project_lane<-paste(data@meta.data$lane,data@meta.data$donor_final, sep="_")
data<-subset(data, subset=pt_project_lane!="220504_lane04_UMCGDDtx00005r")
dim(data)
# reduction to 85750 cells

# exclude epithelial cells from recipients: this is spillover
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Enterocyte"),"compartment_final"]<-"Epithelial"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Paneth"),"compartment_final"]<-"Epithelial"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="TA"),"compartment_final"]<-"Epithelial"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Tuft"),"compartment_final"]<-"Epithelial"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Stem cells"),"compartment_final"]<-"Epithelial"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="BEST4+ epithelial"),"compartment_final"]<-"Epithelial"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Goblet cell"),"compartment_final"]<-"Epithelial"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="I cells (CCK+)"),"compartment_final"]<-"Epithelial"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="K cells (GIP+)"),"compartment_final"]<-"Epithelial"


data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Activated CD8 T"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Activated CD4 T"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Memory B"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="IgA plasma cell"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="cDC2"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="LYVE1+ Macrophage"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Macrophages"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="IgG plasma cell"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Monocytes"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="TRGV2 gdT"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Mast cell"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Cycling B cell"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Naive B"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="ILC3"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="cDC1"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="TRGV5/7 gdT"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="gdT"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="CX3CR1+ CD8 Tmem"),"compartment_final"]<-"Immune"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="TRGV4 gdT"),"compartment_final"]<-"Immune"


data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Stromal 1 (ADAMDEC1+)"),"compartment_final"]<-"Stromal"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="arterial capillary"),"compartment_final"]<-"Stromal"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Stromal 1 (CCL11+)"),"compartment_final"]<-"Stromal"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Stromal 2 (NPY+)"),"compartment_final"]<-"Stromal"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Mature venous EC"),"compartment_final"]<-"Stromal"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Contractile pericyte (PLN+)"),"compartment_final"]<-"Stromal"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Adult Glia"),"compartment_final"]<-"Stromal"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Mature arterial EC"),"compartment_final"]<-"Stromal"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Transitional Stromal 3 (C3+)"),"compartment_final"]<-"Stromal"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="LEC3 (ADGRG3+)"),"compartment_final"]<-"Stromal"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="myofibroblast (RSPO2+)"),"compartment_final"]<-"Stromal"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="LEC1 (ACKR4+)"),"compartment_final"]<-"Stromal"
data@meta.data[which(data@meta.data$predicted.celltype.elmentaiteadultileum=="Stromal 3 (C7+)"),"compartment_final"]<-"Stromal"

data@meta.data$patient_compartment <- paste(data@meta.data$donor_final,data@meta.data$compartment_final, sep="_")

# subset the data
unique(data@meta.data$patient_compartment)
data_2 <- subset(data, patient_compartment != "UMCGDDtx00003r_Epithelial" & patient_compartment != "UMCGDDtx00004r_Epithelial" & patient_compartment !="UMCGDDtx00005r_Epithelial")
unique(data_2@meta.data$patient_compartment)
dim(data_2)
#reduction to 85554 cells

# save
saveRDS(data_2, "/source/ddtx_merged_demultiplexed_clustered_compartment_azi_elmentaiteadultileum_below60pctmito_withoutEpithelialR_demuxlet.rds")
