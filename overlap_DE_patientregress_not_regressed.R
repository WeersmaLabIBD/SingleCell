## single cell sequencing in Crohn's 

# Author: Werna Uniken Venema
# year: 2018

## looking into overlap DE analyses for different patients using two methods:
# 1. comparing DE with and without patient regressed
# 2. comparing DE per patient

#### method 1 - patient regressed

# load DE file, not regressed for patient
DE_allpts<-read.csv("~/Dropbox/scRNAseq paper/Results_paper/markergenes/markers/allcells/eight_celltypes_markers.csv", row.names=1)
row.names(DE_allpts)<-NULL
dim(DE_allpts)
DE_allpts$gene_cluster<-paste0(DE_allpts$gene, sep="_", DE_allpts$cluster)

# use DE file, regressed for patient
# DE with genes >1%, MAST
seuratfile_allcells_patientregr<-SetAllIdent(seuratfile_allcells_patientregr, "eight_cell_types")
allcells_patientregr_DE_markers = FindAllMarkers(seuratfile_allcells_patientregr, min.pct = 0.01, only.pos = T, test.use = "MAST")
row.names(allcells_patientregr_DE_markers)<-NULL
allcells_patientregr_DE_markers<-allcells_patientregr_DE_markers[allcells_patientregr_DE_markers$p_val_adj < 0.05,]
dim(allcells_patientregr_DE_markers)
colnames(allcells_patientregr_DE_markers)[1]<-"p_val_patientregr"
colnames(allcells_patientregr_DE_markers)[2]<-"Av_logFC_patientregr"
colnames(allcells_patientregr_DE_markers)[3]<-"pct.1_patientregr"
colnames(allcells_patientregr_DE_markers)[4]<-"pct.2_patientregr"
colnames(allcells_patientregr_DE_markers)[5]<-"p_val_adj_patientregr"
colnames(allcells_patientregr_DE_markers)[6]<-"cluster_patientregr"
allcells_patientregr_DE_markers$gene_cluster<-paste0(allcells_patientregr_DE_markers$gene, sep="_", allcells_patientregr_DE_markers$cluster_patientregr)
overlap_patient_w_o_regress<-merge(allcells_patientregr_DE_markers, DE_allpts, by="gene_cluster", all=F)

write.csv(overlap_patient_w_o_regress, file="~/Desktop/Single_cell/final_data_paper/DE/overlap_patient_w_o_regress.csv")

####method 2 - overlap patients

# use DE file, regressed for patient
# DE with genes >1%, MAST
seuratfile_pt1cells<-SetAllIdent(seuratfile_pt1cells, "eight_cell_types")
pt1cells_DE_markers = FindAllMarkers(seuratfile_pt1cells, min.pct = 0.01, only.pos = T, test.use = "MAST")
row.names(pt1cells_DE_markers)<-NULL
pt1cells_DE_markers<-pt1cells_DE_markers[pt1cells_DE_markers$p_val_adj < 0.05,]
dim(pt1cells_DE_markers)
colnames(pt1cells_DE_markers)[1]<-"p_val_pt1cells"
colnames(pt1cells_DE_markers)[2]<-"Av_logFC_pt1cells"
colnames(pt1cells_DE_markers)[3]<-"pct.1_pt1cells"
colnames(pt1cells_DE_markers)[4]<-"pct.2_pt1cells"
colnames(pt1cells_DE_markers)[5]<-"p_val_adj_pt1cells"
colnames(pt1cells_DE_markers)[6]<-"cluster_pt1cells"
pt1cells_DE_markers$gene_cluster<-paste0(pt1cells_DE_markers$gene, sep="_", pt1cells_DE_markers$cluster_pt1cells)

# use DE file, regressed for patient
# DE with genes >1%, MAST
seuratfile_pt2cells<-SetAllIdent(seuratfile_pt2cells, "eight_cell_types")
pt2cells_DE_markers = FindAllMarkers(seuratfile_pt2cells, min.pct = 0.01, only.pos = T, test.use = "MAST")
row.names(pt2cells_DE_markers)<-NULL
pt2cells_DE_markers<-pt2cells_DE_markers[pt2cells_DE_markers$p_val_adj < 0.05,]
dim(pt2cells_DE_markers)
colnames(pt2cells_DE_markers)[1]<-"p_val_pt2cells"
colnames(pt2cells_DE_markers)[2]<-"Av_logFC_pt2cells"
colnames(pt2cells_DE_markers)[3]<-"pct.1_pt2cells"
colnames(pt2cells_DE_markers)[4]<-"pct.2_pt2cells"
colnames(pt2cells_DE_markers)[5]<-"p_val_adj_pt2cells"
colnames(pt2cells_DE_markers)[6]<-"cluster_pt2cells"
pt2cells_DE_markers$gene_cluster<-paste0(pt2cells_DE_markers$gene, sep="_", pt2cells_DE_markers$cluster_pt2cells)

## calculate overlap
overlap_patient1_2<-merge(pt1cells_DE_markers, pt2cells_DE_markers, by="gene_cluster", all=F)

write.csv(overlap_patient1_2, file="~/Desktop/Single_cell/final_data_paper/DE/overlap_patient1_2.csv")

# use DE file, regressed for patient
# DE with genes >1%, MAST
seuratfile_pt3cells<-SetAllIdent(seuratfile_pt3cells, "eight_cell_types")
pt3cells_DE_markers = FindAllMarkers(seuratfile_pt3cells, min.pct = 0.01, only.pos = T, test.use = "MAST")
row.names(pt3cells_DE_markers)<-NULL
pt3cells_DE_markers<-pt3cells_DE_markers[pt3cells_DE_markers$p_val_adj < 0.05,]
dim(pt3cells_DE_markers)
colnames(pt3cells_DE_markers)[1]<-"p_val_pt3cells"
colnames(pt3cells_DE_markers)[2]<-"Av_logFC_pt3cells"
colnames(pt3cells_DE_markers)[3]<-"pct.1_pt3cells"
colnames(pt3cells_DE_markers)[4]<-"pct.2_pt3cells"
colnames(pt3cells_DE_markers)[5]<-"p_val_adj_pt3cells"
colnames(pt3cells_DE_markers)[6]<-"cluster_pt3cells"
pt3cells_DE_markers$gene_cluster<-paste0(pt3cells_DE_markers$gene, sep="_", pt3cells_DE_markers$cluster_pt3cells)

## calculate overlap
overlap_patient1_2_3<-merge(overlap_patient1_2, pt3cells_DE_markers, by="gene_cluster", all=F)
write.csv(overlap_patient1_2_3, file="~/Desktop/Single_cell/final_data_paper/DE/overlap_patient1_2_3.csv")

