## extract CD risk genes from DE genes table
# open DE .txt file
DE_mucosa_CD4_CD8<-read.table("..")

#create column 'Gene'
DE_mucosa_CD4_CD8$Gene<-sapply(strsplit(row.names(DE_mucosa_CD4_CD8), split='-', fixed=TRUE), function(x) (x[2]))

#load genes of interest file
GOI<-read.table("..", sep=";", header=TRUE)

# match DE list and GOI
GOI_pos<-merge(GOI, DE_mucosa_CD4_CD8, by="Gene", all=FALSE)

# subset DE_mucosa_CD4_CD8 with only genes above DE value for two-sided testing (<-1.96 and >1.96)
GOI_DE_CD4pos<-subset(GOI_pos, cZ>1.96)
write.csv(GOI_DE_CD4pos, "..")
GOI_DE_CD8pos<-subset(GOI_pos, cZ<(-1.96))
write.csv(GOI_DE_CD8pos, "..")

#load genes of interest file
DRUGTARGETGENES<-read.table("..", sep=";", header=TRUE)

# match DE list and DRUGTARGETGENES
DRUGTARGETGENES_pos<-merge(DRUGTARGETGENES, DE_blood_CD4_CD8, by="Gene", all=FALSE)
# save matched list
write.csv(DRUGTARGETGENES_pos, "..")

# subset DE_mucosa_CD4_CD8 with only genes above DE value for two-sided testing (<-1.96 and >1.96)
DRUGTARGETGENES_CD4pos<-subset(DRUGTARGETGENES_pos, cZ>1.95)
write.csv(DRUGTARGETGENES_CD4pos, "..")
DRUGTARGETGENES_CD8pos<-subset(DRUGTARGETGENES_pos, cZ<(-1.95))
write.csv(DRUGTARGETGENES_CD8pos, "..")
