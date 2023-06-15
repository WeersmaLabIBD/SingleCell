lanes<-c("lane1", "lane2", "lane3", "lane4", "lane5", "lane6", "lane7", "lane8", "lane9", "lane10")

x<-NULL
for(i in lanes){

y<-read.table(paste0("220504_", i, ".best"), header=T)
y$barcode_lane = substr(y$BARCODE,1,nchar(y$BARCODE)-2)
y$barcode_lane<-paste0(y$barcode_lane, "_220504_", i)
y<-y[c(5,13,15,21)]


metadata<-read.table("/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/seurat_preprocess_samples/ddtx_merged_demultiplexed_clustered_compartment_azi_elmentaiteadultileum_below60pctmito_epifilter_meta.csv")

merge<-merge(y, metadata, by="barcode_lane", all=F)

table(y$SNG.BEST.GUESS, y$DROPLET.TYPE)
table(metadata$donor_final, metadata$lane)


merge_single<-merge[merge$DROPLET.TYPE=="SNG",]
table(y$SNG.BEST.GUESS, y$DROPLET.TYPE)
table(merge_single$SNG.BEST.GUESS, merge_single$DROPLET.TYPE)
table(merge_single$donor_final, merge_single$SNG.BEST.GUESS)

x<-rbind(x, y)}

data<-readDRS("/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/seurat_preprocess_samples/objects/ddtx_merged_demultiplexed.rds")

#add demuxlet output to seuratfile
metadata<-data@meta.data
metadata<-merge(metadata, x,by="barcode_lane")
rownames(metadata)<-metadata$barcode_lane
metadata<-metadata[c(35,36,37)]
data <- AddMetaData(data, metadata)

# add clinical metadata
reference<-read.csv("/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/seurat_preprocess_samples/reference_2.csv", sep=";")
reference<-reference[-c(1,10)]
reference<-reference[-18,]

meta<-data@meta.data
dim(meta)
meta$pt_project_lane<-paste(meta$lane, meta$SNG.BEST.GUESS, sep="_")
meta<-merge(meta, reference, by="pt_project_lane", all=T)
dim(meta)
meta<-meta[c(7,39:47)]
head(meta)
rownames(meta)<-meta$barcode_lane
meta<-meta[-1]
data<-AddMetaData(data, meta)

# save
saveRDS(data, "/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/seurat_preprocess_samples/ddtx_merged_demultiplexed_demuxlet.rds")


