## comparison demuxlet cell and doublet assignment results to souporcell results
# WTC june 2023

# prep demuxlet results for comparing to souporcell results

lane1<-read.table("/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/demultiplexing/demuxlet/output/220504_lane01.best", header=T)
lane1$barcode_lane = substr(lane1$BARCODE,1,nchar(lane1$BARCODE)-2)
lane1$barcode_lane = substr(lane1$BARCODE,1,nchar(lane1$BARCODE)-2)
lane1$barcode_lane<-paste0(lane1$barcode_lane, "_220504_lane01")

lane2<-read.table("/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/demultiplexing/demuxlet/output/220504_lane02.best", header=T)
lane2$barcode_lane = substr(lane2$BARCODE,1,nchar(lane2$BARCODE)-2)
lane2$barcode_lane = substr(lane2$BARCODE,1,nchar(lane2$BARCODE)-2)
lane2$barcode_lane<-paste0(lane2$barcode_lane, "_220504_lane02")


lane3<-read.table("/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/demultiplexing/demuxlet/output/220504_lane03.best", header=T)
lane3$barcode_lane = substr(lane3$BARCODE,1,nchar(lane3$BARCODE)-2)
lane3$barcode_lane = substr(lane3$BARCODE,1,nchar(lane3$BARCODE)-2)
lane3$barcode_lane<-paste0(lane3$barcode_lane, "_220504_lane03")

lane4<-read.table("/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/demultiplexing/demuxlet/output/220504_lane04.best", header=T)
lane4$barcode_lane = substr(lane4$BARCODE,1,nchar(lane4$BARCODE)-2)
lane4$barcode_lane = substr(lane4$BARCODE,1,nchar(lane4$BARCODE)-2)
lane4$barcode_lane<-paste0(lane4$barcode_lane, "_220504_lane04")

lane5<-read.table("/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/demultiplexing/demuxlet/output/220504_lane05.best", header=T)
lane5$barcode_lane = substr(lane5$BARCODE,1,nchar(lane5$BARCODE)-2)
lane5$barcode_lane = substr(lane5$BARCODE,1,nchar(lane5$BARCODE)-2)
lane5$barcode_lane<-paste0(lane5$barcode_lane, "_220504_lane05")

lane6<-read.table("/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/demultiplexing/demuxlet/output/220504_lane06.best", header=T)
lane6$barcode_lane = substr(lane6$BARCODE,1,nchar(lane6$BARCODE)-2)
lane6$barcode_lane = substr(lane6$BARCODE,1,nchar(lane6$BARCODE)-2)
lane6$barcode_lane<-paste0(lane6$barcode_lane, "_220504_lane06")

lane7<-read.table("/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/demultiplexing/demuxlet/output/220504_lane07.best", header=T)
lane7$barcode_lane = substr(lane7$BARCODE,1,nchar(lane7$BARCODE)-2)
lane7$barcode_lane = substr(lane7$BARCODE,1,nchar(lane7$BARCODE)-2)
lane7$barcode_lane<-paste0(lane7$barcode_lane, "_220504_lane07")

lane8<-read.table("/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/demultiplexing/demuxlet/output/220504_lane08.best", header=T)
lane8$barcode_lane = substr(lane8$BARCODE,1,nchar(lane8$BARCODE)-2)
lane8$barcode_lane = substr(lane8$BARCODE,1,nchar(lane8$BARCODE)-2)
lane8$barcode_lane<-paste0(lane8$barcode_lane, "_220504_lane08")

lane9<-read.table("/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/demultiplexing/demuxlet/output/220504_lane09.best", header=T)
lane9$barcode_lane = substr(lane9$BARCODE,1,nchar(lane9$BARCODE)-2)
lane9$barcode_lane = substr(lane9$BARCODE,1,nchar(lane9$BARCODE)-2)
lane9$barcode_lane<-paste0(lane9$barcode_lane, "_220504_lane09")

lane10<-read.table("/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/demultiplexing/demuxlet/output/220504_lane10.best", header=T)
lane10$barcode_lane = substr(lane10$BARCODE,1,nchar(lane10$BARCODE)-2)
lane10$barcode_lane = substr(lane10$BARCODE,1,nchar(lane10$BARCODE)-2)
lane10$barcode_lane<-paste0(lane10$barcode_lane, "_220504_lane10")

x<-rbind(lane1, lane2)
x<-rbind(x, lane3)
x<-rbind(x, lane4)
x<-rbind(x, lane5)
x<-rbind(x, lane6)
x<-rbind(x, lane7)
x<-rbind(x, lane8)
x<-rbind(x, lane9)
x<-rbind(x, lane10)


x<-x[c(5,13,15,21)]
dim(x)

# load mito, epi and mislabeled sample filtered souporcell results
metadata<-read.table("/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/seurat_preprocess_samples/ddtx_merged_demultiplexed_clustered_compartment_azi_elmentaiteadultileum_below60pctmito_epifilter_meta.csv")
dim(metadata)

## demuxlet assigned singlets largely overlap with souporcell assigned singlets 
merge<-merge(x, metadata, by="barcode_lane", all=F)
dim(merge)
merge_single<-merge[merge$DROPLET.TYPE=="SNG",]
table(merge_single$donor_final, merge_single$SNG.BEST.GUESS, merge_single$DROPLET.TYPE)

#read demultiplexed Seuratobject, not doublet filtered
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


