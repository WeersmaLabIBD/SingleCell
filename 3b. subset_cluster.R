## Subset cluster

# for example cluster 0
cluster0<-SubsetData(datafile, ident=0)

# Elbowplot
PCElbowPlot(cluster0.1, num.pc = 40)

# find clusters
cluster0 <- FindClusters(cluster0, pc.use = 1:30, resolution = 0.6, print.output = 0, save.SNN = T)

# Run/Plot TSNE
cluster0.1=RunTSNE(cluster0.1, dims.use=1:30)
TSNEPlot(cluster0.1)
table(cluster0.1@data.info$res.1, cluster0.1@data.info$patient)
table(cluster0.1@data.info$res.1, cluster0.1@data.info$tissue)

dev.copy(pdf,paste("TSNEPlot.pdf"))
dev.off()

# define RNA markers
markers_cluster0.1<-FindAllMarkers(cluster0.1)
markers_cluster0.1 %>% group_by(cluster) %>% top_n(10, avg_diff) -> top10
write.csv(markers_cluster0.1, "markers_cluster0.1.csv")

