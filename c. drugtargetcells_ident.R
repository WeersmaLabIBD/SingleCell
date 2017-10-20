## identify drugtargetcells

# CD8+ cells, DE in mucosa vs blood
a=c("C1QC","PPIF","TNF")
b=allptsmucosa_regtissue_CD8
DotPlot(b,a, cex.use=20)
FeaturePlot(b,a,cols.use = c( "azure2","deeppink1"))
TSNEPlot(b)

# CD8+ cells, DE in blood vs mucosa
a=c("ALOX5","FCGR3A") # FCGR3A marks cluster 2
b=allptsblood_regtissue_CD8
DotPlot(b,a, cex.use=20)
FeaturePlot(b,a,cols.use = c( "azure2","deeppink1"))
TSNEPlot(b)

# CD8- cells, DE in blood vs mucosa 
a=c("ALOX5","CRBN","FCGR1A","FCGR2A","FCGR3A","IMPDH1") # cluster6
b=allptsblood_regtissue_CD4
DotPlot(b,a, cex.use=20)
FeaturePlot(b,a,cols.use = c( "azure2","deeppink1"))
TSNEPlot(b)

# CD8- cells, DE in mucosa vs blood
a=c("C1S","IKBKB","NR3C1")
b=allptsmucosa_regtissue_CD4
DotPlot(b,a, cex.use=20)
FeaturePlot(b,a,cols.use = c( "azure2","deeppink1"))
TSNEPlot(b)

# mucosacells, DE in CD8- population
a=c("ATG9B","IFNG","MKI67","PLCG2")
b=allptsmucosa_regtissue_CD4
DotPlot(b,a, cex.use=20)
FeaturePlot(b,a,cols.use = c( "azure2","deeppink1"))
TSNEPlot(b)

# bloodcells, DE in CD8- population
a=c("ATG2A","ATG4A","ATG4C","FASLG","NOTCH1","OPN3")
b=allptsblood_regtissue_CD4
DotPlot(b,a, cex.use=20)
FeaturePlot(b,a,cols.use = c( "azure2","deeppink1"))
TSNEPlot(b)

# allcells, DE in CD8- cells
a=c("ATG4C","CSF2","PLCG2","RORC","SPP1")
b=..
DotPlot(b,a, cex.use=20)
FeaturePlot(b,a,cols.use = c( "azure2","deeppink1"))
TSNEPlot(b)