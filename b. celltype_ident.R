##CD8 mucosa	(18 res 0.6)
#cytotoxic		(cluster ?)
b=datafile
DotPlot(b, c("CD160", "TGFBR1", "SLAMF7", "CAPG"), dot.scale=15)	
FeaturePlot(b, c("CD160", "TGFBR1", "SLAMF7", "CAPG"), cols.use = c( "azure2","deeppink1"), pt.size=2, no.axes = T)	
VlnPlot(b, c("CD160", "TGFBR1", "SLAMF7", "CAPG"))
b=RenameIdent(b, 0, "CT-cells")
# central memory/naïve T (cluster ?)
DotPlot(b, c("CCR7",	"SELL", "CD27", "TCF7", "ID3"), dot.scale=15)
FeaturePlot(b, c("CCR7",	"SELL", "CD27", "TCF7", "ID3"), cols.use = c( "azure2","deeppink1"), pt.size=2, no.axes = T)
VlnPlot(b, c("CCR7",	"SELL", "CD27", "TCF7", "ID3"))
b=RenameIdent(b, , "CM/N T-cells")
# tissue resident effector memory	 (cluster 1)
DotPlot(b, c("IFNG",	"IL17A", "RGS1"), dot.scale=15)	
FeaturePlot(b, c("IFNG",	"IL17A", "RGS1"), cols.use = c( "azure2","deeppink1"), pt.size=2, no.axes = T)	
VlnPlot(b, c("IFNG",	"IL17A", "RGS1"))
b=RenameIdent(b, 1, "EM T-cells")

## CD4 mucosa	(17pcs, res 0.8)
#NK T cell	(cluster 2)
b=allptsmucosa_cd4
DotPlot(b, c("CAPG", "CD160"	, "GZMA", "NKG7"), dot.scale=15)		
FeaturePlot(b, c("CAPG", "CD160","GZMA", "NKG7"),cols.use = c( "azure2","deeppink1"), pt.size=2)
VlnPlot(b, c("CAPG", "CD160","GZMA", "NKG7"))
TSNEPlot(b)
FeaturePlot(b,"RGS1")
b=RenameIdent(b, c(2), "NK T-cells")

# memory effector cells, tissue resident effector cells Th1, Th2, Treg, Th17, Tfh (cluster 0)
FeaturePlot(b,c("IL2", "IFNG", "TNF", "LTA", "STAT4", "TBX21", "HLX"), cols.use = c( "azure2","deeppink1"), pt.size=2)
FeaturePlot(b,c("GATA3", "MAF", "IRF4", "GFI1", "STAT6", "IL4"), cols.use = c( "azure2","deeppink1"), pt.size=2)
# low in SELL, high in RGS1 - tissue resident effector cells
FeaturePlot(b, c("SELL", "RGS1"),cols.use = c( "azure2","deeppink1"), pt.size=2)
b=RenameIdent(b, "TRE-cells", "EM T-cells")

# central memory (cluster 3)
DotPlot(b, c("CCR7", "SELL"), dot.scale=15)
FeaturePlot(b, c("CCR7", "SELL"),cols.use = c( "azure2","deeppink1"), pt.size=2)
VlnPlot(b, c("CCR7", "SELL"))
b=RenameIdent(b, c(3), "CM-cells")

# Th17 (within cluster 0)
DotPlot(b, c("IFNG", "CSF2", "IL17A"), dot.scale=15)
FeaturePlot(b, c("IFNG", "CSF2", "IL17A"),cols.use = c( "azure2","deeppink1"), pt.size=2)
VlnPlot(b, c("IFNG", "CSF2", "IL17A"))

#Treg (cluster 1)
DotPlot(b, c("CTLA4", "FOXP3", "TNFRSF18", "TIGIT", "CCR7"), dot.scale=15)
FeaturePlot(b, c("CTLA4", "FOXP3", "TNFRSF18", "TNFRSF4","TIGIT", "CCR7", "ID2", "ID3"),cols.use = c( "azure2","deeppink1"), pt.size=1)
VlnPlot(b, c("CTLA4", "FOXP3", "TNFRSF18","TNFRSF4", "TIGIT", "CCR7", "SELL"))
b=RenameIdent(b,1, "CMTreg-cells")



## CD4 blood	(20)
#NK T cell (cluster 4)
b=blood_cd4_no_outliers
a=c("GZMA", "TBX21", "FCGR3A", "NCR1", "CD160", "NKG7", "SLAMF6")
DotPlot(b, a, dot.scale=15)
FeaturePlot(b, a)
VlnPlot(b,a)			
# central memory/naïve T (cluster 3)
a=c("SELL", "CCR7", "LEF1", "FCGR3A")
FeaturePlot(object = b, features.plot = a, cols.use = c("grey", "yellow", "orange", "red"), overlay = TRUE, no.legend = FALSE, pt.size=2)
DotPlot(b, a, dot.scale=15)
FeaturePlot(b, a)
VlnPlot(b,a)
# Th17 (cluster )
a=c("KLRB1", "TGFBR1", "IL12RB1", "ENTPD1", "IFNG", "ABCB1", "CCR6", "CCL20", "RORA", "CSF2", "RORC")
DotPlot(b, a, dot.scale=15)
FeaturePlot(b, a)
VlnPlot(b,a)
# Th22 - cluster 0
a=c("CCR10", "AHR")

# Treg # cluster 1
a=c("IL2RA", "TNFRSF18", "TNFRSF4", "CTLA4", "CAPG", "AHR", "FOXP3", "TIGIT")
DotPlot(allptsblood_regtissue_CD4, a, dot.scale=15)
FeaturePlot(b, a)


## CD8 blood	(15)
# cytotoxic	(cluster 2)
a=c("NKG7", "GNLY")
b=allptsblood_regtissue_CD8
DotPlot(b, a, dot.scale = 20)
FeaturePlot(b, a)
VlnPlot(b,a)
# central memory/naïve T		(cluster (0),1)
a=c("CCR7","SELL", "CCL5")
DotPlot(b, a, dot.scale = 20)
FeaturePlot(b, a)
VlnPlot(b,a)
# Treg?? (cluster )
a=c("TNFRSF4", "TNFRSF18")
DotPlot(b, a, dot.scale = 20)
FeaturePlot(b, a)
VlnPlot(b,a)
TSNEPlot(b)


## CD8 celltypes liver cancer article
a=c("LEF1",  "CCR7") # teruggevonden
FeaturePlot(b, a)
c=blood_cd8_no_outliers
FeaturePlot(c, a)
a=c("CX3CR1", "FGFBP2", "FCGR3A") # teruggevonden

## CD4 celltypes liver cancer article
a=c("CCR7", "TCF7") ## ook TCF7-only cluster
b=allptsblood_regtissue_CD4
c=blood_cd4_no_outliers
FeaturePlot(b, a)
FeaturePlot(c, a)

a=c("SELL", "IL2RA", "FOXP3" ) # not FOXP3, anders cluster 5 en 0
a=c("TIGIT", "CTLA4") # no co=clusters
a=c("CCR8", "TNFRSF9", "LAYN" ) # no co-clusters
a=c("TNFRSF4", "TNFRSF18") # cluster 5 and 0
a=c("SELL", "IL2RA","TIGIT", "CTLA4","TNFRSF4", "TNFRSF18")
FeaturePlot(object = c, features.plot = a, cols.use = c("grey", "yellow", "orange", "red"), overlay = TRUE, no.legend = FALSE, pt.size=2)

b=allptsblood_regtissue_CD4
FeaturePlot(b, a) 
a=c("THSD4", "RP11.510N19.5", "DNAH9", "MTND2P5", "FABP1", "TRERNA1") # high expression marking cd4_blood_cluster_5 
a=c("HNRNPL", "EIF5A", "PCBP2", "CTNNB1", "CALR", "IFI6", "ACOT13", "RAB12", "NUTF2") # cluster 2 low expression of these genes
a=c("GNLY", "CCL5", "NKG7", "GZMH", "PRF1", "GZMB", "GPR56", "CST7", "GZMH") # high expression cluster 4
a=c("RPS19", "TMSB10") # cluster0 low expression blood cd4, hgh expression "HNRNPL" DCAKD ADAM19 RNF2 SEC22B
a=c("NOG", "CCR7") # high expression cluster3, low expression GZMH
a=c("S100A9", "LYZ", "CYBB", "S100A8", "CSF3R") # cluster 6
a=c("TNFRSF4", "PDCD5", "TAGLN2", "TNFRSF18") # cluster 1

## CD8 blood	(15)
# cytotoxic	(cluster 2)
a=c("NKG7", "GNLY", "CCR7","SELL", "CCL5","TNFRSF4", "TNFRSF18" )
b=allptsmucosa_cd8
DotPlot(b, a, dot.scale = 20)
FeaturePlot(b, a, cols.use = c( "azure2","deeppink1"), pt.size=2)
VlnPlot(b,a)
# central memory/naïve T		(cluster 0)
a=c("CCR7","SELL", "CCL5")
DotPlot(b, a, dot.scale = 20)
FeaturePlot(b, a)
VlnPlot(b,a)
# Treg?? (cluster 1 )
a=c("TNFRSF4", "TNFRSF18") 
DotPlot(b, a, dot.scale = 20)
FeaturePlot(b, a)
VlnPlot(b,a)
TSNEPlot(b)


## CD8 celltypes liver cancer article
a=c("LEF1",  "CCR7") # teruggevonden, cluster 0,1
FeaturePlot(b, a)
a=c("CX3CR1", "FGFBP2", "FCGR3A") # teruggevonden in cluster 2 subset!
FeaturePlot(b, a)
allptsblood_regtissue_CD8<-UpdateSeuratObject(allptsblood_regtissue_CD8)
b<-FindClusters(b, dims.use = 3:35)
b <- RunTSNE(object = b, dims.use = 3:35, do.fast = TRUE)
TSNEPlot(b)
FeaturePlot(allptsmucosa_regtissue_CD4, c("MKI67", "CTLA4", "FOXP3", "TNFRSF18", "TNFRSF4","TIGIT", "CCR7", "ID2", "ID3"),cols.use = c( "azure2","deeppink1"), pt.size=1)

FeaturePlot(allptsmucosa_regtissue_CD4, c("CXCR5", "PRDM1", "ICOS", "CXCR3", "CXCR4", "CTLA4", "SELL", "ID2", "ID3"),cols.use = c( "azure2","deeppink1"), pt.size=1)
TSNEPlot(b)

FeaturePlot(b,c("IL7R","IL2RA", "FOXP3"), cols.use = c( "azure2","deeppink1"), pt.size=2)
FeaturePlot(LPL_IEL_sig,c("CTLA4", "FOXP3", "TNFRSF18", "TNFRSF4","TIGIT", "CCR7", "SELL"),cols.use = c( "azure2","deeppink1"), pt.size=1)
FeaturePlot(b,c("TGIF1"))
FeaturePlot(b,c("LAG3", "ENTPD1", "SELL", "IL7R", "ITGAE", "TNFRSF18"),cols.use = c( "azure2","deeppink1"), pt.size=1)

FeaturePlot(LPL_IEL_sig,c("CD44", "SELL", "CD69", "ITGAE"),cols.use = c( "azure2","deeppink1"), pt.size=1) #TRM
FeaturePlot(LPL_IEL_sig,c("CD44", "SELL", "CD69", "ITGAE", "IL7R"),cols.use = c( "azure2","deeppink1"), pt.size=1)


TRM: CD44hiCD62L−CD69+CD103+, TCM: CD44hiCD62L+ CD69−CD103−CD127+, and TEM: CD44hiCD62L−CD69-CD103−CD127+