## change rownames
CellsMeta = datafile@data
head(CellsMeta)
#create df with genes
gene<-data.frame(sapply(strsplit(rownames(datafile@data), split='-', fixed=TRUE), function(x) (x[2])))

#make genes unique rownames of seurat file
rownames(gene) = make.names(gene[,1], unique=TRUE)
rownames(CellsMeta) <- rownames(gene)
head(CellsMeta)
datafile@data = CellsMeta
head(datafile@data)


