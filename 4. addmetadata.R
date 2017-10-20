CellsMeta = datafile@meta.data
head(CellsMeta)

patient<-data.frame(substr(colnames(datafile@data), 1, 1))
row.names(patient)=colnames(datafile@data)

CellsMeta["patient"] <- patient[,1]
head(CellsMeta)

CellsMetaTrim <- subset(CellsMeta, select = c("patient"))
head(CellsMetaTrim)

file_per_patient <- AddMetaData(allpts_perplexity500, CellsMetaTrim)

head(file_per_patient@data.info)

VlnPlot(file_per_patient, c("nGene", "nUMI", "patient"), nCol = 3)

## add location
#extract file that is a copy from @data.info
CellsMeta = datafile@data.info
head(CellsMeta)

#make a dataframe containing all elements you want to add to metadata
tissue<-data.frame(sapply(strsplit(colnames(datafile@data), split='_', fixed=TRUE), function(x) (x[2])))
#patient<-data.frame(substr(colnames(datafile@data), 1, 1))
#make sure this file has the colnames of your CellsMeta file as rownames
row.names(tissue)=colnames(datafile@data)

#add column 1 of "patient" as a column to CellsMeta file
CellsMeta["tissue"] <- tissue[,1]
head(CellsMeta)

#subset column to add
CellsMetaTrim <- subset(CellsMeta, select = c("tissue"))
head(CellsMetaTrim)

#add this column as metadata
datafile <- AddMetaData(datafile, CellsMetaTrim)

head(datafile@data.info)

colnames(datafile@data.info)

#check
VlnPlot(datafile, c("nGene", "nUMI", "tissue"), nCol = 3)



#### add metadata
CellsMeta = datafile@data.info
head(CellsMeta)

patient<-data.frame(sapply(strsplit(colnames(datafile@data), split='_', fixed=TRUE), function(x) (x[1])))

#patient<-data.frame(substr(colnames(datafile@data), 1, 2))
row.names(patient)=colnames(datafile@data)

CellsMeta["patient"] <- patient[,1]
head(CellsMeta)

CellsMetaTrim <- subset(CellsMeta, select = c("patient"))
head(CellsMetaTrim)

datafile <- AddMetaData(datafile, CellsMetaTrim)

head(datafile@data.info)
