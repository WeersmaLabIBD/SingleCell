**extract CD8 pos and neg bloodcell names from old seurat file and make new seurat file with raw data and these cells**
---

Author: WTC
Date: 20171020

**take filtered dataset containing all cells raw data**
```
datafile<-SetAllIdent(datafile, "tissue")
```
identify in 'cells_blood' which cells are the bloodcells that are left after filtering of the allcells-dataset
```
cells_blood<-WhichCells(datafile, "BLOOD")
```
subset data with only bloodcells
```
BLOODCELLS<-SubsetData(datafile, cells.use = cells_blood)
```
identify in 'CD8*_BLOODCELLS' which cells are CD8-positive and which are not
```
BLOODCELLS<-SetAllIdent(BLOODCELLS, "CD8ab.positive")
CD8pos_BLOODCELLS<-WhichCells(BLOODCELLS, "1")
CD8neg_BLOODCELLS<-WhichCells(BLOODCELLS, "0")
```
create dataframes containing only raw data of CD8*_mucosacells
```
raw_data<-data.frame(datafile@raw.data, check.names = F)
raw_data_CD8pos_blood<-raw_data[,CD8pos_BLOODCELLS]
raw_data_CD8neg_blood<-raw_data[,CD8neg_BLOODCELLS]
```

identify in 'cells_mucosa' which cells are the mucosacells that are left after filtering of the allcells-dataset
```
cells_mucosa<-WhichCells(datafile, c("IEL", "LPL"))
```
subset data with only mucosacells
```
MUCOSACELLS<-SubsetData(datafile, cells.use = cells_mucosa)
```
identify in 'CD8*_MUCOSACELLS' which cells are CD8-positive and which are not
```
MUCOSACELLS<-SetAllIdent(MUCOSACELLS, "CD8ab.positive")
CD8pos_MUCOSACELLS<-WhichCells(MUCOSACELLS, "1")
CD8neg_MUCOSACELLS<-WhichCells(MUCOSACELLS, "0")
```
create dataframes containing only war data of CD8*_mucosacells
```
raw_data_CD8pos_mucosa<-raw_data[,CD8pos_MUCOSACELLS]
raw_data_CD8neg_mucosa<-raw_data[,CD8neg_MUCOSACELLS]
```

**in the same way, you can subset IEL, LPL and BLOOD cells**
