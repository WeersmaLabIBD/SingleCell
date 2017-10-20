
**blood vs mucosa differential expression analysis using SCDE on cluster**
====
Author: WTC
Date: 20171020

*takes a lot of time, so performing on cluster is advised*

```
library("methods")
library(devtools)
library(ggplot2)
library(Matrix) # Sparse matrices
library(dplyr) # Dataframe manipulation
library(flexmix)
library(scde) # Differential Expression
```
**load data**
```
cd<-read.csv("raw_data_seurat.csv", header=TRUE, row.names=1, check.names=FALSE)
metadata<-read.csv("data_info_seurat.csv", header=TRUE, row.names=1)
```

filter raw data for cells<200 genes (expressed = >1 count) and genes<3 cells
```
cd2 <- clean.counts(cd, min.lib.size=200, min.reads = 1, min.detected = 3)
```
change metadatafile so that it contains only th filtered cells from datafile 'cd2'
```
dim(metadata)
dim(cd2)
xx=colnames(cd2)
metadata2<-metadata[xx,]
dim(metadata2)
dim(cd2)
```

make celltype file which will be used to define which cell belongs to which tissue category, using metadata2 file and converting IEL and LPL to MUCOSA
```
celltype = as.factor(gsub("IEL|LPL","MUCOSA",metadata2$tissue))
```
giving each value in the factor celltype the name of the corresponding cell
```
names(celltype) <- colnames(cd2)
```
calculate error models
---

```
err.mod <- scde.error.models(counts = cd2, groups = celltype, n.cores = 15,
           threshold.segmentation = TRUE, save.crossfit.plots = FALSE,
           save.model.plots = FALSE, verbose = 1)
```
estimate gene expression prior
---

```
gene.ex.prior <- scde.expression.prior(models = err.mod, counts = cd2, length.out = 400,
                                       show.plot = F)
```
Testing for differential expression
---
```
ediff <- scde.expression.difference(err.mod, cd2, gene.ex.prior, groups = celltype,
                    n.randomizations = 100, n.cores = 15, verbose = 1)
```
write out a table with all the results, showing most significantly different genes
```
write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ], file = “..”, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
```
