library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)
library(tidyr)
library(MAST)
library(readxl)
library(openxlsx)
library(reshape2)

# extract counts to matrix
counts=as.data.frame.matrix(table(data$Final_HTO, data$celltypes))

# calculate proportions
total_row = apply(counts, 1, sum)
pcts = lapply(counts, function(x) {
  x / total_row
})
frequencies = as.data.frame(pcts)
rownames(counts) = rownames(frequencies)

# add state
frequencies$disease =  sapply(strsplit(rownames(frequencies),"-"), `[`, 1)
frequencies$inflammation =  sapply(strsplit(rownames(frequencies),"-"), `[`, 2)
frequencies$state =  apply( frequencies[ ,c(43,44)] , 1 , paste , collapse = "-" )
frequencies = frequencies[,-c(43,44)]


# remove doublets and MT-Hi
frequencies=frequencies[,-c(6,10,15,19,25,32)]

# plot
df = melt(frequencies)
df$state = factor(df$state, levels = c("HC-NI", "UC-NI", "UC-I", "PSC-NI", "PSC-I"))
df$variable = factor(df$state, levels = c())
ggplot(df, aes(x=state, y=value)) + 
  geom_boxplot(outlier.size = 1) + 
  facet_wrap( ~ variable, scales="free") +
  theme_bw(base_size = 1)

ggsave("Results/Figures/frequencies.pdf", width = 33, height = 19)
