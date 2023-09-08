
###################################################
# Title: basic DE analysis
# Date: 28-08-2023
###################################################


######################
# Library
######################

library(ggplot2) 
library(dplyr) 
library(tidyr) 
library(stringr)
library(patchwork)
library(grid)
library(gridExtra) 
library(ggpubr)
library(RColorBrewer) 
library(Seurat)
library(readr)
library(MAST)
library(readxl)
library(openxlsx)
library(enrichR)
require(tidyverse)
library(packcircles)

######################
# Main codes
######################

# read in dataframes
epi <- readRDS("Nieuw/epi_azimuth_duox2.rds")
imm <- readRDS("Nieuw/imm_azimuth_with_plasma.rds")
str <- readRDS("Nieuw/stro_azimuth.rds")

DefaultAssay(epi) = "RNA"
DefaultAssay(imm) = "RNA"
DefaultAssay(str) = "RNA"
Idents(epi) <- "celltype.final"
Idents(imm) <- "predicted.cell_type.pred"
Idents(str) <- "predicted.cell_type.pred"

# Create a table of frequencies
epi_frequencies <- as.data.frame.matrix(table(epi$celltype.final, epi$disease))
str_frequencies <- as.data.frame.matrix(table(str$predicted.cell_type.pred, str$disease))
imm_frequencies <- as.data.frame.matrix(table(imm$predicted.cell_type.pred, imm$disease))
frequencies <- rbind(epi_frequencies, str_frequencies, imm_frequencies)

#write_csv(frequencies, "Results/DE_2023/frequencies.csv")

# Create a table with rows celltype and column disease/analysis, with number of DE genes per analysis for each celltype for UC and PSC
files <- list.files(path="Results/DE_2023/DE_PSCINI_all/", pattern="*.csv", full.names=TRUE, recursive=FALSE)
total_PSC <- data.frame()
lapply(files, function(x) {
  t <- read.csv(x, header=TRUE) # load file
  celltype <- tools::file_path_sans_ext(basename(x))
  t$celltype <- celltype
  out <- t[t$p_val_adj < 0.05,]
    # write to file
  total_PSC <<- rbind(total_PSC, out)
})
total_PSC <- total_PSC[,c(6,8,9)]
DEpercelltype_PSC <- total_PSC %>%
  group_by(celltype) %>%
  summarize(num_genes = n())

files <- list.files(path="Results/DE_2023/DE_UCINI_all/", pattern="*.csv", full.names=TRUE, recursive=FALSE)
total_UC <- data.frame()
lapply(files, function(x) {
  t <- read.csv(x, header=TRUE) # load file
  celltype <- tools::file_path_sans_ext(basename(x))
  t$celltype <- celltype
  out <- t[t$p_val_adj < 0.05,]
  # write to file
  total_UC <<- rbind(total_UC, out)
})
total_UC <- total_UC[,c(6,8,9)]
DEpercelltype_UC <- total_UC %>%
  group_by(celltype) %>%
  summarize(num_genes = n())

# Filter celltypes with n genes > 9
DEpercelltype_PSC_10 <- DEpercelltype_PSC[DEpercelltype_PSC$num_genes > 9,]
DEpercelltype_UC_10 <- DEpercelltype_UC[DEpercelltype_UC$num_genes > 9,]

# Make figure for PSC
  # Generate the layout. This function return a dataframe with one line per bubble. 
  # It gives its center (x and y) and its radius, proportional of the value
packing <- circleProgressiveLayout(DEpercelltype_PSC_10$num_genes, sizetype='area')

  # We can add these packing information to the initial data frame
data <- cbind(DEpercelltype_PSC_10, packing)

  # Check that radius is proportional to value. We don't want a linear relationship, since it is the AREA that must be proportional to the value
plot(data$radius, data$value)

  # The next step is to go from one center + a radius to the coordinates of a circle that
  # is drawn by a multitude of straight lines.
dat.gg <- circleLayoutVertices(packing, npoints=50)

ggplot() + 
  
  # Make the bubbles
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill=as.factor(id)), colour = "black", alpha = 0.6) +
  
  # Add text in the center of each bubble + control its size
  geom_text(data = data, aes(x, y, size=num_genes, label = celltype)) +
  scale_size_continuous(range = c(1,4)) +
  
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()

ggsave("/Users/amberbangma/Documents/R/PSC/Results/Circle_PSC.pdf")

# Make figure for UC
# Generate the layout. This function return a dataframe with one line per bubble. 
# It gives its center (x and y) and its radius, proportional of the value
packing <- circleProgressiveLayout(DEpercelltype_UC_10$num_genes, sizetype='area')

# We can add these packing information to the initial data frame
data <- cbind(DEpercelltype_UC_10, packing)

# Check that radius is proportional to value. We don't want a linear relationship, since it is the AREA that must be proportional to the value
plot(data$radius, data$value)

# The next step is to go from one center + a radius to the coordinates of a circle that
# is drawn by a multitude of straight lines.
dat.gg <- circleLayoutVertices(packing, npoints=50)


ggplot() + 
  
  # Make the bubbles
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill=as.factor(id)), colour = "black", alpha = 0.6) +
  
  # Add text in the center of each bubble + control its size
  geom_text(data = data, aes(x, y, size=num_genes, label = celltype)) +
  scale_size_continuous(range = c(1,4)) +
  
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()

ggsave("/Users/amberbangma/Documents/R/PSC/Results/Circle_UC.pdf",  width=6, height=6)
