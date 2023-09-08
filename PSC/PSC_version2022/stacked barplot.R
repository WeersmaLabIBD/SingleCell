

library(Seurat)
library(ggplot2)
library(cowplot)

# read PBMC Seurat object (V2)
# PBMC_object <- "/Users/s.qs/Documents/Chapters/2 Chapter 2 - vedo2/2 Vedo2_PBMC_dataset/batch 1+2/20221027_ECCO/vedo2_pbmc_T0T1T4T5_no_doublet_Eryth_HSPC_Platelet.rds"
PBMC_object <- "/Users/s.qs/Documents/Chapters/2 Chapter 2 - vedo2/2 Vedo2_PBMC_dataset/batch 1+2/20221207_stacked_bar_plot/vedo2_pbmc_T0T1T4T5_no_doublet_Eryth_HSPC_Platelet.rds"
vedo2_PBMC <- readRDS(PBMC_object)

plot_loc <- "/Users/s.qs/Documents/Chapters/2 Chapter 2 - vedo2/2 Vedo2_PBMC_dataset/batch 1+2/20221207_stacked_bar_plot/output/"

get_nr_cells_per_group <- function(seurat_object, groups){
  # get the numbers table
  numbers_table <- data.frame(table(seurat_object@meta.data[, groups]))
  return(numbers_table)
}

create_barplot_frame <- function(numbers_table, number_column, stack_column, group_columns){
  # add a new column that has the groups
  numbers_table[['group']] <- apply(numbers_table, 1, function(x){
    # get the value of each column
    group_parts <- rep(NA, times = length(group_columns))
    for(i in 1:length(group_columns)){
      group_parts[i] <- x[group_columns[i]]
    }
    # and paste it together
    return(paste(group_parts, collapse = '_'))
  })
  # add a percentage
  numbers_table[['group_pct']] <- apply(numbers_table, 1, function(x){
    # get the group for this row
    group <- x['group']
    # and the number
    number <- as.numeric(x[number_column])
    # get the total of the group
    group_total <- sum(as.numeric(numbers_table[numbers_table[['group']] == group, number_column]))
    # calculate what percentage this entry makes up of the total
    pct <- number / group_total
    return(pct)
  })
  return(numbers_table)
}

create_barplot_percentages <- function(numbers_table, number_column, stack_column, group_columns, cell_types_to_plot=NULL, group_renames=NULL){
  numbers_table <- create_barplot_frame(numbers_table, number_column, stack_column, group_columns)
  # subset to specific cell types
  if (!is.null(cell_types_to_plot)) {
    numbers_table <- numbers_table[numbers_table[[stack_column]] %in% cell_types_to_plot, ]
  }
  if (!is.null(group_renames)) {
    numbers_table[['group']] <- as.vector(unlist(group_renames[numbers_table[['group']]]))
  }
  # grab the colours
  cc <- get_color_coding_dict()
  # subset to what we have
  colour_map <- cc[intersect(names(cc), numbers_table[[stack_column]])]
  # make the fill
  fillScale <- scale_fill_manual(name = "cell type",values = colour_map)
  # create the plot
  p <- ggplot(data = numbers_table, mapping = aes(x = numbers_table[['group']], y = numbers_table[['group_pct']], fill = numbers_table[[stack_column]])) +
    geom_bar(stat = 'identity', position = 'stack') +
    fillScale +
    xlab('group') +
    ylab('fraction') +
    theme(panel.border = element_rect(color='black', fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  return(p)
}
get_color_coding_dict <- function(){
  color_coding <- list()
  #ecco box
  color_coding[["NKs"]] <- "#ADD8E6"
  color_coding[["ILCs"]] <- "#8AACEB"
  color_coding[["DCs"]] <- "#6781F0"
  color_coding[["Monocytes"]] <- "#0B048C"
  color_coding[["B_cells"]] <- "#EC98CA"
  color_coding[["Plasma_cells"]] <- "#FFB6C1"
  color_coding[["Tregs"]] <- "#D979D3"
  color_coding[["CD8_T_cells"]] <- "#C55BDD"
  color_coding[["CD4_T_cells"]] <- "#A020F0"
  color_coding[["T_cells_other"]] <- "#7A04C4"
  return(color_coding)
}
#####get the number per status
cell_numbers <- get_nr_cells_per_group(vedo2_PBMC, c('ecco_plot_celltype', 'PGA_response', 'new_timepints'))
#
celltype_ratio <- create_barplot_frame(cell_numbers, 'Freq', 'ecco_plot_celltype', c('PGA_response', 'new_timepints'))

# all cell types
levels(cell_numbers$ecco_plot_celltype)
cell_numbers$ecco_plot_celltype <- factor(cell_numbers$ecco_plot_celltype, levels = levels(cell_numbers$ecco_plot_celltype)[c(4,6,5,7,1,2,3,8,10,9)])
celltype_ratio$group <- as.factor(celltype_ratio$group)
levels(celltype_ratio$group)
celltype_ratio$group <- factor(celltype_ratio$group, levels = levels(celltype_ratio$group)[c(4,3,2,1)])


all_celltypes <- create_barplot_percentages(cell_numbers, 'Freq', 'ecco_plot_celltype', c('PGA_response', 'new_timepints'), 
                                            cell_types_to_plot = c('NKs', 'ILCs', 'DCs', 'Monocytes', 'Plasma_cells', 'B_cells', 'CD8_T_cells', 'CD4_T_cells', 'Tregs', "T_cells_other"),
                                            group_renames = list('no_before' = 'NR T0', 'no_after' = 'NR T4', 'yes_before' = 'R T0', 'yes_after' = 'R T4'))

ggsave("1_all_cell_types_stacked_bar_plot.tiff", width = 4, height = 6, path = plot_loc)


#
innate_plot <- create_barplot_percentages(cell_numbers, 'Freq', 'ecco_plot_celltype', c('PGA_response', 'new_timepints'), cell_types_to_plot = c('NKs', 'ILCs', 'DCs', 'Monocytes'), group_renames = list('no_before' = 'NR T0', 'no_after' = 'NR T4', 'yes_before' = 'R T0', 'yes_after' = 'R T4'))
ggsave("2_innate_stacked_bar_plot.tiff", width = 10, height = 10, path = plot_loc)

adaptive_plot <- create_barplot_percentages(cell_numbers, 'Freq', 'ecco_plot_celltype', c('PGA_response', 'new_timepints'), cell_types_to_plot = c('Plasma_cells', 'B_cells', 'CD8_T_cells', 'CD4_T_cells', 'Tregs', "T_cells_other"), group_renames = list('no_before' = 'NR T0', 'no_after' = 'NR T4', 'yes_before' = 'R T0', 'yes_after' = 'R T4'))
ggsave("3_adaptive_stacked_bar_plot.tiff", width = 10, height = 10, path = plot_loc)






