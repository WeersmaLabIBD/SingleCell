# SingleCell

Repository of the code used in the LPMC v2 / CHANGE project

The steps for replicating the results of this study, are grouped in categories

# Contents

-   Alignment of the sequence data to the GRCh38 reference Human Genome
-   Demultiplexing_and_doublet_detection: demultiplexing and assignment of individuals to cells
-   General quality control and filtering
-   Normalization and data integration
-   Dimensional reduction and clustering
-   Cell type classification
-   Expression of integrins
-   Cell abundances analysis
-   Differential gene expression analysis and pathway analysis
-   Cell-cell interaction analysis

# Workflow

steps should be done in this order, to replicate the results

## alignment

1.  '*pre_preprocessing/alignment/lpmcv2_upload_badly_named_lanes.sh*' This will upload the BGI lanes to the cluster for the 2021 data

2.  '*pre_preprocessing/alignment/lpmcv2_rename_lanes.sh*' This will rename the BGI lanes to the format that CellRanger understands for the 2021 data

3.  '*pre_preprocessing/lpmcv2_move_batch2_fastqs*' this will move the BGI lanes to the correct directories for the 2022 data

4.  '*pre_preprocessing/alignment/lpmcv2_create_cellranger_jobs_introns*' This will create jobs to run each 10x lane through CellRanger. These jobs are then to be submitted

5.  '*pre_preprocessing/alignment/lpmcv2_download_websummaries*' This will download the web summary file outputs of Cellranger

## genotyping

6.  '*pre_preprocessing/genotyping/lpmcv2_prepare_genotypes_for_imputation_2021.R*' This will prepare the PSAM files to be in the format that the genotype imputation pipeline accepts them (2021 data)

7.  '*pre_preprocessing/genotyping/lpmcv2_prepare_genotypes_for_imputation_2022.R' This will prepare the PSAM files to be in the format that the genotype imputation pipeline accepts them (2022 data)

8.  '*pre_preprocessing/genotyping/lpmcv2_subset_genotype_data.sh*' This will subset the genotype data for each 10x lane, with only SNPs with a MAF\>0.05 for that lane

## demultiplexing

10. '*demultiplexing_and_doublet_detection/lpmcv2_create_souping_samples_jobs.sh*' This will create Souporcell demultiplexing jobs. These jobs are then to be submitted. Be sure that there is an unzipped version of the barcodes in the counts directory

## merging

11. '*preprocessing/lpmcv2_lane_to_seurat.R*' This will take a 10x lane that is the output of cellranger, and turn it into a Seurat object

12. '*preprocessing/lpmcv2_create_lane_to_seurat_jobs.sh*' This script will create jobs that read a 10x lane and turn it into Seurat objects using the lane_to_seurat.R. These jobs then need to be submitted.

13. '*preprocessing/lpmcv2_filter_and_merge_lanes.R*' This will merge the Seurat objects for each lane, and combine them into one object. Then the Souporcell output is read, and doublets are removed. Finally the cells with a mitochondrial gene percentage of 70 percent or higher are removed

## normalization

14.  '*preprocessing/lpmcv2_sct_normalize_and_cluster.R*' perform SCT normalization, clustering and UMAP dimensional reduction

## cell type assignment

15. '*celltype_assigning/split_compartments.R*' This will do a compartment assignment of immune, stromal and epipithelial to the Seurat object, based on marker gene expression

16. '*celltype_assigning/lpmcv2_create_elmentaite2021_objects.R*' create a Seurat object from the Elmentaite 2021 study

17. '*celltype_assigning/lpmcv2_create_martin_object.R*' create a Seurat object from the Martin 2019 study

18. '*celltype_assigning/lpmcv2_add_azimuth_classifications.R' add the cell type assignments from Azimuth to the object

19. '*celltype_assigning/lpmcv2_add_lower_celltypes.R' add lower resolution cell types
