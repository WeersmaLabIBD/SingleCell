singularity exec --bind /groups/umcg-weersma/tmp01/ /groups/umcg-weersma/tmp01/tools/sc-eqtlgen-consortium-pipeline/wg1/wg1-pipeline-20230308_2.simg popscle demuxlet \
        --sam /groups/umcg-weersma/tmp01/projects/ddtx/ongoing/demultiplexing/popscle_tools_filtered/alignment/filtered_alignment//ddtx_220504_laneX_genofiltered_possorted_genome_bam.bam \
        --tag-group CB \
        --tag-UMI UB \
        --field GT \
        --vcf /groups/umcg-weersma/tmp01/projects/ddtx/ongoing/demultiplexing/popscle_tools_filtered/genotype/sorted_genotype/ddtx_bamsorted_220504_laneX_maf005.vcf.gz \
        --out /groups/umcg-weersma/tmp01/projects/ddtx/ongoing/demultiplexing/demuxlet/output/220504_laneX \
        --group-list /groups/umcg-weersma/tmp01/projects/ddtx/processed/alignment/cellranger_output/GRCh38/220504_laneX/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
