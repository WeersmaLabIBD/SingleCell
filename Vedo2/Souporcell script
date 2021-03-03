#!/bin/bash
#SBATCH --job-name=souporcell_200612_lane1
#SBATCH --output=/groups/umcg-weersma/tmp01/Emilia/souporcell/souporcell_output/souporcell_200612_lane1.out
#SBATCH --error=/groups/umcg-weersma/tmp01/Emilia/souporcell/souporcell_output/souporcell_200612_lane1.err
#SBATCH --time=6-23:59:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=100gb
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --tmp=1000gb

ml PythonPlus
python /groups/umcg-weersma/tmp01/Amber/souporcell/souporcell/renamer_v2.py \
--bam /groups/umcg-weersma/tmp01/Emilia/final_output/200612_lane1/cellranger_200612_lane1/outs/possorted_genome_bam.bam \
--barcodes /groups/umcg-weersma/tmp01/Emilia/souporcell/barcodes/barcodes_200612_lane1.tsv \
--out /groups/umcg-weersma/tmp01/Emilia/souporcell/remap/200612_lane1_remap.fq.fq

/groups/umcg-weersma/tmp01/Amber/souporcell/minimap2-2.17_x64-linux/minimap2 -ax splice -t 8 -G50k -k 21 -w 11 --sr -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N2
0 -f1000,5000 -n2 -m20 -s40 -g2000 -2K50m \
--secondary=no /groups/umcg-weersma/tmp01/Amber/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa \
/groups/umcg-weersma/tmp01/Emilia/souporcell/remap/200612_lane1_remap.fq.fq \
> /groups/umcg-weersma/tmp01/Emilia/souporcell/remap/minimap_200612_lane1.sam

python /groups/umcg-weersma/tmp01/Amber/souporcell/souporcell/retag_v2.py \
--sam /groups/umcg-weersma/tmp01/Emilia/souporcell/remap/minimap_200612_lane1.sam \
--out /groups/umcg-weersma/tmp01/Emilia/souporcell/remap/minitagged_200612_lane1.bam

ml SAMtools/1.9-foss-2018b
samtools sort /groups/umcg-weersma/tmp01/Emilia/souporcell/remap/minitagged_200612_lane1.bam \
> /groups/umcg-weersma/tmp01/Emilia/souporcell/remap/minitagged_sorted_200612_lane1.bam
samtools index /groups/umcg-weersma/tmp01/Emilia/souporcell/remap/minitagged_sorted_200612_lane1.bam

ml freebayes
freebayes -f /groups/umcg-weersma/tmp01/Amber/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa -iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6 -g 100000 \
/groups/umcg-weersma/tmp01/Emilia/souporcell/remap/minitagged_sorted_200612_lane1.bam \
> /groups/umcg-weersma/tmp01/Emilia/souporcell/variants/200612_lane1.vcf

/groups/umcg-weersma/tmp01/Amber/souporcell/vartrix-v1.1.4-x86_64-linux/vartrix \
--umi \
--mapq 30 \
-b /groups/umcg-weersma/tmp01/Emilia/final_output/200612_lane1/cellranger_200612_lane1/outs/possorted_genome_bam.bam \
-c /groups/umcg-weersma/tmp01/Emilia/souporcell/barcodes/barcodes_200612_lane1.tsv \
--scoring-method coverage \
--threads 8 \
--ref-matrix /groups/umcg-weersma/tmp01/Emilia/souporcell/allele_counting/200612_lane1_wes_ref.mtx \
--out-matrix /groups/umcg-weersma/tmp01/Emilia/souporcell/allele_counting/200612_lane1_wes_alt.mtx \
-v /groups/umcg-weersma/tmp01/Emilia/souporcell/variants/200612_lane1.vcf \
--fasta /groups/umcg-weersma/tmp01/Amber/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa

/groups/umcg-weersma/tmp01/Amber/souporcell/souporcell/souporcell/target/release/souporcell \
--alt_matrix /groups/umcg-weersma/tmp01/Emilia/souporcell/allele_counting/200612_lane1_wes_alt.mtx \
--ref_matrix /groups/umcg-weersma/tmp01/Emilia/souporcell/allele_counting/200612_lane1_wes_ref.mtx \
--num_clusters 4 \
--threads 8 \
--barcodes /groups/umcg-weersma/tmp01/Emilia/souporcell/barcodes/barcodes_200612_lane1.tsv \
> /groups/umcg-weersma/tmp01/Emilia/souporcell/final_output/clusters_200612_lane1_tmp.tsv

/groups/umcg-weersma/tmp01/Amber/souporcell/souporcell/troublet/target/release/troublet \
-a /groups/umcg-weersma/tmp01/Emilia/souporcell/allele_counting/200612_lane1_wes_alt.mtx \
-r /groups/umcg-weersma/tmp01/Emilia/souporcell/allele_counting/200612_lane1_wes_ref.mtx \
--clusters /groups/umcg-weersma/tmp01/Emilia/souporcell/final_output/clusters_200612_lane1_tmp.tsv \
> /groups/umcg-weersma/tmp01/Emilia/souporcell/final_output/clusters_200612_lane1.tsv
