#!/bin/bash
#SBATCH --time=0:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task 1
#SBATCH --job-name=souts/li_crc_adobo_seurat
#SBATCH --output=%x-%j.out

module load python/3.8.2

source ~/ENV/bin/activate

#python normalise_fqn.py ./jerby_arnon/filtered/counts.tsv.xz ./jerby_arnon/filtered/fqn.tsv.xz

(time python runs_adobo.py ./li_crc/counts_no_ensembl_int.tsv ./li_crc/crc_clusters.csv seurat ./output/li_crc/seurat_results.tsv) 2>times/crc.txt

