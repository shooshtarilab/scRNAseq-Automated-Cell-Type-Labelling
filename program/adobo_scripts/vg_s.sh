#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task 1
#SBATCH --job-name=souts/vg_adobo
#SBATCH --output=%x-%j.out

module load python/3.8.2

source ~/ENV/bin/activate


(time python runs_adobo.py ./vg/counts.tsv.xz ./vg/van_galen_clusters.csv seurat ./output/vg/seurat_results.tsv) 2>times/vg.txt

