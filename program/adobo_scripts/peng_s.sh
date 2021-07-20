#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=116G
#SBATCH --cpus-per-task 1
#SBATCH --job-name=souts/peng_adobo
#SBATCH --output=%x-%j.out

module load python/3.8.2

source ~/ENV/bin/activate


(time python runs_adobo.py ./peng/counts.tsv ./peng/peng_clusters.csv seurat ./output/peng/seurat_results.tsv) 2>times/peng.txt

