#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task 1
#SBATCH --job-name=souts/lam_lung_adobo
#SBATCH --output=%x-%j.out

module load python/3.8.2

source ~/ENV/bin/activate

#python normalise_fqn.py ./jerby_arnon/filtered/counts.tsv.xz ./jerby_arnon/filtered/fqn.tsv.xz

(time python runs_adobo.py ./llc/counts.tsv.xz ./llc/lambrechts_lung_clusters.csv seurat ./output/llc/seurat_results.tsv) 2>times/llc.sh

