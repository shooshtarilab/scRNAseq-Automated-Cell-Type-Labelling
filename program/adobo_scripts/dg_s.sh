#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task 1
#SBATCH --job-name=souts/dg_adobo_seurat
#SBATCH --output=%x-%j.out

module load python/3.8.2

source ~/ENV/bin/activate

#python normalise_fqn.py ./jerby_arnon/filtered/counts.tsv.xz ./jerby_arnon/filtered/fqn.tsv.xz

(time python runs_adobo.py ./dg/counts.tsv.xz ./dg/darmanis_clusters.csv seurat ./output/dg/seurat_results.tsv) 2>times/dg.txt

