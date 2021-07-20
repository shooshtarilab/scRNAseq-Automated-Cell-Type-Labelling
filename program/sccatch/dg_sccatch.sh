#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task 1
#SBATCH --job-name=dg_sccatch
#SBATCH --output=%x-%j.out

module load nixpkgs/16.09
module load gcc/7.3.0 
module load r/4.0.2

(time Rscript run_sccatch.R ./data/dg_counts.tsv.xz ./data/darmanis_clusters.csv 'seurat' 'Glioblastoma' 'Blood' 'Brain') 2> times/dg_sccatch.txt

