#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task 1
#SBATCH --job-name=tm_sccatch
#SBATCH --output=%x-%j.out

module load nixpkgs/16.09
module load gcc/7.3.0 
module load r/4.0.2

(time Rscript run_sccatch.R ./data/tm_counts.tsv.xz ./data/tirosh_melanoma_clusters.csv 'seurat' 'Melanoma' Blood 'Peripheral blood' Skin) 2> times/tm_sccatch.txt

