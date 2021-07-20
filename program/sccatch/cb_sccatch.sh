#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task 1
#SBATCH --job-name=cb_sccatch
#SBATCH --output=%x-%j.out

module load nixpkgs/16.09
module load gcc/7.3.0 
module load r/4.0.2

(time Rscript run_sccatch.R ./data/cb_counts.tsv.xz ./data/chung_breast_clusters.csv 'seurat' 'Breast Cancer' 'Blood' 'Breast' 'Mammary gland') 2> times/cb_sccatch.txt

