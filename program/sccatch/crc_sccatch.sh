#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task 1
#SBATCH --job-name=crc_sccatch
#SBATCH --output=%x-%j.out

module load nixpkgs/16.09
module load gcc/7.3.0 
module load r/4.0.2

(time Rscript run_crc.R ./data/crc_counts.tsv.xz ./data/crc_clusters.csv 'seurat' 'Colorectal Cancer' 'Blood' 'Colon' Colorectum 'Gastrointestinal tract' Intestine Liver Lung 'Venous blood') 2> times/crc_sccatch.txt

