#!/bin/bash
#SBATCH --time=0:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task 1
#SBATCH --job-name=souts/tm_adobo_seurat
#SBATCH --output=%x-%j.out

module load python/3.8.2

source ~/ENV/bin/activate

(time python run_tm.py) 2>times/tm.txt

