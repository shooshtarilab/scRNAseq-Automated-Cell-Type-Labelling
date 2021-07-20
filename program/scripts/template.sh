#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task 4
#SBATCH --job-name=
#SBATCH --output=%x-%j.out
module load gcc
module load r/3.5.2
module load perl/5.22.2
module load java/1.8

export PERL5LIB=./bin/perl_modules/

