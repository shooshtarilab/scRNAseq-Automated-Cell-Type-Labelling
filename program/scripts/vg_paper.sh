#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task 4
#SBATCH --job-name=vg_paper
#SBATCH --output=%x-%j.out
#module load gcc
#module load r/3.5.2
#module load perl/5.22.2
#module load java/1.8

export PERL5LIB=./bin/perl_modules/

#Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/van_galen/paper_counts_E_xy_matrix.tsv -c tme_inputs/van_galen/gene_signatures.gmt -g tme_inputs/van_galen/paper_gold_std.tsv -o output/vg_results/ -p vg_paper_counts -t y -m CIBERSORT_BINARY -n CIBERSORT_BINARY -v y -a 0.4,1 -b 0,1 -s 100,100,0 -r 1-1 -y ranks -z 2000

Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/van_galen/paper_fqn_E_xy_matrix.tsv -c tme_inputs/van_galen/gene_signatures.gmt -g tme_inputs/van_galen/paper_gold_std.tsv -o output/vg_results/ -p vg_paper_fqn -t y -m GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -n GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -v y -a 0.4,1 -b 0,1 -s 100,100,0 -r 1-1 -y ranks -z 2000
##Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/van_galen/filtered_counts_E_xy_matrix.tsv -c tme_inputs/van_galen/gene_signatures.gmt -g tme_inputs/van_galen/filtered_gold_std.tsv -o output/vg_results/ -p vg_paper_counts -t y -m GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -n GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -v y -a 0.4,1 -b 0,1 -s 100,100,0 -r 1-1 -y ranks -z 2000
