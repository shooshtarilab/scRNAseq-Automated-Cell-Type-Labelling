#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task 4
#SBATCH --job-name=cb_paper
#SBATCH --output=%x-%j.out
#module load gcc
#module load r/3.5.2
#module load perl/5.22.2
#module load java/1.8

export PERL5LIB=./bin/perl_modules/
cd ~/Desktop/program/


Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/chung_breast/seurat_fqn_E_xy_matrix.tsv -c tme_inputs/chung_breast/gene_signatures.gmt -g tme_inputs/chung_breast/seurat_gold_std.tsv -o seurat_output/cb_results/ -p cb_seurat_fqn -t y -m GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -n GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -v y -a 0.4,1 -b 0,1  -s 100,100,0 -r 1-1 -y ranks -z 2000

#/usr/bin/time -o output/cbgsva.txt Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/chung_breast/paper_fqn_E_xy_matrix.tsv -c tme_inputs/chung_breast/gene_signatures.gmt -g tme_inputs/chung_breast/paper_gold_std.tsv -o output/null/ -p cb_paper_fqn -t y -m GSVA -n GSVA -s 100,100,0 -r 1-1 -y ranks -z 2000

#/usr/bin/time -o output/cbmeta.txt Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/chung_breast/paper_fqn_E_xy_matrix.tsv -c tme_inputs/chung_breast/gene_signatures.gmt -g tme_inputs/chung_breast/paper_gold_std.tsv -o output/null/ -p cb_paper_fqn -t y -m METANEIGHBOR_BINARY -n METANEIGHBOR_BINARY -s 100,100,0 -r 1-1 -y ranks -z 2000

#/usr/bin/time -o output/cbora.txt Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/chung_breast/paper_fqn_E_xy_matrix.tsv -c tme_inputs/chung_breast/gene_signatures.gmt -g tme_inputs/chung_breast/paper_gold_std.tsv -o output/null/ -p cb_paper_fqn -t y -m ORA -n ORA -s 100,100,0 -r 1-1 -y ranks -z 2000

#/usr/bin/time -o output/cbciber.txt Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/chung_breast/paper_fqn_E_xy_matrix.tsv -c tme_inputs/chung_breast/gene_signatures.gmt -g tme_inputs/chung_breast/paper_gold_std.tsv -o output/null/ -p cb_paper_fqn -t y -m CIBERSORT_BINARY -n CIBERSORT_BINARY -s 100,100,0 -r 1-1 -y ranks -z 2000

