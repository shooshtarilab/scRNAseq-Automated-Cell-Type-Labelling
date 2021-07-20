#!/bin/bash

#module load gcc
#module load r/3.5.2
#module load perl/5.22.2
#module load java/1.8

cd ./program/
pwd

export PERL5LIB=./bin/perl_modules/

Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/chung_breast/seurat_fqn_E_xy_matrix.tsv -c tme_inputs/chung_breast/gene_signatures.gmt -g tme_inputs/chung_breast/seurat_gold_std.tsv -o seurat_output/cb_results/ -p cb_seurat_fqn -t y -m GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -n GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -v y -a 0.4,1 -b 0,1  -s 100,100,0 -r 1-1 -y ranks -z 2000

Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/darmanis_glioblastoma/seurat_fqn_E_xy_matrix.tsv -c tme_inputs/darmanis_glioblastoma/gene_signatures.gmt -g tme_inputs/darmanis_glioblastoma/seurat_gold_std.tsv -o output/dg_results/ -p dg_seurat_fqn -t y -m GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY  -n GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -v y -a 0.4,1 -b 0,1 -s 100,100,0 -r 1-1 -y ranks -z 2000

Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/jerby_arnon_melanoma/seurat_fqn_E_xy_matrix.tsv -c tme_inputs/jerby_arnon_melanoma/gene_signatures.gmt -g tme_inputs/jerby_arnon_melanoma/seurat_gold_std.tsv -o seurat_output/jam_results/ -p jam_seurat_fqn -t y -m GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -n GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -v y -a 0.4,1 -b 0,1 -s 100,100,0 -r 1-1 -y ranks -z 2000

Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/li_crc/seurat_fqn_E_xy_matrix.tsv -c tme_inputs/li_crc/gene_signatures.gmt -g tme_inputs/li_crc/seurat_gold_std.tsv -o seurat_output/li_crc_results/ -p li_crc_seurat_fqn -t y -m GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY  -n GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY  -v y -a 0.4,1 -b 0,1 -s 100,100,0 -r 1-1 -y ranks -z 2000

Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/lambrechts_lung_carcinoma/seurat_fqn_E_xy_matrix.tsv -c tme_inputs/lambrechts_lung_carcinoma/gene_signatures.gmt -g tme_inputs/lambrechts_lung_carcinoma/seurat_gold_std.tsv -o output/llc_results/ -p llc_seurat_fqn -t y -m GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -n GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -v y -a 0.4,1 -b 0,1 -s 100,100,0 -r 1-1 -y ranks -z 2000

Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/peng_pancreatic/seurat_fqn_matrix.tsv -c tme_inputs/peng_pancreatic/gene_signatures.gmt -g tme_inputs/peng_pancreatic/seurat_gold_std.tsv -o seurat_output/peng_results/ -p peng_seurat_fqn -t y -m GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -n GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -v y -a 0.4,1 -b 0,1 -s 100,100,0 -r 1-1 -y ranks -z 2000

Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/tirosh_melanoma/seurat_fqn_E_xy_matrix.tsv -c tme_inputs/tirosh_melanoma/gene_signatures.gmt -g tme_inputs/tirosh_melanoma/seurat_gold_std.tsv -o seurat_output/tm_results/ -p tm_seurat_fqn -t y -m GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -n GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -v y -a 0.4,1 -b 0,1 -s 100,100,0 -r 1-1 -y ranks -z 2000
#TODO update paths to the right files for van galen tme inputs
Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/van_galen/combined_paper_combined_fqn_E_xy_matrix.tsv -c tme_inputs/van_galen/combined_gene_sigs.gmt -g tme_inputs/van_galen/combined_paper_combined_gold_std.tsv -o output/vg_results_combined/ -p vg_paper_fqn -t y -m GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -n GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -v y -a 0.4,1 -b 0,1 -s 100,100,0 -r 1-1 -y ranks -z 2000