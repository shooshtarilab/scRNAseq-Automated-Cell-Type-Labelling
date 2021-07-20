#!/bin/bash

export PERL5LIB=./bin/perl_modules/
cd ~/Desktop/program/
#run on paper

#Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/li_crc/paper_counts_E_xy_matrix.tsv -c tme_inputs/li_crc/gene_signatures.gmt -g tme_inputs/li_crc/paper_gold_std.tsv -o output/li_crc_results/ -p li_crc_paper_counts -t y -m CIBERSORT_BINARY -n CIBERSORT_BINARY -v y -a 0.4,1 -b 0,1 -s 100,100,0 -r 1-1 -y ranks -z 2000

#Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/li_crc/paper_cqn_E_xy_matrix.tsv -c tme_inputs/li_crc/gene_signatures.gmt -g tme_inputs/li_crc/paper_gold_std.tsv -o output/li_crc_results/ -p li_crc_paper_cqn -t y -m GSEA,GSVA,METANEIGHBOR_BINARY,ORA -n GSEA,GSVA,METANEIGHBOR_BINARY,ORA -v y -a 0.4,1 -b 0,1 -s 100,100,0 -r 1-1 -y ranks -z 2000
Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/li_crc/seurat_fqn_E_xy_matrix.tsv -c tme_inputs/li_crc/gene_signatures.gmt -g tme_inputs/li_crc/seurat_gold_std.tsv -o seurat_output/li_crc_results/ -p li_crc_seurat_fqn -t y -m GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY  -n GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY  -v y -a 0.4,1 -b 0,1 -s 100,100,0 -r 1-1 -y ranks -z 2000
