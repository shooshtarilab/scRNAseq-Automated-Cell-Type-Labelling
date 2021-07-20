#!/bin/bash

#NOTE: this "short" example could still take hours to run

(time python program/runs_adobo.py ./data/chung_breast/counts.tsv.xz ./data/chung_breast/clusters.csv ./program/output/cb_adobo.tsv) 2> ./times/cb_adobo.txt

program/tme_inputs/avg_expr.ipynb


cd ./program/
export PERL5LIB=./bin/perl_modules/
Rscript bin/main_wrapper/subsamples_gene_classes_and_runs_enrichment_scripts.R -i tme_inputs/chung_breast/seurat_fqn_E_xy_matrix.tsv -c tme_inputs/chung_breast/gene_signatures.gmt -g tme_inputs/chung_breast/seurat_gold_std.tsv -o seurat_output/cb_results/ -p cb_seurat_fqn -t y -m GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -n GSEA,GSVA,METANEIGHBOR_BINARY,ORA,CIBERSORT_BINARY -v y -a 0.4,1 -b 0,1  -s 100,100,0 -r 1-1 -y ranks -z 2000
cd ..

(time Rscript ./program/sccatch/run_sccatch.R ./data/chung_breast/counts.tsv.xz ./data/chung_breast/clusters.csv 'seurat' 'Breast Cancer' 'Blood' 'Breast' 'Mammary gland') 2> ./times/cb_sccatch.txt

cd ./cell_based_program/
Rscript CV.R ../data/chung_breast/counts.tsv.xz ../data/chung_breast/ ../data/chung_breast/clusters.csv
python CV_r2py.py ../data/chung_breast/
cd ..
Rscript ./cell_based_program/R_methods.R chung_breast ./data/
python ./cell_based_program/Python_methods.py ./data/ chung_breast
python ./cell_based_program/run_LAmbDA.py ./data/ chung_breast
python ./cell_based_program/run_scVItool.py ./data/ chung_breast

#TODO how much of the scripts after this line require all datasets
# to be run?
cell_based_program/other_scripts/results_table.ipynb
cell_based_program/other_scripts/time_sim.py
cell_based_program/other_scripts/time_bar_plot.ipynb
cell_based_program/other_scripts/time_coefficient_plot.ipynb
cell_based_program/other_scripts/#TODO heatmap.ipynb

result_gathering.ipynb

#TODO bootstrap

#TODO i know these are hardcoded to require all datasets
sensitivity_plots.ipynb
main_figures.ipynb
#TODO don't have ping's scripts for subsampling
subsampling_all_cells/looking_to_automate_singletons.ipynb
rare_cell_types.ipynb
#TODO don't have ping's scripts for patient training
patient_data/predictions_results/score_patients.ipynb
supplementary_figures.ipynb
