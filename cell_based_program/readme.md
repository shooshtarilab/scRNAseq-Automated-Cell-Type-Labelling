General Usage
1. Generate input data files based on the following format
Data.csv.gz : cells X genes count matrix, barcodes as row names and gene names as column names
Labels.csv : Annotations file which contains cell type labels correspond to the cells in the count matrix. First row of the file contains the header.

2. Generate cross-validation file used in R scripts (CV_folds.RData) and Python scripts (CV_folds.pkl). CV.R script can be used to produce training and test indices (saved in .RData format) for cross validation in the R-based methods. For Python-based methods, CV_r2py.py can convert the .RData file to pickle file, which can be used in python.

3. Run R_methods.R for evaluating R-based methods and Python_methods.py for the evaluation of Python-based methods. Users have to revise the 'path_to_data_folder' and 'path_to_save_output' in these two files so that data files and prediction results can be correctly loaded. The input data files (Data.csv.gz, Labels.csv and CV_folds.RData/CV_folds.pkl) are assumed to be placed in the same folder (path_to_data_folder). Users need to revise the corresponding path if these files are stored in different locations.

Some other scripts
results_table.ipynb can combine all prediction result files to a single tsv table.
time_sim.py combines all prediction time files together and generates files for subsequent analysis
time_bar_plot.ipynb generates bar plots of running time showed in the supplementary
time_coefficient_plot.ipynb generates the heatmap of coefficients which reveals the importance of dataset dimensions


Notice:
1. Comment out the algorithms if you don't want to evaluate them.
2. LAmbDA and scANVI are implemented in different files.
3. Users need to install TensorFlow 1.8 to successfully run Cell_BLAST and LAmbDA. For scANVI (run_scVItool.py), the latest version of TensorFlow is preferred.
4. All scripts were implemented for cross-validation. Users have to revise them if they want to apply the algorithm to predict labels for a new dataset.
