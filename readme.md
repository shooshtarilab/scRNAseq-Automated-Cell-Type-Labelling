# Cell Type Labelling Analysis
These are scripts accompanying our paper on automated methods for labelling cell
types in scRNAseq data from the tumour microenvironment.

#TODO doi pending.

## Disclaimer
Many of these scripts require the entire analysis to be run first, since they 
were written not for a standalone pipeline, but to process the data for our 
paper. They can be re-run to replicate the results, but this can take an 
extremely long time.

As such, we have included the inputs and outputs for the most relevant stages of
the pipeline, hopefully facilitating the ease of review. To generate the figures
in our paper, simply clone the repo, install python, jupyter, pandas, numpy,
scipy, scikit-learn, matplotlib, and seaborn, then run main_figures.ipynb,
singletons.ipynb, underrepresented_cell_types.ipynb, score_patients.ipynb, and
supplemental_figures.ipynb 

Most of the jupyter notebooks have their last run saved, so the outputs of the 
code can be evaluated without needing to run the cells again and also target 
relative file paths, allowing users to clone the repo and run the notebooks 
without having to modify file paths within the notebook (provided they are run 
on Linux or macOS since all file paths use "/" as a directory separator.)

run_order.sh contains more information on what files are input to each script,
and run_short_example.sh has examples on how to call some scripts for one 
dataset.

### Preprocessing Scripts
* [./program/runs_adobo.py](./program/runs_adobo.py)
    * Runs adobo predictions and normalises the data using FQN
    * Inputs: `cluster.csv` and `counts.tsv.xz` file for each dataset 
    * Outputs: Adobo predictions file and `fqn.tsv.xz` (fqn normalised expression)
    * Arguments: The path to the counts file, the path to the clusters file, and
     the output path for the predictions
* [./program/tme_inputs/avg_expr.ipynb](./program/tme_inputs/avg_expr.ipynb)
    * Reads in the counts and fqn gene x cell matrices and computes the average 
    expression per cluster
    * Inputs: `counts.tsv.xz`, `fqn.tsv.xz` and `clusters.csv` for each dataset
    * Outputs: Average expression per cluster in raw counts and fqn, gold 
    standards for the clusters
    * Arguments: See python notebook - paths to the input and output files

### Cluster labelling Algorithms
* [./program/run.sh](./program/run.sh)
    * Runs GSEA, GSVA, ORA, MetaNeighbor, and CIBERSORT
    * Inputs: Average expression per cluster for all datasets, gold standards 
    for clusters, and cell type signature gene sets
    * Outputs: Run data from each algorithm, including prediction files
        * These outputs are meant to be post-processed by a later script
    * The scripts called in here are not my own, and can be found 
    [here](https://github.com/jdime/scRNAseq_cell_cluster_labeling), along with
    more extensive documentation.
    * Note that I cannot distribute GSEA and CIBERSORT, so you must download 
    the jar files for each method and place them in `./program/`
        * [CIBERSORT](https://cibersort.stanford.edu/download.php)
        * [GSEA](http://software.broadinstitute.org/gsea/downloads.jsp)
* [./program/sccatch/run_sccatch.R](./program/sccatch/run_sccatch.R)
    * Runs scCATCH predictions
    * Inputs: `counts.tsv.xz` and `clusters.csv` for each dataset
    * Outputs: Prints scCATCH predictions to console
    * Arguments: 
        * Path to the genes x cell counts file 
        * Path to the file containing seurat clusters for each cell
        * Column name to read clusters from in the clusters file ('seurat')
        * The cancer type of the dataset
        * A list of tissues available for the dataset as documented by scCATCH
    * Examples for some of the datasets are present in `./program/sccatch/*.sh`

### Cell Based Labelling Algorithms
* [./cell_based_program/CV.R](./cell_based_program/CV.R)
    * Generates the 5-fold cross validation groups as .Rdata files
    * Inputs: A file containing labels for each cell
    * Outputs: An Rdata file containing cross validation groups
    * Arguments
        * Path to the genes x cells counts file
        * Path to the folder in which to save the Rdata file
        * Path to the folder containing the cell labels
* [./cell_based_program/CV_r2py.py](./cell_based_program/CV_r2py.py)
    * Converts the Rdata files from CV.R to pkl files for the python methods
    * Inputs: The Rdata file
    * Outputs: The Pkl file, in the same directory as the Rdata file
    * Arguments: Path to the folder containing the Rdata files  
* [./cell_based_program/R_methods.R](./cell_based_program/R_methods.R)
    * Runs all R based methods on a single dataset
    * Inputs: A counts matrix, labels file, and Rdata CV file for the dataset
    * Outputs: Predictions in ./cell_based_program/output/\<dataset\>
    * Arguments:
        * The name of the dataset
        * The path to the folder containing the counts matrix, labels, and Rdata
* [./cell_based_program/Python_methods.py](./cell_based_program/Python_methods.py)
    * Runs all python based methods on a single dataset
        * Lambda and scVI have their own scripts
    * Inputs: A counts matrix, labels file, and pkl CV file for the dataset
    * Outputs: Predictions in ./cell_based_program/output/\<dataset\>
    * Arguments:
        * The path to the folder containing the counts, labels, and pkl
        * The name of the dataset
    * **DEPENDENCIES IN ./cell_based_program/OtherMethods_packages.txt**
* [./cell_based_program/run_LAmbDA.py](./cell_based_program/run_LAmbDA.py)
    * Runs LAmbDA on a single dataset
    * Inputs: A counts matrix, labels file, and pkl CV file for the dataset
    * Outputs: Predictions in ./cell_based_program/output/\<dataset\>
    * Arguments:
        * The path to the folder containing the counts, labels, and pkl
        * The name of the dataset
* [./cell_based_program/run_scVItool.py](./cell_based_program/run_scVItool.py)
    * Runs scVI on a single dataset
    * Inputs: A counts matrix, labels file, and pkl CV file for the dataset
    * Outputs: Predictions in ./cell_based_program/output/\<dataset\>
    * Arguments:
        * The path to the folder containing the counts, labels, and pkl
        * The name of the dataset
    * **DEPENDENCIES IN ./cell_based_program/scANVI_packages.txt**
        * This requires a different version of tensorflow, it is recommended
        to use a virtual environment to manage it.
* [./cell_based_program/other_scripts/results_table.ipynb](./cell_based_program/other_scripts/results_table.ipynb)
    * Gathers the predictions from all methods into a single file
* [./cell_based_program/other_scripts/time_sim.py](./cell_based_program/other_scripts/time_sim.py)
    * Gathers the running times into a single file
* [./cell_based_program/other_scripts/time_bar_plot.ipynb](./cell_based_program/other_scripts/time_bar_plot.ipynb)
    * Makes bar plots of the running time as seen in our supplemental figures
* [./cell_based_program/other_scripts/time_coefficient_plot.ipynb](./cell_based_program/other_scripts/time_coefficient_plot.ipynb)
    * Makes the dataframe needed for the time coefficient plot generated in 
    main_figures.ipynb
* [./cell_based_program/other_scripts/heatmap.ipynb](./cell_based_program/other_scripts/heatmap.ipynb)
    * Makes the dataframe needed for the running time heatmap generated in 
    main_figures.ipynb

## **NOTE THAT ALL FOLLOWING SCRIPTS ASSUME THAT ALL ALGORITHMS HAVE BEEN RUN FOR ALL DATASETS AND MAY NOT WORK WITHOUT ALL DATA PRESENT** 
### Post-processing
* [result_gathering.ipynb](result_gathering.ipynb)
    * Gathers the predictions from the cluster labelling methods and combines 
    them with the predictions from the cell based methods
    * **NOTE** Because adobo and scCATCH are semi-supervised and use their own 
    gene signature database, their predictions may not exactly match what is 
    present in the gold standards (for instance we may only know for sure that a
    cell is a T Cell, but it may be predicted to be a memory T Cell or helper T 
    Cell). I could not find a way to make these match automatically, so the 
    labels are renamed manually prior to running this script
* [fmeasure_scripts/F-Measure-Analysis-Bootstrap.R](fmeasure_scripts/F-Measure-Analysis-Bootstrap.R)
    * Obtains bootstrapped F1 Scores for each algorithm on each dataset
    * **Disclaimer** This takes a _long_ time to run. Multiple hours, especially
    if it is not modified to exlude some algorithms or the largest datasets.
    * Inputs: Prediction files for each dataset
    * Ouputs: Large dataframe containing bootstrapped F1 scores for each dataset
* [score_all_methods.ipynb](score_all_methods.ipynb)
    * Obtains per-class F1 score among other metrics for all datasets and 
    algorithms, gathered into a single large dataframe
    * Inputs: Prediction files for each dataset
    * Ouputs: Large dataframe containing per-class as well as macro and micro 
    averaged scores for each dataset-algorithm pair
* [main_figures.ipynb](main_figures.ipynb)
    * Generates most of the main figures for our paper
    * Inputs:
        * data_sizes.tsv: A hardcoded dataframe outlining aspects of the 
        datasets for figure 1
        * Rdata/F-Measure-Bootstrap-Ensemble.tsv: A file containing bootstrapped
        F1 scores from all methods. Output from F-Measure-Analysis-Bootstrap.R
        * times/df_for_heatmap.tsv: A dataframe with running times output from 
        an earlier script
        * times/df_coef.tsv: A dataframe with regression coefficients from
        ./cell_based_program/other_scripts/time_coefficient_plot.ipynb
        * performance/seurat/bigdf.tsv: The large dataframe output 
        from score_all_methods.ipynb
    * Outputs: plots for the paper under ./plots
* [underrepresented_cell_types.ipynb](underrepresented_cell_types.ipynb)
    * Makes the main and supplemental underrepresented cell types plots
    * Inputs: performance/seurat/bigdf.tsv: The large dataframe output from 
    score_all_methods.ipynb
    * Outputs: The two underrepresented cell types plots

### Subsampling Experiment
These scripts subsample the datasets as outlined in our paper
* [subsampling_all_cells/sub_sample.py](subsampling_all_cells/sub_sample.py)
    * Creates subsampled datasets to use in the cell-based algorithms
    * After running this script, regenerate the CV folds from the subsampled 
    data and run the cell based algorithms again
* [subsampling_all_cells/looking_to_automate_singletons.ipynb](subsampling_all_cells/looking_to_automate_singletons.ipynb)
    * Evaluates the performance of all algorithms on the subsampled datasets vs 
    the original performance
    * Creates all plots related to this analysis seen in our paper

### Patient training experiment
These scripts create new cross validation groups with cells from each patient 
always grouped together and trains/tests the algorithms on these groups
* [patient_data/Patient_CV.ipynb](./patient_data/Patient_CV.ipynb)
    * Generate the CV folds with patients grouped together
* [patient_data/Patient_CV_py2R.R](./patient_data/Patient_CV_py2R.R)
    * Convert the pkl folds to Rdata
    * After running this script, re-run all cell-based algorithms
* [patient_data/predictions_results/score_patients.ipynb](patient_data/predictions_results/score_patients.ipynb)
    * Reads the predictions from the patient data as well as the patient groups
    and creates the plots for our paper

### Supplemental figures
* supplementary_figures.ipynb generates all supplemental figures that are not 
generated by other scripts.
