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