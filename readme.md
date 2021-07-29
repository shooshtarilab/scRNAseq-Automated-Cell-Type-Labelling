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
the pipeline, hopefully facilitating the ease of review.

Most of the jupyter notebooks have their last run saved, so the outputs of the 
code can be evaluated without needing to run the cells again and also target 
relative file paths, allowing users to clone the repo and run the notebooks 
without having to modify file paths within the notebook (provided they are run 
on Linux or macOS since all file paths use "/" as a directory separator.)

run_order.sh contains more information on what files are input to each script,
and run_short_example.sh has examples on how to call some scripts for one 
dataset.

Each section listed below is modular, and the inputs needed are provided so
that they can be run independent of the others. You could for instance run
the pre-processing and data science/plot generation modules without re-running
the entire predictions pipeline (as the prediction pipeline takes a very long
time)

## Replicating All Data Science
Since running the complete prediction pipeline for this project takes a large
amount of compute resources (weeks of running time and hundreds of gigabytes of
memory), we have provided as many of our results and inputs as possible.

In order to re-run only the post-processing and data-science results, users can
simply clone the repo, install the listed pacakges, and run the following
notebooks/scripts. (To run on windows you will need to replace "/" with "\" in
file paths)
* [main_figures.ipynb](main_figures.ipynb)
* [singletons.ipynb](./subsampling_all_cells/singletons.ipynb)
* [underrepresented_cell_types.ipynb](underrepresented_cell_types.ipynb)
* [score_patients.ipynb](./patient_data/predictions_results/score_patients.ipynb)
* [supplemental_figures.ipynb](supplemental_figures.ipynb)
* Dependencies (I've listed the version originally used, but newer releases
should work just as well)
    * Python 3.8
    * jupyter 1.0.0
    * pandas 1.1.3
    * numpy 1.18.5
    * scipy 1.5.3
    * scikit-learn 0.23.2
    * matplotlib 3.3.2
    * seaborn 0.11.0

--------------------------------------------------------------------------------
### [Preprocessing Scripts](./preprocessing.md)
* These scripts prepare inputs for the cluster labelling algorithms, namely by
normalising the data and getting average expression for each cluster

### [Cluster labelling Algorithms](./cluster_labelling.md)
* These scripts run cluster labelling algorithms in order to predict cell types
in each cluster
* You must run the preprocessing scripts first

### [Cell Based Labelling Algorithms](./cell_labelling.md)
* These scripts prepare and run the cell based labelling algorithms
* No scripts must be run first, but there are additional dependencies outlined
in the linked readme. You must use a virtual environment since two different
tensorflow versions are required

## **NOTE THAT ALL FOLLOWING SCRIPTS ASSUME THAT ALL ALGORITHMS HAVE BEEN RUN FOR ALL DATASETS AND MAY NOT WORK WITHOUT ALL DATA PRESENT** 
### [Post-processing](./post-processing.md)
* These scripts gather all of the predictions into a single file, score them
against the gold standards, and plot the main results
* These scripts rely on ALL previous scripts having been run.

### [Subsampling Experiment](./subsampling.md)
* These scripts subsample the datasets and analyse the results as outlined in
our paper
* There are no precursor scripts besides what is outlined in the readme, but all
cell based algorithms must be re-run multiple times for each dataset as part of
this analysis

### [Patient training experiment](./patients.md)
* These scripts create new cross validation groups with cells from each patient 
always grouped together and trains/tests the algorithms on these groups
* There are no precursor scripts besides what is outlined in the readme, but all
cell based algorithms must be re-run for each dataset as part of this analysis


### [Supplemental figures](./supplemental.md)
* These scripts generate the supplemental figures included in our paper.
* This relies on all previous analysis having been run.
    * We have provided inputs from our runs so that this can be re-run easily

## Package Versions
* Python Version 3.8
    * jupyter 1.0.0
    * pandas 1.1.3
    * numpy 1.18.5
    * scipy 1.5.3
    * scikit-learn 0.23.2
    * matplotlib 3.3.2
    * seaborn 0.11.0
    * adobo
* R Version 3.5.2
    * optparse
    * vioplot
    * GSA
    * data.table
    * precrec
    * ROCR
    * Seurat
    * dplyr
    * Rserve
    * e1071
    * colorRamps
    * stats
    * GSVA (BioConductor)
    * qvalue (BioConductor)
    * preprocessCore (BioConductor)
* Perl Version 5
    * Date::Calc
* Java
    * [CIBERSORT](https://cibersort.stanford.edu/download.php)
    * [GSEA 3.0](http://software.broadinstitute.org/gsea/downloads.jsp)
