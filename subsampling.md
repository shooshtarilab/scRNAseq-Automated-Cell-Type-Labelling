### Subsampling Experiment
These scripts subsample the datasets as outlined in our paper
* [subsampling_all_cells/sub_sample.py](subsampling_all_cells/sub_sample.py)
    * Creates subsampled datasets to use in the cell-based algorithms
    * After running this script, regenerate the CV folds from the subsampled 
    data and run the cell based algorithms again
* [subsampling_all_cells/singletons.ipynb](subsampling_all_cells/singletons.ipynb)
    * Evaluates the performance of all algorithms on the subsampled datasets vs 
    the original performance
    * Creates all plots related to this analysis seen in our paper