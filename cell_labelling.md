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
    * **DEPENDENCIES IN [./cell_based_program/OtherMethods_packages.txt](./cell_based_program/OtherMethods_packages.txt)**
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
    * **DEPENDENCIES IN [./cell_based_program/scANVI_packages.txt](./cell_based_program/scANVI_packages.txt)**
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