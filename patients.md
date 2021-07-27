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