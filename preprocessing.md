
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