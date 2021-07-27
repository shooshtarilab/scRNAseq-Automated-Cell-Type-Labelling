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