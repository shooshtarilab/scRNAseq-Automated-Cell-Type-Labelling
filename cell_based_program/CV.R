Cross_Validation <- function(LabelsPath, col_Index = 1,OutputDir){
  "
  Cross_Validation
  Function returns train and test indices for 5 folds stratified across unique cell populations,
  It return a 'CV_folds.RData' file which then used as input to classifiers wrappers.

  Parameters
  ----------
  LabelsPath : Cell population annotations file path (.csv).
  col_Index : column index (integer) defining which level of annotation to use,
  in case of multiple cell type annotations (default is 1)
  OutputDir : Output directory defining the path of the exported file.
  "

  Labels <- as.matrix(read.csv(LabelsPath))
  Labels <- as.vector(Labels[,col_Index])

  # Getting training and testing folds indices
  library(rBayesianOptimization)
  n_folds = 5
  Folds <- KFold(Labels,nfolds = n_folds, stratified = TRUE)
  Test_Folds <- c(n_folds:1)
  Train_Idx <- list()
  Test_Idx <- list()
  for (i in c(1:length(Folds))){
    Temp_Folds <- Folds
    Temp_Folds[Test_Folds[i]] <- NULL
    Train_Idx[i] <- list(unlist(Temp_Folds))
    Test_Idx[i] <- Folds[Test_Folds[i]]
  }
  remove(Temp_Folds,i,Folds)
  #setwd(OutputDir)
  save(n_folds,Train_Idx,Test_Idx,col_Index,file = file.path(OutputDir,'CV_folds.RData'))
}

args <- commandArgs(trailingOnly=TRUE)
main_dir <- args[1] #'path_to_dataset'
output_dir <- args[2] #'path_to_output_folder'
labels_path <- args[3] #path to labels file

#Cross_Validation(file.path(main_dir, 'Labels.csv'), 1, output_dir)
Cross_Validation(labels_path, 1, output_dir)
