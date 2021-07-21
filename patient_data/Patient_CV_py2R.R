# Convert the pickle file to RData format

library(reticulate)

main_dir <- 'path_to_pickle_file'
dataset <- 'path_to_output'

pd <- import("pandas")
pickle_data <- pd$read_pickle(paste(main_dir, '/', dataset, '/CV_data.pkl', sep=""))

n_folds <- pickle_data[['nfolds']]

col_Index <- 1

Train_Idx <- list()
Test_Idx <- list()

for (i in 0:(n_folds-1)) {
  Train_Idx[[i+1]] <- pickle_data[['train_ind']][[toString(i)]]
  Test_Idx[[i+1]] <- pickle_data[['test_ind']][[toString(i)]]
}

OutputDir <- paste(main_dir, '/', dataset, sep="")
save(n_folds,Train_Idx,Test_Idx,col_Index,file = file.path(OutputDir,'CV_data.RData'))
