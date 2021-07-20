run_CHETAH<-function(Data,Labels,CV_RDataPath,OutputDir){

  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])

  library(CHETAH)
  library(SingleCellExperiment)
  True_Labels_CHETAH <- list()
  Pred_Labels_CHETAH <- list()
  Total_Time_CHETAH <- list()
  Data = t(as.matrix(Data))

  for (i in c(1:n_folds)){
    sce <- SingleCellExperiment(assays = list(counts = Data[,Train_Idx[[i]]]),
                                colData = data.frame(celltypes = Labels[Train_Idx[[i]]]))

    sce_test <- SingleCellExperiment(assays = list(counts = Data[,Test_Idx[[i]]]),
                                     colData = data.frame(celltypes = Labels[Test_Idx[[i]]]))
    start_time <- Sys.time()
    sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce)
    end_time <- Sys.time()

    Total_Time_CHETAH[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    True_Labels_CHETAH[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_CHETAH[i] <- list(sce_test$celltype_CHETAH)
  }
  True_Labels_CHETAH <- as.vector(unlist(True_Labels_CHETAH))
  Pred_Labels_CHETAH <- as.vector(unlist(Pred_Labels_CHETAH))
  Total_Time_CHETAH <- as.vector(unlist(Total_Time_CHETAH))

  setwd(OutputDir)

  write.csv(True_Labels_CHETAH,'CHETAH_True_Labels.csv',row.names = FALSE)
  write.csv(Pred_Labels_CHETAH,'CHETAH_Pred_Labels.csv',row.names = FALSE)
  write.csv(Total_Time_CHETAH,'CHETAH_Total_Time.csv',row.names = FALSE)
}
