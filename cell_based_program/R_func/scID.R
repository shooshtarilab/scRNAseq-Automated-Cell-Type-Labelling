run_scID<-function(Data,Labels,CV_RDataPath,OutputDir){

  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])

  library(scID)
  library(Seurat)
  True_Labels_scID <- list()
  Pred_Labels_scID <- list()
  Total_Time_scID <- list()
  Data = t(as.matrix(Data))

  for (i in c(1:n_folds)){
    Train_Labels <- list(Labels[Train_Idx[[i]]])
    names(Train_Labels[[1]]) <- colnames(Data[,Train_Idx[[i]]])
    start_time <- Sys.time()
    scID_output <- scid_multiclass(Data[,Test_Idx[[i]]], Data[,Train_Idx[[i]]], Train_Labels[[1]])
    end_time <- Sys.time()

    Total_Time_scID[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    True_Labels_scID[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_scID[i] <- list(as.vector(scID_output$labels))
  }
  True_Labels_scID <- as.vector(unlist(True_Labels_scID))
  Pred_Labels_scID <- as.vector(unlist(Pred_Labels_scID))
  Total_Time_scID <- as.vector(unlist(Total_Time_scID))

  setwd(OutputDir)

  write.csv(True_Labels_scID,'scID_True_Labels.csv',row.names = FALSE)
  write.csv(Pred_Labels_scID,'scID_Pred_Labels.csv',row.names = FALSE)
  write.csv(Total_Time_scID,'scID_Total_Time.csv',row.names = FALSE)
}
