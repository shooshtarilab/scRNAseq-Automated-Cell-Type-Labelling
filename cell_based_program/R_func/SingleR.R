run_SingleR<-function(Data,Labels,CV_RDataPath,OutputDir){
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])

  library(SingleR)
  library(Seurat)
  True_Labels_SingleR <- list()
  Pred_Labels_SingleR <- list()
  Total_Time_SingleR <- list()
  Data = t(as.matrix(Data))

  for (i in c(1:n_folds)){
    start_time <- Sys.time()
    singler = SingleR(method = "single", Data[,Test_Idx[[i]]], Data[,Train_Idx[[i]]], Labels[Train_Idx[[i]]]))
    end_time <- Sys.time()

    Total_Time_SingleR[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    True_Labels_SingleR[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_SingleR[i] <- list(as.vector(singler$labels))
  }
  True_Labels_SingleR <- as.vector(unlist(True_Labels_SingleR))
  Pred_Labels_SingleR <- as.vector(unlist(Pred_Labels_SingleR))
  Total_Time_SingleR <- as.vector(unlist(Total_Time_SingleR))

  setwd(OutputDir)

  write.csv(True_Labels_SingleR,'SingleR_True_Labels.csv',row.names = FALSE)
  write.csv(Pred_Labels_SingleR,'SingleR_Pred_Labels.csv',row.names = FALSE)
  write.csv(Total_Time_SingleR,'SingleR_Total_Time.csv',row.names = FALSE)
}
