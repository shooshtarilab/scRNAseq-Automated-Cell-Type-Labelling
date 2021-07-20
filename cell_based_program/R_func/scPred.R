run_scPred<-function(Data,Labels,CV_RDataPath, OutputDir){

  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])

  library(scPred)
  library(Seurat)
  library(magrittr)

  True_Labels_scPred <- list()
  Pred_Labels_scPred <- list()
  Training_Time_scPred <- list()
  Testing_Time_scPred <- list()
  Data = t(as.matrix(Data))

  for (i in c(1:n_folds)){
    obj.train <- CreateSeuratObject(counts=Data[,Train_Idx[[i]]])
    df <- as.data.frame(Labels[Train_Idx[[i]]], row.names=colnames(obj.train))
    obj.train <- AddMetaData(obj.train, metadata = df, col.name = 'cell_type')

    obj.test <- CreateSeuratObject(counts=Data[,Test_Idx[[i]]])

    obj.train <- obj.train %>%
      NormalizeData(verbose=F) %>%
      FindVariableFeatures(verbose=F) %>%
      ScaleData(verbose=F) %>%
      RunPCA(verbose=F)

    obj.test <- NormalizeData(obj.test, verbose=F)

    # Training
    start_time <- Sys.time()

    obj.train <- getFeatureSpace(obj.train, "cell_type")
    obj.train <- trainModel(obj.train)

    end_time <- Sys.time()
    Training_Time_scPred[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    # Prediction
    start_time <- Sys.time()
    obj.test <- scPredict(obj.test, obj.train)
    end_time <- Sys.time()

    Testing_Time_scPred[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    True_Labels_scPred[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_scPred[i] <- list(obj.test@meta.data$scpred_prediction)
  }

  True_Labels_scPred <- as.vector(unlist(True_Labels_scPred))
  Pred_Labels_scPred <- as.vector(unlist(Pred_Labels_scPred))
  Training_Time_scPred <- as.vector(unlist(Training_Time_scPred))
  Testing_Time_scPred <- as.vector(unlist(Testing_Time_scPred))

  setwd(OutputDir)

  write.csv(True_Labels_scPred,'scPred_True_Labels.csv',row.names = FALSE)
  write.csv(Pred_Labels_scPred,'scPred_Pred_Labels.csv',row.names = FALSE)
  write.csv(Training_Time_scPred,'scPred_Training_Time.csv',row.names = FALSE)
  write.csv(Testing_Time_scPred,'scPred_Testing_Time.csv',row.names = FALSE)
}
