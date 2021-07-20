run_singleCellNet<-function(Data,Labels,CV_RDataPath,OutputDir){

  colnames(Data) <- gsub('_','.',colnames(Data), fixed = TRUE)
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])

  library(singleCellNet)
  library(dplyr)
  True_Labels_singleCellNet <- list()
  Pred_Labels_singleCellNet <- list()
  Training_Time_singleCellNet <- list()
  Testing_Time_singleCellNet <- list()
  Data = t(as.matrix(Data))

  for(i in c(1:n_folds)){
    DataTrain <- Data[,Train_Idx[[i]]]
    DataTest <- Data[,Test_Idx[[i]]]

    start_time <- Sys.time()
    cgenes2<-findClassyGenes(DataTrain, data.frame(Annotation = Labels[Train_Idx[[i]]]), "Annotation")
    cgenesA<-cgenes2[['cgenes']]
    grps<-cgenes2[['grps']]
    DataTrain<-as.matrix(DataTrain[cgenesA,])
    xpairs<-ptGetTop(DataTrain, grps, quickPairs = FALSE)#, ncores = 1)
    
    pdTrain<-query_transform(DataTrain[cgenesA, ], xpairs)
    rf<-sc_makeClassifier(pdTrain[xpairs,], genes=xpairs, groups=grps)
    end_time <- Sys.time()
    Training_Time_singleCellNet[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    start_time <- Sys.time()
    DataTest<-query_transform(DataTest[cgenesA,], xpairs)
    classRes <-rf_classPredict(rf, DataTest)
    end_time <- Sys.time()
    Testing_Time_singleCellNet[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    True_Labels_singleCellNet[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_singleCellNet[i] <- list((rownames(classRes)[apply(classRes,2,which.max)])[1:length(Test_Idx[[i]])])
  }

  True_Labels_singleCellNet <- as.vector(unlist(True_Labels_singleCellNet))
  Pred_Labels_singleCellNet <- as.vector(unlist(Pred_Labels_singleCellNet))
  Training_Time_singleCellNet <- as.vector(unlist(Training_Time_singleCellNet))
  Testing_Time_singleCellNet <- as.vector(unlist(Testing_Time_singleCellNet))

  setwd(OutputDir)

  write.csv(True_Labels_singleCellNet,'singleCellNet_True_Labels.csv',row.names = FALSE)
  write.csv(Pred_Labels_singleCellNet,'singleCellNet_Pred_Labels.csv',row.names = FALSE)
  write.csv(Training_Time_singleCellNet,'singleCellNet_Training_Time.csv',row.names = FALSE)
  write.csv(Testing_Time_singleCellNet,'singleCellNet_Testing_Time.csv',row.names = FALSE)
}
