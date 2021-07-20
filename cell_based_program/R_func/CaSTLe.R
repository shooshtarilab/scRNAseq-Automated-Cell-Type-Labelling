run_CaSTLe<-function(Data,Labels,CV_RDataPath, OutputDir){

  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])

  library(igraph)
  library(xgboost)
  True_Labels_Castle <- list()
  Pred_Labels_Castle <- list()
  Training_Time_Castle <- list()
  Testing_Time_Castle <- list()

  BREAKS=c(-1, 0, 1, 6, Inf)
  nFeatures = 100

  for(i in c(1:n_folds)){
    ds1 = Data[Train_Idx[[i]],]
    ds2 = Data[Test_Idx[[i]],]

    sourceCellTypes = as.factor(Labels[Train_Idx[[i]]])
    targetCellTypes = as.factor(Labels[Test_Idx[[i]]])

    start_time <- Sys.time()

    source_n_cells_counts = apply(ds1, 2, function(x) { sum(x > 0) } )
    target_n_cells_counts = apply(ds2, 2, function(x) { sum(x > 0) } )
    common_genes = intersect( colnames(ds1)[source_n_cells_counts>10],
                              colnames(ds2)[target_n_cells_counts>10])
    remove(source_n_cells_counts, target_n_cells_counts)
    ds1 = ds1[, colnames(ds1) %in% common_genes]
    ds2 = ds2[, colnames(ds2) %in% common_genes]
    ds = rbind(ds1[,common_genes], ds2[,common_genes])
    isSource = c(rep(TRUE,nrow(ds1)), rep(FALSE,nrow(ds2)))
    remove(ds1, ds2)

    topFeaturesAvg = colnames(ds)[order(apply(ds, 2, mean), decreasing = T)]
    end_time <- Sys.time()
    Training_Time_Castle[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    start_time <- Sys.time()

    L = length(levels(sourceCellTypes))
    targetClassification = as.data.frame(matrix(rep(0,L*sum(!isSource)), nrow=L), row.names = levels(sourceCellTypes))

    for (cellType in levels(sourceCellTypes)) {
      inSourceCellType = as.factor(ifelse(sourceCellTypes == cellType, cellType, paste0("NOT",cellType)))

      topFeaturesMi = names(sort(apply(ds[isSource,],2,function(x) { compare(cut(x,breaks=BREAKS),inSourceCellType,method = "nmi") }), decreasing = T))

      selectedFeatures = union(head(topFeaturesAvg, nFeatures) , head(topFeaturesMi, nFeatures) )

      tmp = cor(ds[,selectedFeatures], method = "pearson")
      tmp[!lower.tri(tmp)] = 0
      selectedFeatures = selectedFeatures[apply(tmp,2,function(x) any(x < 0.9))]
      remove(tmp)

      dsBins = apply(ds[, selectedFeatures], 2, cut, breaks= BREAKS)

      nUniq = apply(dsBins, 2, function(x) { length(unique(x)) })

      ds0 = model.matrix(~ . , as.data.frame(dsBins[,nUniq>1]))
      remove(dsBins, nUniq)

      cat(paste0("<h2>Classifier for ",cellType,"</h2>"))

      inTypeSource = sourceCellTypes == cellType

      xg=xgboost(data=ds0[isSource,] ,
                 label=inTypeSource,
                 objective="binary:logistic",
                 eta=0.7 , nthread=1, nround=20, verbose=0,
                 gamma=0.001, max_depth=5, min_child_weight=10)

      inTypeProb = predict(xg, ds0[!isSource, ])

      targetClassification[cellType,] = inTypeProb
    }
    end_time <- Sys.time()
    Testing_Time_Castle[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    True_Labels_Castle[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_Castle[i] <- list(rownames(targetClassification)[apply(targetClassification,2,which.max)])
  }
  True_Labels_Castle <- as.vector(unlist(True_Labels_Castle))
  Pred_Labels_Castle <- as.vector(unlist(Pred_Labels_Castle))
  Training_Time_Castle <- as.vector(unlist(Training_Time_Castle))
  Testing_Time_Castle <- as.vector(unlist(Testing_Time_Castle))

  setwd(OutputDir)

  write.csv(True_Labels_Castle,'CaSTLe_True_Labels.csv',row.names = FALSE)
  write.csv(Pred_Labels_Castle,'CaSTLe_Pred_Labels.csv',row.names = FALSE)
  write.csv(Training_Time_Castle,'CaSTLe_Training_Time.csv',row.names = FALSE)
  write.csv(Testing_Time_Castle,'CaSTLe_Testing_Time.csv',row.names = FALSE)
}
