.libPaths("/hpf/largeprojects/ccmbio/amahalanabis/tools/R/library")
suppressPackageStartupMessages(library(SingleCellExperiment))

F.Measure.Func <- function(C, K) {
  
  if (length(C) != length(K)) {
    print("Error. The length of vectors should be the same.")
    F.measure <- NA
  } else {
    C <- as.numeric(as.factor(C))
    K <- as.numeric(as.factor(K))
    
    n <- max(C)
    m <- max(K)
    
    a <- Pr <- Re <- Fm <- matrix(0, nrow=n, ncol=m)
    
    F.measure <- 0
    for (i in 1:n) {
      for (j in 1:m) {
        a[i,j] <- length(intersect(which(C == i), which(K == j)))
        
        if (a[i,j] == 0) { Fm[i,j]} else {
          Pr[i,j] <- a[i,j]/length(which(K == j))
          Re[i,j] <- a[i,j]/length(which(C == i))
          Fm[i,j] <- (2 * Pr[i,j] * Re[i,j])/(Pr[i,j] + Re[i,j])
        }
      }
      w <- length(which(C == i))/length(C)
      F.measure <- F.measure + w * max(Fm[i,])
    }
    F.vector <- apply(Fm, 1, max)
  }
  return(list(F.measure=F.measure, F.vector=F.vector))
}


ARI.Func <- function(C, K) {
  
  library("mclust")
  
  if (length(C) != length(K)) {
    print("Error. The length of vectors should be the same.")
    ARI <- NA
  } else {
    ARI <- adjustedRandIndex(C,K)
  }
  return(ARI)
}

VI.Func <- function(C, K) {
  library("mcclust")
  
  if (length(C) != length(K)) {
    print("Error. The length of vectors should be the same.")
    VI.measure <- NA
  } else {
    VI.measure <- vi.dist(C, K)
  }
  return(VI.measure)
}

Majority.Func <- function (cell, truth, alg) {
  
  if (length(C) != length(K)) {
    print("Error. The length of vectors should be the same.")
    F.measure <- NA
  } else {
  
    data <- data.frame(cell=cell, truth=as.numeric(as.factor(truth)), alg=alg)
    
    Majority <- 0
    Majority.vector <- c()
    
    num_clusters <- length(unique(data$alg))
    
    for (cluster in unique(data$alg)){
      tmp <- data[data$alg == cluster,]
      freqs <- prop.table(as.table(table(tmp$truth)))
      majority_pct = max(freqs)
      Majority.vector[cluster] = majority_pct
    }
  }
  
  Majority.vector <- na.omit(Majority.vector)
  Majority = median(Majority.vector)
  
  return(list(Majority=Majority, Majority.vector=Majority.vector))
}


AMI.Func <- function(C,K) {
  
  library("aricode")
  
  if (length(C) != length(K)) {
    print("Error. The length of vectors should be the same.")
    AMI <- NA
  } else {
    AMI <- AMI(C,K)
  }
  return(AMI)
}

F.Bootstrap.Func <- function(F.vector) {
  
  out <- unlist(lapply(as.list(1:10000), function(x) {
    mean(sample(F.vector, size=length(F.vector), replace=T))
  }))

  return(out)
}

Homogeneity.Func <- function(C, K) {
  library("infotheo")
  mi <- mutinformation(C, K)
  entropy.a <- entropy(C)
  entropy.b <- entropy(K)
  if (entropy.a == 0.0) {
    homogeneity <- 1.0
  } else {
    homogeneity <- mi / entropy.a
  }
  return(homogeneity)
}
