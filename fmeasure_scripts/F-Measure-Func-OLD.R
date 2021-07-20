F.Measure.Func <- function(C, K) {
  
  if (length(C) != length(K)) {
    print("Error. The length of vectors should be the same.")
    F.measure <- NA
  } else {
    # Here is the problem, they're indepentently converted to numerics
    # Which assumes that K is a subset of C.
    #C <- as.numeric(as.factor(C))
    #K <- as.numeric(as.factor(K))
    # c is the truth labels, k is predictions

    #TODO try using setdiff to remove labels in K that are not present in C
      # can't do this because the weights get messed up
      # need to keep the cells in the truth set but not the prediction set
    
    #n <- max(C)
    #m <- max(K)
    n <- length(unique(C))
    ni <- unique(C)
    m <- length(unique(K))
    mi <- unique(K) # check and see if i need this
    
    a <- Pr <- Re <- Fm <- Sp <- matrix(0, nrow=n, ncol=m)
    
    F.measure <- Sensitivity <- Specificity <- Precision <- 0
    for (i in 1:n) #check each truth label
    {
      for (j in 1:m) { # check each prediction
        #a[i,j] <- length(intersect(which(C == i), which(K == j)))
        a[i,j] <- length(intersect(which(C == ni[i]), which(K == ni[j]))) #true positives
          # select values where where the prediction is the same as the truth label
            # which(c==n) subsets the cells to only those that're actually a truth label
            # which(k==n) subsets the cell to only those that're predicted to be a truth label
            # which(k==n) subsets the cells to be those that're predicted to be a potentially different thing

          #TODO do the arrays n and m need to be the same?
            # do we even need a nested loop? we're only checking for cells 
            # types in the set of truth labels, other labels are unknowns.
        
        if (a[i,j] == 0) { Fm[i,j]} else {
          #TODO check the k==n... comparison, 
            # if we compare k==m, then we aren't getting positives, 
            # k==n is all values predicted to be a class
          Pr[i,j] <- a[i,j]/length(which(K == ni[j])) 
          Re[i,j] <- a[i,j]/length(which(C == ni[i]))
          Fm[i,j] <- (2 * Pr[i,j] * Re[i,j])/(Pr[i,j] + Re[i,j])
        }
      }
      w <- length(which(C == ni[i]))/length(C)
      F.measure <- F.measure + w * max(Fm[i,]) # weighted f-measure, w is proportion of total cells that is a particular cells by the truth labels
      Sensitivity <- Sensitivity + w * max(Re[i,]) # weighted f-measure, w is proportion of total cells that is a particular cells by the truth labels
      Precision <- Precision + w * max(Pr[i,]) # weighted f-measure, w is proportion of total cells that is a particular cells by the truth labels
    }
    F.vector <- apply(Fm, 1, max) # this is a per-cell type f-measure
  }
  return(list(F.measure=F.measure, F.vector=F.vector, sens=Sensitivity, prec=Precision))
}


# ARI.Measure.Func <- function(C, K) {
#   library("mclust")
#   
#   if (length(C) != length(K)) {
#     print("Error. The length of vectors should be the same.")
#     ARI.measure <- NA
#   } else {
#     ARI.measure <- adjustedRandIndex(C, K)
#   }
#   return(ARI.measure)
# }



F.Bootstrap.Func <- function(F.vector) {
  
  out <- unlist(lapply(as.list(1:10000), function(x) {
    mean(sample(F.vector, size=length(F.vector), replace=T))
  }))
  return(out)
}
