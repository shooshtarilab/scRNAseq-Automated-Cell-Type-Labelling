F.Measure.Func <- function(C, K) {
  
  if (length(C) != length(K)) {
    print("Error. The length of vectors should be the same.")
    F.measure <- NA
  } else {
    # c is the truth labels, k is predictions
    n <- length(unique(C))
    ni <- unique(C)
    m <- length(unique(K))
    
    a <- Pr <- Re <- Fm <- Sp <- matrix(0, nrow=n, ncol=n)
    
    F.measure <- Sensitivity <- Specificity <- Precision <- 0
    #TODO mention that changing this to a single loop nearly exactly matches sklearn's classification_report weighted avg
    for (i in 1:n) #check each truth label
    {
      j <- i
      a[i,j] <- length(intersect(which(C == ni[i]), which(K == ni[i]))) #true positives
        # select values where where the prediction is the same as the truth label
          # which(c==n) subsets the cells to only those that're actually a truth label
          # which(k==n) subsets the cell to only those that're predicted to be the same truth label
      
      if (a[i,j] == 0) { Fm[i,j]} else {
        Pr[i,j] <- a[i,j]/length(which(K == ni[i])) 
        Re[i,j] <- a[i,j]/length(which(C == ni[i]))
        Fm[i,j] <- (2 * Pr[i,j] * Re[i,j])/(Pr[i,j] + Re[i,j])
      }
      w <- length(which(C == ni[i]))/length(C)
      F.measure <- F.measure + w * max(Fm[i,]) # weighted f-measure, w is proportion of total cells that is a particular cells by the truth labels
      Sensitivity <- Sensitivity + w * max(Re[i,]) # weighted f-measure, w is proportion of total cells that is a particular cells by the truth labels
      Precision <- Precision + w * max(Pr[i,]) # weighted f-measure, w is proportion of total cells that is a particular cells by the truth labels
    }
    F.vector <- apply(Fm, 1, max) # this is a per-cell type f-measure
  }
  return(list(F.measure=F.measure, F.vector=F.vector, sens=Sensitivity, prec=Precision))
  #TODO this does not always give f.measure between sens and prec, but it also seems dependent on result
  # when stepping through the code, PER CLASS f measure seems to always be between recall and precision
  # as the values are weighted, f-measure ends up being unbounded by the two
    # i think this means we had the right idea (i.e. the reason it is unbounded 
    #  is because sometimes prec<rec sometimes rec<prec) but that was occuring at a different step.
    #  Some classes prec<rec, others rec<prec.
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
