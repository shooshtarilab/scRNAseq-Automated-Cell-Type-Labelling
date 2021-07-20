
adjustedRandIndex <- function (x, y) 
{
  x <- as.vector(x) # the truth set
  y <- as.vector(y) # the test set
  if(length(x) != length(y)) 
    stop("arguments must be vectors of the same length")
  tab <- table(x,y)
  if(all(dim(tab)==c(1,1))) return(1)
  a <- sum(choose(tab, 2))              # number of pairs of elements in the same subset in x and y
  b <- sum(choose(rowSums(tab), 2)) - a # number of pairs of elements in different subsets in x and y
  c <- sum(choose(colSums(tab), 2)) - a # number of pairs of elements in the same subset in x and different in y
  d <- choose(sum(tab), 2) - a - b - c  # number of pairs of elements in different subsets in x and same in y

  #TODO how should this change to reflect that items must be given a particular label to be correct?
  # i think i just need to make sure these are calculated correctly.
    # a = TP
    # b = TN
    # c = FN
    # d = FP

  ARI <- (a - (a + b) * (a + c)/(a + b + c + d)) /
    ((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}
