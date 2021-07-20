suppressPackageStartupMessages(library(SingleCellExperiment))
source("Functions.R")

library(dplyr)
library(tidyr)
library(ggplot2)
library(tidytext)
library(janeaustenr)
library(forcats)

##################################################################
#######                 Setting Parameters                ########
##################################################################
ref <- "truth"
all.datasets <- c("cb","dg","jam","li_crc","tm","llc","peng","vg")
# Additional algorithms - To be added later: "bigScale", "raceid", "simlr"
all.algorithms <- c("cibersort",
                    "gsea",
                    "gsva",
                    "metaneighbor",
                    "ora",
                    "adobo",
                    "sccatch",
                    "SVM",
                    "SVMrej",	
                    "RF",	
                    "LDA",	
                    "LDArej",	
                    "NMC",	
                    "kNN9",	
                    "ACTINN",	
                    "scVI",	
                    "Cell_BLAST",	
                    "SingleCellNet",
                    "LAmbDA",	
                    "scPred",	
                    "CaSTLe",	
                    "CHETAH",	
                    "scID",	
                    "scmapcell",	
                    "scmapcluster",	
                    "singleR"
                   )
all.algorithms <- c('cibersort', 'SVM')
##################################################################
#######               Calculating F Measure               ########
##################################################################
ARI <- low.CI <- high.CI <- matrix(0, nrow=length(all.algorithms), ncol=length(all.datasets), dimnames=list(all.algorithms, all.datasets))

df <- c()
for (dataset in all.datasets) {
  print(dataset)
  data.path <- paste0("../cellres/", dataset, "_paper_predictions.tsv")
  data <- read.csv(data.path, header=T,sep='\t')
  
  for (alg in intersect(all.algorithms, colnames(data))) {
    print(alg)
    C <- data[, ref]
    K <- data[, alg]
    
    N <- nrow(data)
    bootstrap.func <- function(x) {
      indx <- sample(1:N, size=N, replace=T)
      C1 <- C[indx]
      K1 <- K[indx]
      return(ARI.Func(C1, K1))
    }
    
    # Identify ARI measures
    rslt <- unlist(lapply(1:10000, bootstrap.func))
    CI <- quantile(rslt, probs = c(0.025, 0.5, 0.975))
    low.CI[alg, dataset] <- round(CI[1], 2)
    high.CI[alg, dataset] <- round(CI[3], 2)
    ARI[alg, dataset] <- round(CI[2], 2)
    x <- c(dataset, alg, ARI[alg, dataset], low.CI[alg, dataset],  high.CI[alg, dataset])
    df <- rbind(df, x)
  }
  colnames(df) <- c("dataset", "algorithm", "ARI", "ARI_low_CI", "ARI_high_CI")
  write.table(df, file="../Rdata/ARI_bootstrapping.tsv", row.names = FALSE, sep='\t') 
}

df <- as.data.frame(df)

##################################################################
#######         Calculating Scores and Ranks              ########
##################################################################

#CALCULATE SCORE WHICH IS RANK - NUMBER OF ALGORITHMS
df$score = 1
for (dataset in all.datasets) {
  df[df$dataset == dataset,"score"] <- rank(df[df$dataset == dataset,"ARI"])
}

for (i in 1:nrow(df)) {
  current_dataset <- df[i,1]
  num_algorithms <- nrow(df[df$dataset == current_dataset,])
  df[i, "score"] <- 1 + num_algorithms - df[i, "score"]
}

#ASSIGN EACH ALGORITHM A RANK WHICH IS THE MEDIAN SCORE
df$rank = 1
for (alg in all.algorithms) {
  median.score <- median(df[df$algorithm == alg, "score"])
  df[df$algorithm == alg, "rank"] = median.score
}

write.table(df, file="../Rdata/ARI_bootstrapping.tsv", row.names = FALSE, sep='\t')


