.libPaths("/hpf/largeprojects/ccmbio/amahalanabis/tools/R/library")
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
all.datasets <- c("darmanis_glioblastoma_allCells", "JA_melanoma_allCells", "li_crc_allCells", 
                  "tirosh_melanoma_allCells", "chung_breast_allCells",
                  "darmanis_glioblastoma_nonTumour", "JA_melanoma_nonTumour", "li_crc_nonTumour", 
                  "tirosh_melanoma_nonTumour", "chung_breast_nonTumour",
                  "peng_pancreatic_allCells", "peng_pancreatic_nonTumour",
                  "vangalen_AML_allCells", "vangalen_AML_nonTumour", 
                  "lambrechts_lung_allCells", "lambrechts_lung_nonTumour")

# Additional algorithms - To be added later: "bigScale", "raceid", "simlr"
all.algorithms <- c("altAnalyze",	"ascend", "backspin", "bigScale",
                    "cellranger",	"cidr",	"countClust",	"monocle",	"pcaReduce",
                    "phenograph",	"raceid", "rca",	"sc3",	"scran",	"seurat",
                    "simlr", "sincera",	"tscan")
##################################################################
#######               Calculating F Measure               ########
##################################################################
homogeneity <- low.CI <- high.CI <- matrix(0, nrow=length(all.algorithms), ncol=length(all.datasets), dimnames=list(all.algorithms, all.datasets))

df <- c()
for (dataset in all.datasets) {
  print(dataset)
  data.path <- paste0("./Clustering_Results/", dataset, ".csv")
  data <- read.csv(data.path, header=T)
  
  for (alg in intersect(all.algorithms, colnames(data))) {
    print(alg)
    C <- data[, ref]
    K <- data[, alg]
    
    N <- nrow(data)
    cell <- data[,"cell"]
    bootstrap.func <- function(x) {
      indx <- sample(1:N, size=N, replace=T)
      C1 <- C[indx]
      K1 <- K[indx]
      cell1 <- cell[indx]
      return(Homogeneity.Func(C1,K1))
    }
    
    # Identify Homogeneity measures
    rslt <- unlist(lapply(1:10000, bootstrap.func))
    
    # Measure 95% confidence intervals using bootstrapping
    CI <- quantile(rslt, probs = c(0.025, 0.5, 0.975))
    low.CI[alg, dataset] <- round(CI[1], 2)
    high.CI[alg, dataset] <- round(CI[3], 2)
    print (CI)
    print(low.CI[alg, dataset])
    print(high.CI[alg,dataset])
    homogeneity[alg, dataset] <- round(CI[2], 2)
    
    x <- c(dataset, alg, homogeneity[alg, dataset], low.CI[alg, dataset], high.CI[alg, dataset])
    df <- rbind(df, x)
  }
  write.csv(df, file="dataframe_results/Homogeneity_bootstrapping.csv", row.names = FALSE) 
}

colnames(df) <- c("dataset", "algorithm", "metric_mean", "metric_low_CI", "metric_high_CI")
df$metric = "Homogeneity"
write.csv(df, file="Homogeneity_bootstrapping.csv", row.names = FALSE)
df <- read.csv("Homogeneity_bootstrapping.csv", stringsAsFactors = FALSE)
colnames(df) <- c("dataset", "algorithm", "metric_mean", "metric_low_CI", "metric_high_CI")
df$metric = "homogeneity"

for (i in 1:nrow(df)) {
  split_name = strsplit(df[i, "dataset"], "_")
  df[i, "dataset_name"] <- paste0(split_name[[1]][1], "_", split_name[[1]][2])
  df[i, "cell_type"] <- split_name[[1]][3]
}

##################################################################
#######         Calculating Scores and Ranks              ########
##################################################################

#CALCULATE SCORE WHICH IS RANK - NUMBER OF ALGORITHMS
df$score = 1
for (dataset in all.datasets) {
  df[df$dataset == dataset,"score"] <- rank(df[df$dataset == dataset,"metric_mean"])
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

##################################################################
#######  Generate Score Plots coloured by algorithm rank  ########
##################################################################

for (i in 1:nrow(df)) {
  split_name = strsplit(df[i, "dataset"], "_")
  df[i, "dataset_name"] <- paste0(split_name[[1]][1], "_", split_name[[1]][2])
  df[i, "cell_type"] <- split_name[[1]][3]
}
df <- cbind(name=paste0(df[,"metric"], "-", df[,"cell_type"]), df)
write.csv(df, file="dataframe_results/homogeneity_bootstrapping.csv", row.names = FALSE)

# 
# for (i in 1:nrow(df)) {
#   split_name = strsplit(df[i, "dataset"], "_")
#   df[i, "dataset_name"] <- paste0(split_name[[1]][1], "_", split_name[[1]][2])
#   df[i, "cell_type"] <- split_name[[1]][3]
# }
# 
# ##################################################################
# #######         Calculating Scores and Ranks              ########
# ##################################################################
# 
# #CALCULATE SCORE WHICH IS RANK - NUMBER OF ALGORITHMS
# df$score = 1
# for (dataset in all.datasets) {
#   df[df$dataset == dataset,"score"] <- rank(df[df$dataset == dataset,"ARI"])
# }
# 
# for (i in 1:nrow(df)) {
#   current_dataset <- df[i,1]
#   num_algorithms <- nrow(df[df$dataset == current_dataset,])
#   df[i, "score"] <- 1 + num_algorithms - df[i, "score"]
# }
# 
# #ASSIGN EACH ALGORITHM A RANK WHICH IS THE MEDIAN SCORE
# df$rank = 1
# for (alg in all.algorithms) {
#   median.score <- median(df[df$algorithm == alg, "score"])
#   df[df$algorithm == alg, "rank"] = median.score
# }
# 
# tmp <- read.csv("dataframe_results/AMI_bootstrapping.csv", stringsAsFactors = FALSE)
# df <- rbind(df, tmp)
# df <- cbind(name=paste0(df[,"metric"], "-", df[,"cell_type"]), df)
