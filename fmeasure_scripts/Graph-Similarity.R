rm(list = ls())
Work.Dir <- "/home/erik/Documents/school/thesis/results/"
source(paste0(Work.Dir, "/fmeasure_scripts/F-Measure-Func.R"))

all.datasets <- c("cb","dg","jam","li_crc","llc","peng","tm","vg")

# Additional algorithms - To be added later: "bigScale", "raceid", "simlr"
#all.algorithms <- c("cibersort","gsea","gsva","metaneighbor","ora","adobo","sccatch")
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


low.CI <- high.CI <- f.measure <- cov.data <- cor.data <- f.measure.avg <- list()
for (dataset in all.datasets) {
  low.CI[[dataset]] <- high.CI[[dataset]] <- f.measure[[dataset]] <- cor.data[[dataset]] <- cov.data[[dataset]] <- f.measure.avg[[dataset]] <- matrix(0, nrow=length(all.algorithms), ncol=length(all.algorithms), dimnames=list(all.algorithms, all.algorithms))
  
  print(dataset)
  data.path <- paste0(Work.Dir,"/cellres/", dataset, "_paper_predictions.tsv")
  data <- read.csv(data.path, header=T, sep='\t')
  
  for (alg1 in all.algorithms) {
    for (alg2 in all.algorithms) {
      C <- data[, alg1]
      K <- data[, alg2]
      
      # Identify Mean F measures
      F.Measure.List <- F.Measure.Func(C, K)
      
      # Measure 95% confidence intervals using bootstrapping 
      F.vector <- F.Measure.List$F.vector
      CI <- quantile(F.Bootstrap.Func(F.vector), probs = c(0.025, 0.5, 0.975))
      low.CI[[dataset]][alg1, alg2] <- round(CI[1], 2)
      high.CI[[dataset]][alg1, alg2] <- round(CI[3], 2)
      f.measure[[dataset]][alg1, alg2] <- F.Measure.List$F.measure # round(CI[2], 2) #
      
      # cov.data[[dataset]][alg1, alg2] <- cov(C,K)
      # cor.data[[dataset]][alg1, alg2] <- cor(C,K, method = "spearman")
    }
  }
  
  # Measuring average F-Measure
  for (alg1 in all.algorithms) {
    for (alg2 in all.algorithms) {
      f.measure.avg[[dataset]][alg1, alg2] <- round((f.measure[[dataset]][alg1, alg2] + f.measure[[dataset]][alg2, alg1])/2, 2)
    }
  }
}


f.measure.total <- (f.measure.avg[[1]][all.algorithms, all.algorithms] +
                      f.measure.avg[[2]][all.algorithms, all.algorithms] +
                      f.measure.avg[[3]][all.algorithms, all.algorithms] +
                      f.measure.avg[[4]][all.algorithms, all.algorithms])/4


pdf(file=paste0(Work.Dir, "/Figures/Similar-Algorithms-Heatmap.pdf"), width=8, height=8)
heatmap(f.measure.total, keep.dendro=T, scale="none", margins = c(7, 7))
dev.off()


f.data <- 1-f.measure.total
hc <- hclust(d=dist(f.data), method = 'complete')

pdf(file=paste0(Work.Dir, "/Figures/Similar-Algorithms-Hierarchical.pdf"), width=8, height=8)
plot(hc)
dev.off()
memb <- cutree(hc, k = 8)