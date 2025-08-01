rm(list = ls())
Work.Dir <- "~/Documents/school/research/labelling_project_git/"
source(paste0(Work.Dir, "/fmeasure_scripts/F-Measure-Func.R"))

library(dplyr)
library(tidyr)
library(ggplot2)
library(tidytext)
library(janeaustenr)
library(forcats)

status <- "Bootstrap-Ensemble" #  "Actual" # 

##################################################################
#######                 Setting Parameters                ########
##################################################################
ref <- "truth"
all.datasets <- c("cb","dg","jam","li_crc","tm","llc","peng","vg")
##all.datasets <- c("cb","dg","jam","li_crc","llc","peng","tm")
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
#all.algorithms <- c("cibersort","LAmbDA")


##################################################################
#######               Calculating F Measure               ########
##################################################################
f.measure <- low.CI <- high.CI <- sens <- low.CI.sens <- high.CI.sens <- prec <- low.CI.prec <- high.CI.prec <- matrix(0, nrow=length(all.algorithms), ncol=length(all.datasets), dimnames=list(all.algorithms, all.datasets))

df <- c()
for (dataset in all.datasets) {
  print(dataset)
  #data.path <- paste0(Work.Dir,"/cellres/", dataset, "_seurat_predictions.tsv")
  data.path <- paste0(Work.Dir,"/predictions/", dataset, "_predictions.tsv")
  data <- read.csv(data.path, header=T, sep='\t')
  
  for (alg in all.algorithms) {
    print(alg)
    C <- data[, ref] #truth labels
    K <- data[, alg] #algorithm being tested
    
    N <- nrow(data)
    bootstrap.func <- function(x) {
      indx <- sample(1:N, size=N, replace=T)
      C1 <- C[indx]
      K1 <- K[indx]
      f <- F.Measure.Func(C1, K1)
      return(c(f$F.measure,f$sens, f$prec)) #This returns the F.measure not the vector from the function
      #return(F.Measure.Func(C1, K1)$F.measure) #This returns the F.measure not the vector from the function
    }
    
    # Identify F measures
    ##rslt <- unlist(lapply(1:10000, bootstrap.func))
      #lapply should return a list of three-tuple
      # do.call() will turn the list of three-tuples into the matrix like the for loop
    #TODO try do.call(rbind(lapply(),boots))
    rslt <- c()
    for (i in 1:10000){
      rslt <- rbind(rslt,bootstrap.func())
    }
    colnames(rslt) <- c('fm', 'sens','prec')

    # Measure 95% confidence intervals using bootstrapping
    CI <- quantile(rslt[,'fm'], probs = c(0.025, 0.5, 0.975))
    low.CI[alg, dataset] <- round(CI[1], 2)
    high.CI[alg, dataset] <- round(CI[3], 2)
    f.measure[alg, dataset] <- round(CI[2], 2)
    
    
    # Measure 95% confidence intervals using bootstrapping
    CI.prec <- quantile(rslt[,'prec'], probs = c(0.025, 0.5, 0.975))
    low.CI.prec[alg, dataset] <- round(CI.prec[1], 2)
    high.CI.prec[alg, dataset] <- round(CI.prec[3], 2)
    prec[alg, dataset] <- round(CI.prec[2], 2)
    
    # Measure 95% confidence intervals using bootstrapping
    CI.sens <- quantile(rslt[,'sens'], probs = c(0.025, 0.5, 0.975))
    low.CI.sens[alg, dataset] <- round(CI.sens[1], 2)
    high.CI.sens[alg, dataset] <- round(CI.sens[3], 2)
    sens[alg, dataset] <- round(CI.sens[2], 2)

    #x <- c(dataset, alg, f.measure[alg, dataset], low.CI[alg, dataset], high.CI[alg, dataset])
    x <- c(dataset, alg, 
           f.measure[alg, dataset], low.CI[alg, dataset], high.CI[alg, dataset],
           sens[alg, dataset], low.CI.sens[alg, dataset], high.CI.sens[alg, dataset],
           prec[alg, dataset], low.CI.prec[alg, dataset], high.CI.prec[alg, dataset]
          )
    df <- rbind(df, x)
  }
}

##################################################################
#######                 Ranking Algorithms                ########
##################################################################
score <- matrix(0, nrow=length(all.algorithms), ncol=length(all.datasets), dimnames=list(all.algorithms, all.datasets))
for (dataset in all.datasets) {
  score[, dataset] <- rank(f.measure[, dataset])
}

# avg.score <- sort(rowMeans(score), decreasing = T)
avg.score <- sort(apply(score, 1, median), decreasing = T)
rnk.score <- 1:length(all.algorithms)
names(rnk.score) <- names(avg.score)


##################################################################
#######             Preparing the Data Frames             ########
##################################################################
df <- as.data.frame(df)
#colnames(df) <- c("dataset", "algorithm", "F.Mean", "F.Low", "F.High")
colnames(df) <- c("dataset", "algorithm", 
                  "F.Mean", "F.Low", "F.High",
                  "Rec.Mean", "Rec.Low", "Rec.High",
                  "Prec.Mean", "Prec.Low", "Prec.High" 
                  #"Spec.Mean", "Spec.Low", "Spec.High"
                 )
df <- cbind(name=paste0(df[,"dataset"], "-", df[,"algorithm"]), df)
# df <- cbind(df, col=as.numeric(as.factor(df$algorithm)))
df <- cbind(df, Rank=rnk.score[df$algorithm])


df <- cbind(df, score=rep(0, nrow(df)))
for (i in 1:nrow(df)) {
  df[i, "score"] <- 1+length(all.algorithms) - score[df[i, "algorithm"], df[i, "dataset"]]
}

df$F.Mean <- as.numeric(df$F.Mean)
df$F.Low <- as.numeric(df$F.Low)
df$F.High <- as.numeric(df$F.High)
df$score <- as.numeric(df$score)


##################################################################
############             Generating Plots             ############
##################################################################
plt <- df %>%
  group_by(dataset) %>%
  top_n(length(all.algorithms)) %>%
  ungroup %>%
  mutate(dataset = as.factor(dataset),
         algorithm = reorder_within(algorithm, F.Mean, dataset)) %>%
  ggplot(aes(algorithm, F.Mean, fill = Rank)) +
  geom_errorbar( aes(x=algorithm, ymin=F.Low, ymax=F.High), width=0.4, colour="grey", alpha=0.9, size=1) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~dataset, scales = "free_y") +
  coord_flip() +
  scale_x_reordered() +
  scale_y_continuous(expand = c(0,0)) +
  labs(y = "F Measure",
       x = "Algorithms",
       title = "What are the top algorithms for each dataset?")


save(df, file=paste0(Work.Dir, "/Rdata/F-Measure-", status, "-vg.Rdata"))
write.table(df, file=paste0(Work.Dir, "/Rdata/F-Measure-", status, "-vg.tsv"), sep='\t')

pdf(paste0(Work.Dir, "/Figures/F-Measures-", status, "-vg.pdf"), width=12, height=6)
print(plt)
dev.off()

plt2 <- ggplot(df, aes(x = reorder(algorithm, -Rank), y = score, fill = Rank)) +
  geom_boxplot() +
  coord_flip() +
  labs(y = "Rank of Algorithm",
       x = "Algorithms")

pdf(paste0(Work.Dir, "/Figures/Ranks-", status, "-vg.pdf"), width=8, height=4)
print(plt2)
dev.off()

# n.cl<- sort(apply(data[,all.algorithms], 2, function(x) {length(unique(x))}))
# r <- score[all.algorithms, all.datasets[[2]]]
# pdf(paste0(Work.Dir, "/Figures/Fmeasure-NumCl-", status, ".pdf"), width=4, height=4)
# plot(n.cl, r)
# dev.off()
