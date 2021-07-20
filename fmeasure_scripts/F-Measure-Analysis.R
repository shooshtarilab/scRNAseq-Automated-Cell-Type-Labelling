rm(list = ls())
Work.Dir <- "/home/erik/Documents/school/thesis/results/"
source(paste0(Work.Dir, "/fmeasure_scripts/F-Measure-Func.R"))

library(dplyr)
library(tidyr)
library(ggplot2)
library(tidytext)
library(janeaustenr)
library(forcats)

status <- "Actual-Ensemble" # "Median" # 

##################################################################
#######                 Setting Parameters                ########
##################################################################
ref <- "truth"
#all.datasets <- c("darmanis_glioblastoma", "JA_melanoma", "li_crc", "tirosh_melanoma")
all.datasets <- c("cb","dg","jam","li_crc","llc","peng","tm","vg")
# Additional algorithms - To be added later: "bigScale", "raceid", "simlr"
#all.algorithms <- c("altAnalyze",	"ascend",	"cellRanger",	"cidr",	"countClust",	"monocle",	"pcaReduce",	"phenograph",	"rca",	"sc3",	"scran",	"seurat",	"sincera",	"tscan", "clue", "SAME_AIC", "SAME_BIC")
#all.algorithms <- c("cibersort","gsea","gsva","metaneighbor","ora","adobo", "sccatch")
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


##################################################################
#######               Calculating F Measure               ########
##################################################################
f.measure <- low.CI <- high.CI <- matrix(0, nrow=length(all.algorithms), ncol=length(all.datasets), dimnames=list(all.algorithms, all.datasets))

df <- c()
for (dataset in all.datasets) {
  print(dataset)
  #data.path <- paste0(Work.Dir, "/Data/Full-Results/", dataset, ".csv")
  data.path <- paste0(Work.Dir,"/cellres/", dataset, "_paper_predictions.tsv")
  data <- read.csv(data.path, header=T, sep='\t')

  for (alg in intersect(all.algorithms, colnames(data))) {
    C <- data[, ref]
    K <- data[, alg]
    
    # Identify Mean F measures
    F.Measure.List <- F.Measure.Func(C, K)
    
    # Measure 95% confidence intervals using bootstrapping
    F.vector <- F.Measure.List$F.vector
    CI <- quantile(F.Bootstrap.Func(F.vector), probs = c(0.025, 0.5, 0.975))
    low.CI[alg, dataset] <- round(CI[1], 2)
    high.CI[alg, dataset] <- round(CI[3], 2)
    f.measure[alg, dataset] <- F.Measure.List$F.measure # round(CI[2], 2) # 
    
    x <- c(dataset, alg, f.measure[alg, dataset], low.CI[alg, dataset], high.CI[alg, dataset])
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
colnames(df) <- c("dataset", "algorithm", "F.Mean", "F.Low", "F.High")
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

pdf(paste0(Work.Dir, "/Figures/F-Measures", status, ".pdf"), width=12, height=6)
print(plt)
dev.off()

plt2 <- ggplot(df, aes(x = reorder(algorithm, -Rank), y = score, fill = Rank)) +
  geom_boxplot() +
  coord_flip() +
  labs(y = "Rank of Algorithm",
       x = "Algorithms")

pdf(paste0(Work.Dir, "/Figures/Ranks", status, ".pdf"), width=8, height=4)
print(plt2)
dev.off()

n.cl<- sort(apply(data[,all.algorithms], 2, function(x) {length(unique(x))}))
r <- score[all.algorithms, all.datasets[[2]]]
pdf(paste0(Work.Dir, "/Figures/Fmeasure-NumCl", status, ".pdf"), width=4, height=4)
plot(n.cl, r)
dev.off()


# plt2 <- ggplot(df, aes(x=algorithm, y=score, fill = col)) +
#  geom_boxplot() +
#  coord_flip()


# ggplot(df, aes(x=algorithm, y=F.Mean, group=dataset)) +
#   # ggplot(df, aes(group=dataset)) +
#   geom_bar( aes(x=reorder(algorithm, F.Mean), y=F.Mean), stat="identity", fill="blue", alpha=0.9) +
#   # geom_bar( aes(x=algorithm, y=F.Mean), stat="identity", fill="blue", alpha=0.9) +
#   geom_errorbar( aes(x=algorithm, ymin=F.Low, ymax=F.High), width=0.4, colour="orange", alpha=0.9, size=1) +
#   facet_wrap(~ dataset, ncol=2) +
#   theme(axis.text.x=element_text(angle =+45, vjust = 0.5))
# # facet_grid(. ~ dataset)