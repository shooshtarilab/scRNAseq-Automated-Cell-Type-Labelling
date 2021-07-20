.libPaths("/hpf/largeprojects/ccmbio/amahalanabis/tools/R/library")
suppressPackageStartupMessages(library(SingleCellExperiment))

all.datasets <- c("li_crc_allCells", "darmanis_glioblastoma_allCells", "JA_melanoma_allCells",
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

silhouette_matrix <- matrix(0, nrow=length(all.datasets), ncol=length(all.algorithms), dimnames=list(all.datasets, all.algorithms))

for (dataset in all.datasets) {
  
  filename = paste0(dataset, "_TPM.RData")
  print(filename)
  load(filename)
  data.path <- paste0("Clustering_Results/", dataset, ".csv")
  data <- read.csv(data.path, header=T, stringsAsFactors = FALSE)
  
  mat <- counts(sce)[,colnames(counts(sce)) %in% data[,"cell"]]
  rownames(mat)=NULL
  data <- data[match(colnames(mat), data$cell),]
  mat <- log2(mat + 1)
  mat <- t(mat)
  mat <- mat[ , which(apply(mat, 2, var) != 0)]
  matrix.pca <- prcomp(mat, center = TRUE,scale. = TRUE) 
  D = dist(matrix.pca$x[,1:10])
  
  for (alg in intersect(all.algorithms, colnames(data))) {
    print (alg)
    sil <- cluster::silhouette(data[,alg], D)
    silhouette_matrix[dataset, alg] <- mean(sil[,3])
    write.csv(silhouette_matrix, file="dataframe_results/silhouette_matrix.csv")
  }
}

write.csv(silhouette_matrix, file="silhouette_matrix.csv")
