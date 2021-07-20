library(Seurat)
library(scCATCH)
library(data.table)

run_sccatch <- function(data, meta, cluster_method, cancer, tissue){
    seurat.data <- CreateSeuratObject(counts=data,project='sccatch', min.cell=0, min.features=0)
    seurat.data <- NormalizeData(seurat.data)
    for (i in 1:max(meta[,cluster_method])){ 
        #get cells in cluster and set idents in seurat
        use_cells = meta[meta[,cluster_method] == i,]$cell 
        Idents(seurat.data, cells=use_cells) <-i
    }
    cluster_markers <- findmarkergenes(seurat.data,
                                    species = 'Human',
                                    cluster = 'All',
                                    match_CellMatch = T,
                                    cancer = cancer,
                                    tissue = tissue)
    cluster_annos <- scCATCH(cluster_markers,
                        species = 'Human',
                        cancer = cancer,
                        tissue = tissue)
    
    #TODO add the truth column from metadata to cluster_annos
    return(cluster_annos)
}
run_sccatch_nocm <- function(data, meta, cluster_method, cancer, tissue){
    seurat.data <- CreateSeuratObject(counts=data,project='sccatch', min.cell=0, min.features=0)
    seurat.data <- NormalizeData(seurat.data)
    for (i in 1:max(meta[,cluster_method])){ 
        #get cells in cluster and set idents in seurat
        use_cells = meta[meta[,cluster_method] == i,]$cell 
        Idents(seurat.data, cells=use_cells) <-i
    }
    cluster_markers <- findmarkergenes(seurat.data,
                                    species = 'Human',
                                    cluster = 'All',
                                    match_CellMatch = F
                                      )
    cluster_annos <- scCATCH(cluster_markers,
                        species = 'Human',
                        cancer = cancer,
                        tissue = tissue)
    
    #TODO add the truth column from metadata to cluster_annos
    return(cluster_annos)
}

args <- commandArgs(trailingOnly = T)

data <- read.csv(args[1],sep="\t",row.names=1,check.names=F)
data.meta <- read.csv(args[2])

cluster_annos <- run_sccatch(data, data.meta, args[3], args[4], args[-1:-4])
cluster_annos
