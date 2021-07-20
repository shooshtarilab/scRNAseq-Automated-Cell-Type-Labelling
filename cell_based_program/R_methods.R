rm(list = ls())

args <- commandArgs(trailingOnly=TRUE)

main_dir <- "cell_based_program"
out_dir <- file.path(main_dir,"output")

dataset <- args[1] #the dataset to run
input_dir <- args[2] #path to the folder with inputs 

source(file.path(main_dir, 'R_func','CaSTLe.R'))
source(file.path(main_dir, 'R_func','CHETAH.R'))
source(file.path(main_dir, 'R_func','scID.R'))
source(file.path(main_dir, 'R_func','scmap.R'))
source(file.path(main_dir, 'R_func','singleCellNet.R'))
source(file.path(main_dir, 'R_func','SingleR.R'))
source(file.path(main_dir, 'R_func','scPred.R'))

DataPath <- file.path(input_dir,dataset,'counts.tsv.xz')
LabelsPath <- file.path(input_dir,dataset,'clusters.csv')
CV_RDataPath <- file.path(input_dir,dataset,'CV_folds.RData')
OutputPath <- file.path(out_dir,dataset,)


# Data.csv.gz : cells X genes count matrix, barcodes as row names and
# gene names as column names
# Labels.csv : Annotations file which contains cell type labels correspond to
# the cells in the count matrix
# CV_folds.RData : file that contains index and parameters used in the cross-validation


Data <- read.table(DataPath,row.names = 1,header=T,sep='\t')
Data <- t(Data)
Labels <- as.matrix(read.csv(LabelsPath))

# Comment out the algorithms that you don't need
CaSTLe(Data, Labels, CV_RDataPath, OutputPath)
scID(Data, Labels, CV_RDataPath, OutputPath)
CHETAH(Data, Labels, CV_RDataPath, OutputPath)
scmap(Data, Labels, CV_RDataPath, OutputPath)
SingleR(Data, Labels, CV_RDataPath, OutputPath)
singleCellNet(Data, Labels, CV_RDataPath, OutputPath)
scPred(Data, Labels, CV_RDataPath, OutputPath)
