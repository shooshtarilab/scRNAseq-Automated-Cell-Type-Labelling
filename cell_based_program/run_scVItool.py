import os
import numpy as np
import pandas as pd
import time as tm
import pickle

import scvi
import scanpy as sc

from sys import argv

main_dir = argv[1] #'path_to_data_folder'
out_dir = "cell_based_program/output"#'path_to_save_output'

dataset = argv[2] #'GSE72056'

DataPath = main_dir+dataset+'/counts.tsv.xz'
LabelsPath = main_dir+dataset+'/clusters.csv'
CV_Path = main_dir+dataset+'/CV_folds.pkl'
OutputDir = out_dir+dataset+'/' # Output directory

"""
Data.csv.gz : cells X genes count matrix, barcodes as row names and
gene names as column names
Labels.csv : Annotations file which contains cell type labels correspond to
the cells in the count matrix
CV_folds.pkl : pickle file that contains index used in the cross-validation
"""

#os.chdir(OutputDir)

Data = pd.read_csv(DataPath, index_col=0)
Label = pd.read_csv(LabelsPath, header=0, index_col=None)

Label.index = Data.index.values

with open(CV_Path, 'rb') as f:
    CV_data = pickle.load(f)

nfolds = CV_data['nfolds']
test_ind = CV_data['test_ind']
train_ind = CV_data['train_ind']

pred = []
tot=[]

for i in range(np.squeeze(nfolds)):
    test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
    train_ind_i = np.array(train_ind[i], dtype = 'int') - 1

    start = tm.time()

    adata = sc.AnnData(Data)
    labels = Label.copy()
    labels.iloc[test_ind_i,] = "Unknown"
    adata.obs["celltype_scanvi"] = labels

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=2000,
        layer="counts",
        subset=True
    )

    scvi.data.setup_anndata(
        adata,
        layer="counts",
        batch_key=None,
        labels_key="celltype_scanvi",
    )

    lvae = scvi.model.SCANVI(adata, "Unknown", use_cuda=True, n_latent=30, n_layers=2)
    lvae.train(n_epochs_semisupervised = 50)

    adata.obs["C_scANVI"] = lvae.predict(adata)

    tot.append(tm.time()-start)

    pred.extend(adata.obs["C_scANVI"].iloc[test_ind_i,])

pred = pd.DataFrame(pred)
tot_time = pd.DataFrame(tot)
os.chdir(OutputDir)
pred.to_csv("scVI_Pred_Labels.csv", index = False)
tot_time.to_csv("scVI_Total_Time.csv", index = False)
