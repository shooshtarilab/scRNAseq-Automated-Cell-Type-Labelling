# -*- coding: utf-8 -*-
"""
Created on Mon May  4 17:45:17 2020

@author: luoping

Convert the .RData file to pickle file that can be used in python
"""

import rpy2.robjects as robjects
import numpy as np
import pickle
from sys import argv

args = argv

#main_dir = argv[1] #'C:/Users/luopi/scRNAseq_docker/'

#dataset = argv[2] #'GSE70630'

#CV_RDataPath = main_dir+dataset+'/CV_folds.RData'
CV_RDataPath = args[1] +'/CV_folds.RData'

robjects.r['load'](CV_RDataPath)

CV_data = {}

CV_data['nfolds'] = np.array(robjects.r['n_folds'], dtype = 'int')
CV_data['tokeep'] = np.array(robjects.r['Cells_to_Keep'], dtype = 'bool')
col = np.array(robjects.r['col_Index'], dtype = 'int')
CV_data['col'] = col - 1
test_ind = np.array(robjects.r['Test_Idx'])
train_ind = np.array(robjects.r['Train_Idx'])

CV_data['test_ind'] = {}
CV_data['train_ind'] = {}

for i in range(np.squeeze(CV_data['nfolds'])):
    CV_data['test_ind'][i] = np.array(test_ind[i], dtype = 'int')
    CV_data['train_ind'][i] = np.array(train_ind[i], dtype = 'int')


#with open(main_dir+dataset+'/CV_folds.pkl', 'wb') as f:
with open(args[1]+'/CV_folds.pkl', 'wb') as f:
    pickle.dump(CV_data, f)
