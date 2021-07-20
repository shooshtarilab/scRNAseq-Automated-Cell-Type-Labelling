import os
import numpy as np
from pathlib import Path
import pandas as pd
import time as tm
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neighbors import NearestCentroid
from sklearn.ensemble import RandomForestClassifier
import pickle

import warnings
warnings.filterwarnings("ignore")

main_dir = 'path_to_data_folder'
out_dir = 'path_to_save_output'

dataset = 'GSE72056'

def main():
    Spath = out_dir+dataset+'/'
    Path(Spath).mkdir(parents=True, exist_ok=True)
    os.chdir(Spath) # Output directory

    data, labels, CV_data = ReadData(main_dir+dataset+'/Data.csv.gz',
                                     main_dir+dataset+'/Labels.csv',
                                     main_dir+dataset+'/CV_folds.pkl')
                                     
    """
    Data.csv.gz : cells X genes count matrix, barcodes as row names and
    gene names as column names
    Labels.csv : Annotations file which contains cell type labels correspond to
    the cells in the count matrix
    CV_folds.pkl : pickle file that contains index used in the cross-validation
    """

    SVM(data, labels, CV_data)
    SVMrej(data, labels, CV_data)
    kNN9(data, labels, CV_data)
    LDA(data, labels, CV_data)
    LDArej(data, labels, CV_data)
    NMC(data, labels, CV_data)
    RF(data, labels, CV_data)
    ACTINN(data, labels, CV_data)
    Cell_BLAST(main_dir+dataset+'/Data.csv.gz',
               main_dir+dataset+'/Labels.csv',
               CV_data)

def ReadData(DataPath, LabelsPath, CV_Path):
    # read the CV file
    with open(CV_Path, 'rb') as f:
        CV_data = pickle.load(f)
    col = CV_data['col']

    # read the data
    data = pd.read_csv(DataPath,index_col=0,sep=',')
    labels = pd.read_csv(LabelsPath, header=0,index_col=None, sep=',', usecols = col)

    return data, labels, CV_data

def SVM(data, labels, CV_data):
    #Get the CV info
    nfolds = CV_data['nfolds']
    tokeep = CV_data['tokeep']

    test_ind = CV_data['test_ind']
    train_ind = CV_data['train_ind']

    labels = labels.iloc[tokeep]
    data = data.iloc[tokeep]

    # normalize data
    data = np.log10(data+1)

    Classifier = LinearSVC()

    tr_time=[]
    ts_time=[]
    truelab = []
    pred = []

    for i in range(np.squeeze(nfolds)):
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1

        train=data.iloc[train_ind_i]
        test=data.iloc[test_ind_i]
        y_train=labels.iloc[train_ind_i]
        y_test=labels.iloc[test_ind_i]

        start=tm.time()
        Classifier.fit(train, y_train)
        tr_time.append(tm.time()-start)

        start=tm.time()
        predicted = Classifier.predict(test)
        ts_time.append(tm.time()-start)

        truelab.extend(y_test.values)
        pred.extend(predicted)

    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)

    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)

    truelab.to_csv("SVM_True_Labels.csv", index = False)
    pred.to_csv("SVM_Pred_Labels.csv", index = False)
    tr_time.to_csv("SVM_Training_Time.csv", index = False)
    ts_time.to_csv("SVM_Testing_Time.csv", index = False)

def SVMrej(data, labels, CV_data, Threshold = 0.7):
    nfolds = CV_data['nfolds']
    test_ind = CV_data['test_ind']
    train_ind = CV_data['train_ind']

    # normalize data
    data = np.log10(data+1)

    Classifier = LinearSVC()
    clf = CalibratedClassifierCV(Classifier)

    tr_time=[]
    ts_time=[]
    truelab = []
    pred = []

    for i in range(np.squeeze(nfolds)):
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1

        train=data.iloc[train_ind_i]
        test=data.iloc[test_ind_i]
        y_train=labels.iloc[train_ind_i]
        y_test=labels.iloc[test_ind_i]

        start=tm.time()
        clf.fit(train, y_train)
        tr_time.append(tm.time()-start)

        start=tm.time()
        predicted = clf.predict(test)
        prob = np.max(clf.predict_proba(test), axis = 1)
        unlabeled = np.where(prob < Threshold)
        predicted[unlabeled] = 'Unknown'
        ts_time.append(tm.time()-start)

        truelab.extend(y_test.values)
        pred.extend(predicted)

    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)

    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)

    truelab.to_csv("SVMrej_True_Labels.csv", index = False)
    pred.to_csv("SVMrej_Pred_Labels.csv", index = False)
    tr_time.to_csv("SVMrej_Training_Time.csv", index = False)
    ts_time.to_csv("SVMrej_Testing_Time.csv", index = False)

def kNN9(data, labels, CV_data):
    nfolds = CV_data['nfolds']
    test_ind = CV_data['test_ind']
    train_ind = CV_data['train_ind']

    # normalize data
    data = np.log10(data+1)

    Classifier = KNeighborsClassifier(n_neighbors=9)

    tr_time=[]
    ts_time=[]
    truelab = []
    pred = []

    for i in range(np.squeeze(nfolds)):
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1

        train=data.iloc[train_ind_i]
        test=data.iloc[test_ind_i]
        y_train=labels.iloc[train_ind_i]
        y_test=labels.iloc[test_ind_i]

        start=tm.time()
        Classifier.fit(train, y_train)
        tr_time.append(tm.time()-start)

        start=tm.time()
        predicted = Classifier.predict(test)
        ts_time.append(tm.time()-start)

        truelab.extend(y_test.values)
        pred.extend(predicted)

    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)

    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)

    truelab.to_csv("kNN9_True_Labels.csv", index = False)
    pred.to_csv("kNN9_Pred_Labels.csv", index = False)
    tr_time.to_csv("kNN9_Training_Time.csv", index = False)
    ts_time.to_csv("kNN9_Testing_Time.csv", index = False)

def LDA(data, labels, CV_data):
    nfolds = CV_data['nfolds']
    test_ind = CV_data['test_ind']
    train_ind = CV_data['train_ind']

    # normalize data
    data = np.log10(data+1)

    Classifier = LinearDiscriminantAnalysis()

    tr_time=[]
    ts_time=[]
    truelab = []
    pred = []

    for i in range(np.squeeze(nfolds)):
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1

        train=data.iloc[train_ind_i]
        test=data.iloc[test_ind_i]
        y_train=labels.iloc[train_ind_i]
        y_test=labels.iloc[test_ind_i]

        start=tm.time()
        Classifier.fit(train, y_train)
        tr_time.append(tm.time()-start)

        start=tm.time()
        predicted = Classifier.predict(test)
        ts_time.append(tm.time()-start)

        truelab.extend(y_test.values)
        pred.extend(predicted)

    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)
    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)

    truelab.to_csv("LDA_True_Labels.csv", index = False)
    pred.to_csv("LDA_Pred_Labels.csv", index = False)
    tr_time.to_csv("LDA_Training_Time.csv", index = False)
    ts_time.to_csv("LDA_Testing_Time.csv", index = False)

def LDArej(data, labels, CV_data, Threshold = 0.7):
    nfolds = CV_data['nfolds']
    test_ind = CV_data['test_ind']
    train_ind = CV_data['train_ind']

    # normalize data
    data = np.log10(data+1)

    Classifier = LinearDiscriminantAnalysis()

    tr_time=[]
    ts_time=[]
    truelab = []
    pred = []

    for i in range(np.squeeze(nfolds)):
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1

        train=data.iloc[train_ind_i]
        test=data.iloc[test_ind_i]
        y_train=labels.iloc[train_ind_i]
        y_test=labels.iloc[test_ind_i]

        start=tm.time()
        Classifier.fit(train, y_train)
        tr_time.append(tm.time()-start)

        start=tm.time()
        predicted = Classifier.predict(test)
        prob = np.max(Classifier.predict_proba(test), axis = 1)
        unlabeled = np.where(prob < Threshold)
        predicted[unlabeled] = 'Unknown'
        ts_time.append(tm.time()-start)

        truelab.extend(y_test.values)
        pred.extend(predicted)

    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)

    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)

    truelab.to_csv("LDArej_True_Labels.csv", index = False)
    pred.to_csv("LDArej_Pred_Labels.csv", index = False)
    tr_time.to_csv("LDArej_Training_Time.csv", index = False)
    ts_time.to_csv("LDArej_Testing_Time.csv", index = False)

def NMC(data, labels, CV_data):
    nfolds = CV_data['nfolds']
    test_ind = CV_data['test_ind']
    train_ind = CV_data['train_ind']

    # normalize data
    data = np.log10(data+1)

    Classifier = NearestCentroid()

    tr_time=[]
    ts_time=[]
    truelab = []
    pred = []

    for i in range(np.squeeze(nfolds)):
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1

        train=data.iloc[train_ind_i]
        test=data.iloc[test_ind_i]
        y_train=labels.iloc[train_ind_i]
        y_test=labels.iloc[test_ind_i]

        start=tm.time()
        Classifier.fit(train, y_train)
        tr_time.append(tm.time()-start)

        start=tm.time()
        predicted = Classifier.predict(test)
        ts_time.append(tm.time()-start)

        truelab.extend(y_test.values)
        pred.extend(predicted)

    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)

    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)

    truelab.to_csv("NMC_True_Labels.csv", index = False)
    pred.to_csv("NMC_Pred_Labels.csv", index = False)
    tr_time.to_csv("NMC_Training_Time.csv", index = False)
    ts_time.to_csv("NMC_Testing_Time.csv", index = False)

def RF(data, labels, CV_data):
    nfolds = CV_data['nfolds']
    test_ind = CV_data['test_ind']
    train_ind = CV_data['train_ind']

    # normalize data
    data = np.log10(data+1)

    Classifier = RandomForestClassifier(n_estimators = 50)

    tr_time=[]
    ts_time=[]
    truelab = []
    pred = []

    for i in range(np.squeeze(nfolds)):
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1

        train=data.iloc[train_ind_i]
        test=data.iloc[test_ind_i]
        y_train=labels.iloc[train_ind_i]
        y_test=labels.iloc[test_ind_i]

        start=tm.time()
        Classifier.fit(train, y_train)
        tr_time.append(tm.time()-start)

        start=tm.time()
        predicted = Classifier.predict(test)
        ts_time.append(tm.time()-start)

        truelab.extend(y_test.values)
        pred.extend(predicted)

    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)

    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)

    truelab.to_csv("RF_True_Labels.csv", index = False)
    pred.to_csv("RF_Pred_Labels.csv", index = False)
    tr_time.to_csv("RF_Training_Time.csv", index = False)
    ts_time.to_csv("RF_Testing_Time.csv", index = False)

def ACTINN(data, labels, CV_data):
    nfolds = CV_data['nfolds']
    test_ind = CV_data['test_ind']
    train_ind = CV_data['train_ind']

    tot=[]
    truelab = []
    pred = []

    for i in range(np.squeeze(nfolds)):
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1

        train=data.iloc[train_ind_i]
        test=data.iloc[test_ind_i]
        y_train=labels.iloc[train_ind_i]
        y_test=labels.iloc[test_ind_i]

        train = train.transpose()
        test = test.transpose()

        train.to_csv("train.csv")
        test.to_csv("test.csv")
        y_train.to_csv("train_lab.csv", header = False, index = True, sep = '\t')
        y_test.to_csv("test_lab.csv", header = False, index = True, sep = '\t')

        tm.sleep(60)

        os.system("python /Users/pluo/ACTINN-master/actinn_format.py -i train.csv -o train -f csv")
        os.system("python /Users/pluo/ACTINN-master/actinn_format.py -i test.csv -o test -f csv")

        start = tm.time()
        print('start to predict')
        os.system("/Users/pluo/ACTINN-master/actinn_predict.py -trs train.h5 -trl train_lab.csv -ts test.h5")
        tot.append(tm.time()-start)

        tm.sleep(60)

        truelab.extend(y_test.values)
        predlabels = pd.read_csv('predicted_label.txt',header=0,index_col=None, sep='\t', usecols = [1])
        pred.extend(predlabels.values)

    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)
    tot_time = pd.DataFrame(tot)

    truelab.to_csv("ACTINN_True_Labels.csv", index = False)
    pred.to_csv("ACTINN_Pred_Labels.csv", index = False)
    tot_time.to_csv("ACTINN_Total_Time.csv", index = False)

def Cell_BLAST(DataPath, LabelsPath, CV_data):
    import tensorflow as tf
    tf.logging.set_verbosity(0)

    import Cell_BLAST as cb
    from numpy import genfromtxt as gft

    nfolds = CV_data['nfolds']
    col = CV_data['col']
    test_ind = CV_data['test_ind']
    train_ind = CV_data['train_ind']

    # read the data and labels
    data_old = cb.data.ExprDataSet.read_table(DataPath,orientation="cg", sep=",", index_col = 0, header = 0, sparsify = True).normalize()
    labels = pd.read_csv(LabelsPath, header=0,index_col=None, sep=',', usecols = col)

    data = cb.data.ExprDataSet(data_old.exprs[tokeep],data_old.obs.iloc[tokeep],data_old.var,data_old.uns)

    labels = gft(LabelsPath, dtype = "str", skip_header = 1, delimiter = ",", usecols = col)

    truelab = []
    pred = []
    tr_time = []
    ts_time = []

    for i in range(np.squeeze(nfolds)):
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1

        train=data[train_ind_i,:]
        test=data[test_ind_i,:]
        y_train = labels[train_ind_i]
        y_test = labels[test_ind_i]

        train.obs['cell_type'] = y_train

        start = tm.time()

        # reduce dimensions
        num_epoch = 30
        models = []

        for j in range(4):
            models.append(cb.directi.fit_DIRECTi(train, epoch=num_epoch, patience=10, random_seed = j, path="%d" % j))

        # train model
        blast = cb.blast.BLAST(models, train)
        tr_time.append(tm.time()-start)

        # predict labels
        start = tm.time()
        test_pred = blast.query(test).annotate('cell_type')
        ts_time.append(tm.time()-start)

        truelab.extend(y_test)
        pred.extend(test_pred.values)

    #write results
    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)

    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)

    truelab.to_csv("Cell_BLAST_True_Labels.csv", index = False)
    pred.to_csv("Cell_BLAST_Pred_Labels.csv", index = False)
    tr_time.to_csv("Cell_BLAST_Training_Time.csv", index = False)
    ts_time.to_csv("Cell_BLAST_Testing_Time.csv", index = False)

main()
