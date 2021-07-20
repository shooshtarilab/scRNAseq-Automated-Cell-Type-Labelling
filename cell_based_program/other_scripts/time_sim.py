#Time plot

import pandas as pd
import numpy as np
import os

dataset_list = ['GSE72056','GSE75688','GSE81861','GSE115978','GSE84465','GSE116256',
                'Lambrechts','Peng']

column_name = ['Tirosh_Melanoma','Chung_BC','Li_CRC','JA_Melanoma','Darmanis_GBM',
               'vanGalan_AML','Lambrechts_LC','Peng_PC']

time_dir = 'path_to_prediction_folder'

tt_list = ["CaSTLe","Cell_BLAST","kNN9","LAmbDA","LDA","LDArej","NMC","RF",
           "scmapcell","scmapcluster","scPred","SingleCellNet","SVM","SVMrej"]
total_list = ["scVI","ACTINN","scID","CHETAH","singleR"]

mtx = np.zeros(((len(tt_list) + len(total_list)), len(dataset_list)))
Std = np.zeros(mtx.shape)

trainM = np.zeros((len(tt_list), len(dataset_list)))
testM = np.zeros(trainM.shape)

trainStd = np.zeros(trainM.shape)
testStd = np.zeros(testM.shape)

for i in range(len(dataset_list)):
    dataset = dataset_list[i]

    for j in range(len(tt_list)):
        item = tt_list[j]
        train = time_dir + dataset + '/' + item + '_Training_Time.csv'
        test = time_dir + dataset + '/' + item + '_Testing_Time.csv'

        if os.path.isfile(train):
            train = pd.read_csv(train).iloc[:,0].tolist()
            test = pd.read_csv(test).iloc[:,0].tolist()

            train = np.log1p(np.array(train))#log here#
            test = np.log1p(np.array(test))#log here#

            mtx[j,i] = np.mean(train+test)
            Std[j,i] = np.std(train+test)

            trainM[j,i] = np.mean(train)
            testM[j,i] = np.mean(test)
            trainStd[j,i] = np.std(train)
            testStd[j,i] = np.std(test)

    for j in range(len(total_list)):
        item = total_list[j]
        total = time_dir + dataset + '/' + item + '_Total_Time.csv'
        total = pd.read_csv(total).iloc[:,0].tolist()
        total = np.log1p(np.array(total)) #log here#

        mtx[j+len(tt_list),i] = np.mean(total)
        Std[j+len(tt_list),i] = np.std(total)

# result = pd.DataFrame(mtx, index=tt_list+total_list, columns=column_name)
# result.to_csv('Ping_total_mean.tsv', sep='\t')
# result = pd.DataFrame(Std, index=tt_list+total_list, columns=column_name)
# result.to_csv('Ping_total_std.tsv', sep='\t')

os.chdir('./log_time/')

np.save("mean_total.npy",mtx)
np.save("std_total.npy",Std)

np.save("mean_train.npy",trainM)
np.save("mean_test.npy",testM)
np.save("std_train.npy",trainStd)
np.save("std_test.npy",testStd)
