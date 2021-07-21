import pandas as pd
from random import choices
from random import sample

df = pd.read_csv("Labels.csv") #path to dataset Labels.csv file

Labels = list(df.iloc[:,0])

# Select the cell types used in the subsampling experiment
type_list = ['Dendritic cells', 'Fibroblasts', 'Mast cells', 'T cells',
             'Langerhans cells', 'Endothelial cell', 'Cancer cells',
             'Alveolar cell', 'Natural killer cells', 'Macrophages', 'B cells',
             'Basal cells', 'Granulocytes', 'Lymphatic EC', 'Erythroblasts',
             'Secretory club cells', 'Epithelial cells']

typeDic = {}

for item in type_list:
    typeDic[item] = []

for i in range(len(Labels)):
    if Labels[i] in type_list:
        typeDic[Labels[i]].append(i)

N = 800 #number of cells for each cell type

Chosen = []

for key in type_list:
    if len(typeDic[key]) < N:
        Chosen += choices(typeDic[key], k=N)
    else:
        Chosen += sample(typeDic[key], k=N)

Cells = list(pd.read_csv("celllist.tsv", header=None).iloc[:,0]) # row names corresponds to the Labels file.

cLabels = [Labels[i].upper() for i in Chosen]
cLabels = [item.replace(' ', '_') for item in cLabels]
cCells = [Cells[i] for i in Chosen]

# save the Labels and their corresponding cell names
pd.DataFrame(cCells).to_csv('./sub_sample/'+ str(N)+'/celllist_sub.tsv', index=False, header=False)
pd.DataFrame(cLabels).to_csv('./sub_sample/'+str(N)+'/Labels.csv', index=False)

# Read the count matrix
data = pd.read_csv(main_dir+"Lambrechts"+'/Data.csv.gz',index_col=0,sep=',')

data.index = Cells
data = data.loc[cCells]
print(data.shape)

#revise rownames
data.index = ['Sample_'+str(i) for i in range(data.shape[0])]
data.to_csv(main_dir+"Lambrechts/sub_sample/"+"Data.csv")
