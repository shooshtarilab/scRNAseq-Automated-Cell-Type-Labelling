import pandas as pd 
import scipy.stats as st
import numpy as np

datasets = ['cb', 'dg', 'jam', 'li_crc', 'llc', 'peng', 'tm', 'vg']

pca = pd.DataFrame(columns=['dataset', 'algorithm','percentage','p_high', 'p_low'])

for dataset in datasets:
    print(f"-------------------------------{dataset}-------------------------------------")
    preds = pd.read_csv(f'../cellres/{dataset}_paper_predictions.tsv', sep='\t')
    for col in preds.columns:
        if col in ['cell','truth','seurat','dataset', 'cluster','cluster_truth', 'mergecell', 'paper']:
            continue
        print(f"{col}")
        vals = []
        for i in range(10):
            indx = preds.sample(n=len(preds['truth']), replace=True, axis=0)
            vals.append(len(indx[indx['truth']==indx[col]])/len(indx['truth']))
        ci = st.norm.interval(alpha=0.95, loc=np.mean(vals), scale=st.sem(vals))
        pca = pca.append({'dataset':dataset,
                          'algorithm':col,
                          'percentage':round(np.mean(vals),2),
                          'p_high': round(ci[1],2),
                          'p_low': round(ci[0],2)
                          },
                         ignore_index=True)

print(pca.head())

pca.to_csv('../Rdata_paper/percentage_correctly_assigned_bootstrap.tsv',sep='\t',index=False)