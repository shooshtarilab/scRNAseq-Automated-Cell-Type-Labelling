import sys
import adobo as ad
import pandas as pd

def run_adobo(df, clusterfile, methodname,outfile, symbol_switch=False, norm_tech='fqn'):
    
    # for some reason passing the dataframe in doesn't work for the fqn normalisation
    #if df.isnull().values.any():
    #    exp = ad.dataset(df.fillna(0))
    #else:
    #    exp = ad.dataset(df)
    cbdf = pd.read_csv(df, 
                   sep='\t',
                   index_col=0)
    exp = ad.IO.load_from_file(df,
                               desc='dataset',
                               sep='\t',
                               header=True,
                               sparse=False)
    
    
    clust = pd.read_csv(clusterfile) # read from the file, make sure you've renamed the two columns
    clust = clust.set_index('cell') # set the index to the first column
    clust = pd.Series(clust[methodname]) # make a series with the cluster numbers
                                            # need to subtract one so the cluster numbers start at zero
                                            # this may fix my missing clusters/weird labels
    clust = clust[cbdf.columns] #remove cells that were filtered
    
    #display(clust)

    #cells were removed from clustering, remove them from the input data
    diff = [item for item in cbdf.columns.values if item not in clust.index.values]
    print(f"there were {len(diff)} cells removed from the raw inputs")
    #exp.norm_data['standard'] = {}
    #if len(diff)>0:
    #    new_df = df.drop(columns=diff, axis=1)
    #    exp.norm_data['standard']['data'] = new_df #store the normalised data (in this case with missing values removed)
    #else:
    #    exp.norm_data['standard']['data'] = df #store the normalised data
    
    #normalisation
    ad.normalize.norm(exp, method = norm_tech)
    ndf = exp.norm_data[norm_tech]['data']
    ndf.columns = cbdf.columns
    ndf.to_csv(df.rsplit('/',1)[0]+"fqn.tsv.xz",sep="\t")
    
    #make the dict structure
    #exp.norm_data[norm_tech]['clusters'] = {}
    exp.norm_data[norm_tech]['clusters']['custom'] = {}
    exp.norm_data[norm_tech]['clusters']['custom']['membership'] = clust # store the cluster membership
    
    if symbol_switch:
        ad.preproc.symbol_switch(exp, species='human')
    
    ad.bio.cell_type_predict(exp, clustering='custom',verbose=True) # predict on custom clusters
    #display(exp.norm_data[norm_tech]['clusters']['custom']['cell_type_prediction']) # view the predictions
    print(exp.norm_data[norm_tech]['clusters']['custom']['cell_type_prediction']) # view the predictions
    
    exp.norm_data[norm_tech]['clusters']['custom']['cell_type_prediction'].to_csv(outfile,
                                                                                   sep='\t') # view the predictions 
    return


if len(sys.argv) == 6:
    run_adobo(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
else:
    run_adobo(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

