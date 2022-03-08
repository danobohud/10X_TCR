import pandas as pd
import numpy as np

def noise_correct(df,neg,pos):
    '''Subtract maximum negative binding signal from each positive dextramer per cell'''
    print('Correcting for negative dextramer binding')
    out = dict()
    for i in range(len(df)):
        cell = df.iloc[i]
        ctrl = np.max([cell[c] for c in neg])
        positives=[c for c in pos if (cell[c]-ctrl)>0]
        if len(positives)>0:
            out[cell['barcode']]={p: cell[p]-ctrl for p in positives}
    df=df[df['barcode'].isin(out.keys())]

    return df, out
    
def batch_correct(vals):
    print('Normalising binding signals per cell')
    out = dict()
    for barcode in list(vals.keys()):
        denominator = np.sum([vals[barcode][key] for key in vals[barcode].keys()])
        out[barcode] = {k:v/denominator for k,v in vals[barcode].items()}
    return out

def read_clonotypes(contigs,dataframe):
    print('Reading in clonotype data')
    contigs=pd.read_csv(contigs)
    clondict={contigs.iloc[i]['barcode']:contigs.iloc[i]['raw_clonotype_id'] for i in range(len(contigs))}
    dataframe['clonotype']=[clondict[barcode] for barcode in dataframe['barcode'].unique()]
    return dataframe, clondict

def clonotype_correct(dataframe,original_counts,positive_binders,clonotypes):
    print('Normalising binding signals across clonotypes')
    denoms=dict()
    for clonotype in dataframe['clonotype'].unique():
        sub=dataframe[dataframe['clonotype']==clonotype]
        denoms[clonotype]= {p: np.sum(sub[p].values) for p in positive_binders}

    normdict={}
    for barcode in list(original_counts.keys()):
        clonotype=clonotypes[barcode]
        normdict[barcode]= {binder: original_counts[barcode][binder]/denoms[clonotype][binder] for binder in original_counts[barcode].keys()}

    return normdict

def combine_batch_clonotype(vals,batch,clonotype_norm):
    print('Combining batch and clonotype normalisation')
    finaldict=dict()
    for barcode in vals.keys():
        finaldict[barcode]={binder: (batch[barcode][binder]**2)*clonotype_norm[barcode][binder] for binder in vals[barcode].keys()}
    finaldict2=dict()
    for barcode in finaldict.keys():
        finaldict2[barcode]=max(finaldict[barcode],key=finaldict[barcode].get)
    
    return finaldict2

def dextramer_normalise(matrix_path,contig_path):
    
    print('Reading in dextramer binding data')
    df=pd.read_csv(matrix_path)
    neg = [c for c in df.columns[18:] if '_NC' in c and 'binder' not in c]
    pos = [c for c in df.columns[18:] if c not in neg and 'binder' not in c]
    df2, vals= noise_correct(df,neg,pos)

    batch = batch_correct(vals)
    df2, clonotypes = read_clonotypes(contig_path,df2)
    clonotype_norm = clonotype_correct(df2,vals,pos,clonotypes)
    out = combine_batch_clonotype(vals,batch,clonotype_norm)

    return out