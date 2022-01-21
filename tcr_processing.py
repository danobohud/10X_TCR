import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
from matplotlib import pyplot as plt, cm as mpl_cm
from cycler import cycler
import pickle

sc.set_figure_params(figsize=(4, 4))
sc.settings.verbosity = 2  # verbosity: errors (0), warnings (1), info (2), hints (3)


def save_pickle(data, savefile):
    assert type(savefile)==str
    print('Saving file to ',savefile)
    with open(savefile, 'wb') as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

def load_pickle(loadfile):
    assert type(loadfile)==str
    print('Loading files from', loadfile)
    with open(loadfile,'rb') as f:
        data = pickle.load(f)
    return data

def get_clondict(consensus_annotations):
    clon = pd.read_csv(consensus_annotations)
    clondict = {clon.iloc[x]['cdr3']:{'chain': clon.iloc[x]['chain'], 'clonotype_id': clon.iloc[x]['clonotype_id'][9:]} for x in range(len(clon))}
    CDR3=pd.DataFrame.from_dict(clondict,orient='index')
    return CDR3