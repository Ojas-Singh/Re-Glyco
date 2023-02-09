import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from lib import pdb,graph,dihedral,clustering,tfindr
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
from mpl_toolkits.mplot3d import Axes3D
from random import randint



name="bisecting"
f="data/"+name+"/"+name+".dry.pdb"

pdbdata, frames = pdb.multi(f)
df = pdb.to_DF(pdbdata)


# idx_noH=df.loc[(df['Element']!="H"),['Number']].iloc[:]['Number']-1
# pca2,pca20 = clustering.pcawithG(frames,idx_noH)
# pca2.to_csv("data/"+name+"/"+name+'_G_pca.csv',index_label="i")  
# pca20.to_csv("data/"+name+"/"+name+'_G_pca_full.csv',index_label="i") 


pairs = tfindr.savetotorparts(name)
torsion_names = dihedral.pairtoname(pairs,df)

torsiondataDF= dihedral.pairstotorsion(pairs,frames,torsion_names)
torsiondataDF.to_csv("data/"+name+"/"+name+'_torsions.csv')


#plots all diherdal KDEs
# dihedral.kdemax(torsiondataDF)

