import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from lib import pdb,graph,dihedral,clustering,tfindr
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
from mpl_toolkits.mplot3d import Axes3D
from random import randint



# name="A3G3S1-F3"
name =  "bisecting"
f="data/"+name+"/"+name+".dry.pdb"

pdbdata, frames = pdb.multi(f)
df = pdb.to_DF(pdbdata)


# idx_noH=df.loc[(df['Element']!="H"),['Number']].iloc[:]['Number']-1
# pca2,pca20 = clustering.pcawithG(frames,idx_noH)
# pca2.to_csv("data/"+name+"/"+name+'_G_pca.csv',index_label="i")  
# pca20.to_csv("data/"+name+"/"+name+'_G_pca_full.csv',index_label="i") 


pairs,external,internal = tfindr.torsionspairs(name)
pairs = np.asarray(pairs)

torsion_names = dihedral.pairtoname(external,df)
ext_DF = dihedral.pairstotorsion(external,frames,torsion_names)
for i in range(len(internal)):
    torsion_names.append("internal")
torsiondataDF= dihedral.pairstotorsion(pairs,frames,torsion_names)
torsiondataDF.to_csv("data/"+name+"/"+name+'_torsions.csv',index_label="i")

tor=clustering.normalizetorsion(ext_DF)
pca2,pca20 = clustering.pcawithT(tor)
pca2.to_csv("data/"+name+"/"+name+'_T_pca.csv',index_label="i")  
# pca20.to_csv("data/"+name+"/"+name+'_T_pca_full.csv',index_label="i") 

#plots all diherdal KDEs
# dihedral.kdemax(torsiondataDF)

