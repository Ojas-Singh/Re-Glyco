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


idx_noH=df.loc[(df['Element']!="H"),['Number']].iloc[:]['Number']-1
pcaG= clustering.pcawithG(frames,idx_noH,3)
pcaG.to_csv("data/"+name+"/"+name+'_G_pca.csv',index_label="i")

tsneG = clustering.tsnewithG(frames,idx_noH,3)
tsneG.to_csv("data/"+name+"/"+name+'_G_tsne.csv',index_label="i")  

pairs,external,internal = tfindr.torsionspairs(name)
pairs = np.asarray(pairs)

torsion_names = dihedral.pairtoname(external,df)
ext_DF = dihedral.pairstotorsion(external,frames,torsion_names)
for i in range(len(internal)):
    torsion_names.append("internal")
torsiondataDF= dihedral.pairstotorsion(pairs,frames,torsion_names)
torsiondataDF.to_csv("data/"+name+"/"+name+'_torsions.csv',index_label="i")

tor=clustering.normalizetorsion(ext_DF)
pcaT = clustering.pcawithT(tor,3)
pcaT.to_csv("data/"+name+"/"+name+'_T_pca.csv',index_label="i")  

tsneT = clustering.tsnewithT(tor,3)
tsneT.to_csv("data/"+name+"/"+name+'_T_tsne.csv',index_label="i") 


#plots all diherdal KDEs
# dihedral.kdemax(torsiondataDF)

