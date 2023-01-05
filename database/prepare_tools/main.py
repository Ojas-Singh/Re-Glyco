import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from lib import pdb,graph,dihedral


name="GlycanId"
f="../../data/prodfull.dry.pdb"
pdbdata, frames = pdb.multi(f)
pdbdataDF = pdb.to_DF(pdbdata)

tormeta=[[2,3,4,4],[3,4,4,4],[4,5,3,3],[4,7,4,4],[4,8,6,6],[8,9,2,2],[5,6,2,2]]
torsionmeta = []
for i in tormeta:
    torsionmeta.append(dihedral.res2input(i[0],i[1],i[2],i[3],pdbdataDF))
torsiondataDF= dihedral.torsions(torsionmeta,pdbdata,frames,f)
torsions=torsiondataDF.to_numpy()


b = np.asanyarray(torsionmeta,dtype=object)
a = np.asanyarray(pdbdataDF,dtype=object)
np.savez_compressed(name, a=a, b=b,allow_pickle=True)
torsiondataDF.to_csv('output/'+name+'torsions.csv')











# G = np.zeros((len(frames),int(len(pdbdata[0])*(len(pdbdata[0])+1)/2)))
# for i in range(len(frames)):
#     G[i]= graph.G_flatten(frames[i])

# pca = PCA(n_components=20)
# t = pca.fit_transform(G)
# x=[]
# y=[]
# for i in t:
#     x.append(i[0])
#     y.append(i[1])
# df = pd.DataFrame(data=np.column_stack((x,y)), columns = ['X','Y'])
# df.to_csv('output/'+name+'G_pca.csv',index_label="i")  
# PCA_components = pd.DataFrame(t)
# PCA_components.to_csv('output/'+name+'G_pca_full.csv',index_label="i") 



