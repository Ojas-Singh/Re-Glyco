import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from lib import pdb,graph,dihedral,clustering,tfindr
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
from mpl_toolkits.mplot3d import Axes3D
from random import randint
import plotly.express as px

name = "glycan"
# name =  "A3g3S1-F3"
# name = "bisecting"
f="data/"+name+"/"+name


pca = pd.read_csv(f+"_G_pca_NB.csv")
tsne = pd.read_csv(f+"_G_tsne.csv")
clustering.plot3dkde(pca)
clustering.plot3dkde(tsne)

Tdata = pd.read_csv(f+"_torsions.csv")
idx,popp = clustering.filterlow(pca,.1)
print(popp)
numofcluster = 7
data,label = clustering.cluster(pca,numofcluster,idx)

data.to_csv("data/"+name+"/"+name+'_clustered.csv',index=False) 
pca["cluster"] = label
tsne["cluster"] = label
Tdata["cluster"] = label
fig=px.scatter(pca, x="0", y="1",
              color="cluster")
fig.show()

# df1 = pd.read_csv(f+"_clustered.csv")

# fig=px.scatter(Tdata, x="4_8_phi", y="4_8_psi",
#               color="cluster")
# fig.show()

# df = pd.read_csv("data/bisecting/bisecting_torsions.csv")
# print(np.shape(clustering.normalizetorsion(df)))