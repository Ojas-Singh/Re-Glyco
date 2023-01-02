from cProfile import label
import pandas as pd
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import matplotlib
import plotly.express as px
import plotly.graph_objects as go
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn.cluster import BisectingKMeans
from sklearn.cluster import SpectralClustering
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import MeanShift
import hdbscan
from lib import pdb,graph,dihedral

# df = pd.read_csv('output/hep.fulltorsion.csv') 

# l=[[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10],[10,11]] #organic0
# l=[[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10],[10,12]] #organic1
# l=[[2,3],[3,4],[4,6],[6,7],[7,8],[8,9],[9,10],[10,11],[11,12]] #organic2
# col=list(df.columns)
# print(len(col))
# dihedral.plot(df,l)


# df = pd.read_csv('output/selected.csv') 
# fig = px.scatter(df, x="4ZB_UYA_psi", y="4ZB_UYA_phi",color="i")

# fig.show()





# import seaborn as sns
# sns.set_theme(style="darkgrid")
# # Set up the figure
# f, ax = plt.subplots(figsize=(8, 8))
# ax.set_aspect("equal")
# # Draw a contour plot to represent each bivariate density
# a="4ZB_4YA_psi"
# b="4ZB_4YA_phi"
# sns.kdeplot(
#     data=df,
#     x=a,
#     y=b,
#     thresh=.1,color="orange",label=a
# )
# c="4ZB_UYA_psi"
# d="4ZB_UYA_phi"
# sns.kdeplot(
#     data=df,
#     x=c,
#     y=d,thresh=.1,color="blue",label=c
# )
# plt.show()

# df = pd.read_csv('output/hep.fullG_pca2.csv') 
# fig = px.scatter(df, x="X", y="Y",color="i")
# fig.show()


# df = pd.read_csv('output/hep.fulltorsion_pca.csv') 
# data= df[['4ZB_UYA_psi', '4ZB_UYA_phi']].copy()
# x=[]
# y=[]

# for i in range(len(data)):
#     if i%10==0:
#         x.append(data.iloc[i,0])
#         y.append(data.iloc[i,1])
# s=pd.DataFrame([x,y])
# s=s.transpose()

# clustering = hdbscan.HDBSCAN(min_cluster_size=2,max_cluster_size=5)
# clustering_labels = clustering.fit_predict(s)
# print(clustering_labels)

# clustering = DBSCAN(eps=3.2, min_samples=40).fit(s)
# clustering= AffinityPropagation(damping=0.5, max_iter=100, convergence_iter=15, copy=True, preference=None, affinity='euclidean', verbose=False, random_state=None).fit(s)
# clustering = SpectralClustering(n_clusters=6,
# assign_labels='discretize',
#  random_state=0).fit(s)
# print(s)
# clustering = KMeans(n_clusters=3).fit(s)
# clustering =MeanShift().fit(s)
# clustering_labels= clustering.labels_

# fig=px.scatter(s, x=0, y=1,
#               color=clustering_labels.astype(float))


df=pd.read_csv('output\GlycanIdG_pca.csv')
fig=px.scatter(df, x="X", y="Y",
              color="i")
fig.show()




