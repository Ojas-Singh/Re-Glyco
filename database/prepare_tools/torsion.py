import numpy as np
import pandas as pd
from lib import pdb,graph,dihedral
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import plotly.express as px


name="GlycanId"
df = pd.read_csv('output/'+name+'torsions.csv')

#plots all diherdal KDEs
# dihedral.kdemax(df)


# corr = df.corr().abs()
# print(corr.values[np.triu_indices_from(corr.values,1)].mean())

# plt.matshow(corr, cmap='magma')
# plt.colorbar()
# plt.show()

pca_df=pd.read_csv('output\GlycanIdG_pca_full.csv')

# clustering = KMeans(n_clusters=10).fit(pca_df[['0','1','2','3','4','5']])
# clustering_labels= clustering.labels_
# pca_df.insert(1,"cluster",clustering_labels,False)
# df1 = pca_df.loc[pca_df["cluster"] ==0]
# print(df1)


# fig=px.scatter(pca_df, x='0', y='1',
#               color=clustering_labels.astype(float))

# fig.show()

for i in range(1,20):
    clustering = KMeans(n_clusters=i).fit(pca_df[['0','1']])
    clustering_labels= clustering.labels_
    
    p=pca_df.copy()
    p.insert(1,"cluster",clustering_labels,False)
    kk=0
    for j in range(0,i):
        df1 = p.loc[p["cluster"] ==j]
        df1=df1.loc[:, df1.columns != 'cluster']
        df1=df1.loc[:, df1.columns != 'i']
        corr = df1.corr().abs()
        # corr = df1.cov().abs()
        # print(i,j)
        # plt.matshow(corr, cmap='magma')
        # plt.colorbar()
        # plt.show()
        # print(corr.unstack().mean())
        kk+=corr.unstack().mean()
        # print(corr.values[np.triu_indices_from(corr.values,1)].mean())
    print(i,kk/i)
    
