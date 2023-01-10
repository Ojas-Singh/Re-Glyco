import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from lib import pdb,graph,dihedral


name="bisecting"
f="data/prodfull.dry.pdb"
pdbdata, frames = pdb.multi(f)
pdbdataDF = pdb.to_DF(pdbdata)



tormeta=[[2,3,4,4],[3,4,4,4],[4,5,3,3],[4,7,4,4],[4,8,6,6],[8,9,2,2],[5,6,2,2]]
torsionmeta = []
for i in tormeta:
    torsionmeta.append(dihedral.res2input(i[0],i[1],i[2],i[3],pdbdataDF))


adata = np.asanyarray(pdbdataDF,dtype=object)
bdata = np.asanyarray(torsionmeta,dtype=object)
torsiondataDF= dihedral.torsions(torsionmeta,pdbdata,frames,f)
torsions=torsiondataDF.to_numpy()  #c
torsiondataDF.to_csv('output/'+name+'torsions.csv')


import networkx as nx
connect =[]
f2 = "data/bi.txt"
with open(f2, 'r') as f:
            lines = f.readlines()
            i=1
            for line in lines:
                a = int((line[8:12]).strip(" "))
                b = int((line[13:17]).strip(" "))
                connect.append((a,b))

anotherlist=[]
anotherlist2=[]
torsionpoints=[]
for i in torsionmeta:
    torsionpoints.append(i[1])
    torsionpoints.append(i[2])
    if len(i)==4:
        torsionpoints.append(i[3])
for i in torsionpoints:
    # Create a sample graph
    G = nx.Graph()
    G.add_edges_from(connect)

    # Choose the nodes to split the graph
    node1 = i[1]+1
    node2 = i[2]+1

    # Remove the edge between the two nodes
    G.remove_edge(node1, node2)
    arr = np.ones(len(pdbdata[0]),dtype=bool) 
    # Split the graph into two parts
    for i in nx.node_connected_component(G, node1):
        arr[i-1]=False
    anotherlist.append([[x-1 for x in list(nx.node_connected_component(G, node1))],[x-1 for x in list(nx.node_connected_component(G, node2))]])
    anotherlist2.append(arr)
print(anotherlist[0])
anotherlist = np.asanyarray(anotherlist,dtype=object)
np.savez_compressed("data/"+name, a = adata, b = bdata, c = torsions, d = torsionpoints, e = anotherlist,f = anotherlist2,allow_pickle=True)





#plots all diherdal KDEs
# dihedral.kdemax(torsiondataDF)







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



