import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from lib import pdb,graph,dihedral
import matplotlib.pyplot as plt
import networkx as nx
from mpl_toolkits.mplot3d import Axes3D
from random import randint
import networkx as nx
import igraph as ig


name="A3G3S1-F3"
f="data/"+name+"/"+name+".pdb"
f2 = "data/"+name+"/"+"graph"


pdbdata = pdb.parse(f)
df = pdb.to_DF(pdbdata)



connect =[]


with open(f2, 'r') as f:
            lines = f.readlines()
            i=1
            for line in lines:
                a = int((line[8:13]).strip(" "))
                b = int((line[13:19]).strip(" "))
                connect.append((a,b))

G = nx.Graph()
G.add_edges_from(connect)
H_list = df.loc[(df['Element']=="H"),['Number']].iloc[:]['Number']
# print(H_list)
for i in H_list:
    G.remove_node(i)

# print(list(nx.find_cycle(G, orientation="ignore")))
cycle=[]
for i in nx.cycle_basis(G):
    for j in i:
        cycle.append(j)
color_map = []
pairs=[]
for node in G:
    if G.degree(node) ==2 and node not in cycle:
        color_map.append('red')
        j = list(G.neighbors(node))[0]
        k = list(G.neighbors(node))[1]
        if G.degree(j)>1:
            pairs.append([list(G.neighbors(j))[0],j,node,k])
        if G.degree(k)>1:
            pairs.append([list(G.neighbors(k))[0],k,node,j])
    elif df.loc[(df['Number']==node),['Element']].iloc[0]['Element']=="O":
        color_map.append('green')
    else: 
        color_map.append('blue') 

nx.draw(G, pos=nx.kamada_kawai_layout(G),node_color=color_map, with_labels = True)
plt.show()

anotherlist=[]


for i in pairs:
    G = nx.Graph()
    G.add_edges_from(connect)
    nx.draw(G)
    node1 = i[1]
    node2 = i[2]
    G.remove_edge(node1, node2)
    arr = np.ones(len(pdbdata[0]),dtype=bool) 
    for i in nx.node_connected_component(G, node1):
        arr[i-1]=False
    # anotherlist.append([[x-1 for x in list(nx.node_connected_component(G, node1))],[x-1 for x in list(nx.node_connected_component(G, node2))]])
    anotherlist.append(arr)
# print(anotherlist[0])
# anotherlist = np.asanyarray(anotherlist,dtype=object)
np.savez_compressed("data/"+name+"/"+name+"_torparts", a = pairs,b = anotherlist, allow_pickle=True)


