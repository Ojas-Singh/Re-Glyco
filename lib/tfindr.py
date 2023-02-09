import numpy as np
from lib import pdb,graph,dihedral
import matplotlib.pyplot as plt
import networkx as nx



def savetotorparts(name):
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
    for i in H_list:
        G.remove_node(i)
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
                l = list(G.neighbors(j))
                l = [x for x in m if x != node]
                pairs.append([l[0],j,node,k])
            if G.degree(k)>1:
                m = list(G.neighbors(k))
                m = [x for x in m if x != node]
                pairs.append([m[0],k,node,j])
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
        anotherlist.append(arr)
    np.savez_compressed("data/"+name+"/"+name+"_torparts", a = pairs,b = anotherlist, allow_pickle=True)
    return pairs

