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
                l = [x for x in l if x != node]
                pairs.append([l[0],j,node,k])
            if G.degree(k)>1:
                m = list(G.neighbors(k))
                m = [x for x in m if x != node]
                pairs.append([m[0],k,node,j])
        # elif df.loc[(df['Number']==node),['Element']].iloc[0]['Element']=="O":
        #     color_map.append('green')
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

import numpy as np
from lib import pdb,graph,dihedral
import matplotlib.pyplot as plt
import networkx as nx



def torsionspairs(name):
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
    external=[]
    red=[]
    internal=[]
    for node in G:
        if G.degree(node) ==2 and node not in cycle:
            color_map.append('red')
            red.append(node)
        # elif df.loc[(df['Number']==node),['Element']].iloc[0]['Element']=="O":
        #     color_map.append('green')
        else: 
            color_map.append('blue') 
    for node in red:
        j = list(G.neighbors(node))[0]
        k = list(G.neighbors(node))[1]
        if (j in red) or (k in red):
            if  df.loc[(df['Number']==node),['Name']].iloc[0]['Name']=="O6":
                if j in red:
                    l = list(G.neighbors(j))
                    l = [x for x in l if x != node]
                    m = list(G.neighbors(k))
                    m = [x for x in m if x != node]
                    aa = df.loc[(df['Number']==m[0]),['Name']].iloc[0]['Name']
                    external.append([l[0],j,node,k])
                    if (df.loc[(df['Number']==k),['Name']].iloc[0]['Name'] =="C1" ) or (df.loc[(df['Number']==k),['Name']].iloc[0]['Name'] =="C2"):
                        if aa=="O5" or aa=="C1":
                            external.append([m[0],k,node,j])
                        else:
                            external.append([m[1],k,node,j])
                    external.append([df.loc[(df['ResId']==df.loc[(df['Number']==l[0]),['ResId']].iloc[0]['ResId']) & (df['Name']== "C4"),['Number']].iloc[0]['Number'],l[0],j,node])
                    
                if k in red:
                    l = list(G.neighbors(k))
                    l = [x for x in l if x != node]
                    m = list(G.neighbors(j))
                    m = [x for x in m if x != node]
                    aa = df.loc[(df['Number']==m[0]),['Name']].iloc[0]['Name']
                    external.append([l[0],k,node,j])
                    if (df.loc[(df['Number']==j),['Name']].iloc[0]['Name'] =="C1") or (df.loc[(df['Number']==j),['Name']].iloc[0]['Name'] =="C2"):
                        if aa=="O5" or aa=="C1":
                            external.append([m[0],j,node,k])
                        else:
                            external.append([m[1],j,node,k])
                    external.append([df.loc[(df['ResId']==df.loc[(df['Number']==l[0]),['ResId']].iloc[0]['ResId']) & (df['Name']== "C4"),['Number']].iloc[0]['Number'],l[0],k,node])
        else:
            if G.degree(j)>1:
                l = list(G.neighbors(j))
                l = [x for x in l if x != node]
                if (l[1]in cycle or l[0] in cycle) and k in cycle:
                    if df.loc[(df['Number']==j),['Name']].iloc[0]['Name']=="C1":
                        if df.loc[(df['Number']==l[0]),['Name']].iloc[0]['Name']=="O5":
                            external.append([l[0],j,node,k])
                        else:
                            external.append([l[1],j,node,k])
                    elif G.degree(j)==4:
                        external.append([df.loc[(df['ResId']==df.loc[(df['Number']==j),['ResId']].iloc[0]['ResId']) & (df['Name']== "C1"),['Number']].iloc[0]['Number'],j,node,k])
                    else:
                        jj = df.loc[(df['Number']==j),['Name']].iloc[0]['Name']
                        external.append([k,node,j,df.loc[(df['ResId']==df.loc[(df['Number']==j),['ResId']].iloc[0]['ResId']) & (df['Name']== "C"+str(int(jj.strip("C"))+1)),['Number']].iloc[0]['Number']])    
                else:
                    internal.append([l[0],j,node,k])
            if G.degree(k)>1:
                m = list(G.neighbors(k))
                m = [x for x in m if x != node]
                if (m[1] in cycle or m[0] in cycle) and j in cycle:
                    if df.loc[(df['Number']==k),['Name']].iloc[0]['Name']=="C1":
                        if df.loc[(df['Number']==m[0]),['Name']].iloc[0]['Name']=="O5":
                            external.append([m[0],k,node,j])
                        else:
                            external.append([m[1],k,node,j])
                    elif G.degree(k)==4:
                        external.append([df.loc[(df['ResId']==df.loc[(df['Number']==k),['ResId']].iloc[0]['ResId']) & (df['Name']== "C1"),['Number']].iloc[0]['Number'],k,node,j])
                    else:
                        kk= df.loc[(df['Number']==k),['Name']].iloc[0]['Name']
                        external.append([j,node,k,df.loc[(df['ResId']==df.loc[(df['Number']==k),['ResId']].iloc[0]['ResId']) & (df['Name']== "C"+str(int(kk.strip("C"))+1)),['Number']].iloc[0]['Number']])    
                
                else:
                    internal.append([m[0],k,node,j])
    nx.draw(G, pos=nx.kamada_kawai_layout(G),node_color=color_map, with_labels = True)
    plt.show()

    anotherlist=[]
    pairs = external + internal

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

    
    return pairs,external,internal



