import numpy as np
import pandas as pd
from lib import pdb,graph,dihedral
import pickle
import matplotlib.pyplot as plt
loaded = np.load('data/bisecting.npz',allow_pickle=True)
# print(loaded['a'])
a=loaded['f']
print(np.where(a[0])[0])

# import networkx as nx
# connect =[]
# f2 = "data/bi.txt"
# with open(f2, 'r') as f:
#             lines = f.readlines()
#             i=1
#             for line in lines:
#                 a = int((line[8:12]).strip(" "))
#                 b = int((line[13:17]).strip(" "))
#                 connect.append((a,b))

# for i in loaded['d']:
#     G = nx.Graph()
#     G.add_edges_from(connect)
#     print(nx.is_connected(G))

#     G.remove_edge(i[1]+1, i[2]+1)
#     print(nx.is_connected(G))
#     print(nx.node_connected_component(G, i[1]+1))
#     nx.draw_spring(G,with_labels=True)
#     plt.draw()
#     plt.show()