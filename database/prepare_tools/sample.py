import numpy as np
import pandas as pd
from lib import pdb,graph,dihedral
import pickle


G = pdb.parse("../../data/bisecting.pdb")
G = pdb.to_DF(G)
G = np.asanyarray(G,dtype=object)
df = pd.DataFrame(G, columns = ['Number','Name','ResName','Chain','ResId','X','Y','Z','Element'])

print(df)
loaded = np.load('file.npz',allow_pickle=True)
print(loaded['a'][0])