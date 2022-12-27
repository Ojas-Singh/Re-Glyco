import numpy as np
import pandas as pd
from lib import pdb,graph,dihedral
import matplotlib.pyplot as plt

name="GlycanId"
df = pd.read_csv('output/'+name+'torsions.csv')

#plots all diherdal KDEs
# dihedral.kdemax(df)
corr = df.corr().abs()
# corr.style.background_gradient(cmap='magma').set_precision(2)

plt.matshow(corr, cmap='magma')
plt.colorbar()
plt.show()

