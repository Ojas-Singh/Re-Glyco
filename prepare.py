import numpy as np
import pandas as pd
from lib import pdb,graph,dihedral,clustering,tfindr
import config
from scipy import stats



name="glycan"
f=config.data_dir+name+"/"+name+".dry.pdb"

pdbdata, frames = pdb.multi(f)
df = pdb.to_DF(pdbdata)


# idx_noH=df.loc[(df['Element']!="H"),['Number']].iloc[:]['Number']-1
# pcaG= clustering.pcawithG(frames,idx_noH,config.number_of_dimensions)
# pcaG.to_csv(config.data_dir+name+"/"+name+'_G_pca_NB.csv',index_label="i")

# tsneG = clustering.tsnewithG(frames,idx_noH,config.number_of_dimensions)
# tsneG.to_csv(config.data_dir+name+"/"+name+'_G_tsne.csv',index_label="i")  

pairs,external,internal = tfindr.torsionspairs(name)
pairs = np.asarray(pairs)

torsion_names = dihedral.pairtoname(external,df)
ext_DF = dihedral.pairstotorsion(external,frames,torsion_names)
for i in range(len(internal)):
    torsion_names.append("internal")

torsiondataDF= dihedral.pairstotorsion(pairs,frames,torsion_names)
torsiondataDF.to_csv(config.data_dir+name+"/"+name+'_torsions_full.csv',index_label="i")

# df = df[(np.abs(stats.zscore(df.loc[:, df.columns!='i'])) < 3)]
k=0
outliers=[]
ok=[]
for i in np.asarray(np.abs(stats.zscore(ext_DF)) < 2,dtype=bool):
    p=True
    for j in i:
        if not j:
            p=False
    if p :
        ok.append(k)
    else:
        outliers.append(k)
    k+=1
print("Outliers : ",len(outliers), "Total points : ",len(ext_DF.iloc[:]))
ext_DF["cluster"]=["0" for x in range(len(ext_DF.iloc[:]))]
ext_DF["cluster"][outliers]=str(-1)

exDF_clean = ext_DF.loc[ext_DF["cluster"]!="-1"]
exDF = exDF_clean.loc[:, exDF_clean .columns!='cluster']
print(len(ok))


# tor=clustering.normalizetorsion(exDF)
# pcaT = clustering.pcawithT(tor,config.number_of_dimensions)
# pcaT["i"]=ok
# pcaT.to_csv(config.data_dir+name+"/"+name+'_T_pca.csv',index_label=False)  

# tsneT = clustering.tsnewithT(tor,config.number_of_dimensions)
# tsneT["i"]=ok
# tsneT.to_csv(config.data_dir+name+"/"+name+'_T_tsne.csv',index_label=False) 


exDF.to_csv(config.data_dir+name+"/"+name+'_torsions.csv',index_label="i")



#plots all diherdal KDEs
# dihedral.kdemax(torsiondataDF)

