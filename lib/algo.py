import random
import pandas as pd
import numpy as np
import pandas as pd
from lib import pdb 
import copy

# Chain is P for protien And G for Test Glycan
def steric_score_function(glycoprotein):
    r=0
    return r

def attach(protein,glycans,glycosylation_locations):
    protein_df= pdb.to_DF(protein)
    Parr=protein_df[['X','Y','Z']].to_numpy()
    glycoprotein_final = copy.deepcopy(protein_df)
    gly=[]
    ChainId= ["B","C","D","E","F","G","H","I"]
    k=0
    for i in glycosylation_locations:
        target_ResId= int(i["begin"])
        OD1 = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'OD1'),['Number']].iloc[0]['Number'] -1
        CG = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'CG'),['Number']].iloc[0]['Number'] -1
        ND2 = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'ND2'),['Number']].iloc[0]['Number'] -1
        G = sampling(glycans[0])
        C1 = G.loc[(G['ResId']==2) & (G['Name']== 'C1'),['Number']].iloc[0]['Number'] -1
        O5 = G.loc[(G['ResId']==2) & (G['Name']== 'O5'),['Number']].iloc[0]['Number'] -1
        Garr = G[['X','Y','Z']].to_numpy()
        Garr = Garr-Garr[C1]
        Garr = Garr + Parr[ND2]
        Gn =  pd.DataFrame(Garr, columns = ['X','Y','Z'])
        G.update(Gn)
        G = G.drop([0,1])
        G["Number"] = glycoprotein_final["Number"].iloc[-1] + G["Number"] 
        G["Chain"] = ChainId[k]
        k+=1
        glycoprotein_final= pd.concat([glycoprotein_final,G])
    return glycoprotein_final

def sampling(Glycanid):
    G = pdb.parse("data/bisecting.pdb")
    return pdb.to_DF(G)

# def shift_glycan(glycan,,theta)
#     return glycan_new