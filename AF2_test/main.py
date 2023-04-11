import sys
sys.path.append('../')
from lib import prep,algo,pdb
import time,os,config
from io import StringIO
import pandas as pd
import numpy as np
import copy

def append_to_results(filename, uni_id, target_ResId, timer,r,glycan):
    with open(filename, "a") as f:
        tt = uni_id.strip('\n')
        f.write(f"uniprot_id: {tt} , spot: {target_ResId}, time: {timer}, steric: {r}, Glycan: {glycan} \n")





phisd=(-130,-63)
psisd=(152,205)

test_glycans=[
    "DGalpb1-4DGlcpNAcb1-2DManpa1-3[DGalpb1-4DGlcpNAc[6S]b1-2DManpa1-6]DManpb1-4DGlcpNAcb1-4[LFucpa1-6]DGlcpNAca1-OH",

]

def attach_single(protein,glycan,target_ResId,phisd,psisd):
    protein_df= pdb.to_DF(protein)
    Parr=protein_df[['X','Y','Z']].to_numpy(dtype=float)
    glycoprotein_final = copy.deepcopy(protein_df)
    timer=[]
    ChainId= ["B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
    k=0
    CB = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'CB'),['Number']].iloc[0]['Number'] -1
    CG = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'CG'),['Number']].iloc[0]['Number'] -1
    ND2 = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'ND2'),['Number']].iloc[0]['Number'] -1
    all_G,loaded = algo.sampling(glycan)
    G = all_G[0]
    C1 = G.loc[(G['ResId']==2) & (G['Name']== 'C1'),['Number']].iloc[0]['Number'] -1
    O5 = G.loc[(G['ResId']==2) & (G['Name']== 'O5'),['Number']].iloc[0]['Number'] -1
    O1 = G.loc[(G['ResId']==1) & (G['Name']== 'O1'),['Number']].iloc[0]['Number'] -1
    Garr = G[['X','Y','Z']].to_numpy(dtype=float)
    torsionpoints = loaded["a"]
    torsionparts  = loaded["b"]
    clash= False 
    s=time.time()
    FGarr= Garr
    rf=100000000
    for Gi in all_G:
        Garr = Gi[['X','Y','Z']].to_numpy(dtype=float)
        Garr,r = algo.optwithwiggle(Garr,O1,CB,CG,ND2,C1,O5,Parr,torsionpoints,torsionparts,phisd,psisd)
        if r<rf:
            FGarr=Garr
            rf=r
        if r==1.0:
            print("wow!")
            break
        else:
            print()
    timer.append(time.time()-s)
    print("Spot : ",target_ResId," Psi : ",int(algo.fastest_dihedral(Parr[CB],Parr[CG],Parr[ND2],FGarr[C1]))," Phi : ",int(algo.fastest_dihedral(Parr[CG],Parr[ND2],FGarr[C1],FGarr[O5])))
    Gn =  pd.DataFrame(FGarr, columns = ['X','Y','Z'])
    G.update(Gn)
    G = G.drop([0,1])
    G["Number"] = glycoprotein_final["Number"].iloc[-1] + G["Number"] 
    G["Chain"] = ChainId[k]
    k+=1
    glycoprotein_final= pd.concat([glycoprotein_final,G])
    Parr=glycoprotein_final[['X','Y','Z']].to_numpy(dtype=float)

    return glycoprotein_final,rf,timer




def fetch(uni_id):
    fold = prep.download_and_prepare_alphafoldDB_model(uni_id,"temp/")
    out= prep.query_uniprot_for_glycosylation_locations(uni_id)
    return fold,out

with open('important.txt', "r") as f:
    content = f.readlines()
    print(len(content))
    for uni_id in content:
        fold, out = fetch(uni_id.strip("\n"))
        with open(fold) as ifile:
            system = "".join([x for x in ifile])
            protein = pdb.parse(fold)
            confidence= []
            p=1
            lines = system.split("\n")
            for x in lines:
                if x.startswith("ATOM"):
                    if int((x[23:27]).strip(" "))==p:
                        confidence.append(float((x[61:67]).strip(" ")))
                        p+=1
            glycosylation_locations = out["glycosylations"]
            glycosylation_locations_N=[]
            for i in range(len(glycosylation_locations)):
                if not glycosylation_locations[i]["description"].startswith('N-linked'):
                    pass
                else:
                    glycosylation_locations_N.append(glycosylation_locations[i])
            dirlist = [ item for item in os.listdir(config.data_dir) if os.path.isdir(os.path.join(config.data_dir, item)) ]
            for target in glycosylation_locations_N:
                target_ResId= int(target["begin"])
                for k in test_glycans:
                    g,clash,timer = attach_single(protein,k,target_ResId,phisd,psisd)
                    append_to_results("result.txt", uni_id, target_ResId, timer,clash,k)