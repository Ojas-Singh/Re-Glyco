import random
import pandas as pd
import numpy as np
from lib import pdb 
import copy
from numba import njit
import streamlit as st
import os,time
import config
import glycors


@njit(fastmath=True)
def fastest_angle(p0,p1,p2):
    cosine_angle = np.dot(np.subtract(p0,p1),np.subtract(p2,p1)) / (np.linalg.norm(np.subtract(p0,p1)) * np.linalg.norm(np.subtract(p2,p1)))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)

def Garrfromtorsion(Garr,torsionpoints,torsionparts):
    for i in range(len(torsionpoints)):
            M1 = glycors.rotation_mat(Garr[torsionpoints[i][2]]-Garr[torsionpoints[i][1]],np.radians(random.uniform(-10, 10)))
            Garr = Garr-Garr[torsionpoints[i][1]]
            for j in np.where(torsionparts[i])[0]:
                Garr[j] = np.dot(M1,Garr[j])
            Garr = Garr+Garr[torsionpoints[i][1]]
    return Garr


def optwithwiggle(GarrM,O1,CB,CG,ND2,C1,O5,Parr,torsionpoints,torsionparts,phisd,psisd):
        r = 100000000
        GarrF= GarrM
        phiF=0
        psiF=0
        for i in range(20):
            Garr = Garrfromtorsion(GarrM,torsionpoints,torsionparts)
            Garr = Garr-Garr[O1]
            Garr = Garr + Parr[ND2]
            axis = np.cross(Parr[CG]-Parr[ND2],Garr[C1]-Parr[ND2])
            an=fastest_angle(Parr[CG],Parr[ND2],Garr[C1])
            theta = np.radians(random.uniform(130, 110) - an)
            Garr = Garr - Parr[ND2]
            M0 = glycors.rotation_mat(axis, theta)
            for i in range(len(Garr)):
                Garr[i] = np.dot(M0,Garr[i])
            Garr = Garr + Parr[ND2]
            phi,psi,ri = glycors.opt_genetic(CB,CG,ND2,C1,O5,Garr,Parr,psisd,phisd)
            if ri<r:
                GarrF= Garr
                phiF=phi
                psiF=psi
                r=ri
            if ri==1.0:
                break
        return glycors.adjust_dihedrals(phiF,psiF,CB,CG,ND2,C1,O5,GarrF,Parr),r,phiF,psiF


def get_number_after_underscore(filename):
    return float(filename.split("_")[1].split(".")[0])

def sampling(Glycanid,linkage):
    folder_path = ""
    if linkage == "alpha":
        folder_path = config.data_dir+Glycanid+"/clusters/alpha/"
    elif linkage == "beta":
        folder_path = config.data_dir+Glycanid+"/clusters/beta/"
    filenames = os.listdir(folder_path)
    pdb_files = [filename for filename in filenames if filename.endswith(".pdb")]
    sorted_pdb_files = sorted(pdb_files, key=get_number_after_underscore,reverse=True)
    all_pdbs=[]
    for i in sorted_pdb_files:
        all_pdbs.append(pdb.to_DF(pdb.parse(folder_path+i)))
    loaded = np.load(config.data_dir+Glycanid+"/output/torparts.npz",allow_pickle=True)
    return all_pdbs,loaded

def attach(protein,glycans,glycosylation_locations):
    protein_df= pdb.to_DF(protein)
    Parr=protein_df[['X','Y','Z']].to_numpy(dtype=float)
    glycoprotein_final = copy.deepcopy(protein_df)
    timer=[]
    ChainId= ["B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
    k=0
    for i in range(len(glycosylation_locations)):
        try: 
            target_ResId= int(glycosylation_locations[i]["begin"])
        except:
            target_ResId=int(glycosylation_locations[i])
        if protein_df.loc[(protein_df['ResId'] == target_ResId), 'ResName'].iloc[0] in config.O_linked["Res"]:
            psisd = config.O_linked["psi"]
            phisd = config.O_linked["phi"]
            CB = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'CB'),['Number']].iloc[0]['Number'] -1
            CG = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'CA'),['Number']].iloc[0]['Number'] -1
            ND2 = protein_df.loc[(protein_df['ResId']==target_ResId) & ((protein_df['Name']== 'OG1')|(protein_df['Name']== 'OG')),['Number']].iloc[0]['Number'] -1
            if glycans[i] == "None":
                continue
            else:
                all_G,loaded = sampling(glycans[i],"alpha")
                G = all_G[0]
                C1 = G.loc[(G['ResId']==2) & (G['Name']== 'C1'),['Number']].iloc[0]['Number'] -1
                O5 = G.loc[(G['ResId']==2) & (G['Name']== 'O5'),['Number']].iloc[0]['Number'] -1
                O1 = G.loc[(G['ResId']==1) & (G['Name']== 'O1'),['Number']].iloc[0]['Number'] -1
                Garr = G[['X','Y','Z']].to_numpy(dtype=float)
                torsionpoints = loaded["a"]
                torsionparts  = loaded["b"]
                clash= True 
                s=time.time()
                box = st.empty()
                FGarr= Garr
                rf=100000000
                box.info(f"Trying cluster 1/{len(all_G)}...")
                for (clus_num,Gi) in enumerate(all_G):
                    Garr = Gi[['X','Y','Z']].to_numpy(dtype=float)
                    Garr,r,phi,psi = optwithwiggle(Garr,O1,CB,CG,ND2,C1,O5,Parr,torsionpoints,torsionparts,phisd,psisd)
                    if r<rf:
                        FGarr=Garr
                        rf=r
                    if r==1.0:
                        clash=False
                        box.info(f"Clash Solved for {glycans[i]} at ASN :{target_ResId} with phi : {int(phi)} and psi : {int(psi)}")
                        break
                    else:
                        box.info(f"Clash detected, trying cluster {clus_num+2}/{len(all_G)}...")
                if clash:
                    box.warning(f"Clash exist for {glycans[i]} at ASN :{target_ResId} with phi : {int(phi)} and psi : {int(psi)}")

                timer.append(time.time()-s)
                Gn =  pd.DataFrame(FGarr, columns = ['X','Y','Z'])
                G.update(Gn)
                G = G.drop([0,1])
                G["Number"] = glycoprotein_final["Number"].iloc[-1] + G["Number"] 
                G["Chain"] = ChainId[k]
                k+=1
                glycoprotein_final= pd.concat([glycoprotein_final,G])
                Parr=glycoprotein_final[['X','Y','Z']].to_numpy(dtype=float)
        if protein_df.loc[(protein_df['ResId'] == target_ResId), 'ResName'].iloc[0] in config.N_linked["Res"]:
            psisd = config.N_linked["psi"]
            phisd = config.N_linked["phi"]
            CB = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'CB'),['Number']].iloc[0]['Number'] -1
            CG = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'CG'),['Number']].iloc[0]['Number'] -1
            ND2 = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'ND2'),['Number']].iloc[0]['Number'] -1
            if glycans[i] == "None":
                continue
            else:
                all_G,loaded = sampling(glycans[i],"beta")
                G = all_G[0]
                C1 = G.loc[(G['ResId']==2) & (G['Name']== 'C1'),['Number']].iloc[0]['Number'] -1
                O5 = G.loc[(G['ResId']==2) & (G['Name']== 'O5'),['Number']].iloc[0]['Number'] -1
                O1 = G.loc[(G['ResId']==1) & (G['Name']== 'O1'),['Number']].iloc[0]['Number'] -1
                Garr = G[['X','Y','Z']].to_numpy(dtype=float)
                torsionpoints = loaded["a"]
                torsionparts  = loaded["b"]
                clash= True 
                s=time.time()
                box = st.empty()
                FGarr= Garr
                rf=100000000
                box.info(f"Trying cluster 1/{len(all_G)}...")
                for (clus_num,Gi) in enumerate(all_G):
                    Garr = Gi[['X','Y','Z']].to_numpy(dtype=float)
                    Garr,r,phi,psi = optwithwiggle(Garr,O1,CB,CG,ND2,C1,O5,Parr,torsionpoints,torsionparts,phisd,psisd)
                    if r<rf:
                        FGarr=Garr
                        rf=r
                    if r==1.0:
                        clash=False
                        box.info(f"Clash Solved for {glycans[i]} at ASN :{target_ResId} with phi : {int(phi)} and psi : {int(psi)}")
                        break
                    else:
                        box.info(f"Clash detected, trying cluster {clus_num+2}/{len(all_G)}...")
                if clash:
                    box.warning(f"Clash exist for {glycans[i]} at ASN :{target_ResId} with phi : {int(phi)} and psi : {int(psi)}")

                timer.append(time.time()-s)
                Gn =  pd.DataFrame(FGarr, columns = ['X','Y','Z'])
                G.update(Gn)
                G = G.drop([0,1])
                G["Number"] = glycoprotein_final["Number"].iloc[-1] + G["Number"] 
                G["Chain"] = ChainId[k]
                k+=1
                glycoprotein_final= pd.concat([glycoprotein_final,G])
                Parr=glycoprotein_final[['X','Y','Z']].to_numpy(dtype=float)

    return glycoprotein_final,clash



