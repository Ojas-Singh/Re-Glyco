import random
import pandas as pd
import numpy as np
from lib import pdb 
import copy
from numba import njit
import os,time
import config
from config import RESIDUE_MAP
import glycors


@njit(fastmath=True)
def fastest_angle(p0,p1,p2):
    cosine_angle = np.dot(np.subtract(p0,p1),np.subtract(p2,p1)) / (np.linalg.norm(np.subtract(p0,p1)) * np.linalg.norm(np.subtract(p2,p1)))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)

def Garrfromtorsion(Garr,torsionpoints,torsionparts,factor):
    for i in range(len(torsionpoints)):
            M1 = glycors.rotation_mat(Garr[torsionpoints[i][2]]-Garr[torsionpoints[i][1]],np.radians(random.uniform(-10*factor, 10*factor)))
            Garr = Garr-Garr[torsionpoints[i][1]]
            for j in np.where(torsionparts[i])[0]:
                Garr[j] = np.dot(M1,Garr[j])
            Garr = Garr+Garr[torsionpoints[i][1]]
    return Garr


def optwithwiggle(GarrM,O1,CB,CG,ND2,C1,O5,Parr,torsionpoints,torsionparts,phisd,psisd,wiggle_tries,wiggle_method):
        r = 100000000
        GarrF= GarrM
        phiF=0
        psiF=0
        if wiggle_method =="none":
            for i in range(wiggle_tries):
                Garr = Garrfromtorsion(GarrM,torsionpoints,torsionparts,factor=(i/wiggle_tries))
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
        elif wiggle_method=="progressive_random":
            for i in range(wiggle_tries):
                Garr = Garrfromtorsion(GarrM,torsionpoints,torsionparts,factor=(i/wiggle_tries))
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
        elif wiggle_method=="random":
            for i in range(wiggle_tries):
                Garr = Garrfromtorsion(GarrM,torsionpoints,torsionparts,factor=1)
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
    return int(filename.split("_")[1].strip("cluster")[0])

def sampling(Glycanid,linkage):
    folder_path = config.data_dir+Glycanid+"/"
    ending =""
    if linkage == "alpha":
        ending="alpha.pdb"
    elif linkage == "beta":
        ending="beta.pdb"
    filenames = os.listdir(folder_path)
    pdb_files = [filename for filename in filenames if filename.endswith(ending)]
    sorted_pdb_files = sorted(pdb_files, key=get_number_after_underscore,reverse=True)
    all_pdbs=[]
    for i in sorted_pdb_files:
        all_pdbs.append(pdb.to_DF(pdb.parse(folder_path+i)))
    loaded = np.load(config.data_dir+Glycanid+"/output/torparts.npz",allow_pickle=True)
    return all_pdbs,loaded

def attach(protein,glycans,glycosylation_locations):
    clash=True
    protein_df= pdb.to_DF(protein)
    glycoprotein_final = copy.deepcopy(protein_df)
    for i in range(len(glycosylation_locations)):
        try: 
            target_ResId= int(glycosylation_locations[i]["begin"])
        except:
            target_ResId=int(glycosylation_locations[i])
        target_Chain = "A"
        resname =protein_df.loc[(protein_df['ResId'] == target_ResId), 'ResName'].iloc[0] 
        residue_data = RESIDUE_MAP.get(resname)
        glycoprotein_final, Parr, box = attach_glycan(glycoprotein_final, glycan= glycans[i], linkage= residue_data, target_ResId = target_ResId,target_Chain = target_Chain)
    return glycoprotein_final,clash



def attach_glycan(glycoprotein_final, glycan, linkage, target_ResId,target_Chain):
    box = ""
    protein_df = glycoprotein_final
    Parr=glycoprotein_final[['X','Y','Z']].to_numpy(dtype=float)
    psisd = linkage["psi"]
    phisd = linkage["phi"]
    link  = linkage["link"]

    CB  = protein_df.loc[(protein_df['Chain']==target_Chain) & (protein_df['ResId']==target_ResId) & (protein_df['Name'].isin(linkage["A"])),['Number']].iloc[0]['Number'] -1
    CG  = protein_df.loc[(protein_df['Chain']==target_Chain) & (protein_df['ResId']==target_ResId) & (protein_df['Name'].isin(linkage["B"])),['Number']].iloc[0]['Number'] -1
    ND2 = protein_df.loc[(protein_df['Chain']==target_Chain) & (protein_df['ResId']==target_ResId) & (protein_df['Name'].isin(linkage["C"])),['Number']].iloc[0]['Number'] -1
    
    if glycan == "None":
        return glycoprotein_final, Parr, box
    else:
        all_G,loaded = sampling(glycan,link)
        G = all_G[0]
        C1 = G.loc[(G['ResId']==2) & (G['Name']== 'C1'),['Number']].iloc[0]['Number'] -1
        O5 = G.loc[(G['ResId']==2) & (G['Name']== 'O5'),['Number']].iloc[0]['Number'] -1
        O1 = G.loc[(G['ResId']==1) & (G['Name']== 'O1'),['Number']].iloc[0]['Number'] -1
        Garr = G[['X','Y','Z']].to_numpy(dtype=float)
        torsionpoints = loaded["a"]
        torsionparts  = loaded["b"]
        clash= True 
        FGarr= Garr
        cluster_best = 0
        rf=100000000
        box+= f"Trying cluster 1/{len(all_G)}..."
        for (clus_num,Gi) in enumerate(all_G):
            Garr = Gi[['X','Y','Z']].to_numpy(dtype=float)
            Garr,r,phi,psi = optwithwiggle(Garr,O1,CB,CG,ND2,C1,O5,Parr,torsionpoints,torsionparts,phisd,psisd,wiggle_tries=1,wiggle_method="none")
            if r<rf:
                FGarr=Garr
                rf=r
                cluster_best = clus_num
            if r==1.0:
                clash=False
                box+= f"Clash Solved for {glycan} at ASN :{target_ResId} with phi : {int(phi)} and psi : {int(psi)}"
                break
            else:
                box+= f"Clash detected, trying cluster {clus_num+2}/{len(all_G)}..."
        if rf != 1.0:
            box+= f"Clash detected, trying cluster {cluster_best+1}/{len(all_G)} with wiggle..."
            Gi = all_G[cluster_best]
            Garr = Gi[['X','Y','Z']].to_numpy(dtype=float)
            Garr,r,phi,psi = optwithwiggle(Garr,O1,CB,CG,ND2,C1,O5,Parr,torsionpoints,torsionparts,phisd,psisd,wiggle_tries=40,wiggle_method="random")
            if r<rf:
                FGarr=Garr
                rf=r
            if r==1.0:
                clash=False
                box+= f"Clash Solved for {glycan} at ASN :{target_ResId} with phi : {int(phi)} and psi : {int(psi)}"
                pass
            else:
                box+= f"Clash detected, trying cluster {clus_num+2}/{len(all_G)} with wiggle..."
        if clash:
            box+= f"Clash exist for {glycan} at ASN :{target_ResId} with phi : {int(phi)} and psi : {int(psi)}"
        Gn =  pd.DataFrame(FGarr, columns = ['X','Y','Z'])
        G.update(Gn)
        G = G.drop([0,1])
        G["Number"] = glycoprotein_final["Number"].iloc[-1] + G["Number"] 
        G["Chain"] = "G"
        glycoprotein_final= pd.concat([glycoprotein_final,G])
        Parr=glycoprotein_final[['X','Y','Z']].to_numpy(dtype=float)

        return glycoprotein_final, Parr, box