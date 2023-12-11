import random
import pandas as pd
import numpy as np
from lib import pdb 
import copy
import re
# from numba import njit
import os,time
import config
from config import RESIDUE_MAP
import glycors


# @njit(fastmath=True)
# def glycors.fastest_angle(p0,p1,p2):
#     cosine_angle = np.dot(np.subtract(p0,p1),np.subtract(p2,p1)) / (np.linalg.norm(np.subtract(p0,p1)) * np.linalg.norm(np.subtract(p2,p1)))
#     angle = np.arccos(cosine_angle)
#     return np.degrees(angle)

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
                Garr = Garrfromtorsion(GarrM,torsionpoints,torsionparts,factor=0.1)
                Garr = Garr-Garr[O1]
                Garr = Garr + Parr[ND2]
                axis = np.cross(Parr[CG]-Parr[ND2],Garr[C1]-Parr[ND2])
                an=glycors.fastest_angle(Parr[CG],Parr[ND2],Garr[C1])
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
                an=glycors.fastest_angle(Parr[CG],Parr[ND2],Garr[C1])
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
                an=glycors.fastest_angle(Parr[CG],Parr[ND2],Garr[C1])
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
    folder_path = config.data_dir+Glycanid+"/PDB_format/"
    ending =""
    if linkage == "alpha":
        ending="alpha.PDB.pdb"
    elif linkage == "beta":
        ending="beta.PDB.pdb"
    filenames = os.listdir(folder_path)
    pdb_files = [filename for filename in filenames if filename.endswith(ending)]
    sorted_pdb_files = sorted(pdb_files, key=get_number_after_underscore,reverse=True)
    all_pdbs=[]
    for i in sorted_pdb_files:
        all_pdbs.append(pdb.to_DF(pdb.parse_gly(folder_path+i)))
    loaded = np.load(config.data_dir+Glycanid+"/output/torparts.npz",allow_pickle=True)
    return all_pdbs,loaded

def attach_skip(protein,glycans,glycosylation_locations):
    link_pairs =[]
    clashes=[]
    protein_df= pdb.to_DF(protein)
    glycoprotein_final = copy.deepcopy(protein_df)
    s = time.time()
    Box = f"Calculation started\n"
    last_chain = protein_df['Chain'].iloc[-1]
    currentGlycanChain = chr(ord(last_chain) + 1)
    for i in range(len(glycosylation_locations)):
        try: 
            target_ResId= int(glycosylation_locations[i]["begin"])
            target_Chain = "A"
        except:
            target_ResId=int(glycosylation_locations[i].split("_")[0])
            target_Chain = str(glycosylation_locations[i].split("_")[1])
        resname =protein_df.loc[(protein_df['ResId'] == target_ResId) & (protein_df['Chain'] == target_Chain), 'ResName'].iloc[0] 
        residue_data = RESIDUE_MAP.get(resname)
        glycoprotein_final_, Parr, box ,clash ,link_pair = attach_glycan(glycoprotein_final, glycan= glycans[i], linkage= residue_data, target_ResId = target_ResId,target_Chain = target_Chain,glycan_chain=currentGlycanChain)
        Box += box
        clashes.append(clash)
        if not clash:
            glycoprotein_final = glycoprotein_final_
            currentGlycanChain = chr(ord(currentGlycanChain) + 1) 
            link_pairs.append(link_pair)
        else:
            pass
    Box += f"Finished in {time.time() -s } seconds"
    return glycoprotein_final,any(clashes), Box,link_pairs




def attach(protein,glycans,glycosylation_locations):
    link_pairs=[]
    clashes=[]
    protein_df= pdb.to_DF(protein)
    glycoprotein_final = copy.deepcopy(protein_df)
    s = time.time()
    Box = f"Calculation started\n"
    last_chain = protein_df['Chain'].iloc[-1]
    currentGlycanChain = chr(ord(last_chain) + 1)
    for i in range(len(glycosylation_locations)):
        try: 
            target_ResId= int(glycosylation_locations[i]["begin"])
            target_Chain = "A"
        except:
            target_ResId=int(glycosylation_locations[i].split("_")[0])
            target_Chain = str(glycosylation_locations[i].split("_")[1])
        resname =protein_df.loc[(protein_df['ResId'] == target_ResId) & (protein_df['Chain'] == target_Chain), 'ResName'].iloc[0] 
        residue_data = RESIDUE_MAP.get(resname)
        glycoprotein_final, Parr, box ,clash, link_pair = attach_glycan(glycoprotein_final, glycan= glycans[i], linkage= residue_data, target_ResId = target_ResId,target_Chain = target_Chain,glycan_chain=currentGlycanChain)
        link_pairs.append(link_pair)
        Box += box
        clashes.append(clash)
        currentGlycanChain = chr(ord(currentGlycanChain) + 1) 
    Box += f"Finished in {time.time() -s } seconds"
    return glycoprotein_final,any(clashes), Box ,link_pairs



def attach_glycan(glycoprotein_final, glycan, linkage, target_ResId,target_Chain,glycan_chain):
    
    protein_df = glycoprotein_final
    
    Parr=glycoprotein_final[['X','Y','Z']].to_numpy(dtype=float)
    protein_df.to_csv('df.csv')
    # print(protein_df.to_string())
    
    last_sugar_residue =  re.split(r'\)|\]', glycan)[-1]
    psisd = linkage["sugars"][last_sugar_residue]["psi"]
    phisd = linkage["sugars"][last_sugar_residue]["phi"]
    link  = linkage["sugars"][last_sugar_residue]["link"]

    # CB  = protein_df.loc[(protein_df['Chain']==target_Chain) & (protein_df['ResId']==target_ResId) & (protein_df['Name'].isin(linkage["A"])),['Number']].iloc[0]['Number'] -1
    # CG  = protein_df.loc[(protein_df['Chain']==target_Chain) & (protein_df['ResId']==target_ResId) & (protein_df['Name'].isin(linkage["B"])),['Number']].iloc[0]['Number'] -1
    # ND2 = protein_df.loc[(protein_df['Chain']==target_Chain) & (protein_df['ResId']==target_ResId) & (protein_df['Name'].isin(linkage["C"])),['Number']].iloc[0]['Number'] -1
    
    CB = protein_df.index[(protein_df['Chain'] == target_Chain) & (protein_df['ResId'] == target_ResId) & (protein_df['Name'].isin(linkage["A"]))][0]
    CG = protein_df.index[(protein_df['Chain'] == target_Chain) & (protein_df['ResId'] == target_ResId) & (protein_df['Name'].isin(linkage["B"]))][0]
    ND2 = protein_df.index[(protein_df['Chain'] == target_Chain) & (protein_df['ResId'] == target_ResId) & (protein_df['Name'].isin(linkage["C"]))][0]

    link_pair=(ND2,len(Parr))

    if glycan == "None":
        return glycoprotein_final, Parr, box
    else:
        all_G,loaded = sampling(glycan,link)
        box = f"Residue : {target_ResId}{target_Chain}\n {glycan} has {len(all_G)} Clusters\n"
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
        box+= f"Trying cluster 0...\n"
        for (clus_num,Gi) in enumerate(all_G):
            Garr = Gi[['X','Y','Z']].to_numpy(dtype=float)
            Garr,r,phi,psi = optwithwiggle(Garr,O1,CB,CG,ND2,C1,O5,Parr,torsionpoints,torsionparts,phisd,psisd,wiggle_tries=1,wiggle_method="none")
            if r<rf:
                FGarr=Garr
                rf=r
                cluster_best = clus_num
            if r==1.0:
                clash=False
                box+= f"Clash Solved for {glycan} at residue :{target_ResId}{target_Chain} with phi : {int(phi)} and psi : {int(psi)}, cluster {clus_num}\n"
                break
            else:
                box+= f"Clash detected, trying cluster {clus_num+1}...\n"
        if rf != 1.0:
            box+= f"Clash detected, trying cluster {cluster_best} with wiggle...\n"
            Gi = all_G[cluster_best]
            Garr = Gi[['X','Y','Z']].to_numpy(dtype=float)
            Garr,r,phi,psi = optwithwiggle(Garr,O1,CB,CG,ND2,C1,O5,Parr,torsionpoints,torsionparts,phisd,psisd,wiggle_tries=40,wiggle_method="random")
            if r<rf:
                FGarr=Garr
                rf=r
            if r==1.0:
                clash=False
                box+= f"Clash Solved for {glycan} at residue :{target_ResId}{target_Chain} with phi : {int(phi)} and psi : {int(psi)}, cluster {cluster_best}\n"
                pass
            else:
                box+= f"\n"
        if clash:
            box+= f"Clash exist for {glycan} at residue :{target_ResId}{target_Chain} with phi : {int(phi)} and psi : {int(psi)}, cluster {cluster_best}\n"
        Gn =  pd.DataFrame(FGarr, columns = ['X','Y','Z'])
        G.update(Gn)
        G = G.drop([0,1])
        G["Number"] = glycoprotein_final["Number"].iloc[-1] + G["Number"] 
        G["Chain"] = glycan_chain
        # G["Chain"] = target_Chain
        glycoprotein_final= pd.concat([glycoprotein_final,G])
        Parr=glycoprotein_final[['X','Y','Z']].to_numpy(dtype=float)

        return glycoprotein_final, Parr, box , clash ,link_pair