import random
import pandas as pd
import numpy as np
from lib import pdb , sasa
import os,time,copy,re
import config
from config import RESIDUE_MAP
import glycors
import json
import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter
from MDAnalysis.coordinates.PDB import PDBWriter


def sample(Glycanid,linkage,frame_number,cluster):
    folder_path = config.data_dir+Glycanid+"/output/"
    ending =""
    if linkage == "alpha":
        ending="alpha.npy"
    elif linkage == "beta":
        ending="beta.npy"
    frames = np.load(folder_path+ending)
    if cluster == 'all':
        # Randomly select frame_number frames
        if len(frames) < frame_number:
            # If there are fewer frames than requested, return all frames
            sampled_frames = frames
        else:
            # Randomly sample frame_number frames
            sampled_indices = np.random.choice(len(frames), frame_number, replace=False)
            sampled_frames = frames[sampled_indices]
    path = config.data_dir+Glycanid+"/PDB_format_HETATM/"+Glycanid+"_cluster0_alpha.PDB.pdb"
    pdbdf = pdb.to_DF(pdb.parse_gly(path))
    return pdbdf,sampled_frames






def attach(protein,glycans,glycosylation_locations):
    lowest_frame=100000
    total_frames_with_data =[]
    link_pairs=[]
    clashes=[]
    protein_df= pdb.to_DF(protein)
    glycoprotein_final = copy.deepcopy(protein_df)
    s = time.time()
    Box = f"Calculation started\n"
    Box += f"Will try 200 conformations first!\n"
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
        glycoprotein_final, Parr, box ,clash, link_pair, outframes, df = attach_glycan(glycoprotein_final, glycan= glycans[i], linkage= residue_data, target_ResId = target_ResId,target_Chain = target_Chain,glycan_chain=currentGlycanChain)
        frame_number = len(outframes)
        print(frame_number)
        
        if frame_number==0:
            Box += f"all conforamtion failed\n"
            
        else:
            total_frames_with_data.append((outframes,currentGlycanChain,df))
            if frame_number < lowest_frame:
                lowest_frame = frame_number
        link_pairs.append(link_pair)
        Box += box
        clashes.append(clash)
        currentGlycanChain = chr(ord(currentGlycanChain) + 1) 
    Box += f"Showing results with frames {lowest_frame}\n"
    Box += f"Finished in {time.time() -s } seconds"
    print(len(total_frames_with_data))
        
    
    return glycoprotein_final,any(clashes), Box ,link_pairs,total_frames_with_data,lowest_frame



def attach_glycan(glycoprotein_final, glycan, linkage, target_ResId,target_Chain,glycan_chain):
    
    protein_df = glycoprotein_final
    
    Parr=glycoprotein_final[['X','Y','Z']].to_numpy(dtype=float)
    
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
        target_frame_count = 50 
        total_output_frames = []
        total_attempt =0
        while len(total_output_frames) < target_frame_count and total_attempt<5:
            # Sample frames
            all_G, frames = sample(glycan, link, 200, 'all')

            G = all_G
            C1 = G.loc[(G['ResId']==2) & (G['Name']== 'C1'),['Number']].iloc[0]['Number'] -1
            O5 = G.loc[(G['ResId']==2) & (G['Name']== 'O5'),['Number']].iloc[0]['Number'] -1
            O1 = G.loc[(G['ResId']==1) & (G['Name']== 'O1'),['Number']].iloc[0]['Number'] -1
            # Find steric interactions
            outframes =glycors.find_steric_interactions(CB,CG,ND2,C1,O5,O1,frames,Parr,psisd,phisd)
            # Append new output frames to the total output frames
            total_output_frames.extend(outframes)
            total_attempt+=1

        # all_G,frames = sample(glycan,link,1000,'all')
        
        # G = all_G
        # C1 = G.loc[(G['ResId']==2) & (G['Name']== 'C1'),['Number']].iloc[0]['Number'] -1
        # O5 = G.loc[(G['ResId']==2) & (G['Name']== 'O5'),['Number']].iloc[0]['Number'] -1
        # O1 = G.loc[(G['ResId']==1) & (G['Name']== 'O1'),['Number']].iloc[0]['Number'] -1

        # outframes =glycors.find_steric_interactions(CB,CG,ND2,C1,O5,O1,frames,Parr,psisd,phisd)
            
            
        box = f"Residue : {target_ResId}{target_Chain}\n {glycan} passed {len(outframes)} ensemble\n"
        # Gn =  pd.DataFrame(outframes[0], columns = ['X','Y','Z'])
        # G.update(Gn)
        # G = G.drop([0,1])
        # G["Number"] = glycoprotein_final["Number"].iloc[-1] + G["Number"] 
        # G["Chain"] = glycan_chain
        # # G["Chain"] = target_Chain
        # glycoprotein_final= pd.concat([glycoprotein_final,G])
        # Parr=glycoprotein_final[['X','Y','Z']].to_numpy(dtype=float)
        clash = False
        
        return glycoprotein_final, Parr, box , clash ,link_pair, total_output_frames, all_G
    

def runSASA(pdb_file, xtc_file, first_frame_pdb_file,outfilename):
    # Load the multi-frame PDB file
    u = mda.Universe(pdb_file)

    # Create an XTC writer with the same number of atoms as the universe
    with XTCWriter(xtc_file, n_atoms=u.atoms.n_atoms) as W:
        # Iterate through all frames and write each frame to the XTC file
        for ts in u.trajectory:
            W.write(u)

    # Export the first frame to a separate PDB file
    with PDBWriter(first_frame_pdb_file, n_atoms=u.atoms.n_atoms) as W:
        u.trajectory[0]  # Go to the first frame
        W.write(u)

    sasa.run([first_frame_pdb_file],[xtc_file],probe=0.25,ndots=15,endtime=0,outfile=outfilename)




def sasa_molj(filename,output_sasa):
    # Attempting to parse the content as JSON
    # Read the file content again
    with open('output.molj', 'r') as file:
        content = file.read()

    # New URL to be used
    new_url = "https://glycoshape.io/output/"+filename


    # Replace the old URL with the new URL in the file content
    updated_content = content.replace("https://glycoshape.io/output/output2.pdb", new_url)

    # Write the updated content back to the file
    with open(output_sasa, 'w') as file:
        file.write(updated_content)


    
