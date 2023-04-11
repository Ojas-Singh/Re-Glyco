import random
import pandas as pd
import numpy as np
from lib import pdb 
import copy
from numba import njit
import numba as nb
import streamlit as st
import math
import os,time
import config
import glycors

@njit(fastmath=True,parallel=True)
def eucl_opt(A, B):
    assert A.shape[1]==B.shape[1]
    C=np.empty((A.shape[0],B.shape[0]),A.dtype)
    I_BLK=32
    J_BLK=32
    
    #workaround to get the right datatype for acc
    init_val_arr=np.zeros(1,A.dtype)
    init_val=init_val_arr[0]
    
    #Blocking and partial unrolling
    #Beneficial if the second dimension is large -> computationally bound problem 
    # 
    for ii in nb.prange(A.shape[0]//I_BLK):
        for jj in range(B.shape[0]//J_BLK):
            for i in range(I_BLK//4):
                for j in range(J_BLK//2):
                    acc_0=init_val
                    acc_1=init_val
                    acc_2=init_val
                    acc_3=init_val
                    acc_4=init_val
                    acc_5=init_val
                    acc_6=init_val
                    acc_7=init_val
                    for k in range(A.shape[1]):
                        acc_0+=(A[ii*I_BLK+i*4+0,k] - B[jj*J_BLK+j*2+0,k])**2
                        acc_1+=(A[ii*I_BLK+i*4+0,k] - B[jj*J_BLK+j*2+1,k])**2
                        acc_2+=(A[ii*I_BLK+i*4+1,k] - B[jj*J_BLK+j*2+0,k])**2
                        acc_3+=(A[ii*I_BLK+i*4+1,k] - B[jj*J_BLK+j*2+1,k])**2
                        acc_4+=(A[ii*I_BLK+i*4+2,k] - B[jj*J_BLK+j*2+0,k])**2
                        acc_5+=(A[ii*I_BLK+i*4+2,k] - B[jj*J_BLK+j*2+1,k])**2
                        acc_6+=(A[ii*I_BLK+i*4+3,k] - B[jj*J_BLK+j*2+0,k])**2
                        acc_7+=(A[ii*I_BLK+i*4+3,k] - B[jj*J_BLK+j*2+1,k])**2
                    C[ii*I_BLK+i*4+0,jj*J_BLK+j*2+0]=np.sqrt(acc_0)
                    C[ii*I_BLK+i*4+0,jj*J_BLK+j*2+1]=np.sqrt(acc_1)
                    C[ii*I_BLK+i*4+1,jj*J_BLK+j*2+0]=np.sqrt(acc_2)
                    C[ii*I_BLK+i*4+1,jj*J_BLK+j*2+1]=np.sqrt(acc_3)
                    C[ii*I_BLK+i*4+2,jj*J_BLK+j*2+0]=np.sqrt(acc_4)
                    C[ii*I_BLK+i*4+2,jj*J_BLK+j*2+1]=np.sqrt(acc_5)
                    C[ii*I_BLK+i*4+3,jj*J_BLK+j*2+0]=np.sqrt(acc_6)
                    C[ii*I_BLK+i*4+3,jj*J_BLK+j*2+1]=np.sqrt(acc_7)
        #Remainder j
        for i in range(I_BLK):
            for j in range((B.shape[0]//J_BLK)*J_BLK,B.shape[0]):
                acc_0=init_val
                for k in range(A.shape[1]):
                    acc_0+=(A[ii*I_BLK+i,k] - B[j,k])**2
                C[ii*I_BLK+i,j]=np.sqrt(acc_0)
    
    #Remainder i
    for i in range((A.shape[0]//I_BLK)*I_BLK,A.shape[0]):
        for j in range(B.shape[0]):
            acc_0=init_val
            for k in range(A.shape[1]):
                acc_0+=(A[i,k] - B[j,k])**2
            C[i,j]=np.sqrt(acc_0)
            
    return C

@njit(fastmath=True)
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

@njit(fastmath=True)
def fastest_angle(p0,p1,p2):
    cosine_angle = np.dot(np.subtract(p0,p1),np.subtract(p2,p1)) / (np.linalg.norm(np.subtract(p0,p1)) * np.linalg.norm(np.subtract(p2,p1)))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)

@njit(fastmath=True)
def fastest_dihedral(p0,p1,p2,p3):
    b1 = p2 - p1
    b0, b1, b2 = -(p1 - p0), b1 / np.sqrt((b1 * b1).sum()), p3 - p2
    v = b0 - (b0[0] * b1[0] + b0[1] * b1[1] + b0[2] * b1[2]) * b1
    w = b2 - (b2[0] * b1[0] + b2[1] * b1[1] + b2[2] * b1[2]) * b1
    x = v[0] * w[0] + v[1] * w[1] + v[2] * w[2]
    y = (b1[1]*v[2] - b1[2]*v[1]) * w[0] + \
        (b1[2]*v[0] - b1[0]*v[2]) * w[1] + \
        (b1[0]*v[1] - b1[1]*v[0]) * w[2]
    return 180 * np.arctan2(y, x) / np.pi

@njit(fastmath=True)
def steric_fast(Garr,Parr):
    r=0
    c = eucl_opt(Garr, Parr)
    c = c[3:][:]
    con1 = c<1.7
    c2 = np.extract(con1, c)
    for i in c2:
        r+=200*np.exp(-((i)**2))
    return r

@njit(fastmath=True)
def Garrfromtorsion(Garr,torsionpoints,torsionparts):
    for i in range(len(torsionpoints)):
            M1 = rotation_matrix(Garr[torsionpoints[i][2]]-Garr[torsionpoints[i][1]],np.radians(random.uniform(-10, 10)))
            Garr = Garr-Garr[torsionpoints[i][1]]
            for j in np.where(torsionparts[i])[0]:
                Garr[j] = np.dot(M1,Garr[j])
            Garr = Garr+Garr[torsionpoints[i][1]]
    return Garr

@njit(fastmath=True)
def Garrfromtorsiondemo(Garr,torsionpoints,torsionrange,torsionparts):
    for i in range(len(torsionpoints)):
            M1 = rotation_matrix(Garr[torsionpoints[i][2]]-Garr[torsionpoints[i][1]],np.radians(random.uniform(-1*torsionrange, torsionrange)))
            Garr = Garr-Garr[torsionpoints[i][1]]
            for j in np.where(torsionparts[i])[0]:
                Garr[j] = np.dot(M1,Garr[j])
            Garr = Garr+Garr[torsionpoints[i][1]]
    return Garr


@njit(fastmath=True)
def rr(phi,psi,CB,CG,ND2,C1,O5,Garr,Parr):
    M1 = rotation_matrix(Parr[ND2]-Parr[CG],np.radians(phi-fastest_dihedral(Parr[CB],Parr[CG],Parr[ND2],Garr[C1])))
    Garr = Garr - Parr[ND2]
    for i in range(len(Garr)):
        Garr[i] = np.dot(M1,Garr[i])
    Garr = Garr + Parr[ND2]
    M2 = rotation_matrix(Garr[C1]-Parr[ND2],np.radians(psi-fastest_dihedral(Parr[CG],Parr[ND2],Garr[C1],Garr[O5])))
    Garr = Garr - Parr[ND2]
    for i in range(len(Garr)):
        Garr[i] = np.dot(M2,Garr[i])
    Garr = Garr + Parr[ND2]
    return Garr


def get_number_after_underscore(filename):
    return float(filename.split("_")[1].split(".")[0])

def sampling(Glycanid):
    folder_path = config.data_dir+Glycanid+"/clusters/beta/"
    filenames = os.listdir(folder_path)
    pdb_files = [filename for filename in filenames if filename.endswith(".pdb")]
    sorted_pdb_files = sorted(pdb_files, key=get_number_after_underscore,reverse=True)
    all_pdbs=[]
    for i in sorted_pdb_files:
        all_pdbs.append(pdb.to_DF(pdb.parse(folder_path+i)))
    loaded = np.load(config.data_dir+Glycanid+"/output/torparts.npz",allow_pickle=True)
    return all_pdbs,loaded

def attach(protein,glycans,glycosylation_locations,phisd,psisd):
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
        CB = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'CB'),['Number']].iloc[0]['Number'] -1
        CG = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'CG'),['Number']].iloc[0]['Number'] -1
        ND2 = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'ND2'),['Number']].iloc[0]['Number'] -1
        all_G,loaded = sampling(glycans[i])
        G = all_G[0]
        C1 = G.loc[(G['ResId']==2) & (G['Name']== 'C1'),['Number']].iloc[0]['Number'] -1
        O5 = G.loc[(G['ResId']==2) & (G['Name']== 'O5'),['Number']].iloc[0]['Number'] -1
        O1 = G.loc[(G['ResId']==1) & (G['Name']== 'O1'),['Number']].iloc[0]['Number'] -1
        Garr = G[['X','Y','Z']].to_numpy(dtype=float)
        torsionpoints = loaded["a"]
        torsionparts  = loaded["b"]
        clash= False
        # Garr = optwithwiggle(Garr,O1,CB,CG,ND2,C1,O5,Parr,torsionpoints,torsionparts,phisd,psisd)
        # c = eucl_opt(Garr, Parr)
        # c = c[3:][:]
        # con1 = c<1.6
        # c2 = np.extract(con1, c)
        # if np.sum(c2)> 0:
        #     print(np.sum(c2))
        #     clash=True     
        s=time.time()

        FGarr= Garr
        rf=100000000
        for Gi in all_G:
            Garr = Gi[['X','Y','Z']].to_numpy(dtype=float)
            Garr,r = optwithwiggle(Garr,O1,CB,CG,ND2,C1,O5,Parr,torsionpoints,torsionparts,phisd,psisd)
            if r<rf:
                FGarr=Garr
                rf=r
            if r==1.0:
                print("wow!")
                st.write("Solved!")
                break
            else:
                st.write("Trying different cluster...")
        timer.append(time.time()-s)
        st.write("Spot : ",target_ResId," Psi : ",int(fastest_dihedral(Parr[CB],Parr[CG],Parr[ND2],FGarr[C1]))," Phi : ",int(fastest_dihedral(Parr[CG],Parr[ND2],FGarr[C1],FGarr[O5])))
        Gn =  pd.DataFrame(FGarr, columns = ['X','Y','Z'])
        G.update(Gn)
        G = G.drop([0,1])
        G["Number"] = glycoprotein_final["Number"].iloc[-1] + G["Number"] 
        G["Chain"] = ChainId[k]
        k+=1
        glycoprotein_final= pd.concat([glycoprotein_final,G])
        Parr=glycoprotein_final[['X','Y','Z']].to_numpy(dtype=float)

    return glycoprotein_final,clash


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
            M0 = rotation_matrix(axis, theta)
            for i in range(len(Garr)):
                Garr[i] = np.dot(M0,Garr[i])
            Garr = Garr + Parr[ND2]
            # phi,psi,ri = glycors.opt(CB,CG,ND2,C1,O5,Garr,Parr,tuple(phisd),tuple(psisd))
            phi,psi,ri = glycors.opt_genetic(CB,CG,ND2,C1,O5,Garr,Parr,psisd,phisd)
            if ri<r:
                GarrF= Garr
                phiF=phi
                psiF=psi
                r=ri
            if ri==1.0:
                break
        return rr(phiF,psiF,CB,CG,ND2,C1,O5,GarrF,Parr),r

