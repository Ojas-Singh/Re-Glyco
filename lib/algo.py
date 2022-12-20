import random
import pandas as pd
import numpy as np
import pandas as pd
from lib import pdb 
import copy
from numba import njit,jit
import numba as nb
from numpy import cross, eye, dot
from scipy.linalg import expm, norm
import time
from scipy.spatial import distance
import numba_scipy
import streamlit as st
import math

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
    #p=[p0,p1,p2,p3]   all xyz numpy array of 4 points.
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
    
    # c = distance.cdist(Garr, Parr, 'euclidean')
    c = eucl_opt(Garr, Parr)
    # c = np.linalg.norm(Garr[:,None,:] - Parr, axis=2)
    con = c<2.5
    c2 = np.extract(con, c)
    for i in c2:
        r+=1/(i+.001)
    return r

def attach(protein,glycans,glycosylation_locations,attempt):
    protein_df= pdb.to_DF(protein)
    Parr=protein_df[['X','Y','Z']].to_numpy(dtype=float)
    glycoprotein_final = copy.deepcopy(protein_df)
    gly=[]
    ChainId= ["B","C","D","E","F","G","H","I"]
    k=0
    for i in glycosylation_locations:
        s=time.time()
        if not i["description"].startswith('N-linked'):
            st.write("Only N Glycosylation yet! Spot :",i["begin"]," is ",i["description"])
            continue
        target_ResId= int(i["begin"])
        st.write("Glycosylating Spot :",i["begin"])
        OD1 = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'OD1'),['Number']].iloc[0]['Number'] -1
        CG = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'CG'),['Number']].iloc[0]['Number'] -1
        ND2 = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'ND2'),['Number']].iloc[0]['Number'] -1
        G = sampling(glycans)
        C1 = G.loc[(G['ResId']==2) & (G['Name']== 'C1'),['Number']].iloc[0]['Number'] -1
        O5 = G.loc[(G['ResId']==2) & (G['Name']== 'O5'),['Number']].iloc[0]['Number'] -1
        O1 = G.loc[(G['ResId']==1) & (G['Name']== 'O1'),['Number']].iloc[0]['Number'] -1
        Garr = G[['X','Y','Z']].to_numpy(dtype=float)
        Garr = Garr-Garr[O1]
        Garr = Garr + Parr[ND2]
        axis = np.cross(Parr[CG]-Parr[ND2],Garr[C1]-Parr[ND2])
        an=fastest_angle(Parr[CG],Parr[ND2],Garr[C1])
        theta = np.radians(109.5 - an)
        Garr = Garr - Parr[ND2]
        M0 = rotation_matrix(axis, theta)
        for i in range(len(Garr)):
            Garr[i] = np.dot(M0,Garr[i])
        Garr = Garr + Parr[ND2]
        phi,psi =opt(OD1,CG,ND2,C1,O5,Garr,Parr,attempt)
        Garr = rr(phi,psi,OD1,CG,ND2,C1,O5,Garr,Parr)
        print(fastest_dihedral(Parr[OD1],Parr[CG],Parr[ND2],Garr[C1]))
        print(fastest_dihedral(Parr[CG],Parr[ND2],Garr[C1],Garr[O5]))
        Gn =  pd.DataFrame(Garr, columns = ['X','Y','Z'])
        G.update(Gn)
        G = G.drop([0,1])
        G["Number"] = glycoprotein_final["Number"].iloc[-1] + G["Number"] 
        G["Chain"] = ChainId[k]
        k+=1
        glycoprotein_final= pd.concat([glycoprotein_final,G])
        Parr=glycoprotein_final[['X','Y','Z']].to_numpy(dtype=float)
        print("exec time",time.time()-s)
    return glycoprotein_final

def sampling(Glycanid):
    if Glycanid== "bisecting":
        G = pdb.parse("data/bisecting.pdb")
    elif Glycanid== "man":
        G = pdb.parse("data/MAN6.pdb")
    return pdb.to_DF(G)

@njit(fastmath=True)
def rr(phi,psi,OD1,CG,ND2,C1,O5,Garr,Parr):
    M1 = rotation_matrix(Parr[ND2]-Parr[CG],np.radians(phi-fastest_dihedral(Parr[OD1],Parr[CG],Parr[ND2],Garr[C1])))
    Garr = Garr - Parr[ND2]
    # Garr = [np.dot(M1,x) for x in Garr]
    for i in range(len(Garr)):
        Garr[i] = np.dot(M1,Garr[i])
    Garr = Garr + Parr[ND2]
    M2 = rotation_matrix(Garr[C1]-Parr[ND2],np.radians(psi-fastest_dihedral(Parr[CG],Parr[ND2],Garr[C1],Garr[O5])))
    Garr = Garr - Parr[ND2]
    # Garr = [np.dot(M2,x) for x in Garr]
    for i in range(len(Garr)):
        Garr[i] = np.dot(M2,Garr[i])
    Garr = Garr + Parr[ND2]
    return Garr


@njit(fastmath=True)
def opt(OD1,CG,ND2,C1,O5,Garr,Parr,attempt):
    phif=0
    psif=0
    r=100000000
    # for phi in range(-180,180,10):
    #     for psi in range(-180,180,10):
    #         Garr = rr(phi,psi,OD1,CG,ND2,C1,O5,Garr,Parr)
    #         ri= steric_fast(Garr,Parr)
    #         if ri<r:
    #             phif= phi
    #             psif= psi
    #             r=ri
    for pp in range(attempt):
        phi = random.uniform(-180, 180)
        psi = random.uniform(-180, 180)
        Garr = rr(phi,psi,OD1,CG,ND2,C1,O5,Garr,Parr)
        ri= steric_fast(Garr,Parr)
        # print(pp,ri)
        if ri<r:
            phif= phi
            psif= psi
            r=ri
    return phif,psif