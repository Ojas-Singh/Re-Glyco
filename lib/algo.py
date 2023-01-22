import random
import pandas as pd
import numpy as np
import pandas as pd
from lib import pdb 
import copy
from numba import njit,jit
import numba as nb
import time
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

# @njit(fastmath=True)


@njit(fastmath=True)
def steric_fast(Garr,Parr):
    r=0
    
    # c = distance.cdist(Garr, Parr, 'euclidean')
    c = eucl_opt(Garr, Parr)
    

    con1 = c<1.6
    # con2 = c>0.01
    # con = np.logical_and(con1,con2)
    c2 = np.extract(con1, c)
    for i in c2:
        # r+=1/(.001+(i-1.5)**2)
        r+=20000*np.exp(-((i)**2))
    # return np.sum(np.reciprocal(c2))
    return r



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
def opt(OD1,CG,ND2,C1,O5,Garr,Parr):
    phif=0
    psif=0
    r=1000000000000000
    for pp in range(2000):
        psi = np.random.normal(-97.5, 33)
        phi = np.random.normal(178, 26)
        Garr = rr(phi,psi,OD1,CG,ND2,C1,O5,Garr,Parr)
        ri= steric_fast(Garr,Parr)
        if ri<r:
            phif= phi
            psif= psi
            r=ri
    
    # for phi in range(152,204):
    #     for psi in range(-130,-64):
    #         Garr = rr(phi,psi,OD1,CG,ND2,C1,O5,Garr,Parr)
    #         ri= steric_fast(Garr,Parr)
    #         if ri<r:
    #             phif= phi
    #             psif= psi
    #             r=ri
    return phif,psif,r


def attach(protein,glycans,glycosylation_locations):
    protein_df= pdb.to_DF(protein)
    Parr=protein_df[['X','Y','Z']].to_numpy(dtype=float)
    glycoprotein_final = copy.deepcopy(protein_df)
    gly=[]
    ChainId= ["B","C","D","E","F","G","H","I"]
    k=0
    for i in range(len(glycosylation_locations)):
        if not glycosylation_locations[i]["description"].startswith('N-linked'):
            st.write("Only N Glycosylation yet! Spot :",glycosylation_locations[i]["begin"]," is ",glycosylation_locations[i]["description"])
            continue
        target_ResId= int(glycosylation_locations[i]["begin"])
        # st.write("Glycosylating Spot :",glycosylation_locations[i]["begin"])
        OD1 = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'CB'),['Number']].iloc[0]['Number'] -1
        CG = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'CG'),['Number']].iloc[0]['Number'] -1
        ND2 = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'ND2'),['Number']].iloc[0]['Number'] -1
        G,loaded = sampling(glycans[i])

        C1 = G.loc[(G['ResId']==2) & (G['Name']== 'C1'),['Number']].iloc[0]['Number'] -1
        O5 = G.loc[(G['ResId']==2) & (G['Name']== 'C2'),['Number']].iloc[0]['Number'] -1
        O1 = G.loc[(G['ResId']==1) & (G['Name']== 'O1'),['Number']].iloc[0]['Number'] -1
        Garr = G[['X','Y','Z']].to_numpy(dtype=float)
        Garr = Garr-Garr[O1]
        Garr = Garr + Parr[ND2]
        axis = np.cross(Parr[CG]-Parr[ND2],Garr[C1]-Parr[ND2])
        an=fastest_angle(Parr[CG],Parr[ND2],Garr[C1])
        theta = np.radians(125 - an)
        Garr = Garr - Parr[ND2]
        M0 = rotation_matrix(axis, theta)
        for i in range(len(Garr)):
            Garr[i] = np.dot(M0,Garr[i])
        Garr = Garr + Parr[ND2]
        phi,psi,dum =opt(OD1,CG,ND2,C1,O5,Garr,Parr)
        Garr = rr(phi,psi,OD1,CG,ND2,C1,O5,Garr,Parr)
        st.write(fastest_dihedral(Parr[OD1],Parr[CG],Parr[ND2],Garr[C1]),fastest_dihedral(Parr[CG],Parr[ND2],Garr[C1],Garr[O5]))
        Gn =  pd.DataFrame(Garr, columns = ['X','Y','Z'])
        G.update(Gn)
        G = G.drop([0,1])
        G["Number"] = glycoprotein_final["Number"].iloc[-1] + G["Number"] 
        G["Chain"] = ChainId[k]
        k+=1
        glycoprotein_final= pd.concat([glycoprotein_final,G])
        Parr=glycoprotein_final[['X','Y','Z']].to_numpy(dtype=float)
        
    return glycoprotein_final


def sampling(Glycanid):
    if Glycanid== "bisecting":
        G = pdb.parse("data/bisecting.pdb")
        loaded = np.load('data/bisecting.npz',allow_pickle=True)
    elif Glycanid== "A3G3S1-F3":
        G = pdb.parse("data/A3G3S1-F3.pdb")
        loaded = np.load('data/A3G3S1-F3.npz',allow_pickle=True)
    elif Glycanid== "a2":
        G = pdb.parse("data/Complex/a2/Cluster1.pdb")
        loaded = np.load('data/file.npz',allow_pickle=True)
    elif Glycanid== "a2g2":
        G = pdb.parse("data/Complex/a2g2/cluster1.pdb")
        loaded = np.load('data/file.npz',allow_pickle=True) 
    elif Glycanid== "a3g3":
        G = pdb.parse("data/Complex/a3g3/cluster1.pdb")
        loaded = np.load('data/file.npz',allow_pickle=True)
    elif Glycanid== "m5":
        G = pdb.parse("data/Oligomannose/man5/cluster1.pdb")
        loaded = np.load('data/file.npz',allow_pickle=True)
    elif Glycanid== "m6_1":
        G = pdb.parse("data/Oligomannose/man6_1/cluster1.pdb")
        loaded = np.load('data/file.npz',allow_pickle=True)
    elif Glycanid== "m6_2":
        G = pdb.parse("data/Oligomannose/man6_2/cluster1.pdb")
        loaded = np.load('data/file.npz',allow_pickle=True)
    elif Glycanid== "m6_3":
        G = pdb.parse("data/Oligomannose/man6_3/cluster1.pdb")
        loaded = np.load('data/file.npz',allow_pickle=True)
    elif Glycanid== "m7_1":
        G = pdb.parse("data/Oligomannose/man7_1/cluster1.pdb")
        loaded = np.load('data/file.npz',allow_pickle=True)
    elif Glycanid== "m7_2":
        G = pdb.parse("data/Oligomannose/man7_2/cluster1.pdb")
        loaded = np.load('data/file.npz',allow_pickle=True)
    elif Glycanid== "m7_3":
        G = pdb.parse("data/Oligomannose/man7_3/cluster1.pdb")
        loaded = np.load('data/file.npz',allow_pickle=True)
    elif Glycanid== "m7_4":
        G = pdb.parse("data/Oligomannose/man7_4/cluster1.pdb")
        loaded = np.load('data/file.npz',allow_pickle=True)
    elif Glycanid== "m8_1":
        G = pdb.parse("data/Oligomannose/man8_1/cluster1.pdb")
        loaded = np.load('data/file.npz',allow_pickle=True)
    elif Glycanid== "m8_2":
        G = pdb.parse("data/Oligomannose/man8_2/cluster1.pdb")
        loaded = np.load('data/file.npz',allow_pickle=True)
    elif Glycanid== "m8_3":
        G = pdb.parse("data/Oligomannose/man8_3/cluster1.pdb")
        loaded = np.load('data/file.npz',allow_pickle=True)
    elif Glycanid== "m9":
        G = pdb.parse("data/Oligomannose/man9/cluster1.pdb")
        loaded = np.load('data/file.npz',allow_pickle=True)
    return pdb.to_DF(G),loaded


def attachwithwiggle(protein,glycans,glycosylation_locations):
    protein_df= pdb.to_DF(protein)
    Parr=protein_df[['X','Y','Z']].to_numpy(dtype=float)
    glycoprotein_final = copy.deepcopy(protein_df)
    gly=[]
    ChainId= ["B","C","D","E","F","G","H","I"]
    k=0
    for i in range(len(glycosylation_locations)):
        if not glycosylation_locations[i]["description"].startswith('N-linked'):
            st.write("Only N Glycosylation yet! Spot :",glycosylation_locations[i]["begin"]," is ",glycosylation_locations[i]["description"])
            continue
        target_ResId= int(glycosylation_locations[i]["begin"])
        # st.write("Glycosylating Spot :",i["begin"])
        OD1 = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'CB'),['Number']].iloc[0]['Number'] -1
        CG = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'CG'),['Number']].iloc[0]['Number'] -1
        ND2 = protein_df.loc[(protein_df['ResId']==target_ResId) & (protein_df['Name']== 'ND2'),['Number']].iloc[0]['Number'] -1
        G,loaded = sampling(glycans[i])

        C1 = G.loc[(G['ResId']==2) & (G['Name']== 'C1'),['Number']].iloc[0]['Number'] -1
        O5 = G.loc[(G['ResId']==2) & (G['Name']== 'C2'),['Number']].iloc[0]['Number'] -1
        O1 = G.loc[(G['ResId']==1) & (G['Name']== 'O1'),['Number']].iloc[0]['Number'] -1
        
        Garr = G[['X','Y','Z']].to_numpy(dtype=float)
        tormeta = loaded["b"]
        torsions = loaded["c"]
        torsionpoints = loaded["d"]
        torsionparts  = loaded["f"]
        torsionparts = np.asarray(torsionparts)
        torsionpoints= np.asarray(torsionpoints)
        Garr = optwithwiggle(Garr,O1,OD1,CG,ND2,C1,O5,Parr,torsionpoints,torsions,torsionparts)
        st.write("Phi : ",fastest_dihedral(Parr[OD1],Parr[CG],Parr[ND2],Garr[C1]))
        st.write("Psi : ",fastest_dihedral(Parr[CG],Parr[ND2],Garr[C1],Garr[O5]))
        Gn =  pd.DataFrame(Garr, columns = ['X','Y','Z'])
        G.update(Gn)
        G = G.drop([0,1])
        G["Number"] = glycoprotein_final["Number"].iloc[-1] + G["Number"] 
        G["Chain"] = ChainId[k]
        k+=1
        glycoprotein_final= pd.concat([glycoprotein_final,G])
        Parr=glycoprotein_final[['X','Y','Z']].to_numpy(dtype=float)
        
        
    return glycoprotein_final


@njit(fastmath=True)
def optwithwiggle(GarrM,O1,OD1,CG,ND2,C1,O5,Parr,torsionpoints,torsions,torsionparts):
        r = 100000000
            
        GarrF= GarrM
        phiF=0
        psiF=0
        for i in range(100):
            Garr = Garrfromtorsion(GarrM,torsionpoints,torsions,torsionparts)
            Garr = Garr-Garr[O1]
            Garr = Garr + Parr[ND2]
            axis = np.cross(Parr[CG]-Parr[ND2],Garr[C1]-Parr[ND2])
            an=fastest_angle(Parr[CG],Parr[ND2],Garr[C1])
            theta = np.radians(125 - an)
            Garr = Garr - Parr[ND2]
            M0 = rotation_matrix(axis, theta)
            for i in range(len(Garr)):
                Garr[i] = np.dot(M0,Garr[i])
            Garr = Garr + Parr[ND2]
            phi,psi,ri =opt(OD1,CG,ND2,C1,O5,Garr,Parr)
            # Garr = rr(178,-97.5,OD1,CG,ND2,C1,O5,Garr,Parr)
            # ri = steric_fast(Garr,Parr)
            if ri<r:
                GarrF= Garr
                phiF=phi
                psiF=psi
        return rr(phiF,psiF,OD1,CG,ND2,C1,O5,GarrF,Parr)
        # return GarrF

@njit(fastmath=True)
def Garrfromtorsion(Garr,torsionpoints,torsions,torsionparts):
    # randomidx = random.randint(0,len(torsions)-1)
    # torsion = torsions[randomidx]
    for i in range(len(torsionpoints)):
            # M1 = rotation_matrix(Garr[torsionpoints[i][2]]-Garr[torsionpoints[i][1]],np.radians(torsion[i]-fastest_dihedral(Garr[torsionpoints[i][0]],Garr[torsionpoints[i][1]],Garr[torsionpoints[i][2]],Garr[torsionpoints[i][3]])))
            M1 = rotation_matrix(Garr[torsionpoints[i][2]]-Garr[torsionpoints[i][1]],np.radians(random.uniform(-25, 25)))

            Garr = Garr-Garr[torsionpoints[i][1]]
            # for j in torsionparts[i][1]:
            for j in np.where(torsionparts[i])[0]:
                Garr[j] = np.dot(M1,Garr[j])
            Garr = Garr+Garr[torsionpoints[i][1]]
    return Garr


@njit(fastmath=True)
def Garrfromtorsiondemo(Garr,torsionpoints,torsionrange,torsionparts):

    for i in range(len(torsionpoints)):
            M1 = rotation_matrix(Garr[torsionpoints[i][2]]-Garr[torsionpoints[i][1]],np.radians(random.uniform(-1*torsionrange, torsionrange)))
            Garr = Garr-Garr[torsionpoints[i][1]]
            # for j in torsionparts[i][1]:
            for j in np.where(torsionparts[i])[0]:
                Garr[j] = np.dot(M1,Garr[j])
            Garr = Garr+Garr[torsionpoints[i][1]]
    return Garr

