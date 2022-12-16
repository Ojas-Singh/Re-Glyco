import random
import pandas as pd
import numpy as np
import pandas as pd
from lib import pdb 
import copy
from numba import njit
from numpy import cross, eye, dot
from scipy.linalg import expm, norm



import math

@njit()
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





@njit()
def fastest_angle(p0,p1,p2):
    cosine_angle = np.dot(np.subtract(p0,p1),np.subtract(p2,p1)) / (np.linalg.norm(np.subtract(p0,p1)) * np.linalg.norm(np.subtract(p2,p1)))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)


@njit()
def fastest_dihedral(p):
    #p=[p0,p1,p2,p3]   all xyz numpy array of 4 points.
    b1 = p[2] - p[1]
    b0, b1, b2 = -(p[1] - p[0]), b1 / np.sqrt((b1 * b1).sum()), p[3] - p[2]
    v = b0 - (b0[0] * b1[0] + b0[1] * b1[1] + b0[2] * b1[2]) * b1
    w = b2 - (b2[0] * b1[0] + b2[1] * b1[1] + b2[2] * b1[2]) * b1
    x = v[0] * w[0] + v[1] * w[1] + v[2] * w[2]
    y = (b1[1]*v[2] - b1[2]*v[1]) * w[0] + \
        (b1[2]*v[0] - b1[0]*v[2]) * w[1] + \
        (b1[0]*v[1] - b1[1]*v[0]) * w[2]
    return 180 * np.arctan2(y, x) / np.pi


# Chain is P for protien And G for Test Glycan
@njit()
def steric_score_function(Garr,Parr):
    r=0
   
    for i in Garr:
        for j in Parr:
            if 0<np.linalg.norm(i - j)<2:
                r+= 1/np.linalg.norm(i - j)
            elif np.linalg.norm(i - j)==0:
                r+= 1000000
            else:
                r+=0

    return r

def attach(protein,glycans,glycosylation_locations):
    protein_df= pdb.to_DF(protein)
    Parr=protein_df[['X','Y','Z']].to_numpy(dtype=float)
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
        phif=0
        psif=0
        r=100000000
        Darr= copy.deepcopy(Garr)
        for pp in range(100):
            phi = random.uniform(-180, 180)
            psi = random.uniform(-180, 180)
            Darr = rr(phi,psi,OD1,CG,ND2,C1,O5,Darr,Parr)
            Aarr = np.delete(Darr, (0), axis=0)
            Aarr = np.delete(Aarr, (0), axis=0)
            ri= steric_score_function(Aarr,Parr)
            if ri<r:
                phif= phi
                psif= psi
                r=ri
            print(r)
        print("final phi,psi",phif,psif)
        Garr = rr(phif,psif,OD1,CG,ND2,C1,O5,Garr,Parr)
        print(fastest_dihedral([Parr[OD1],Parr[CG],Parr[ND2],Garr[C1]]))
        print(fastest_dihedral([Parr[CG],Parr[ND2],Garr[C1],Garr[O5]]))
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

def rr(phi,psi,OD1,CG,ND2,C1,O5,Garr,Parr):
    M1 = rotation_matrix(Parr[ND2]-Parr[CG],np.radians(phi-fastest_dihedral([Parr[OD1],Parr[CG],Parr[ND2],Garr[C1]])))
    Garr = Garr - Parr[ND2]
    for i in range(len(Garr)):
        Garr[i] = np.dot(M1,Garr[i])
    Garr = Garr + Parr[ND2]
    M2 = rotation_matrix(Garr[C1]-Parr[ND2],np.radians(psi-fastest_dihedral([Parr[CG],Parr[ND2],Garr[C1],Garr[O5]])))
    Garr = Garr - Parr[ND2]
    for i in range(len(Garr)):
        Garr[i] = np.dot(M2,Garr[i])
    Garr = Garr + Parr[ND2]
    return Garr