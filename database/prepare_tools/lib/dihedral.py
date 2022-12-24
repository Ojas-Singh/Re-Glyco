import numpy as np
from numba import njit
import pandas as pd
import scipy.stats as st
import matplotlib.pyplot as plt
import pdb



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

def psiphi(A,B,frames,pdbdata):
    df= pdb.to_DF(pdbdata)
    #psi = C1 O Cx' Cx+1'
    #phi = O5 C1 O Cx'

    # 1-4 aplha
    a=df.loc[(df['ResId']==A) & (df['name']== 'C1'),['Number']].iloc[0]['Number'] -1
    b=df.loc[(df['ResId']==A) & (df['name']== 'O4'),['Number']].iloc[0]['Number'] -1
    c=df.loc[(df['ResId']==B) & (df['name']== 'C4'),['Number']].iloc[0]['Number'] -1
    d=df.loc[(df['ResId']==B) & (df['name']== 'C5'),['Number']].iloc[0]['Number'] -1
    e=df.loc[(df['ResId']==A) & (df['name']== 'O5'),['Number']].iloc[0]['Number'] -1
    f=df.loc[(df['ResId']==A) & (df['name']== 'C1'),['Number']].iloc[0]['Number'] -1
    g=df.loc[(df['ResId']==A) & (df['name']== 'O4'),['Number']].iloc[0]['Number'] -1
    h=df.loc[(df['ResId']==B) & (df['name']== 'C4'),['Number']].iloc[0]['Number'] -1
    psi=[]
    phi=[]
    for i in range(len(frames)):
        psi.append(fastest_dihedral([frames[i][a],frames[i][b],frames[i][c],frames[i][d]]))
        phi.append(fastest_dihedral([frames[i][e],frames[i][f],frames[i][g],frames[i][h]]))
    return psi,phi

# def psiphiomega(A,B,frames,pdbdata):
#     #psi = C1 O Cx' Cx+1'
#     #phi = O5 C1 O Cx'
#     #omega = O C6' C5' C4'
#     psi=[]
#     phi=[]
#     omega=[]
#     for i in range(len(frames)):
#         psi.append(fastest_dihedral(p))
#         phi.append(fastest_dihedral(q))
#         omega.append(fastest_dihedral(r))
#     return


def angle(l,frames,pdbdata):
    rf=np.zeros((len(frames),2*len(l)))
    df= pdb.to_DF(pdbdata)
    S=[]
    k=0
    for i in l:
        #psi = C1 O Cx' Cx+1'
        #phi = O5 C1 O Cx'
        # 1-4 aplha

        a=df.loc[(df['ResId']==i[0]) & (df['name']== 'C1'),['Number']].iloc[0]['Number'] -1
        b=df.loc[(df['ResId']==i[0]) & (df['name']== 'O4'),['Number']].iloc[0]['Number'] -1
        c=df.loc[(df['ResId']==i[1]) & (df['name']== 'C4'),['Number']].iloc[0]['Number'] -1
        d=df.loc[(df['ResId']==i[1]) & (df['name']== 'C5'),['Number']].iloc[0]['Number'] -1
        e=df.loc[(df['ResId']==i[0]) & (df['name']== 'O5'),['Number']].iloc[0]['Number'] -1
        f=df.loc[(df['ResId']==i[0]) & (df['name']== 'C1'),['Number']].iloc[0]['Number'] -1
        g=df.loc[(df['ResId']==i[0]) & (df['name']== 'O4'),['Number']].iloc[0]['Number'] -1
        h=df.loc[(df['ResId']==i[1]) & (df['name']== 'C4'),['Number']].iloc[0]['Number'] -1
        for j in range(len(frames)):
            rf[j][k]=fastest_dihedral([frames[j][a],frames[j][b],frames[j][c],frames[j][d]])
            rf[j][k+1]=fastest_dihedral([frames[j][e],frames[j][f],frames[j][g],frames[j][h]])
        k+=2
        r1=df.loc[(df['ResId']==i[0]),['ResName']].iloc[0]['ResName']
        r2=df.loc[(df['ResId']==i[1]),['ResName']].iloc[0]['ResName']
        S.append(str(r1)+"_"+str(r2)+"_"+"psi")
        S.append(str(r1)+"_"+str(r2)+"_"+"phi")
    return rf,S

def plot(df,name):
    col=list(df.columns)
    p=len(col)
    df.iloc[:, 3]
    for i in range(1,p,2):
        psi = df.iloc[:, i]
        phi = df.iloc[:, i+1]
        x=psi
        y=phi
        xmin, xmax = -180, 180
        ymin, ymax = -180, 180

        # Peform the kernel density estimate
        xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        values = np.vstack([x, y])
        kernel = st.gaussian_kde(values)
        f = np.reshape(kernel(positions).T, xx.shape)

        fig = plt.figure()
        ax = fig.gca()
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        # Contourf plot
        cfset = ax.contourf(xx, yy, f, cmap='Blues')
        ## Or kernel density estimate plot instead of the contourf plot
        #ax.imshow(np.rot90(f), cmap='Blues', extent=[xmin, xmax, ymin, ymax])
        # Contour plot
        cset = ax.contour(xx, yy, f, colors='k')
        # Label plot
        # ax.clabel(cset, inline=1, fontsize=10)
        ax.set_xlabel('psi')
        ax.set_ylabel('phi')
        ax.set_title(str(col[i]).strip("psi"))
        plt.savefig(str(col[i]).strip("psi")+'kde.png')

    # plt.show()

def kdemax(df,name):
    col=list(df.columns)
    p=len(col)
    df.iloc[:, 3]
    for i in range(1,p,2):
        psi = df.iloc[:, i]
        phi = df.iloc[:, i+1]
        x=psi
        y=phi
        xmin, xmax = -180, 180
        ymin, ymax = -180, 180

        # Peform the kernel density estimate
        xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        values = np.vstack([x, y])
        kernel = st.gaussian_kde(values)
        f = np.reshape(kernel(positions).T, xx.shape)

        fig = plt.figure()
        ax = fig.gca()
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        # Contourf plot
        cfset = ax.contourf(xx, yy, f, cmap='Blues')
        ## Or kernel density estimate plot instead of the contourf plot
        #ax.imshow(np.rot90(f), cmap='Blues', extent=[xmin, xmax, ymin, ymax])
        # Contour plot
        cset = ax.contour(xx, yy, f, colors='k')
        # Label plot
        # ax.clabel(cset, inline=1, fontsize=10)
        ax.set_xlabel('psi')
        ax.set_ylabel('phi')
        ax.set_title(str(col[i]).strip("psi"))
        plt.savefig(str(col[i]).strip("psi")+'kde.png')

    # plt.show()