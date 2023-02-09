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

def res2input(B,A,C,O,df):
    # a, b ResId
    # c is carbon name.
    if C == 6 :
        carbon = "C"+ str(C)
        carbonminus = "C" + str(C-1)
        oxygen = "O" + str(O)
        carbonminusminus =  "C" + str(C-2)
        #psi = C1 O Cx' Cx+1'
        #phi = O5 C1 O Cx'

        # 1-4 aplha
        a=df.loc[(df['ResId']==A) & (df['Name']== 'C1'),['Number']].iloc[0]['Number'] -1
        b=df.loc[(df['ResId']==B) & (df['Name']== oxygen),['Number']].iloc[0]['Number'] -1
        c=df.loc[(df['ResId']==B) & (df['Name']== carbon),['Number']].iloc[0]['Number'] -1
        d=df.loc[(df['ResId']==B) & (df['Name']== carbonminus),['Number']].iloc[0]['Number'] -1

        e=df.loc[(df['ResId']==A) & (df['Name']== 'O5'),['Number']].iloc[0]['Number'] -1
        f=df.loc[(df['ResId']==A) & (df['Name']== 'C1'),['Number']].iloc[0]['Number'] -1
        g=df.loc[(df['ResId']==B) & (df['Name']== oxygen),['Number']].iloc[0]['Number'] -1
        h=df.loc[(df['ResId']==B) & (df['Name']== carbon),['Number']].iloc[0]['Number'] -1

        i=df.loc[(df['ResId']==B) & (df['Name']== oxygen),['Number']].iloc[0]['Number'] -1
        j=df.loc[(df['ResId']==B) & (df['Name']== carbon),['Number']].iloc[0]['Number'] -1
        k=df.loc[(df['ResId']==B) & (df['Name']== carbonminus),['Number']].iloc[0]['Number'] -1
        l=df.loc[(df['ResId']==B) & (df['Name']== carbonminusminus),['Number']].iloc[0]['Number'] -1


        # [[ResA,ResB],phi,psi,omega]
        return [[A,B],[e,f,g,h],[a,b,c,d],[i,j,k,l]]
    else:
        carbon = "C"+ str(C)
        carbonplus = "C" + str(C+1)
        oxygen = "O" + str(O)
        #psi = C1 O Cx' Cx+1'
        #phi = O5 C1 O Cx'

        # 1-4 aplha
        a=df.loc[(df['ResId']==A) & (df['Name']== 'C1'),['Number']].iloc[0]['Number'] -1
        b=df.loc[(df['ResId']==B) & (df['Name']== oxygen),['Number']].iloc[0]['Number'] -1
        c=df.loc[(df['ResId']==B) & (df['Name']== carbon),['Number']].iloc[0]['Number'] -1
        d=df.loc[(df['ResId']==B) & (df['Name']== carbonplus),['Number']].iloc[0]['Number'] -1
        
        e=df.loc[(df['ResId']==A) & (df['Name']== 'O5'),['Number']].iloc[0]['Number'] -1
        f=df.loc[(df['ResId']==A) & (df['Name']== 'C1'),['Number']].iloc[0]['Number'] -1
        g=df.loc[(df['ResId']==B) & (df['Name']== oxygen),['Number']].iloc[0]['Number'] -1
        h=df.loc[(df['ResId']==B) & (df['Name']== carbon),['Number']].iloc[0]['Number'] -1

        # [[ResA,ResB],phi,psi]
        return [[A,B],[e,f,g,h],[a,b,c,d]]
    

def torsions(torsionmeta,pdbdata,frames,f):
    torsiondata=[]
    meta=[]
    for j in torsionmeta:
        phi=[]
        psi=[]
        omega=[]
        if len(j)==3:
            for i in range(len(frames)):
                phi.append(fastest_dihedral(frames[i][j[1][0]],frames[i][j[1][1]],frames[i][j[1][2]],frames[i][j[1][3]]))
                psi.append(fastest_dihedral(frames[i][j[2][0]],frames[i][j[2][1]],frames[i][j[2][2]],frames[i][j[2][3]]))
            torsiondata.append(psi)
            meta.append("Res_"+str(j[0][0])+"_"+str(j[0][1])+"_phi")
            torsiondata.append(phi)
            meta.append("Res_"+str(j[0][0])+"_"+str(j[0][1])+"_psi")
        else:
            for i in range(len(frames)):
                phi.append(fastest_dihedral(frames[i][j[1][0]],frames[i][j[1][1]],frames[i][j[1][2]],frames[i][j[1][3]]))
                psi.append(fastest_dihedral(frames[i][j[2][0]],frames[i][j[2][1]],frames[i][j[2][2]],frames[i][j[2][3]]))
                omega.append(fastest_dihedral(frames[i][j[3][0]],frames[i][j[3][1]],frames[i][j[3][2]],frames[i][j[3][3]]))
            torsiondata.append(psi)
            meta.append("Res_"+str(j[0][0])+"_"+str(j[0][1])+"_phi")
            torsiondata.append(phi)
            meta.append("Res_"+str(j[0][0])+"_"+str(j[0][1])+"_psi")
            torsiondata.append(omega)
            meta.append("Res_"+str(j[0][0])+"_"+str(j[0][1])+"_omega")
            
    torsiondataDF = pd.DataFrame(data=torsiondata)
    torsiondataDF = torsiondataDF.transpose()
    torsiondataDF.columns = meta
    return torsiondataDF

def kdemax(df):
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

def pairtoname(pairs,df):
    names=[]
    for i in pairs:
        name=""
        a = df.loc[(df['Number']==i[0]),['ResId']].iloc[0]['ResId']
        b = df.loc[(df['Number']==i[1]),['ResId']].iloc[0]['ResId']
        c = df.loc[(df['Number']==i[2]),['ResId']].iloc[0]['ResId']
        d = df.loc[(df['Number']==i[3]),['ResId']].iloc[0]['ResId']
        
        if not (a==b==c==d):
            name+=str(a)+"_"+str(d)+"_"
        else:
            name+=str(a)
        names.append(name)
    return names

def pairstotorsion(pairs,frames,torsion_names):
    torsiondata=[]
    for i in pairs:
        torsion = []
        for j in frames:
            torsion.append(fastest_dihedral(j[int(i[0]-1)],j[int(i[1]-1)],j[int(i[2]-1)],j[int(i[3]-1)]))
        torsiondata.append(torsion)
    
    torsiondataDF = pd.DataFrame(data=torsiondata)
    torsiondataDF = torsiondataDF.transpose()
    torsiondataDF.columns = torsion_names
    return torsiondataDF