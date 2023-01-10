#PDB Format from https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
import numpy as np
import pandas as pd

def to_DF(pdbddata):
    df = pd.DataFrame(data=pdbddata)
    df = df.transpose()
    df.columns = ['Number','Name','ResName','Chain','ResId','X','Y','Z','Element']
    return df

def to_normal(df):
    Number = df['Number'].tolist()
    Name = df['Name'].tolist()
    ResName = df['ResName'].tolist()
    Chain = df['Chain'].tolist()
    ResId = df['ResId'].tolist()
    X = df['X'].tolist()
    Y = df['Y'].tolist()
    Z = df['Z'].tolist()
    Element = df['Element'].tolist()
    pdbdata=[Number,Name,ResName,Chain,ResId,X,Y,Z,Element]
    return pdbdata

def parse(f):
    Number = []
    Name = []
    ResName = []
    Chain = []
    ResId = []
    X = []
    Y = []
    Z = []
    Element = []
    pdbdata=[Number,Name,ResName,Chain,ResId,X,Y,Z,Element]
    with open(f, 'r') as f:
            lines = f.readlines()
            i=1
            for line in lines:
                if line.startswith("ATOM"):
                    pdbdata[0].append(int((line[7:11]).strip(" ")))
                    pdbdata[1].append((line[12:16]).strip(" "))
                    pdbdata[2].append((line[17:20]).strip(" "))
                    pdbdata[3].append((line[20:22]).strip(" "))
                    pdbdata[4].append(int((line[22:26]).strip(" ")))
                    pdbdata[5].append(float(line[31:38]))
                    pdbdata[6].append(float(line[39:46]))
                    pdbdata[7].append(float(line[47:54]))
                    pdbdata[8].append((line[76:78]).strip(" "))
                    i+=1
                if  line.startswith("END"):
                    break
            o = len(pdbdata[0])
    return pdbdata


def exportPDB(fout,pdbdata):
    fn= open(fout,"w+")
    k=""
    for i in range(len(pdbdata[0])):
        line=list("ATOM".ljust(80))
        line[6:10] = str(pdbdata[0][i]).rjust(5) 
        line[12:15] = str(pdbdata[1][i]).ljust(4) 
        line[17:19] = str(pdbdata[2][i]).rjust(3) 
        line[20:21] = str(pdbdata[3][i]).rjust(2) 
        line[22:25] = str(pdbdata[4][i]).rjust(4) 
        line[30:37] = str('{:0.3f}'.format(pdbdata[5][i])).rjust(8) 
        line[38:45] = str('{:0.3f}'.format(pdbdata[6][i])).rjust(8) 
        line[46:53] = str('{:0.3f}'.format(pdbdata[7][i])).rjust(8) 
        line[75:77] = str(pdbdata[8][i]).rjust(3) 
        line= ''.join(line)
        fn.write(line+"\n")
        k=k+line+"\n"
    return k
                
def multi(f):
    frames=[]
    pdbdata = parse(f)
    with open(f, 'r') as f:
            lines = f.readlines()
            mat = np.zeros((len(pdbdata[0]),3))
            j=1
            i=0
            for line in lines:
                if line.startswith("ATOM"):
                    mat[i,0]=float(line[31:38])
                    mat[i,1]=float(line[39:46])
                    mat[i,2]=float(line[47:54])
                    i+=1
                if line.startswith("ENDMDL"):
                    j+=1
                    i=0
                    frames.append(mat)
                    mat = np.zeros((len(pdbdata[0]),3))

    return pdbdata,frames
