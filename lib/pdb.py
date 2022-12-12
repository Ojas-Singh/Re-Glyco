#PDB Format from https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html




import numpy as np

def parse(f):
    N=[]  #ATOM Number
    X=[]  # X Coordinate
    Y=[]  # Y Coordinate
    Z=[]  # Z Coordiante
    A=[]  # Atom Element
    pdbdata=[N,X,Y,Z,A]
    with open(f, 'r') as f:
            lines = f.readlines()
            l=[]
            i=1
            for line in lines:
                if line.startswith("ATOM"):
                    # pdbdata[0].append(int(line[4:11]))
                    pdbdata[0].append(i)
                    pdbdata[1].append(float(line[31:38]))
                    pdbdata[2].append(float(line[39:46]))
                    pdbdata[3].append(float(line[47:54]))
                    pdbdata[4].append((line[76:78]).strip(" "))
                    i+=1
                if  line.startswith("END"):
                    break
            o = len(pdbdata[0])
    return pdbdata

def toxyz(f):
    pdbdata=parse(f)
    f= open(f+".xyz","w+")
    num=str(len(pdbdata[1]))
    f.write(num+"\n")
    f.write("The Molecule\n")
    for i in range(len(pdbdata[1])):
        A=str(pdbdata[4][i-1])
        X='{: f}'.format(pdbdata[1][i-1])
        Y='{: f}'.format(pdbdata[2][i-1])
        Z='{: f}'.format(pdbdata[3][i-1])
        f.write("{} {} {} {}\n".format(A,X,Y,Z))


def fullparse(f):
    N=[]  #ATOM Name
    X=[]  # X Coordinate
    Y=[]  # Y Coordinate
    Z=[]  # Z Coordiante
    A=[]  # Atom Element
    ResName=[] 
    ResId=[]
    Number=[]


    pdbdata=[N,X,Y,Z,A,ResName,ResId,Number]
    with open(f, 'r') as f:
            lines = f.readlines()
            l=[]
            i=1
            for line in lines:
                if line.startswith("ATOM"):
                    # pdbdata[0].append(int(line[4:11]))
                    pdbdata[0].append((line[12:16]).strip(" "))
                    pdbdata[1].append(float(line[31:38]))
                    pdbdata[2].append(float(line[39:46]))
                    pdbdata[3].append(float(line[47:54]))
                    pdbdata[4].append((line[76:78]).strip(" "))
                    pdbdata[5].append((line[17:20]).strip(" "))
                    pdbdata[6].append(int((line[22:26]).strip(" ")))
                    pdbdata[7].append(int((line[7:11]).strip(" ")))
                    i+=1
                if  line.startswith("END"):
                    break
            o = len(pdbdata[0])
    return pdbdata


def multi(f):
    frames=[]
    pdbdata = fullparse(f)
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
                    # np.array([float(line[31:38]),float(line[39:46]),float(line[47:54])])
                    i+=1
                if line.startswith("ENDMDL"):
                    j+=1
                    i=0
                    frames.append(mat)
                    mat = np.zeros((len(pdbdata[0]),3))

    return pdbdata,frames


def connectivity_OFF(f):
    l=[]
    with open(f, 'r') as f:
            lines = f.readlines()
            i= False
            for line in lines:
                if line.startswith("!entry.CONDENSEDSEQUENCE.unit.hierarchy"):
                    break
                if i:
                    r= (line).strip("1\n").strip(" ")
                    e=r.split(" ")
                    l.append([int(x) for x in e])
                if line.startswith("!entry.CONDENSEDSEQUENCE.unit.connectivity"):
                    i=True
    
    return l
    
def out(f):
    l=[]
    with open(f, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("      A V E R A G E S   O V E R"):
                    break
                if line.startswith(" Etot"):
                    r= (line).strip("1\n").strip(" ")
                    e=r.split(" ")
                    l.append(float(e[-1]))
    
    return l

def exportPDB(f,framesid):
    frames=[]
    framesid.sort()
    for i in framesid:
        frames.append([])
    with open(f, 'r') as f:
            lines = f.readlines()
            k=0
            pp=False
            ot=[]
            
            for line in lines:
                if line.startswith("MODEL") and framesid[k]==int(line[10:-1]) and k<len(framesid):
                    pp=True
                    print("ok")

                if pp:
                    frames[k].append(line)

                if pp== True and line.startswith("ENDMDL"):
                    pp=False
                    k+=1
                if k == len(framesid):
                    break
            for i in range(len(framesid)):
                fn= open("output/"+str(framesid[i])+".pdb","w+")
                for line in frames[i]:
                    fn.write(line)
                
