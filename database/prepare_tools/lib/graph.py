import numpy as np
from numba import njit


def to_G(frame):
    G = np.zeros((len(frame),len(frame)))
    for i in range(len(frame)):
        for j in range(len(frame)):
            G[i,j] = np.linalg.norm(frame[i]-frame[j])
            
    return G

def to_Gframes(frames):
    Gframes=[]
    for i in frames:
        Gframes.append(to_G(i))
    return Gframes

@njit
def to_Gflatten(frame):
    f = np.zeros((int(len(frame)*(len(frame)+1)/2)))
    k=0
    for i in range(len(frame)):
        for j in range(i,len(frame)):
            d = np.linalg.norm(frame[i]-frame[j])
            if d==0:
                f[k] = 1
            else:
                f[k] =(1/d)**(0.8)
            k+=1

    return f