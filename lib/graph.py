from scipy.spatial import distance
import numpy as np


def Gfunction(r):
    return 1/(r+.01)

def G_flatten(frame):
    G = distance.cdist(frame, frame, 'euclidean')
    r = np.arange(len(G))
    mask = r[:,None]<r
    G_flat = np.concatenate((G.T[mask],np.diag(G)))
    funcc = np.vectorize(Gfunction)
    G_flat = funcc(G_flat)
    return G_flat

