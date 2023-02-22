from scipy.spatial import distance
import numpy as np
import networkx as nx


def Gfunction(r):
    # return 1/(r+.01)
    return r

def G_flatten(frame):
    G = distance.cdist(frame, frame, 'euclidean')
    r = np.arange(len(G))
    mask = r[:,None]<r
    G_flat = np.concatenate((G.T[mask],np.diag(G)))
    G_flat[G_flat > 5] = 0
    funcc = np.vectorize(Gfunction)
    G_flat = funcc(G_flat)
    return G_flat


