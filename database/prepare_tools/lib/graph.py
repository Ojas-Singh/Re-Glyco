from scipy.spatial import distance
import numpy as np


def G_flatten(frame):
    G = distance.cdist(frame, frame, 'euclidean')
    r = np.arange(len(G))
    mask = r[:,None]<r
    G_flat = np.concatenate((G.T[mask],np.diag(G)))
    return G_flat
