import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from lib import graph
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as stat
from pylab import cm
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans,SpectralCoclustering,SpectralClustering,DBSCAN,MiniBatchKMeans,OPTICS

def normalizetorsion(df):
    tor = df.loc[:, df.columns!='i'].to_numpy()
    tor=tor.T
    for angle in range(len(tor)): 
        x = 0
        y = 0
        for i in tor[angle]:
            x += np.cos(np.deg2rad(i))
            y += np.sin(np.deg2rad(i))
        average_angle = np.arctan2(y, x)
        # print(average_angle)
        for i in range(len(tor[angle])):
            tor[angle][i]=np.arctan(np.tan(np.deg2rad(tor[angle][i])-average_angle))
    tor=tor.T
    return tor

def pcawithT(tor,dim):
    pca = PCA(n_components=dim)
    t = pca.fit_transform(tor)
    return pd.DataFrame(t)

def tsnewithT(tor,dim):
    t= TSNE(n_components=dim, perplexity=30.0, early_exaggeration=12.0, learning_rate='auto', n_iter=5000,n_iter_without_progress=300, min_grad_norm=1e-07, metric='euclidean', metric_params=None, init='pca', verbose=0, random_state=None, method='barnes_hut', angle=0.5, n_jobs=-1).fit_transform(tor)
    return pd.DataFrame(t)

def pcawithG(frames,idx_noH,dim):
    G = np.zeros((len(frames),int(len(frames[0][np.asarray(idx_noH,dtype=int)])*(len(frames[0][np.asarray(idx_noH,dtype=int)])+1)/2)))
    for i in range(len(frames)):
        G[i]= graph.G_flatten(frames[i][np.asarray(idx_noH,dtype=int)])
    pca = PCA(n_components=dim)
    t = pca.fit_transform(G)
    PCA_components = pd.DataFrame(t)
    return PCA_components

def pcawithGNB(frames,graph,dim):
    G = np.zeros((len(frames),int(len(frames[0][np.asarray(idx_noH,dtype=int)])*(len(frames[0][np.asarray(idx_noH,dtype=int)])+1)/2)))
    for i in range(len(frames)):
        G[i]= graph.G_flatten(frames[i][np.asarray(idx_noH,dtype=int)])
    pca = PCA(n_components=dim)
    t = pca.fit_transform(G)
    PCA_components = pd.DataFrame(t)
    return PCA_components



def tsnewithG(frames,idx_noH,dim):
    G = np.zeros((len(frames),int(len(frames[0][np.asarray(idx_noH,dtype=int)])*(len(frames[0][np.asarray(idx_noH,dtype=int)])+1)/2)))
    for i in range(len(frames)):
        G[i]= graph.G_flatten(frames[i][np.asarray(idx_noH,dtype=int)])
    t= TSNE(n_components=dim,  perplexity=30.0, early_exaggeration=12.0, learning_rate='auto', n_iter=5000,n_iter_without_progress=300, min_grad_norm=1e-07, metric='euclidean', metric_params=None, init='pca', verbose=0, random_state=None, method='barnes_hut', angle=0.5, n_jobs=-1).fit_transform(G)
    return pd.DataFrame(t)



def findmaxima(f):
    f = -1*f
    a = np.pad(f, (1, 1), mode='constant',
            constant_values=(np.amax(f), np.amax(f)))
    loc_min = []
    rows = a.shape[0]
    cols = a.shape[1]
    for ix in range(0, rows - 1):
        for iy in range(0, cols - 1):
                    if (a[ix, iy] < a[ix, iy + 1]
                        and a[ix, iy] < a[ix, iy - 1]
                        and a[ix, iy] < a[ix + 1, iy]
                        and a[ix, iy] < a[ix + 1, iy - 1]
                        and a[ix, iy] < a[ix + 1, iy + 1]
                        and a[ix, iy] < a[ix - 1, iy]
                        and a[ix, iy] < a[ix - 1, iy - 1]
                        and a[ix, iy] < a[ix - 1, iy + 1]):
                        temp_pos = (ix-1, iy-1)
                        loc_min.append(temp_pos)
    return loc_min

def nux(x,d,xmin):
    # return xmin + d*x/100
    return (x-xmin)*100/d

def nuy(y,d,ymin):
    # return ymin + d*y/100
    return (y-ymin)*100/d

def nx(x,d,xmin):
    return xmin + d*x/100
    # return (x-xmin)*100/d

def ny(y,d,ymin):
    return ymin + d*y/100
    # return (y-ymin)*100/d

def filterlow(data,k):
    x=data["0"].to_numpy()
    y=data["1"].to_numpy()
    s=pd.DataFrame([x,y])
    s=s.transpose()
    xmin, xmax = min(x), max(x)
    ymin, ymax = min(y), max(y)
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = stat.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    l=[]
    for i in range(len(x)):
        l.append([f[int(nux(x[i],np.abs(xmax-xmin),xmin))-1][int(nuy(y[i],np.abs(ymax-ymin),ymin))-1],i])
    l.sort()
    idx_top=np.ones(len(x),dtype=bool) 
    idx_bottom = np.zeros(len(x),dtype=bool) 
    for i in range(int(k*len(x))):
        idx_top[l[i][1]] = False
        idx_bottom[l[i][1]] = True
    x= np.asarray(x)
    y= np.asarray(y)
    fig = plt.figure()
    ax = fig.gca()
    cfset = ax.contourf(xx, yy, f, cmap='Blues')
    ax.scatter(x[idx_top],y[idx_top],color="#78517C",s=.2)
    ax.scatter(x[idx_bottom],y[idx_bottom],color="#F65058FF",s=.2)
    ax.set_title("Conformation Filter (>10%)")
    plt.savefig('output/PCA_filter.png',dpi=450)
    # plt.show()
    # plt.clf()
    loc_min = findmaxima(f)
    loc=[]
    for i in loc_min:
        loc.append([f[i[0]][i[1]],(i[0],i[1])])
    loc.sort()
    xw=[]
    yw=[]
    for i in range(len(loc)):
        xw.append(nx(loc[i][1][0],np.abs(xmax-xmin),xmin))
        yw.append(ny(loc[i][1][1],np.abs(ymax-ymin),ymin))
    ini=[]
    numofcluster=len(loc)
    print(numofcluster)
    for i in range(numofcluster):
        ini.append([xw[i],yw[i]])
    popp=[]
    for i in ini:
        o=[]
        for j in range(len(data.iloc[:,0])):
            o.append([np.linalg.norm(np.asarray(i)-[data["0"].iloc[j],data["1"].iloc[j]]),data["i"].iloc[j]])
        o.sort()
        popp.append([i,o[0][1]])
    return idx_top,popp

def cluster(data,numofcluster,idx_top):
    from sklearn.metrics.pairwise import pairwise_distances_argmin
    s = data.loc[:, data.columns!='i'].to_numpy()
    # s = normalizetorsion(data)
    s = s[idx_top]
    # clustering = SpectralClustering(n_clusters=numofcluster).fit(s)
    clustering = MiniBatchKMeans(compute_labels=True,n_clusters=numofcluster,init='k-means++', max_iter=5000,batch_size=4).fit(s)
    # mbk_means_cluster_centers = np.sort(clustering.cluster_centers_, axis = 0)
    # clustering_labels = pairwise_distances_argmin(s, mbk_means_cluster_centers)
    # clustering =  DBSCAN(eps=8, min_samples=10).fit(s)
    # clustering = KMeans(n_clusters=numofcluster).fit(s)
    # clustering =  OPTICS(min_samples=100).fit(s)
    clustering_labels= clustering.labels_
    label = []
    k=0
    for i in idx_top:
        if i:
            label.append(clustering_labels[k])
            k+=1
        else:
            label.append(-1)

    data["cluster"] = label
    return data,label

def plot3dkde(data):
    x=data["0"].to_numpy()
    y=data["1"].to_numpy()
    # for i in range(len(data)):
    #     x.append(data.iloc[i,0])
    #     y.append(data.iloc[i,1])
    s=pd.DataFrame([x,y])
    s=s.transpose()
    xmin, xmax = min(x), max(x)
    ymin, ymax = min(y), max(y)
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = stat.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    fig = plt.figure()
    mpl.rcParams['font.family'] = 'Cambria'
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.linewidth'] = 2
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(xx, yy, f, cmap=plt.cm.YlGnBu_r)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title("PCA Space Density")
    ax.view_init(elev=-15, azim=-59)
    ax.contourf(xx, yy, f, zdir='z', offset=-0.8, cmap=plt.cm.YlGnBu_r)
    plt.show()
    # plt.savefig('output/PCA_KDE.png',dpi=450)